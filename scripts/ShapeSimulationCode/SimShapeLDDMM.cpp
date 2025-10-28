#include "Eigen/Dense"
#include "Msg.hpp"
#include "Node.hpp"
#include "SimShapeLDDMM.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"
#include "Tree.hpp"
#include "WriteTSV.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>

SimShapeLDDMM::SimShapeLDDMM(Eigen::MatrixXd lm, Tree* t) : rootShape(lm), tree(t), alpha(0.01), sigma(0.02), d((int)lm.cols()), n((int)lm.rows()){
    //instantiating kernal parmas alpha and sigma with reasonable values
    //j is the number of background grid items
    //assuming lm is nxd where d is the dimensionality of the landmarks and n is the # of landmarks
    //xmin/max and ymin/max are grid bounds
}

SimShapeLDDMM::SimShapeLDDMM(Eigen::MatrixXd lm) : rootShape(lm), alpha(0.01), sigma(0.02), d((int)lm.cols()), n((int)lm.rows()){
    //instantiating kernal parmas alpha and sigma with reasonable values
    //j is the number of background grid items
    //assuming lm is nxd where d is the dimensionality of the landmarks and n is the # of landmarks
    //xmin/max and ymin/max are grid bounds
}

Eigen::MatrixXd SimShapeLDDMM::calculateDiffusionMatrix(Eigen::MatrixXd* shape){
    
    //total time: ~125 microseconds/iter
    Eigen::MatrixXd r(n, n);
    for (int i = 0; i < n; i++)
     for (int j = 0; j < n; j++)
         r(i, j) = std::sqrt(1e-7 + distance(shape->row(i), shape->row(j))) / sigma; //directly from Jacky

    Eigen::MatrixXd scalarKernal(n,n);
    for(int i = 0; i < n; i++)
     for(int j = 0; j < n; j++){
         scalarKernal(i, j) = alpha * 2 * ( 3 + 3 * r(i, j) + std::pow(r(i, j), 2) ) * std::exp(-1 * r(i,j)) ;
     }

    Eigen::MatrixXd diffusionMatrix = Eigen::MatrixXd::Zero(d * n, d * n); //d is dimensionality of the landmarks, n is number of landmarks
    Eigen::MatrixXd iD = Eigen::MatrixXd::Identity(d, d);
    for (int i = 0; i < n; i++)
     for (int j = 0; j < n; j++)
         diffusionMatrix.block(d*i, d*j, d,d) = scalarKernal(i, j) * iD;
    
    return(diffusionMatrix);
}

void SimShapeLDDMM::calculateDiffusionMatrixEfficient(Eigen::MatrixXd* shape){
    efficientDiffMat = Eigen::MatrixXd::Zero(d * n, d * n);
    // Pre-compute identity matrix once
    Eigen::MatrixXd iD = Eigen::MatrixXd::Identity(d, d);

    // Single loop to compute everything at once
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            // Compute r value directly without storing intermediate matrix
            double r_val = std::sqrt(1e-7 + distance(shape->row(i), shape->row(j))) / sigma;
            
            // Compute kernel value directly
            double kernel_val = alpha * 2 * (3 + 3 * r_val + r_val * r_val) * std::exp(-r_val);
            
            // Set the block directly
            efficientDiffMat.block(d*i, d*j, d, d) = kernel_val * iD;
        }
    }
};

Eigen::MatrixXd SimShapeLDDMM::calculateDiffusionMatrixRootShape(void){
    return(calculateDiffusionMatrix(&rootShape));
}

bool doLinesIntersect(const Eigen::Vector2d &p1, const Eigen::Vector2d &p2,
                      const Eigen::Vector2d &q1, const Eigen::Vector2d &q2) {
    // Helper function to check if line segments p1p2 and q1q2 intersect.

    auto orientation = [](const Eigen::Vector2d &a, const Eigen::Vector2d &b, const Eigen::Vector2d &c) {
        // Cross product sign: >0 ccw, <0 cw, 0 colinear
        double val = (b.y() - a.y()) * (c.x() - b.x()) - (b.x() - a.x()) * (c.y() - b.y());
        if (val > 0) return 1;    // clockwise
        else if (val < 0) return 2; // counterclockwise
        else return 0;             // colinear
    };

    auto onSegment = [](const Eigen::Vector2d &a, const Eigen::Vector2d &b, const Eigen::Vector2d &c) {
        // Check if point b lies on segment ac
        return b.x() <= std::max(a.x(), c.x()) && b.x() >= std::min(a.x(), c.x()) &&
               b.y() <= std::max(a.y(), c.y()) && b.y() >= std::min(a.y(), c.y());
    };

    int o1 = orientation(p1, p2, q1);
    int o2 = orientation(p1, p2, q2);
    int o3 = orientation(q1, q2, p1);
    int o4 = orientation(q1, q2, p2);

    // General case
    if (o1 != o2 && o3 != o4) return true;

    // Special cases
    if (o1 == 0 && onSegment(p1, q1, p2)) return true;
    if (o2 == 0 && onSegment(p1, q2, p2)) return true;
    if (o3 == 0 && onSegment(q1, p1, q2)) return true;
    if (o4 == 0 && onSegment(q1, p2, q2)) return true;

    return false;
}

bool doSegmentsIntersect3D(const Eigen::Vector3d& p1, const Eigen::Vector3d& p2,
                           const Eigen::Vector3d& q1, const Eigen::Vector3d& q2) {
    // Direction vectors
    Eigen::Vector3d u = p2 - p1;
    Eigen::Vector3d v = q2 - q1;
    Eigen::Vector3d w = p1 - q1;

    Eigen::Vector3d n = u.cross(v);  // normal vector to the plane formed by u and v
    double denom = n.squaredNorm();

    if (denom < 1e-10) {
        // Segments are parallel (and possibly colinear)
        Eigen::Vector3d w2 = p2 - q1;
        Eigen::Vector3d v_norm = v.normalized();

        // Project p1 and p2 onto q1-q2 and check for overlap
        double t0 = v_norm.dot(p1 - q1);
        double t1 = v_norm.dot(p2 - q1);
        double t_min = std::min(t0, t1);
        double t_max = std::max(t0, t1);

        double seg_len = (q2 - q1).norm();
        return !(t_max < 0 || t_min > seg_len);
    }

    // Check if lines intersect and compute the closest points
    double a = u.dot(u);
    double b = u.dot(v);
    double c = v.dot(v);
    double d = u.dot(w);
    double e = v.dot(w);
    double D = a * c - b * b;

    if (std::abs(D) < 1e-10) return false;  // lines are almost parallel

    double sc = (b * e - c * d) / D;
    double tc = (a * e - b * d) / D;

    // Closest points
    Eigen::Vector3d pc = p1 + sc * u;
    Eigen::Vector3d qc = q1 + tc * v;

    // Check if closest points are close enough (tolerance)
    if ((pc - qc).norm() > 1e-6) return false;

    // Check if sc, tc are in [0,1] (segment bounds)
    return (sc >= 0.0 && sc <= 1.0 && tc >= 0.0 && tc <= 1.0);
}

bool SimShapeLDDMM::checkLandmarksCrossed(Eigen::MatrixXd s) {
    if(d == 2){
        int n = (int)s.rows();
        if (n < 4) return false; // fewer than 4 points can’t have crossing in a closed polygon

        // Check dimension is 2 (for simplicity)
        if (s.cols() != 2) {
            throw std::invalid_argument("Landmarks dimension must be 2 for intersection checking.");
        }

        // Check all pairs of edges for intersection, skipping adjacent edges
        for (int i = 0; i < n; i++) {
            Eigen::Vector2d p1 = s.row(i);
            Eigen::Vector2d p2 = s.row((i + 1) % n);

            for (int j = i + 1; j < n; j++) {
                // Skip adjacent edges and the same edge
                if (j == i || j == (i + 1) % n || (i == 0 && j == n - 1)) {
                    continue;
                }

                Eigen::Vector2d q1 = s.row(j);
                Eigen::Vector2d q2 = s.row((j + 1) % n);

                if (doLinesIntersect(p1, p2, q1, q2)) {
                    return true;  // found an intersection
                }
            }
        }

        return false; // no intersections found
    }else{
        int n = (int)s.rows();
        int d = (int)s.cols();

        if (n < 4) return false; // fewer than 4 points can't form a self-intersecting polygon

        if (d != 2 && d != 3) {
            throw std::invalid_argument("Landmarks dimension must be 2 or 3 for intersection checking.");
        }

        for (int i = 0; i < n; i++) {
            Eigen::VectorXd p1 = s.row(i);
            Eigen::VectorXd p2 = s.row((i + 1) % n);

            for (int j = i + 1; j < n; j++) {
                if (j == i || j == (i + 1) % n || (i == 0 && j == n - 1)) {
                    continue;  // skip adjacent or same edges
                }

                Eigen::VectorXd q1 = s.row(j);
                Eigen::VectorXd q2 = s.row((j + 1) % n);

                bool intersect = false;
                if (d == 2) {
                    intersect = doLinesIntersect(p1.head<2>(), p2.head<2>(), q1.head<2>(), q2.head<2>());
                } else if (d == 3) {
                    intersect = doSegmentsIntersect3D(p1, p2, q1, q2);
                }

                if (intersect) return true;
            }
        }

        return false;

    }
}

double SimShapeLDDMM::distance(Eigen::MatrixXd x, Eigen::MatrixXd y){
    if(x.cols() > 1 && x.rows() > 1)
        Msg::error("expecting 1xn or nx1");
    if(y.cols() > 1 && y.rows() > 1)
        Msg::error("expecting 1xn or nx1");
    if(x.cols() == 1 && x.rows() > 1)
        x.transposeInPlace();
    if(y.cols() == 1 && y.rows() > 1)
        y.transposeInPlace();
    return (x.row(0)-y.row(0)).norm();
}

std::map<Node*, Eigen::MatrixXd*> SimShapeLDDMM::getTipShapes(void){
    std::vector<Node*> dpseq = tree->getDownPassSequence();
    std::map<Node*, Eigen::MatrixXd*> tipShapes;
    for(Node* n : dpseq)
        if(n->getIsTip() == true){
            tipNames.push_back(n->getName());
            tipShapes.insert({n, &(nodeShapes[n])});
            }
    return(tipShapes);
}

void SimShapeLDDMM::print(void){
    std::vector<Node*> dpseq = tree->getDownPassSequence();
    for(Node* n : dpseq){
        std::cout << "Shape at node: " << n->getIndex() <<std::endl;
        printEigenR(nodeShapes[n]);
    }
}

void SimShapeLDDMM::printEigen(Eigen::MatrixXd e){
    std::cout << std::setprecision(5);
    for(int i = 0; i < e.rows(); i++){
        for(int j = 0; j < e.cols(); j++)
            std::cout << e(i, j) << "\t";
        std::cout << std::endl;
    }
}

void SimShapeLDDMM::printEigenR(Eigen::MatrixXd e){
    std::cout << std::setprecision(5);
    for(int i = 0; i < e.rows(); i++){
        std::cout <<"c(";
        for(int j = 0; j < e.cols(); j++){
            if(j == e.cols()-1 && i < e.rows()-1){
                std::cout << e(i, j) << "),\t";
            }else if (j == e.cols()-1 && i == e.rows()-1){
                std::cout << e(i, j) << ")\t";
            }else{
                std::cout << e(i, j) << ",\t";
            }
        }
        std::cout << std::endl;
    }
}

void SimShapeLDDMM::runSimulation(double a, double s){
    alpha = a;
    sigma = s;
    
    int gotocounter = 0;
    
    simissue:
    gotocounter++;
    if(gotocounter > 100)
        Msg::error("you're really struggling on this tree champ");
    nodeShapes.clear();
    
    std::vector<Node*> dpseq = tree->getDownPassSequence();
    for (auto i = dpseq.rbegin(); i != dpseq.rend(); i++){
        Node* n = *i;
        if(n == tree->getRoot())
            nodeShapes.insert({n, rootShape});
        else{
            Node* nAnc = n->getAncestor();
            Eigen::MatrixXd previousShape = nodeShapes[nAnc];
            double bl = tree->getBranchLength(n, nAnc);
            if(bl == 0){
                nodeShapes.insert({n, previousShape});
            }else{
                try{
                    nodeShapes.insert({n, simulateShape(previousShape, bl)});
                }catch (const std::runtime_error& e) {
                    std::cout << "Caught runtime_error: " << e.what() << std::endl;
                    goto simissue;
                }
            }
            
        }
    }
}

Eigen::MatrixXd SimShapeLDDMM::simulateShape(Eigen::MatrixXd s, double t){
    std::random_device rd;
    auto seed = rd() ^ (std::chrono::high_resolution_clock::now().time_since_epoch().count());
    RandomVariable rng = RandomVariable(static_cast<uint32_t>(seed));

    Eigen::MatrixXd trueRoot = s;
    Eigen::MatrixXd currShape = s;
    Eigen::VectorXd w = Eigen::VectorXd::Zero(d*n);
    Eigen::VectorXd incrementVec = Eigen::VectorXd::Zero(d*n);
    Eigen::MatrixXd scratch(n,d);
    Eigen::MatrixXd scratch2(n,d);
    bool whileLoopExcept = false;
    int bigWhileCnt = 0;
    
    int numSteps = (int)(1000 * t);
    double timeStep = t / numSteps;
    do{
        currShape = trueRoot;
        whileLoopExcept = false;
        for(int i = 0 ; i<numSteps; i++){
            calculateDiffusionMatrixEfficient(&currShape);
            int whileCnt = 0;
            do{
                for(int j = 0; j < d*n; j++)
                    w[j] = Probability::Normal::rv(&rng, 0, std::sqrt(timeStep));
                incrementVec = efficientDiffMat * w;
                int idx = 0;
                for(int x = 0; x < n; x++){
                    for(int y = 0; y < d; y++){
                        scratch(x, y) = incrementVec[idx];
                        idx++;
                    }
                }
                scratch2 = currShape + scratch;
                whileCnt++;
                if(whileCnt > 500){
                    whileLoopExcept = true;
                    break;
                }
            }while(checkLandmarksCrossed(scratch2) == true);
            
            if(whileLoopExcept == true)
                break;
            currShape = scratch2;
        }
        bigWhileCnt++;
        if(bigWhileCnt > 10)
            throw std::runtime_error("big while loop running forever");
        if(bigWhileCnt > 2)
            Msg::warning("struggling to simulate lddmm shapes without crossing: in big while loop");
    }while(whileLoopExcept == true);

    return currShape;
}

void SimShapeLDDMM::writeShapes(std::string fp){
    std::vector<std::string> rn;
    std::vector<Node*> dpseq = tree->getDownPassSequence();

    int numLandmarks = (int)rootShape.rows();
    int numTips = tree->getNumTaxa();
    Eigen::MatrixXd plot(numLandmarks * numTips, d + 1);

    int rowIdx = 0;

    for (Node* n : dpseq) {
        if (n->getIsTip()) {
            Eigen::MatrixXd shape = nodeShapes[n];
            for (int i = 0; i < numLandmarks; i++) {
                plot(rowIdx, 0) = i + 1;           // landmark index (1-based)
                for(int j = 0; j < d; j++){
                    plot(rowIdx, j+1) = shape(i, j);     // x
                }
                rowIdx++;
                rn.push_back(n->getName());
            }
        }
    }

    WriteTSV w = WriteTSV(fp + "nodeShapes.tsv", true);
    w.addRownamesTSV(rn);
    w.appendDataTSV(plot);
    w.closeTSV();
    
}
