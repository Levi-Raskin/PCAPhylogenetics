#ifndef SimShapeLDDMM_hpp
#define SimShapeLDDMM_hpp

#include "Eigen/Dense"
#include <map>
#include <string>

class Node;
class Tree;

class SimShapeLDDMM{
    public:
                                            SimShapeLDDMM(void) = delete;
                                            SimShapeLDDMM(Eigen::MatrixXd lm, Tree* t);
                                            SimShapeLDDMM(Eigen::MatrixXd lm);
        std::vector<std::string>            getTipNames(void) { return tipNames;}
        std::map<Node*, Eigen::MatrixXd*>   getTipShapes(void);
        void                                print(void);
        void                                runSimulation(double a, double s);
        void                                writeShapes(std::string fp);
    private:
        Eigen::MatrixXd                     calculateDiffusionMatrix(Eigen::MatrixXd* shape);
        void                                calculateDiffusionMatrixEfficient(Eigen::MatrixXd* shape);
        Eigen::MatrixXd                     calculateDiffusionMatrixRootShape(void);
        bool                                checkLandmarksCrossed(Eigen::MatrixXd s);
        double                              distance(Eigen::MatrixXd x, Eigen::MatrixXd y);
        void                                printEigen(Eigen::MatrixXd e);
        void                                printEigenR(Eigen::MatrixXd e);
        Eigen::MatrixXd                     simulateShape(Eigen::MatrixXd s, double t);
        Eigen::MatrixXd                     simulateShapeParallel(Eigen::MatrixXd s, double t);
        double                              alpha; //kernal parameter
        int                                 d; //dimension of landmarks
        Eigen::MatrixXd                     efficientDiffMat; // used in calculateDiffusionMatrixEfficient
        int                                 n; //number of landmarks
        Eigen::MatrixXd                     rootShape; //starting shape
        double                              sigma; //kernal parameter
        std::map<Node*, Eigen::MatrixXd>    nodeShapes;
        std::vector<std::string>            tipNames;
        Tree*                               tree;
    
};

#endif /* SimShapeLDDMM_hpp */
