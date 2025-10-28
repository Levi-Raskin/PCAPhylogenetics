#include "Eigen/Dense"
#include "ReadTSV.hpp"
#include "SimShapeLDDMM.hpp"
#include "Tree.hpp"
#include "Utility.hpp"

#include <iostream>
#include <string>
#include <vector>

int main(int argc, const char* argv[]) {
    //Read in 1,000 trees sampled from posterior distribution
    ReadTSV r = ReadTSV("PATH/PCAPhylogenetics/results/Mongle_et_al_2023_RB/sampledTrees.tsv", true, true, true);

    //add each newick string to a vector called trees; used to initialize Tree class in for loop
    std::vector<std::string> trees;
    for(int i = 0; i < r.getReadStringData().size(); i++)
        trees.push_back(r.getReadStringData()[i][0]);
    
    //Parameters
    int ndataset = 100;
    double sigma = 1.0;
    std::vector<double> alphaValues = {0.1, 0.2, 0.3, 0.4};
    std::vector<int> numLM = {10, 25, 50};

    //Simulation for loop
    for(int lm : numLM){
        //Starting shapes
        Eigen::MatrixXd rootShape2D = Utility::Shapes::generateUnitCirclePoints(lm);
        Eigen::MatrixXd rootShape3D = Utility::Shapes::generateUnitSpherePoints(lm);
        
        for(int t = 0; t < trees.size(); t++){
            //create a new tree object from the newick string
            Tree tree  = Tree(trees[t]);
            std::cout << "Simulating on tree: " << std::endl;
            std::cout << tree.getNewickString() << std::endl;
            
            //2D simulations
            #pragma omp parallel for num_threads(10)
            for(int i = 0; i < ndataset; i++){
                for(double a : alphaValues){
                    SimShapeLDDMM twoD = SimShapeLDDMM(rootShape2D, &tree);
                    twoD.runSimulation(a, sigma);
                    twoD.writeShapes("PATH/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/LDDMMSimRes/twoDimensionSimTreeIndex" + std::to_string(t) +"LM" + std::to_string(lm) + "Alpha" + std::to_string(a) + "Dataset" + std::to_string(i));
                }
                 if(i % 10 == 0)
                    std::cout << "Currently on 2D lm " << lm << " iteration " << i << " for tree " << t << std::endl;
            }
            
            //3D simulations
            #pragma omp parallel for num_threads(10)
            for(int i = 0; i < ndataset; i++){
                for(double a : alphaValues){
                    SimShapeLDDMM threeD = SimShapeLDDMM(rootShape3D, &tree);
                    threeD.runSimulation(a, sigma);
                    threeD.writeShapes("PATH/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/LDDMMSimRes/threeDimensionSimTreeIndex" + std::to_string(t) +"LM" + std::to_string(lm) + "Alpha" + std::to_string(a) + "Dataset" + std::to_string(i));
                }
                 if(i % 10 == 0)
                    std::cout << "Currently on 3D lm " << lm << " iteration " << i << " for tree " << t << std::endl;
            }
        }
    }
    
    return 0;
}
