#include "ForOptimizerOneLeaf.h"

///Show example of usage
int main(int argc, char* argv[]){
    SimplifiedTemplateGeometry g;
    g.D = 19, g.H = 11, g.H_f = 2;
    auto res = performOneLeafSimulation(g, TestConstantParameters(), argc, argv);
    std::cout << *reinterpret_cast<EvaluatedValues*>(&res) << std::endl;
    return 0;
}