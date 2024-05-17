//
// Created by Liogky Alexey on 19.05.2022.
//

#include "LinearSolver.h"

LinearSolver &LinearSolver::operator=(const LinearSolver &ls) {
    if (this == &ls) return *this;
    solver = ls.solver->clone();
    return *this;
}

std::vector<std::string> LinearSolver::getAvailableSolvers() {
    std::vector<std::string> res;
    res.push_back("eigen");
    std::vector<std::string> avsolvers = getAvailableInmostSolvers();
    res.resize(res.size()+avsolvers.size());
    std::move_backward(avsolvers.begin(), avsolvers.end(), res.end());
    return res;
}

void LinearSolver::SetSolver(std::string name) {
    if (name == "eigen") {
        solver = std::make_unique<LinearSolverEigen>();
        return;
    }
    std::vector<std::string> avsolvers;
    if (!(avsolvers = getAvailableInmostSolvers()).empty()){
        auto it = std::find(avsolvers.begin(), avsolvers.end(), name);
        if (it != avsolvers.end()){
            solver = std::make_unique<LinearSolverInmost>(name);
            return;
        }
    }
    throw std::runtime_error("Solver \"" + name + "\" is not available!");
}
