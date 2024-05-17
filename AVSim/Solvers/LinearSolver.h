//
// Created by alex on 10.11.2020.
//

#ifndef AORTIC_VALVE_LINEARSOLVER_H
#define AORTIC_VALVE_LINEARSOLVER_H

#include "LinearSolverInterface.h"
#include "LinearSolverEigen.h"
#include "LinearSolverInmost.h"

class LinearSolver{
    std::unique_ptr<LinearSolverBase> solver;
public:
    LinearSolver(const LinearSolver& ls): solver{ls.solver->clone()} {}
    LinearSolver(LinearSolver&& ls) = default;
    LinearSolver& operator=(const LinearSolver& ls);
    LinearSolver& operator=(LinearSolver&& ls) = default;
    static std::vector<std::string> getAvailableSolvers();
    void SetSolver(std::string name);
    LinearSolver(){
        auto solvs = getAvailableSolvers();
        if (!solvs.empty())
            SetSolver(solvs[0]);
    };
    LinearSolver(std::string name){ SetSolver(name); }
    bool SetMatrix(SparseMatrix& sm, bool ModifiedPattern = true, bool OldPreconditioner = false){
        return solver->SetMatrix(sm, ModifiedPattern, OldPreconditioner);
    }
    bool Solve(std::vector<double>& b, std::vector<double>& x){
        return solver->Solve(b, x);
    }
    bool Solve(long bsize, double* b, long xsize, double* x){
        return solver->Solve(bsize, b, xsize, x);
    }
    bool SetParameter(std::string name, std::string value){
        return solver->SetParameter(name, value);
    }
    bool SetParameterReal(std::string name, double value){
        return solver->SetParameterReal(name, value);
    }
    bool GetParameter(std::string name, std::string& value){
        return solver->GetParameter(name, value); }
    bool GetParameterReal(std::string name, double& value){
        return solver->SetParameterReal(name, value);
    }
    double PreconditionerTime() {
        return solver->PreconditionerTime();
    }
    double IterationsTime() {
        return solver->IterationsTime();
    }
    double ResidualNorm() {
        return solver->ResidualNorm();
    }
    long NumIterations() {
        return solver->NumIterations();
    }
    std::string GetReason() {
        return solver->GetReason();
    }
    LinearSolver& setVerbosityLevel(int lvl){
        return solver->setVerbosityLevel(lvl), *this;
    }
};


#endif //AORTIC_VALVE_LINEARSOLVER_H
