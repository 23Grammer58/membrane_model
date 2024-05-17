//
// Created by alex on 04.12.2020.
//

#ifndef AORTIC_VALVE_LINEARSOLVERINMOST_H
#define AORTIC_VALVE_LINEARSOLVERINMOST_H

#ifdef USE_INMOST_SOLVERS
#include "inmost.h"
#include "LinearSolverInterface.h"

class LinearSolverInmost: public LinearSolverBase {
    INMOST::Sparse::Matrix A;
    INMOST::Sparse::Vector _x, _b;
    INMOST::Solver solver;
    bool inited = false;
public:
    LinearSolverInmost(std::string solverName): solver(solverName) {  }
    bool SetMatrix(SparseMatrix& sm, bool ModifiedPattern = true, bool OldPreconditioner = false) override;
    bool Solve(long bsize, double* b, long xsize, double* x) override;
    static std::vector<double> SAXPY(double a1, SparseMatrix& am, std::vector<double>& x, double a2, std::vector<double>& b); //a1 * A * x + a2 * b

    bool SetParameter(std::string name, std::string value) override;
    bool GetParameter(std::string name, std::string& value) override { value = solver.GetParameter(name); return value != ""; }
    bool GetParameterReal(std::string name, double& value) override {
        std::string v = solver.GetParameter(name);
        if (v == "") return false;
        value = atof(v.c_str());
        return true;
    }
    double PreconditionerTime() override { return solver.PreconditionerTime(); }
    double IterationsTime() override { return solver.IterationsTime(); }
    double ResidualNorm() override { return solver.Residual(); }
    long NumIterations() override { return solver.Iterations(); }
    std::string GetReason() override { return solver.GetReason(); }
    bool SetParameterReal(std::string name, double value) override;
    [[nodiscard]] std::unique_ptr<LinearSolverBase> clone() const override;

};
[[nodiscard]] static std::vector<std::string> getAvailableInmostSolvers() { return INMOST::Solver::getAvailableSolvers(); }
#else
class LinearSolverInmost: public LinearSolverBase {
public:
    LinearSolverInmost(std::string solverName){}
    bool SetMatrix(SparseMatrix& sm, bool ModifiedPattern = true, bool OldPreconditioner = false) override { return false; }
    bool Solve(long bsize, double* b, long xsize, double* x) override { return false; }
    [[nodiscard]] std::unique_ptr<LinearSolverBase> clone() const override { return std::make_unique<LinearSolverInmost>(""); }
};
[[nodiscard]] static std::vector<std::string> getAvailableInmostSolvers() { return std::vector<std::string>(); }
#endif

#endif //AORTIC_VALVE_LINEARSOLVERINMOST_H
