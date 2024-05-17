//
// Created by alex on 10.11.2020.
//

#ifndef AORTIC_VALVE_LINEARSOLVEREIGEN_H
#define AORTIC_VALVE_LINEARSOLVEREIGEN_H
#include "LinearSolverInterface.h"
#include <Eigen/IterativeLinearSolvers>
#include <iostream>


class LinearSolverEigen: public LinearSolverBase {
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double, Eigen::RowMajor>, Eigen::IncompleteLUT<double>> solver;
    Eigen::SparseMatrix<double, Eigen::RowMajor> mat;
    bool inited = false;
    double t_precond = 0, t_solve = 0;
    double resid_norm = -1;
public:
    bool SetMatrix(SparseMatrix& sm, bool ModifiedPattern = true, bool OldPreconditioner = false) override;
    bool Solve(long bsize, double* b, long xsize, double* x) override;
    bool SetParameter(std::string name, std::string value) override;
    bool GetParameter(std::string name, std::string& value) override;
    double PreconditionerTime() override { return t_precond; }
    double IterationsTime() override { return t_solve; }
    double ResidualNorm() override { return resid_norm; }
    long NumIterations() override { return solver.iterations(); }
    std::string GetReason() override;
    std::unique_ptr<LinearSolverBase> clone() const override{
        return std::make_unique<LinearSolverEigen>();
    }
};


#endif //AORTIC_VALVE_LINEARSOLVEREIGEN_H
