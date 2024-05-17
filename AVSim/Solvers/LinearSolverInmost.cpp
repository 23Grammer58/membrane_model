//
// Created by alex on 04.12.2020.
//

#include "LinearSolverInmost.h"

#ifdef USE_INMOST_SOLVERS
using namespace INMOST;

bool LinearSolverInmost::SetMatrix(SparseMatrix &sm, bool ModifiedPattern, bool OldPreconditioner) {
    A.Clear();
    A.SetInterval(0, sm.rows());
    for (int i = 0; i < sm.rows(); ++i) {
        auto& row = sm[i];
        auto& rowA = A[i];
        for (auto &j: row)
            rowA[j.first] = j.second;
    }
    if (verb_lvl >= 2)
    std::cout << "\tStart preconding" << std::endl;
    solver.SetMatrix(A, ModifiedPattern, OldPreconditioner);
    if (verb_lvl >= 2)
    std::cout << "\tPreconded" << std::endl;
    inited = true;

    return true;
}

std::vector<double> LinearSolverInmost::SAXPY(double a1, SparseMatrix& sm, std::vector<double>& x, double a2, std::vector<double>& b){
    INMOST::Sparse::Matrix A;
    A.Clear();
    A.SetInterval(0, sm.rows());
    for (int i = 0; i < sm.rows(); ++i) {
        auto& row = sm[i];
        auto& rowA = A[i];
        for (auto &j: row)
            rowA[j.first] = j.second;
    }
    INMOST::Sparse::Vector _x, _b, _r;
    if (_b.Size() != b.size())
    { _b.Clear(); _b.SetInterval(0, b.size()); }
    if (_x.Size() != x.size())
    {_x.Clear(); _x.SetInterval(0, x.size());}
    for (int i = 0; i < b.size(); ++i) _b[i] = b[i];
    for (int i = 0; i < x.size(); ++i) _x[i] = x[i];
    A.MatVec(a1, _x, a2, _b);
    std::vector<double> res(_b.Size());
    for (int i = 0; i < _b.Size(); ++i) res[i] = _b[i];
    return res;
}

bool LinearSolverInmost::Solve(long bsize, double* b, long xsize, double* x) {
    if (!inited) return false;
    if (_b.Size() != bsize)
    { _b.Clear(); _b.SetInterval(0, bsize); }
    if (_x.Size() != xsize)
    {_x.Clear(); _x.SetInterval(0, xsize);}
    for (int i = 0; i < bsize; ++i) _b[i] = b[i];
    for (int i = 0; i < xsize; ++i) _x[i] = x[i];
    if (verb_lvl >= 2) std::cout << "\tStart solving" << std::endl;
    bool res = solver.Solve(_b, _x);
    if (verb_lvl >= 2) std::cout << "\tEnd solving" << std::endl;
    if (verb_lvl >= 1)
        std::cout << "#iterations: " << solver.Iterations() << " residual: " << solver.Residual() <<
                  " time precondition / time iters: " << solver.PreconditionerTime() << " / " << solver.IterationsTime()  << std::endl;
    for (int i = 0; i < xsize; ++i) x[i] = _x[i];

    return res;
}

bool LinearSolverInmost::SetParameter(std::string name, std::string value) {
    solver.SetParameter(name, value);
    std::string res = solver.GetParameter(name);
    return res != "";
}

bool LinearSolverInmost::SetParameterReal(std::string name, double value) {
    solver.SetParameterReal(name, value);
    std::string res = solver.GetParameter(name);
    return res != "";
}

std::unique_ptr<LinearSolverBase> LinearSolverInmost::clone() const {
    auto res =  std::make_unique<LinearSolverInmost>(solver.SolverName());
    res->verb_lvl = verb_lvl;
//    res->solver = solver;
    return res;
}

#else


#endif