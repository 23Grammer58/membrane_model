//
// Created by alex on 10.11.2020.
//

#include "LinearSolverEigen.h"
#include <chrono>

bool LinearSolverEigen::SetMatrix(SparseMatrix &sm, bool ModifiedPattern, bool OldPreconditioner) {
    mat = Eigen::SparseMatrix<double, Eigen::RowMajor>(sm.rows(),sm.cols());
    struct Sizer{
        using value_type = std::map<int, double>::size_type;
        SparseMatrix& sm;
        Sizer(SparseMatrix& sm): sm{sm} {}
        value_type operator[](int i) const { return sm[i].size(); }
    };
    Sizer sizer(sm);
    mat.reserve(sizer);
    for (int i = 0; i < sm.rows(); ++i)
        for (auto& it: sm[i])
            mat.insert(i, it.first) = it.second;
//        mat.makeCompressed();
    std::chrono::time_point<std::chrono::high_resolution_clock> m_beg = std::chrono::high_resolution_clock::now();
    if (ModifiedPattern || !inited)
        solver.analyzePattern(mat);
    solver.factorize(mat);
    t_precond = std::chrono::duration_cast<std::chrono::duration<double, std::ratio<1> >>(std::chrono::high_resolution_clock::now() - m_beg).count();
    inited = true;

    return true;
}

bool LinearSolverEigen::Solve(long bsize, double* b, long xsize, double* x){
    if (!inited) return false;
    Eigen::Map<Eigen::VectorXd, Eigen::Unaligned> bb(b, bsize), xx(x, xsize);
    std::chrono::time_point<std::chrono::high_resolution_clock> m_beg = std::chrono::high_resolution_clock::now();
    xx = solver.solveWithGuess(bb, xx);
    t_solve= std::chrono::duration_cast<std::chrono::duration<double, std::ratio<1> >>(std::chrono::high_resolution_clock::now() - m_beg).count();
    std::cout << "#iterations:     " << solver.iterations() << " estimated error: " << solver.error() << std::endl;
    resid_norm = solver.error() * bb.norm();
    return solver.error() < solver.tolerance();
}

bool LinearSolverEigen::SetParameter(std::string name, std::string value) {
    if (name == "maximum_iterations")
        return solver.setMaxIterations(atoi(value.c_str())), true;
    else if (name == "relative_tolerance")
        return solver.setTolerance(atof(value.c_str())), true;
    else if (name == "drop_tolerance")
        return solver.preconditioner().setDroptol(atof(value.c_str())), true;
    else if (name == "fill_factor")
        return solver.preconditioner().setFillfactor(atof(value.c_str())), true;
    return false;
}

bool LinearSolverEigen::GetParameter(std::string name, std::string& value) {
    if (name == "maximum_iterations")
        return value = std::to_string(solver.maxIterations()), true;
    else if (name == "relative_tolerance")
        return value = std::to_string(solver.tolerance()), true;

    return false;
}

std::string LinearSolverEigen::GetReason() {
    Eigen::ComputationInfo info = solver.info();
    switch (info) {
        case Eigen::ComputationInfo::Success: return "Success";
        case Eigen::ComputationInfo::NoConvergence: return "NoConvergence";
        case Eigen::ComputationInfo::NumericalIssue: return "NumericalIssue";
        case Eigen::ComputationInfo::InvalidInput: return "InvalidInput";
        default: return "Unreachable code";
    }
}
