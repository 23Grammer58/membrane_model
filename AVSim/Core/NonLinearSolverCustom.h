//
// Created by alex on 17.06.2021.
//

#ifndef AORTIC_VALVE_NONLINEARSOLVERCUSTOM_H
#define AORTIC_VALVE_NONLINEARSOLVERCUSTOM_H

#include "NSWorldWrapper.h"
#include "AVSim/Solvers/NonLinearSolverInterface.h"

class NonLinearSolverCustom: public NonLinearSolverBase {
public:
    using NSWorldWrapper = World3d::NSWorldWrapper;

    LinearSolver* m_solver = nullptr;
    SparseMatrix m_sm;
    std::vector<double> m_rhs, m_x, m_dx;
    NSWorldWrapper* m_ns;

    using TauAlgo = std::function<double(int it, double rel_resid, double abs_resid, NonLinearSolverCustom& ns)>;
    using StopCond = std::function<bool(int it, double rel, double abs, int& reason, NonLinearSolverCustom& ns)>;
    TauAlgo tau_algo = nullptr;
    StopCond stop_cond = [](int it, double rel, double abs, int& reason, NonLinearSolverCustom& ns){
        if (it >= 2000 || rel < 1e-7) {
            return reason = 0, true;
        } else if (abs > 1e15 || std::isnan(abs)) {
            return reason = -1, true;
        }
        return reason = 1, false;
    };
    std::function<int(NonLinearSolverCustom& )> newton_algo = [](NonLinearSolverCustom& ns){
        int stat = -1;
        if ((stat = ns.assemble()) < 0) return stat;
        ns.saveResidual();
        if (!ns.solveLinear(ns.m_dx)) return (stat = -4);
        double tau = 1;
        if (ns.tau_algo) tau = ns.tau_algo(ns.m_it, ns.resid_now/ns.resid_init, ns.resid_now, ns);
        if ((stat = ns.applySolution(tau)) < 0) return stat;
        return stat;
    };

    NonLinearSolverCustom& SetInitialGuess(const std::vector<double>& v) { return m_x = v, *this; }
    NonLinearSolverCustom& SetNSWorldWrapper(NSWorldWrapper& ns) { return m_ns = &ns, *this; };
    NonLinearSolverCustom& SetLinearSolver(LinearSolver& s) {return m_solver = &s, *this;};
    int assemble();
    bool solveLinear(std::vector<double>& to_dx);
    void saveResidual();
    int applySolution(double tau);
    bool Solve();

    double computeRhsNorm(std::vector<double>& residual, unsigned int N) { return computeRhsNorm(residual, N, N); }
    double computeRhsNorm(std::vector<double>& residual, unsigned int N1, unsigned int N2);
    std::string GetReason() const override;
    int getReasonFlag() const { return m_reason; }
    long GetNumFuncEvals() const override { return m_fe; }
    long GetNumJacEvals() const override { return m_jace; }
    long GetNumNolinSolvIters() const { return m_lit; }
    double GetResidualNorm() const { return resid_now; }
    long GetNumLinIters() const { return m_lit; }
private:
    int m_fe = 0, m_jace = 0, m_lit = 0;
    double resid_init = -1;
    double resid_now = -1;
    unsigned norm = 2;
    int m_it = 1;
    int m_reason = 0;
};


#endif //AORTIC_VALVE_NONLINEARSOLVERCUSTOM_H
