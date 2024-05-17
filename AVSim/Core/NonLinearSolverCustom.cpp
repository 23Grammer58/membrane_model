//
// Created by alex on 17.06.2021.
//

#include "NonLinearSolverCustom.h"

int NonLinearSolverCustom::assemble() {
    int stat = -1;
    m_fe++; m_jace++;
    if (m_ns){
        if ((stat = m_ns->System(m_x.data(), &m_sm, m_rhs.data())) < 0) return stat;
    } else if (m_prob.m_rhs && m_prob.m_jac){
        if ((stat = m_prob.m_rhs(m_x.data(), m_rhs.data())) < 0) return stat;
        if ((stat = m_prob.m_jac(m_x.data(), &m_sm)) < 0) return stat;
    }
    return stat;
}

bool NonLinearSolverCustom::solveLinear(std::vector<double> &to_dx) {
    if (!m_solver->SetMatrix(m_sm)) return false;
    bool slv = m_solver->Solve(m_rhs, to_dx);
    m_lit += m_solver->NumIterations();
    return slv;
}

void NonLinearSolverCustom::saveResidual() {
    resid_now = computeRhsNorm(m_rhs, norm);
    if (m_it == 1) resid_init = resid_now;
}

bool NonLinearSolverCustom::Solve() {
    m_it = 1;
    resid_init = resid_now = -1;
    m_fe = m_jace = m_lit = 0;
    m_rhs.resize(m_prob.m_dofs);
    m_dx.resize(m_prob.m_dofs);
    m_sm.resize(m_prob.m_dofs);
    if (m_x.empty()) m_x.resize(m_prob.m_dofs, 0.0);
    if (!m_solver) throw std::runtime_error("Lnear solver is not specified");
    int stat = 0;
    do {
        stat = newton_algo(*this);
        m_monitor(m_x.data());
        m_it++;
    } while (stat >= 0 && !stop_cond(m_it, resid_now  / resid_init, resid_now, m_reason, *this));
    if (stat > 0){
        m_fe++;
        int rstat = m_prob.m_rhs(m_x.data(), m_rhs.data());
        if (rstat >= 0) saveResidual();
    }
    if (stat < 0) m_reason = -3;
    if (m_reason >= 0) return true;
    return false;
}

int NonLinearSolverCustom::applySolution(double tau) {
    for (int i = 0; i < m_prob.m_dofs; ++i) m_x[i] -= tau*m_dx[i];
    return 0;
}

double NonLinearSolverCustom::computeRhsNorm(std::vector<double> &residual, unsigned int N1, unsigned int N2) {
    return NSWorldWrapper::VectorNorm(residual.data(), residual.size(), N1, N2);
}

std::string NonLinearSolverCustom::GetReason() const {
    switch (m_reason){
        case 2: return "Reached maximum of iterations";
        case 1: return "Iterations should be continued";
        case 0: return "Converged";
        case -1: return "Diverged";
        case -3: return "Error in assembling of system";
        default: return "Unknown status = " + std::to_string(m_reason);
    }
}
