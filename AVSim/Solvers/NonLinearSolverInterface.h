//
// Created by alex on 08.06.2021.
//

#ifndef AORTIC_VALVE_NONLINEARSOLVERINTERFACE_H
#define AORTIC_VALVE_NONLINEARSOLVERINTERFACE_H
#include "LinearSolver.h"
struct NLProblem{
    using Vector = double*;
    using ConstVector = const double*;
    using Matrix = SparseMatrix*;
    std::function<int(ConstVector x, Vector b)> m_rhs;
    std::function<int(ConstVector x, Matrix A)> m_jac;
    int m_dofs = 0;

    void setResidualFunc(const std::function<int(ConstVector x, Vector b)>& residual) { m_rhs = residual; }
    void setJacobianFunc(const std::function<int(ConstVector x, Matrix A)>& jac) { m_jac = jac; }
    void setProbSize(int n) { m_dofs = n; }
};

class NonLinearSolverBase {
public:
    using Vector = double*;
    using ConstVector = const double*;
    using Matrix = SparseMatrix*;
    NLProblem m_prob;
    std::function<int(ConstVector x)> m_monitor;
    virtual int setProblem(const NLProblem& prob) { return m_prob = prob, 0; }
    virtual int setInterIterationMonitorFunc(const std::function<int(ConstVector x)>& monitor) { return m_monitor = monitor, 0; }

    virtual bool SetParameter(std::string param, std::string val) { return false; }
    virtual bool SetParameterReal(std::string param, double val) { return SetParameter(param, std::to_string(val)); }
    virtual bool SetParameterInt(std::string param, int val) { return SetParameter(param, std::to_string(val)); }
    virtual bool GetParameter(std::string name, std::string& value) const { return false; }
    virtual bool GetParameterReal(std::string name, double& value) const {
        std::string sval = "";
        if (!GetParameter(name, sval)) return false;
        value = stod(sval);
        return true;
    }
    virtual bool GetParameterInt(std::string name, long& value) const {
        std::string sval = "";
        if (!GetParameter(name, sval)) return false;
        value = stol(sval);
        return true;
    }
    virtual std::string GetReason() const { return "Reason is not specified"; };

    virtual bool Solve() = 0;
    virtual long GetNumFuncEvals() const { return 1-LONG_MAX; }
    virtual long GetNumJacEvals() const { return 1-LONG_MAX; }
    virtual long GetNumNolinSolvIters() const { return 1-LONG_MAX; }
    virtual double GetResidualNorm() const { return NAN; }
    virtual long GetNumLinIters() const { return 1-LONG_MAX; }
};


#endif //AORTIC_VALVE_NONLINEARSOLVERINTERFACE_H
