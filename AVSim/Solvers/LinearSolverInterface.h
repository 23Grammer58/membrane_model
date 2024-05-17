//
// Created by alex on 03.11.2020.
//

#ifndef AORTIC_VALVE_LINEARSOLVERINTERFACE_H
#define AORTIC_VALVE_LINEARSOLVERINTERFACE_H
#include <vector>
#include <map>
#include <memory>
#include <cassert>
#include <string>
#include <stdexcept>
//#include "LinearSolverEigen.h"

class SparseMatrix{
    std::vector<std::map<int, double>> m_rows;
public:
    SparseMatrix(int N = 0): m_rows(N) {}
    void clear(){ for (auto& i: m_rows) i.clear(); }
    void resize(int N) { m_rows.resize(N); }
    auto& operator[](int i) { assert(i < rows() && "Wrong dimension"); return m_rows[i]; }
    const auto& operator[](int i) const { assert(i < rows() && "Wrong dimension"); return m_rows[i]; }
//    auto begin() { return m_rows.begin(); }
//    auto end() { return m_rows.end(); }
    double& operator()(int i, int j){ 
        assert(i < rows() && j < cols() && "Wrong dimension");
        auto& row = m_rows[i];
        auto it = row.find(j);
        if (it == row.end()) return row.insert({j, 0.0}).first->second; 
        else return it->second; 
    }
    double operator()(int i, int j) const{
        assert(i < rows() && j < cols() && "Wrong dimension");
        auto& row = m_rows[i];
        auto it = row.find(j);
        if (it == row.end()) return 0.0;
        else return it->second;
    }
    int rows() const { return m_rows.size(); }
    int cols() const { return m_rows.size(); }
    struct IndexFrom{
        int I0 = 0, J0 = 0;
        SparseMatrix* sm;
        IndexFrom(SparseMatrix& sm): sm{&sm} {}
        IndexFrom(SparseMatrix& sm, int I0, int J0): I0{I0}, J0{J0}, sm{&sm} {}
        struct iteratorIndexFrom{
            int J0;
            std::map<int, double>& row;
            iteratorIndexFrom(IndexFrom* ind, int i): row{ind->sm->m_rows[i+ind->I0]}, J0{ind->J0} {}
            auto& operator[](int j) { return row[j + J0]; }
        };
        iteratorIndexFrom operator[](int i) { return iteratorIndexFrom(this, i); }
        auto& operator()(int i, int j){ return sm->operator()(I0 + i, J0 + j); }
    };
    ///res = \sum_i c_i A_i, i = 0, ..., nmtx-1
    ///\note res may be one of matrices B_i
    static void LinSum(SparseMatrix& res, unsigned nmtx, const double* c/*[nmtx]*/, const SparseMatrix** A/*[nmtx]*/);
    static void ScalSum(SparseMatrix& res, double a, const SparseMatrix& A, double b, const SparseMatrix& B){
        const double c[2] = {a, b};
        const SparseMatrix* _A[2] = {&A, &B};
        LinSum(res, 2, c, _A);
    }
    static void Scal(SparseMatrix& res, double a, const SparseMatrix& A){
        const SparseMatrix* pA = &A;
        LinSum(res, 1, &a, &pA);
    }
    static SparseMatrix ScalSum(double a, const SparseMatrix& A, double b, const SparseMatrix& B){
        SparseMatrix res;
        ScalSum(res, a, A, b, B);
        return res;
    }
    static SparseMatrix Identity(int N, double c = 1.0){
        SparseMatrix r(N);
        for (int i = 0; i < N; ++i)
            r[i][i] = c;
        return r;    
    }
};

class LinearSolverBase {
protected:
    int verb_lvl = 0;
public:
    virtual bool SetMatrix(SparseMatrix& sm, bool ModifiedPattern = true, bool OldPreconditioner = false) = 0;
    virtual bool Solve(std::vector<double>& b, std::vector<double>& x){
        return Solve(static_cast<long>(b.size()), b.data(), static_cast<long>(x.size()), x.data());
    };
    virtual bool Solve(long bsize, double* b, long xsize, double* x) = 0;
    /// Parameters:
    /// - "maximum_iterations" - total number of iterations
    /// - "absolute_tolerance" - iterative method will stop on i-th iteration
    ///                          if ||A x(i)-b|| < absolute_tolerance
    /// - "relative_tolerance" - iterative method will stop on i-th iteration
    ///                          if ||A x(i)-b||/||A x(0) - b||
    /// - "divergence_tolerance" - iterative method will fail if
    ///                          ||A x(i) - b|| > divergence_tolerance
    /// - "drop_tolerance"     - tolerance for dropping values during incomplete factorization,
    /// - "fill_factor"        - After the elimination of the row, only the fill largest elements
    ///                          in the L part and the fill largest elements in the U part are kept
    ///                          (in addition to the diagonal element ). Note that fill is computed
    ///                          from the input parameter fill_factor which is used the ratio to control
    ///                          the fill_in relatively to the initial number of nonzero elements.
    virtual bool SetParameter(std::string name, std::string value){ return false; }
    virtual bool GetParameter(std::string name, std::string& value){ return false; }
    virtual bool GetParameterReal(std::string name, double& value){
        std::string sval = "";
        if (!GetParameter(name, sval)) return false;
        value = atof(sval.c_str());
        return true;
    }
    virtual double PreconditionerTime() { return -1; }
    virtual double IterationsTime() { return -1; }
    virtual double ResidualNorm() { return -1; }
    virtual long NumIterations() { return -1; }
    virtual std::string GetReason() { return "Reason not specified"; }
    virtual bool SetParameterReal(std::string name, double value){ return SetParameter(name, std::to_string(value)); }
    virtual std::unique_ptr<LinearSolverBase> clone() const = 0;
    virtual LinearSolverBase& setVerbosityLevel(int lvl) { return verb_lvl = lvl, *this; }
};



#endif //AORTIC_VALVE_LINEARSOLVERINTERFACE_H
