//
// Created by alex on 17.06.2021.
//

#include "NSWorldWrapper.h"

int World3d::NSWorldWrapper::getNodfs() {
    return pw->compStateSize();
}

int World3d::NSWorldWrapper::setX(const double *X) {
    return pw->setStateVector(X);
}

int World3d::NSWorldWrapper::getCurrentX(double *X) {
    return pw->compStateVector(X);
}

std::vector<double> World3d::NSWorldWrapper::getCurrentX() {
    std::vector<double> X(getNodfs());
    getCurrentX(X.data());
    return X;
}

int World3d::NSWorldWrapper::Residual(const double *X, double *R) {
    int stat = 0;
    if ((stat = setX(X)) < 0) return stat;
    if ((stat = compResidual(R)) < 0) return stat;

    return 0;
}

int World3d::NSWorldWrapper::Jacobian(const double *X, SparseMatrix *sm) {
    int stat = 0;
    if ((stat = setX(X)) < 0) return stat;
    if ((stat = compJacobian(sm)) < 0) return stat;

    return 0;
}

int World3d::NSWorldWrapper::System(const double *X, SparseMatrix *sm, double *R) {
    int stat = 0;
    if ((stat = setX(X)) < 0) return stat;
    if ((stat = compResidual(R)) < 0) return stat;
    if ((stat = compJacobian(sm)) < 0) return stat;
    return 0;
}

int World3d::NSWorldWrapper::RenderIteration() {
    if (it_time.elapsed() > render_freq) {
        pw->Render();
        it_time.reset();
    }
    return 0;
}

int World3d::NSWorldWrapper::RenderVector(const double *X) {
    if (it_time.elapsed() > render_freq) {
        for (auto &i: pw->objs())
            if (i.second.second)
                for (auto v: i.second.first.m_mesh.vertices()) {
                    i.second.first.m_next_x[v] = i.second.first.m_x[v];
                    i.second.first.m_x[v] = Point(*(X + 0), *(X + 1), *(X + 2));
                    X += 3;
                }
        pw->Render();
        it_time.reset();
        for (auto &i: pw->objs())
            if (i.second.second)
                for (auto v: i.second.first.m_mesh.vertices()) {
                    i.second.first.m_x[v] = i.second.first.m_next_x[v];
                }
    }
    return 0;
}

int World3d::NSWorldWrapper::compResidual(double *R) {
    return pw->compResidual(R);
}

int World3d::NSWorldWrapper::compJacobian(SparseMatrix *sm) {
    return pw->compJacobian(sm);
}

NLProblem World3d::NSWorldWrapper::makeProblem(World3d::NSWorldWrapper &nsww) {
    NLProblem nlp;
    nlp.m_dofs = nsww.getNodfs();
    nlp.m_rhs = [&nsww](const double *x, double *b) -> int{
        return nsww.Residual(x, b);
    };
    nlp.m_jac = [&nsww](const double *x, SparseMatrix *sm) -> int{
        return nsww.Jacobian(x, sm);
    };
    return nlp;
}

double World3d::NSWorldWrapper::VectorNorm(double *residual, long size, unsigned int N1, unsigned int N2) {
    auto& rhs = residual;
    auto get_loc_norm = [N = N2](auto begin, auto end){
        if (N == 0)
            return *std::max_element(begin, end);
        else if (N == 1){
            double res = 0;
            for (auto i = begin; i != end; ++i) res += fabs(*i);
            return res;
        } else {
            double res = 0;
            for (auto i = begin; i != end; ++i) res += pow(fabs(*i), N);
            return pow(res, 1.0/N);
        }
    };
    if (N1 == N2){
        return get_loc_norm(rhs, rhs + size);
    } else {
        auto get_glob_norm = [N = N1, get_loc_norm](double* v, long size){
            long K = size/3;
            double res = 0;
            if (N == 0){
                for (int i = 0; i < K; ++i) res = std::max(res, get_loc_norm(v+3*i, v + 3*(i+1)));
            } else if (N == 1){
                for (int i = 0; i < K; ++i) res += get_loc_norm(v+3*i, v + 3*(i+1));
            } else {
                for (int i = 0; i < K; ++i) res += pow(get_loc_norm(v+3*i, v + 3*(i+1)), N);
                res = pow(res, 1.0/N);
            }
            return res;
        };
        return get_glob_norm(rhs, size);
    }
}

double World3d::NSWorldWrapper::compCorrentForceNorm(unsigned int Nobj, unsigned int Nmesh, unsigned int Nnode) {
    //int stat = pw->UpdateForces();
    return pw->getWorldNorm("v:force", Nobj, Nmesh, Nnode);
}
