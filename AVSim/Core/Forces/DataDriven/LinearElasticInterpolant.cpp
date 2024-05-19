//
// Created by alex on 01.02.2021.
//

#include "LinearElasticInterpolant.h"
using namespace std;
using namespace World3d;

bool LinearElasticResponse::is_positive(const std::array<double, 6>& C){
    return ((C[0] > 0) && (C[0] * C[2] > C[1] * C[1]) && (C[0]*C[2]*C[4] + 2 * C[1]*C[3]*C[5] > C[0]*C[5]*C[5] + C[2]*C[3]*C[3] + C[4]*C[1]*C[1]));
}

bool LinearElasticResponse::fit(const std::vector<Response>& data){
    least_squares_fit(data);
    if (!is_positive(m_C))
        conditional_optimization(data);
    return true;    
}

static void set_e_to_loc(const array<double, 3>& E, double E_loc[3][6]){
    E_loc[0][0] = E[0], E_loc[0][1] = 2 * E[1], E_loc[0][2] = E[2];
    E_loc[1][1] = E[0], E_loc[1][3] = E[2], E_loc[1][5] = 2 * E[1];
    E_loc[2][2] = E[0], E_loc[2][3] = 2 * E[1], E_loc[2][4] = E[2];
}

void LinearElasticResponse::least_squares_fit(const std::vector<Response>& data){
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> A;
    Eigen::VectorXd b;
    A.resize(3 * data.size(), 6);
    b.resize(3 * data.size());
    auto x = Eigen::Map<Eigen::Matrix<double, 6, 1>>(m_C.data());
    array<double, 3> e{0}, s{0};
    double e_loc[3][6] = {0};
    for (int i = 0; i < data.size(); ++i){
        Response::toStressStrainTensor(data[i].xi, data[i].response, e, s);
        std::swap(e[1], e[2]); std::swap(s[1], s[2]);
        set_e_to_loc(e, e_loc);
        for (int j = 0; j < 3; ++j) {
            for (int l = 0; l < 6; ++l)
                A(3*i+j, l) = e_loc[j][l];
            b[3*i+j] = s[j];
        }
    }
    x = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
}

std::array<double, 3> LinearElasticResponse::evaluate(const double* x){
    double e_k1 = exp(2 * x[0]), e_k2 = exp(2 *x[1]), k3 = x[2];
    std::array<double, 3> E, C, res;
    E[0] = 0.5 * expm1(2*x[0]);
    E[1] = 0.5 * k3 * e_k1;
    E[2] = 0.5 * (expm1(2 * x[1]) + k3 * k3 * e_k1);
    C[0] = m_C[0]*E[0] + 2 * m_C[1] * E[1] + m_C[3] * E[2];
    C[1] = 2 * m_C[1]*E[0] + 4 * m_C[2] * E[1] + 2 * m_C[5] * E[2];
    C[2] = m_C[3]*E[0] + 2 * m_C[5] * E[1] + m_C[4] * E[2];
    res[0] = e_k1 * (C[0] + k3 * (C[1] + k3 * C[2]));
    res[1] = e_k2 * C[2];
    res[2] = e_k1 * (0.5 * C[1] + k3 * C[2]);
    return res;
}

std::array<double, 6> LinearElasticResponse::response_gradient(const std::array<double, 3>& x){
    double k0 = exp(2 * x[0]), k1 = exp(2 *x[1]), k2 = x[2];
    double dk00 = 2*k0, dk11 = 2*k1, dk22 = 1;
    double  E0 = 0.5 * expm1(2*x[0]), 
            E1 = 0.5 * k2 * k0, 
            E2 = 0.5 * (expm1(2 * x[1]) + k2 * k2 * k0);
    double dE00 = k0;
    double dE10 = k2*k0, dE12 = E0;
    double dE20 = k2 * k2 * k0, dE21 = k1, dE22 = k2 * k0;
    
    double C0 = m_C[0]*E0 + 2 * m_C[1] * E1 + m_C[3] * E2;
    double C1 = 2 * m_C[1]*E0 + 4 * m_C[2] * E1 + 2 * m_C[5] * E2;
    double C2 = m_C[3]*E0 + 2 * m_C[5] * E1 + m_C[4] * E2;

    double C00 = m_C[0]*dE00 + 2 * m_C[1] * dE10 + m_C[3] * dE20;
    // double C01 = m_C[3] * dE21;
    double C02 = 2 * m_C[1] * dE12 + m_C[3] * dE22;
    double C10 = 2 * m_C[1]*dE00 + 4 * m_C[2] * dE10 + 2 * m_C[5] * dE20;
    // double C11 = 2 * m_C[5] * dE21;
    double C12 = 4 * m_C[2] * dE12 + 2 * m_C[5] * dE22;
    double C20 = m_C[3]*dE00 + 2 * m_C[5] * dE10 + m_C[4] * dE20;
    double C21 = m_C[4] * dE21;
    double C22 = 2 * m_C[5] * dE12 + m_C[4] * dE22;

    double r00 = k0 * ( 2 * (C0 + k2 * (C1 + k2 * C2)) + (C00 + k2 * (C10 + k2 * C20)) );
    double r02 = k0 * ( C02 + C1 + k2 * (C10 + 2*C2 + k2 * C22));
    double r10 = k1 * C20;
    double r11 = k1 * ( 2*C2 + C21 );
    double r12 = k1 * C22;
    double r22 = k0 * ( 0.5 * C12 + C2 + k2*C22);

    return {r00, r11, r22, r10, r02, r12};
}

std::array<std::array<double, 3>, 3> LinearElasticResponse::gradient(const std::array<double, 3>& xi){
    auto r = response_gradient(xi);
    return {std::array<double, 3>{r[0], r[3], r[4]}, {r[3], r[1], r[5]}, {r[4], r[5], r[2]}};
}

void LinearElasticResponse::conditional_optimization(const std::vector<Response>& data){
    Eigen::Matrix<double, Eigen::Dynamic, 6> A;
    Eigen::VectorXd b;
    A.resize(3 * data.size(), Eigen::NoChange);
    b.resize(3 * data.size());
    array<double, 3> e{0}, s{0};
    double e_loc[3][6] = {0};
    for (int i = 0; i < data.size(); ++i){
        Response::toStressStrainTensor(data[i].xi, data[i].response, e, s);
        std::swap(e[1], e[2]); std::swap(s[1], s[2]);
        set_e_to_loc(e, e_loc);
        for (int j = 0; j < 3; ++j) {
            for (int l = 0; l < 6; ++l)
                A(3*i+j, l) = e_loc[j][l];
            b[3*i+j] = s[j];
        }
    }

    Eigen::Matrix<double, 6, 6> D = A.transpose() * A;
    Eigen::Matrix<double, 6, 1> c = -(A.transpose() * b);
    double c0 = 0;//b.dot(b) / 2;


    Eigen::Matrix<int, 4*3, 6> L;
    for (int i = 0; i < 4; ++i){
        L.row(4*0 + i) << -1, ((i%2) ? 1 : -1), ((i/2) ? 1 : -1), 0, 0, 0;
        L.row(4*1 + i) << 0, ((i%2) ? 1 : -1), -1, 0, 0, ((i/2) ? 1 : -1);
        L.row(4*2 + i) << 0, 0, 0, ((i%2) ? 1 : -1), -1, ((i/2) ? 1 : -1);
    }
    int* Lp[6];
    double *Dp[6];
    for (int i = 0; i < 6; ++i) Lp[i] = L.data() + 4*3*i, Dp[i] = D.data() + 6*i;
    int b_zero[12] = {0};
    bool fl[6] = {false};
    bool *fu = fl;

    typedef CGAL::Quadratic_program_from_iterators
            <int**,                                                // for L
             int*,                                                 // for b
             CGAL::Const_oneset_iterator<CGAL::Comparison_result>, // for r
             bool*,                                                // for fl
             int*,                                                 // for l
             bool*,                                                // for fu
             int*,                                                 // for u
             double**,                                             // for D
             double*>                                              // for c
    Program;
    CGAL::Const_oneset_iterator<CGAL::Comparison_result> r(CGAL::SMALLER);
    Program qp (6, L.rows(), Lp, b_zero, r, fl, b_zero, fu, b_zero, Dp, c.data(), c0);
    auto sol = CGAL::solve_quadratic_program(qp, CGAL::MP_Float());
    if (!sol.is_optimal()){
        throw std::runtime_error(
                std::string("Conditional least-squares didn't reach optimum.\n") +
                "\tunbounded = " + (sol.is_unbounded() ? "true" : "false") +
                "\tinfeasible = " + (sol.is_infeasible() ? "true" : "false")
        );
    }
    double* sit = m_C.data();
    for (auto it = sol.variable_values_begin(); it != sol.variable_values_end(); ++it){
        *(sit++) = CGAL::to_double(*it);
    }
}