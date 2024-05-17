//
// Created by alex on 01.02.2021.
//

#include "LinearZeroSeacher.h"
using namespace std;

bool LinearElastisity::is_positive(const std::array<double, 6>& C){
    return ((C[0] > 0) && (C[0] * C[2] > C[1] * C[1]) && (C[0]*C[2]*C[4] + 2 * C[1]*C[3]*C[5] > C[0]*C[5]*C[5] + C[2]*C[3]*C[3] + C[4]*C[1]*C[1]));
}

void LinearElastisity::init(const cloud_t& data){
    least_squares_fit(data);
    if (!is_positive(m_C))
        conditional_optimization(data);
}

void comp_E_sigma(const ResponseStatistics::Elem & q, array<double, 3>& E, array<double, 3>& S){
    double e_mk1 = exp(-2 * q.ksi[0]), e_mk2 = exp(-2 * q.ksi[1]), k3 = q.ksi[2];
    double e_k1 = 1.0/e_mk1, e_k2 = 1.0 / e_mk2;
    E[0] = 0.5 * (expm1(2 * q.ksi[0]));
    E[1] = 0.5 * k3 * e_k1;
    E[2] = 0.5 * (expm1(2 * q.ksi[1]) + k3 * k3 * e_k1);
    S[0] = e_mk1 * (q.response[0] - 2 * k3 * q.response[2]) + e_mk2 * k3 * k3 * q.response[1];
    S[1] = -e_mk2 * k3 * q.response[1] + e_mk1 * q.response[2];
    S[2] = e_mk2 * q.response[1];
}

void set_e_to_loc(const array<double, 3>& E, double E_loc[3][6]){
    E_loc[0][0] = E[0], E_loc[0][1] = 2 * E[1], E_loc[0][2] = E[2];
    E_loc[1][1] = E[0], E_loc[1][3] = E[2], E_loc[1][5] = 2 * E[1];
    E_loc[2][2] = E[0], E_loc[2][3] = 2 * E[1], E_loc[2][4] = E[2];
}

void LinearElastisity::least_squares_fit(const cloud_t& data){
#ifdef USE_EIGEN_LINSQ
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> A;
    Eigen::VectorXd b;
    A.resize(3 * data.size(), 6);
    b.resize(3 * data.size());
    auto x = Eigen::Map<Eigen::Matrix<double, 6, 1>>(m_C.data());
    array<double, 3> e{0}, s{0};
    double e_loc[3][6] = {0};
    for (int i = 0; i < data.size(); ++i){
        comp_E_sigma(data[i], e, s);
        set_e_to_loc(e, e_loc);
        for (int j = 0; j < 3; ++j) {
            for (int l = 0; l < 6; ++l)
                A(3*i+j, l) = e_loc[j][l];
            b[3*i+j] = s[j];
        }
    }
    x = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
#else
    int n = 3 * data.size(), p = m_C.size();
    double chisq = 0;
    gsl_matrix * E = gsl_matrix_alloc (n, p);
    gsl_vector * sigma = gsl_vector_alloc (n);

    gsl_vector* c = gsl_vector_alloc (p);
    gsl_matrix* cov = gsl_matrix_alloc (p, p);


    array<double, 3> e{0}, s{0};
    double e_loc[3][6] = {0};
    for (int i = 0; i < data.size(); ++i){
        comp_E_sigma(data[i], e, s);
        set_e_to_loc(e, e_loc);
        for (int j = 0; j < 3; ++j) {
            for (int l = 0; l < 6; ++l)
                gsl_matrix_set(E, 3 * i + j, l, e_loc[j][l]);
            gsl_vector_set(sigma, 3*i+j, s[j]);
        }
    }

    {
        gsl_multifit_linear_workspace *work
                = gsl_multifit_linear_alloc(n, p);
        gsl_multifit_linear(E, sigma, c, cov, &chisq, work);
        gsl_multifit_linear_free(work);
    }

    m_C[0] = gsl_vector_get(c, 0);
    m_C[1] = gsl_vector_get(c, 1);
    m_C[2] = gsl_vector_get(c, 5);
    m_C[3] = gsl_vector_get(c, 2);
    m_C[4] = gsl_vector_get(c, 4);
    m_C[5] = gsl_vector_get(c, 3);

    gsl_matrix_free (E);
    gsl_vector_free (sigma);
    gsl_vector_free (c);
    gsl_matrix_free (cov);
#endif
}

std::array<double, 3> LinearElastisity::operator()(const array<double, 3>& x){
    return operator()(x.data());
}

std::array<double, 3> LinearElastisity::operator()(const double* x){
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

void LinearElastisity::conditional_optimization(const cloud_t& data){
    Eigen::Matrix<double, Eigen::Dynamic, 6> A;
    Eigen::VectorXd b;
    A.resize(3 * data.size(), Eigen::NoChange);
    b.resize(3 * data.size());
    array<double, 3> e{0}, s{0};
    double e_loc[3][6] = {0};
    for (int i = 0; i < data.size(); ++i){
        comp_E_sigma(data[i], e, s);
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

LinearZeroSearcher &LinearZeroSearcher::setNearZeroPointChooser(function<Cloud(const Cloud &)> f) {
    nearZeroPointChooser = std::move(f);
    return *this;
}

void LinearZeroSearcher::init(vector<DAT> *cloud) {
    Cloud res = nearZeroPointChooser(*cloud);
    le.init(res);
    InterpolantBase3D::init(cloud);
}

std::array<double, 3> LinearZeroSearcher::operator()(const array<double, 3> &x) {
    m_interpolated = trust_region(x);
    return le(x);
}
