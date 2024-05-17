//
// Created by alex on 01.02.2021.
//

#include "ConjugatePairSearcher.h"

std::vector<double> polynomialfit(const std::set<int> &templ, int data_size, const std::function<double(int)> &xd,
                                  const std::function<double(int)> &yd) {
    gsl_multifit_linear_workspace *ws;
    gsl_matrix *cov, *A;
    gsl_vector *y, *c;
    double chisq;
    int deg = templ.size();

    A = gsl_matrix_alloc(data_size, deg);
    y = gsl_vector_alloc(data_size);
    c = gsl_vector_alloc(deg);
    cov = gsl_matrix_alloc(deg, deg);

    for(int i=0; i < data_size; i++) {
        int j = 0;
        for (const auto& l: templ)
            if (l != 0)
                gsl_matrix_set(A, i, j++, pow(xd(i), l));
            else
                gsl_matrix_set(A, i, j++, 1);

        gsl_vector_set(y, i, yd(i));
    }

    ws = gsl_multifit_linear_alloc(data_size, deg);
    gsl_multifit_linear(A, y, c, cov, &chisq, ws);

    std::vector<double> res(deg);
    for(int i=0; i < deg; i++)
        res[i] = gsl_vector_get(c, i);

    gsl_multifit_linear_free(ws);
    gsl_matrix_free(A);
    gsl_matrix_free(cov);
    gsl_vector_free(y);
    gsl_vector_free(c);

    return res;
}

std::vector<double>
polynomialfit(int deg, int data_size, const std::function<double(int)> &xd, const std::function<double(int)> &yd) {
    std::set<int> templ;
    for (int i = 0; i <= deg; ++i)
        templ.insert(i);
    return polynomialfit(templ, data_size, xd, yd);
}

std::vector<double>
polynomialfit_zero(int deg, int data_size, const std::function<double(int)> &xd, const std::function<double(int)> &yd) {
    std::set<int> templ;
    for (int i = 1; i <= deg; ++i)
        templ.insert(i);
    return polynomialfit(templ, data_size, xd, yd);
}

void ConjugatePairSearcher::init(std::vector<DAT> *cloud) {
    if (cloud->size() >= 5) {
        m_piC = polynomialfit_zero(4, cloud->size(),
                                   [&cloud = *cloud](int i) { return 0.5 * (cloud[i].ksi[0] + cloud[i].ksi[1]); },
                                   [&cloud = *cloud](int i) { return (cloud[i].response[0] + cloud[i].response[1]); });
        m_sigmaC = polynomialfit_zero(5, cloud->size(),
                                      [&cloud = *cloud](int i) { return 0.5 * (cloud[i].ksi[0] - cloud[i].ksi[1]); },
                                      [&cloud = *cloud](int i) { return (cloud[i].response[0] - cloud[i].response[1]); });
        m_tauC = polynomialfit_zero(5, cloud->size(),
                                    [&cloud = *cloud](int i) { return cloud[i].ksi[2]; },
                                    [&cloud = *cloud](int i) { return cloud[i].response[2]; });
    }
    else throw std::runtime_error("Interpolant::ConjugatePairFit: cloud has too few points");
    InterpolantBase3D::init(cloud);
}

std::array<double, 3> ConjugatePairSearcher::operator()(const std::array<double, 3> &x) {
    double delta = 0.5 * (x[0] + x[1]), eps = 0.5 * (x[0] - x[1]), gamma = x[2];
    double pi = _comp(m_piC, delta), sigma = _comp(m_sigmaC, eps), tau = _comp(m_tauC, gamma);
    std::array<double, 3> res;
    res[0] = 0.5 * (pi + sigma);
    res[1] = 0.5 * (pi - sigma);
    res[2] = tau;
    m_interpolated = trust_region(x);
    return res;
}

std::string ConjugatePairSearcher::params_to_str() {
    std::ostringstream ss;
    ss << "pi: f(x) = 0";
    int n;
    n = 0;
    for (auto i: m_piC)
        ss << " + " << i << " * x**" << ++n;
    ss << std::endl;
    ss << "sigma: f(x) = 0";
    n = 0;
    for (auto i: m_sigmaC)
        ss << " + " << i << " * x**" << ++n;
    ss << std::endl;
    ss << "tau: f(x) = 0";
    n = 0;
    for (auto i: m_tauC)
        ss << " + " << i << " * x**" << ++n;
    ss << std::endl;
    return ss.str();
}

std::string ConjugatePairSearcher::pi_points_to_str() {
    std::ostringstream ss;
    ss << "pi points:" << std::endl;
    auto& cloud  = *InterpolantBase3D::m_cloud;
    for (int i = 0; i < cloud.size(); ++i){
        ss  << 0.5 * (cloud[i].ksi[0] + cloud[i].ksi[1]) << ", "
            << (cloud[i].response[0] + cloud[i].response[1]) << std::endl;
    }
    ss << std::endl;
    return ss.str();
}

std::string ConjugatePairSearcher::sigma_points_to_str() {
    std::ostringstream ss;
    auto& cloud  = *InterpolantBase3D::m_cloud;
    ss << "sigma points:" << std::endl;
    for (int i = 0; i < cloud.size(); ++i){
        ss  << 0.5 * (cloud[i].ksi[0] - cloud[i].ksi[1]) << ", "
            << (cloud[i].response[0] - cloud[i].response[1]) << std::endl;
    }
    ss << std::endl;
    return ss.str();
}

std::string ConjugatePairSearcher::tau_points_to_str() {
    std::ostringstream ss;
    auto& cloud  = *InterpolantBase3D::m_cloud;
    ss << "tau points:" << std::endl;
    for (int i = 0; i < cloud.size(); ++i){
        ss  << cloud[i].ksi[2] << ", "
            << cloud[i].response[2]  << std::endl;
    }
    return ss.str();
}

std::string ConjugatePairSearcher::to_str() {
    return params_to_str() + pi_points_to_str() + sigma_points_to_str() + tau_points_to_str();
}

void ConjugatePairSearcher::info(std::string dir) {
    std::ofstream f(dir + "functions.txt", std::ios::trunc);
    f << params_to_str();
    f.close();
    f.open(dir + "pi_pnts.dat");
    f << pi_points_to_str();
    f.close();
    f.open(dir + "sgma_pnts.dat");
    f << sigma_points_to_str();
    f.close();
    f.open(dir + "tau_pnts.dat");
    f << tau_points_to_str();
    f.close();
}
