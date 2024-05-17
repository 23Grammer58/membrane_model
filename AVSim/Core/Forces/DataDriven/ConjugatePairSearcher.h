//
// Created by alex on 01.02.2021.
//

#ifndef AORTIC_VALVE_CONJUGATEPAIRSEARCHER_H
#define AORTIC_VALVE_CONJUGATEPAIRSEARCHER_H

#include <vector>
#include <array>
#include <set>
#include <cfloat>
#include <functional>
#include <fstream>
#include <iostream>
#include <gsl/gsl_multifit.h>
#include "InterpolatorBase.h"
#include "ResponseStatistics.h"

std::vector<double> polynomialfit(const std::set<int>& templ, int data_size, const std::function<double(int)>& xd, const std::function<double(int)>& yd);
std::vector<double> polynomialfit(int deg, int data_size, const std::function<double(int)>& xd, const std::function<double(int)>& yd);
std::vector<double> polynomialfit_zero(int deg, int data_size, const std::function<double(int)>& xd, const std::function<double(int)>& yd);


class ConjugatePairSearcher: public InterpolantBase3D{
public:
    using DAT = ResponseStatistics::Elem;
    using Cloud = std::vector<DAT>;

    std::vector<double> m_piC, m_sigmaC, m_tauC;

    bool m_interpolated = false;
    std::function<bool(const std::array<double, 3>& x)> trust_region = [](const std::array<double, 3>& x) { return true; };
    ConjugatePairSearcher& setTrustRegion(std::function<bool(const std::array<double, 3>& )> f) {
        return trust_region = std::move(f), *this;
    }

    std::string getType() const override { return "ConjugatePairSearcher"; }
    void init(std::vector<DAT>* cloud) override;
    bool interpolated() { return m_interpolated; }
    std::array<double, 3> operator()(const std::array<double, 3>& x) override;
    InterpolantBase* clone() override{ return new ConjugatePairSearcher(*this); }

    std::string params_to_str();
    std::string pi_points_to_str();
    std::string sigma_points_to_str();
    std::string tau_points_to_str();
    std::string to_str();
    void info(std::ostream& out) { out << to_str() << std::endl; }
    void info(std::string dir);
private:
    static double _comp(const std::vector<double>& coefs, double x){
        double res = 0 * x;
        for (int i = coefs.size() - 1; i >= 0; --i) res = (res + coefs[i]) * x;
        return res;
    }
};


#endif //AORTIC_VALVE_CONJUGATEPAIRSEARCHER_H
