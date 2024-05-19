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
#include "MaterialResponse.h"

std::vector<double> polynomialfit(const std::set<int>& templ, int data_size, const std::function<double(int)>& xd, const std::function<double(int)>& yd);
std::vector<double> polynomialfit(int deg, int data_size, const std::function<double(int)>& xd, const std::function<double(int)>& yd);
std::vector<double> polynomialfit_zero(int deg, int data_size, const std::function<double(int)>& xd, const std::function<double(int)>& yd);

namespace World3d{
class ConjugatePairSearcher: public ResponseInterpolant{
public:
    using DAT = Response;
    using Cloud = std::vector<DAT>;

    std::vector<double> m_piC, m_sigmaC, m_tauC;

    std::string getType() const { return "ConjugatePairSearcher"; }
    bool fit(const std::vector<DAT>& cloud) override;
    std::array<double, 3> evaluate(const std::array<double, 3>& x) override;
    std::array<std::array<double, 3>, 3> gradient(const std::array<double, 3>& xi) override { throw std::runtime_error("gradient() for conjugate pairs currently is not implemented"); }
    std::unique_ptr<ResponseInterpolant> clone() const override{ return std::make_unique<ConjugatePairSearcher>(*this); }

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
    const std::vector<DAT>* m_cloud = nullptr;
};

}


#endif //AORTIC_VALVE_CONJUGATEPAIRSEARCHER_H
