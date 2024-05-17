//
// Created by alex on 01.02.2021.
//

#ifndef AORTIC_VALVE_LINEARZEROSEACHER_H
#define AORTIC_VALVE_LINEARZEROSEACHER_H

#define USE_EIGEN_LINSQ

#include <vector>
#include <array>
#include <set>
#include <cfloat>
#include <functional>
#include <fstream>
#include <iostream>

#ifdef USE_EIGEN_LINSQ
#include <Eigen/Dense>
#else
#include <gsl/gsl_multifit.h>
#endif

#include <CGAL/QP_models.h>
#include <CGAL/QP_functions.h>
#include <CGAL/MP_Float.h>

#include "InterpolatorBase.h"
#include "ResponseStatistics.h"

struct LinearElastisity{
    using cloud_t = std::vector<ResponseStatistics::Elem>;

    void init(const cloud_t& data);
    std::array<double, 3> operator()(const std::array<double, 3>& x);
    std::array<double, 3> operator()(const double* x);

private:
    std::array<double, 6> m_C = {0};
    void least_squares_fit(const cloud_t& data);
    static bool is_positive(const std::array<double, 6>& C);
    void conditional_optimization(const cloud_t& data);
};

class LinearZeroSearcher: public InterpolantBase3D{
public:
    using DAT = ResponseStatistics::Elem;
    using Cloud = std::vector<DAT>;
    LinearElastisity le;
    bool m_interpolated = false;

    std::function<Cloud(const Cloud& )> nearZeroPointChooser = [](const Cloud& cloud)-> Cloud {
        Cloud res;
        auto mes = [](std::array<double, 3> p) { return p[0] * p[0] + p[1]*p[1] + p[2]*p[2];};
        for (const auto& i: cloud)
            if (mes(i.ksi) < 0.04 * 0.04)
                res.push_back(i);

        return res;
    };

    std::function<bool(const std::array<double, 3>& x)> trust_region = [](const std::array<double, 3>& x) { return true; };
    LinearZeroSearcher& setTrustRegion(std::function<bool(const std::array<double, 3>& )> f) {
        trust_region = std::move(f);
        return *this;
    }
    LinearZeroSearcher& setNearZeroPointChooser(std::function<Cloud(const Cloud& )> f);
    std::string getType() const override { return "LinearZeroSearcher"; }
    void init(std::vector<DAT>* cloud) override;
    bool interpolated() override { return m_interpolated; }
    std::array<double, 3> operator()(const std::array<double, 3>& x) override;
    InterpolantBase* clone(){ return new LinearZeroSearcher(*this); }
};




#endif //AORTIC_VALVE_LINEARZEROSEACHER_H
