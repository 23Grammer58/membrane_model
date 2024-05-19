//
// Created by alex on 01.02.2021.
//

#ifndef AORTIC_VALVE_LINEARZEROSEACHER_H
#define AORTIC_VALVE_LINEARZEROSEACHER_H

#include <vector>
#include <array>
#include <set>
#include <cfloat>
#include <functional>
#include <fstream>
#include <iostream>

#include <Eigen/Dense>
#include <CGAL/QP_models.h>
#include <CGAL/QP_functions.h>
#include <CGAL/MP_Float.h>

#include "InterpolatorBase.h"
#include "MaterialResponse.h"

namespace World3d{

struct LinearElasticResponse: public ResponseInterpolant{
    /// Fit coeficients for linear elastic potential phi = 1/2 * C_ijkl * E_{ij} E_{kl}
    bool fit(const std::vector<Response>& data) override;

    /// Get response functions d phi / d xi_i
    std::array<double, 3> evaluate(const double* xi /*[3]*/);
    std::array<double, 3> evaluate(const std::array<double, 3>& xi) override { return evaluate(xi.data()); }
    /// Get response function derivatives r_ij = d^2 phi / (d xi_i d xi_j)
    /// return {r_00, r_11, r_22, r_01, r_02, r_12}
    std::array<double, 6> response_gradient(const std::array<double, 3>& xi);

    std::array<std::array<double, 3>, 3> gradient(const std::array<double, 3>& xi) override;
    std::unique_ptr<ResponseInterpolant> clone() const override { return std::make_unique<LinearElasticResponse>(*this); }
private: 
    std::array<double, 6> m_C = {0};
    void least_squares_fit(const std::vector<Response>& data);
    static bool is_positive(const std::array<double, 6>& C);
    void conditional_optimization(const std::vector<Response>& data);    
};

}

#endif //AORTIC_VALVE_LINEARZEROSEACHER_H