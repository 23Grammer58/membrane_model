//
// Created by alex on 29.01.2021.
//

#ifndef AORTIC_VALVE_INTEGRATOR1D_H
#define AORTIC_VALVE_INTEGRATOR1D_H

#include <vector>
#include <cassert>
#include "casadi/casadi.hpp"
#include <functional>

namespace Num_integrator1D{
    std::vector<std::pair<long double, long double>> get_ksi_w(int min_order);
    template<typename RealType, typename FRetType>
    FRetType integrate(const std::function<FRetType(RealType)>& f, RealType a, RealType b, int min_order){
        assert((a.size1() == a.size2()) &&  (b.size1() == b.size2())  && (a.size1() == 1));
        RealType detJ = (b - a) / 2;
        auto getX = [&a, &b](double ksi) { return (b - a) * (ksi + 1) / 2 + a; };
        auto F = [&f, &getX](double ksi) { return f(getX(ksi)); };
        auto ksi_w = get_ksi_w(min_order);
        auto S = ksi_w[0].second * F(ksi_w[0].first);
        for (int i = 1; i < ksi_w.size(); ++i)
            S += ksi_w[i].second * F(ksi_w[i].first);
        return detJ * S;
    }
    using SX = casadi::SX;
    SX integrate(const SX& f, const SX& x, const SX a, const SX b, int min_order);
};


#endif //AORTIC_VALVE_INTEGRATOR1D_H
