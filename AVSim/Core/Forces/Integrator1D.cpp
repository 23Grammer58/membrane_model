//
// Created by alex on 29.01.2021.
//

#include "Integrator1D.h"

std::vector<std::pair<long double, long double>> Num_integrator1D::get_ksi_w(int min_order) {
    std::vector<std::pair<long double, long double>> ksi_w;
    if (min_order <= 1)
        ksi_w.emplace_back(0.0L, 2.0L);
    else if (min_order <= 2){
        ksi_w.emplace_back(-1.0L, 1.0L);
        ksi_w.emplace_back(1.0L, 1.0L);
    }
    else if (min_order <= 3){
        ksi_w.emplace_back(-1.0L, 1.0L / 3);
        ksi_w.emplace_back( 0.0L, 4.0L / 3);
        ksi_w.emplace_back( 1.0L, 1.0L / 3);
    }
    else if (min_order <= 5){
        static long double k = sqrt(0.6L);
        ksi_w.emplace_back(-k, 5.0L / 9);
        ksi_w.emplace_back( 0.0L, 8.0L / 9);
        ksi_w.emplace_back( k, 5.0L / 9);
    }
    else if (min_order <= 7){
        static long double k = sqrt(1.2L);
        static long double w = sqrt(30L) / 36;
        static long double k1 = sqrt((3.0L - 2.0L * k) / 7);
        static long double k2 = sqrt((3.0L + 2.0L * k) / 7);
        ksi_w.emplace_back(-k1, 0.5L + w);
        ksi_w.emplace_back(k1, 0.5L + w);
        ksi_w.emplace_back(-k2, 0.5L - w);
        ksi_w.emplace_back(k2, 0.5L - w);
    }
    else if (min_order <= 9){
        static long double w1 = 13*sqrt(70L) / 900;
        static long double w2 = 322.0L/900;
        static long double k = 2*sqrt(10.0L/7);
        static long double k1 = 1.0L/3 * sqrt(5L - k);
        static long double k2 = 1.0L/3 * sqrt(5L + k);
        ksi_w.emplace_back(0, 128.0L/225);
        ksi_w.emplace_back(-k1, w2+w1);
        ksi_w.emplace_back(k1, w2+w1);
        ksi_w.emplace_back(-k2, w2-w1);
        ksi_w.emplace_back(k2, w2-w1);
    } else if (min_order <= 17) {
        static double k[] = {
                0.968160239507626,
                0.836031107326636,
                0.613371432700590,
                0.324253423403809
        };
        static double w[] = {
                0.081274388361574,
                0.180648160694857,
                0.260610696402935,
                0.312347077040003
        };
        ksi_w.emplace_back(0, 0.330239355001260);
        for (int i = 0; i < 4; ++i){
            ksi_w.emplace_back(k[i], w[i]);
            ksi_w.emplace_back(-k[i], w[i]);
        }
    } else
        throw std::runtime_error("This is not implemented yet");
    return ksi_w;

}

Num_integrator1D::SX
Num_integrator1D::integrate(const Num_integrator1D::SX &f, const Num_integrator1D::SX &x,
                                     const Num_integrator1D::SX a, const Num_integrator1D::SX b,
                                     int min_order) {
    assert((a.size1() == a.size2()) &&  (b.size1() == b.size2())  && (a.size1() == 1));
    SX detJ = (b - a) / 2;
    auto getX = [&a, &b](double ksi) { return (b - a) * (ksi + 1) / 2 + a; };
    auto F = [&f, &x, &getX](double ksi) { return SX::simplify(SX::substitute(f, x, getX(ksi))); };
    SX S(f.size1(), f.size2());
    auto ksi_w = get_ksi_w(min_order);
    for (auto& i: ksi_w)
        S += i.second * F(i.first);
    return detJ * SX::simplify(S);
}