//
// Created by alex on 05.09.2021.
//

#ifndef AORTIC_VALVE_HGO_GSTMODEL_H
#define AORTIC_VALVE_HGO_GSTMODEL_H
#include <cmath>
#include "HyperElastic.h"

namespace World3d{

class HGO_GSTModel: public HyperElasticForceBase {
public:
    DReal _h, _mu, _k1, _k2, _kappa, _q_divg;
    DReal _theta = NAN;
    HGO_GSTModel(double h, double mu, double k1, double k2, double kappa, std::string gen_dir, bool regenerate = true, double q_divg = 50)
    {
        Init(h, mu, k1, k2, kappa, NAN, gen_dir, regenerate, q_divg, "f:aniso_f", "f:aniso_s");
    }
    HGO_GSTModel(double h, double mu, double k1, double k2, double kappa, double theta, std::string gen_dir, bool regenerate = true, double q_divg = 50)
    {
        Init(h, mu, k1, k2, kappa, theta, gen_dir, regenerate, q_divg);
    }
    std::unique_ptr<HyperElasticForceBase> copy() override{
        return std::make_unique<HGO_GSTModel>(*this);
    }
private:
    int Init(double h, double mu, double k1, double k2, double kappa, double theta, std::string gen_dir, bool regenerate, double q_divg,
             std::string aniso_tag_f = "f:aniso_f", std::string aniso_tag_s = "f:aniso_s"){
        _h = h, _mu = mu , _k1 = k1, _k2 = k2, _kappa = kappa, _theta = theta, _q_divg = q_divg;

        type = "HGO_GSTModel";
        SX H  = SX::sym("H" );
        SX Mu = SX::sym("Mu");
        SX K1  = SX::sym("k1" );
        SX K2  = SX::sym("k2" );
        SX kap = SX::sym("kap");
        SX reg = SX::sym("reg");
        auto lv = DefaultHyperElasticForce::getDefaultLocalVars();
        SX Ap = lv["Ap"];
        auto I = DefaultHyperElasticForce::getDefaultInvariants();
        DefaultHyperElasticForce::ArgsDef orth1, orth2;
        if (std::isnan(theta)){
            orth1 = DefaultHyperElasticForce::makeOrthotropicInvariant("I4s", aniso_tag_f);
            orth2 = DefaultHyperElasticForce::makeOrthotropicInvariant("I4f", aniso_tag_s);
        } else {
            std::array<DReal, 3> Av1{cos(theta), sin(theta), 0}, Av2{cos(theta), -sin(theta), 0};
            orth1 = DefaultHyperElasticForce::makeOrthotropicInvariant("I4s", Av1);
            orth2 = DefaultHyperElasticForce::makeOrthotropicInvariant("I4f", Av2);
        }
        SX I1 = I["I1"];
        SX I3 = I["I3"];
        SX I4s = get<2>(orth1)["I4s"], I4f = get<2>(orth2)["I4f"];
        SX I1_3d = I1 + 1.0/I3;
        SX kapI1_3d = kap*I1_3d, onem3kap = 1 - 3*kap;
        SX Es = kapI1_3d + onem3kap*I4s;
        SX Ef = kapI1_3d + onem3kap*I4f;
        SX x = SX::sym("x");
        SX power = SX::if_else_zero(x > 0, K2*SX::sq(x));
        SX expq = SX::exp(reg);
        SX power_Es = SX::substitute(power, x, Es - 1), power_Ef = SX::substitute(power, x, Ef - 1);
        auto res_func = [](SX x) { return 1 + x; };
        SX dphi_s = SX::if_else(power_Es <= reg, SX::exp(power_Es) - 1, expq * res_func(power_Es - reg) - 1);
        SX dphi_f = SX::if_else(power_Ef <= reg, SX::exp(power_Ef) - 1, expq * res_func(power_Ef - reg) - 1);
        //SX dphi_s = SX::if_else_zero(Es > 1, SX::exp(K2*SX::sq(Es - 1)) - 1);
        //SX dphi_f = SX::if_else_zero(Ef > 1, SX::exp(K2*SX::sq(Ef - 1)) - 1);
        SX phi_iso = Mu/2 * (I1_3d - 3);
        SX U = Ap * H * (phi_iso + K1/(2*K2) * (dphi_s + dphi_f));
        auto iv = DefaultHyperElasticForce::getDefaultAvailableInputVars();
        iv.emplace_back("H", H, std::make_unique<ConstDataRequier>(_h));
        iv.emplace_back("Mu", Mu, std::make_unique<ConstDataRequier>(_mu));
        iv.emplace_back("k1", K1, std::make_unique<ConstDataRequier>(_k1));
        iv.emplace_back("k2", K2, std::make_unique<ConstDataRequier>(_k2));
        iv.emplace_back("kap", kap, std::make_unique<ConstDataRequier>(_kappa));
        iv.emplace_back("reg", reg, std::make_unique<ConstDataRequier>(_q_divg));

        HH::ArgsDef args(std::move(iv), std::move(lv), std::move(I));
        args = HH::concat(args, orth1);
        args = HH::concat(args, orth2);
        f = DefaultHyperElasticForce("HGO_GSTModel", std::move(gen_dir), std::move(args), U, regenerate);
        return 0;
    }
};

class MurdockModel: public HyperElasticForceBase {
public:
    DReal _h, _c1, _c2, _k1, _k2, _kappa, _q_divg;
    DReal _theta = NAN;
    MurdockModel(double h, double c1, double c2, double k1, double k2, double kappa, std::string gen_dir, bool regenerate = true, double q_divg = 50)
    {
        Init(h, c1, c2, k1, k2, kappa, NAN, gen_dir, regenerate, q_divg, "f:aniso_f", "f:aniso_s");
    }
    MurdockModel(double h, double c1, double c2, double k1, double k2, double kappa, double theta, std::string gen_dir, bool regenerate = true, double q_divg = 50)
    {
        Init(h, c1, c2, k1, k2, kappa, theta, gen_dir, regenerate, q_divg);
    }
    std::unique_ptr<HyperElasticForceBase> copy() override{
        return std::make_unique<MurdockModel>(*this);
    }
private:
    int Init(double h, double c1, double c2, double k1, double k2, double kappa, double theta, std::string gen_dir, bool regenerate, double q_divg,
             std::string aniso_tag_f = "f:aniso_f", std::string aniso_tag_s = "f:aniso_s"){
        _h = h, _c1 = c1, _c2 = c2, _k1 = k1, _k2 = k2, _kappa = kappa, _theta = theta, _q_divg = q_divg;

        type = "MurdockModel";
        SX H  = SX::sym("H" );
        SX C1 = SX::sym("c1");
        SX C2 = SX::sym("c2");
        SX K1  = SX::sym("k1" );
        SX K2  = SX::sym("k2" );
        SX kap = SX::sym("kap");
        SX reg = SX::sym("reg");
        auto lv = DefaultHyperElasticForce::getDefaultLocalVars();
        SX Ap = lv["Ap"];
        auto I = DefaultHyperElasticForce::getDefaultInvariants();
        DefaultHyperElasticForce::ArgsDef orth1, orth2;
        if (std::isnan(theta)){
            orth1 = DefaultHyperElasticForce::makeOrthotropicInvariant("I4s", aniso_tag_s);
            orth2 = DefaultHyperElasticForce::makeOrthotropicInvariant("I4f", aniso_tag_f);
        } else {
            std::array<DReal, 3> Av1{cos(theta), sin(theta), 0}, Av2{cos(theta), -sin(theta), 0};
            orth1 = DefaultHyperElasticForce::makeOrthotropicInvariant("I4s", Av1);
            orth2 = DefaultHyperElasticForce::makeOrthotropicInvariant("I4f", Av2);
        }
        SX I1 = I["I1"];
        SX I3 = I["I3"];
        SX I4s = get<2>(orth1)["I4s"], I4f = get<2>(orth2)["I4f"];
        SX I1_3d = I1 + 1.0/I3;
        SX kapI1_3d = kap*I1_3d, onem3kap = 1 - 3*kap;
        SX Es = kapI1_3d + onem3kap*I4s;
        SX Ef = kapI1_3d + onem3kap*I4f;
        SX x = SX::sym("x");
        SX power = SX::if_else_zero(x > 0, K2*SX::sq(x));
        SX expq = SX::exp(reg);
        SX power_Es = SX::substitute(power, x, Es - 1), power_Ef = SX::substitute(power, x, Ef - 1);
        auto res_func = [](SX x) { return 1 + x; };
        SX dphi_s = SX::if_else(power_Es <= reg, SX::exp(power_Es) - 1, expq * res_func(power_Es - reg) - 1);
        SX dphi_f = SX::if_else(power_Ef <= reg, SX::exp(power_Ef) - 1, expq * res_func(power_Ef - reg) - 1);
        SX phi_iso = C1 * (exp(C2*(I1_3d - 3)) - 1);
        SX U = Ap * H * (phi_iso + K1/(2*K2) * (dphi_s + dphi_f));
        auto iv = DefaultHyperElasticForce::getDefaultAvailableInputVars();
        iv.emplace_back("H", H, std::make_unique<ConstDataRequier>(_h));
        iv.emplace_back("c1", C1, std::make_unique<ConstDataRequier>(_c1));
        iv.emplace_back("c2", C2, std::make_unique<ConstDataRequier>(_c2));
        iv.emplace_back("k1", K1, std::make_unique<ConstDataRequier>(_k1));
        iv.emplace_back("k2", K2, std::make_unique<ConstDataRequier>(_k2));
        iv.emplace_back("kap", kap, std::make_unique<ConstDataRequier>(_kappa));
        iv.emplace_back("reg", reg, std::make_unique<ConstDataRequier>(_q_divg));

        HH::ArgsDef args(std::move(iv), std::move(lv), std::move(I));
        args = HH::concat(args, orth1);
        args = HH::concat(args, orth2);
        f = DefaultHyperElasticForce("MurdockModel", std::move(gen_dir), std::move(args), U, regenerate);
        return 0;
    }
};

}

#endif //AORTIC_VALVE_HGO_GSTMODEL_H
