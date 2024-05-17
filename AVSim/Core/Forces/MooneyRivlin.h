//
// Created by alex on 29.01.2021.
//

#ifndef AORTIC_VALVE_MOONEYRIVLIN_H
#define AORTIC_VALVE_MOONEYRIVLIN_H

#include "HyperElastic.h"
namespace World3d {
    class MooneyRivlinModel : public HyperElasticForceBase {
    public:
        DReal _mu, _h, _alpha = 0.1, _beta = 0;
        std::array<DReal, 3> _Av = {NAN, NAN, NAN};

        MooneyRivlinModel(double h, double mu, double alpha, std::string gen_dir, bool regenerate = true) :
                _mu{mu}, _h{h}, _alpha{alpha} {
            type = "MooneyRivlinModel";
            SX Mu = SX::sym("mu");
            SX H = SX::sym("H");
            SX A = SX::sym("al");
            auto lv = DefaultHyperElasticForce::getDefaultLocalVars();
            SX Ap = lv["Ap"];
            auto I = DefaultHyperElasticForce::getDefaultInvariants();
            SX I1 = I["I1"];
            SX J = I["J"];
            SX U = Ap * H * Mu / 2 * ((I1 + SX::pow(J, -2) - 3) + A * (SX::pow(J, 2) + I1 * SX::pow(J, -2) - 3));
            auto iv = DefaultHyperElasticForce::getDefaultAvailableInputVars();
            iv.emplace_back("mu", Mu, std::make_unique<ConstDataRequier>(_mu));
            iv.emplace_back("H", H, std::make_unique<ConstDataRequier>(_h));
            iv.emplace_back("al", A, std::make_unique<ConstDataRequier>(_alpha));
            f = DefaultHyperElasticForce("Mooney", gen_dir, std::move(iv), lv, I, U, regenerate);
        }

        MooneyRivlinModel(double h, double mu, double alpha, double beta, std::array<double, 3> Av, std::string gen_dir,
                          bool regenerate = true, std::string aniso_tag_s = "f:aniso_s") :
                _mu{mu}, _h{h}, _alpha{alpha}, _beta{beta}, _Av{Av} {
            type = "MooneyRivlinModel";
            SX Mu = SX::sym("mu");
            SX H = SX::sym("H");
            SX A = SX::sym("al");
            SX B = SX::sym("bet");
            auto lv = DefaultHyperElasticForce::getDefaultLocalVars();
            SX Ap = lv["Ap"];
            auto I = DefaultHyperElasticForce::getDefaultInvariants();
            DefaultHyperElasticForce::ArgsDef orth;
            if (std::isnan(_Av[0])){
                orth = DefaultHyperElasticForce::makeOrthotropicInvariant("I4s", aniso_tag_s);
            }
            else
                orth = DefaultHyperElasticForce::makeOrthotropicInvariant("I4s", _Av);
            SX I1 = I["I1"];
            SX I3 = I["I3"];
            SX I4s = get<2>(orth)["I4s"];
            SX U = Ap * H * Mu * ((I1 + 1.0 / I3 - 3) + A * (I3 + I1 / I3 - 3) + B * /*SX::sq*/(I4s - 1));
            auto iv = DefaultHyperElasticForce::getDefaultAvailableInputVars();
            iv.emplace_back("mu", Mu, std::make_unique<ConstDataRequier>(_mu));
            iv.emplace_back("H", H, std::make_unique<ConstDataRequier>(_h));
            iv.emplace_back("al", A, std::make_unique<ConstDataRequier>(_alpha));
            iv.emplace_back("bet", B, std::make_unique<ConstDataRequier>(_beta));

            HH::ArgsDef args(std::move(iv), std::move(lv), std::move(I));
            args = HH::concat(args, std::move(orth));
            f = DefaultHyperElasticForce("MooneyOrth", gen_dir, std::move(args), U, regenerate);
        }

        std::unique_ptr<HyperElasticForceBase> copy() override {
            if (f.force_name == "MooneyOrth")
                return std::make_unique<MooneyRivlinModel>(_h, _mu, _alpha, _beta, _Av, f.gen_dir, false);
            else if (f.force_name == "Mooney")
                return std::make_unique<MooneyRivlinModel>(_h, _mu, _alpha, f.gen_dir, false);
            else
                throw std::runtime_error("Wrong type");
        }
    };

    class MooneyRivlinModel2 : public HyperElasticForceBase {
    public:
        DReal _mu, _h, _alpha = 0.1, _beta = 0;
        std::array<DReal, 3> _Av = {0, 0, 0};

        MooneyRivlinModel2(double h, double mu, double alpha, std::string gen_dir, bool regenerate = true) :
                _mu{mu}, _h{h}, _alpha{alpha} {
            type = "MooneyRivlinModel2";
            SX Mu = SX::sym("mu");
            SX H = SX::sym("H");
            SX A = SX::sym("al");
            auto lv = DefaultHyperElasticForce::getDefaultLocalVars();
            SX Ap = lv["Ap"];
            auto I = DefaultHyperElasticForce::getDefaultInvariants();
            SX I1 = I["I1"];
            SX J = I["J"];
            SX U = Ap * H * Mu / 2 * ((I1 + SX::pow(J, -2) - 3) + A * (SX::pow(J, 2) + I1 * SX::pow(J, -2) - 3));
            auto iv = DefaultHyperElasticForce::getDefaultAvailableInputVars();
            iv.emplace_back("mu", Mu, std::make_unique<ConstDataRequier>(_mu));
            iv.emplace_back("H", H, std::make_unique<ConstDataRequier>(_h));
            iv.emplace_back("al", A, std::make_unique<ConstDataRequier>(_alpha));
            f = DefaultHyperElasticForce("Mooney2", gen_dir, std::move(iv), lv, I, U, regenerate);
        }

        MooneyRivlinModel2(double h, double mu, double alpha, double beta, std::array<double, 3> Av,
                           std::string gen_dir, bool regenerate = true, std::string aniso_tag_s = "f:aniso_s") :
                _mu{mu}, _h{h}, _alpha{alpha}, _beta{beta}, _Av{Av} {
            type = "MooneyRivlinModel2";
            SX Mu = SX::sym("mu");
            SX H = SX::sym("H");
            SX A = SX::sym("al");
            SX B = SX::sym("bet");
            auto lv = DefaultHyperElasticForce::getDefaultLocalVars();
            SX Ap = lv["Ap"];
            auto I = DefaultHyperElasticForce::getDefaultInvariants();
            DefaultHyperElasticForce::ArgsDef orth;
            if (std::isnan(_Av[0])){
                orth = DefaultHyperElasticForce::makeOrthotropicInvariant("I4s", aniso_tag_s);
            }
            else
                orth = DefaultHyperElasticForce::makeOrthotropicInvariant("I4s", _Av);
            SX I1 = I["I1"];
            SX I3 = I["I3"];
            SX I4s = get<2>(orth)["I4s"];
            SX U = Ap * H * Mu * ((I1 + 1.0 / I3 - 3) + A * (I3 + I1 / I3 - 3) + B * SX::sq(I4s - 1));
            auto iv = DefaultHyperElasticForce::getDefaultAvailableInputVars();
            iv.emplace_back("mu", Mu, std::make_unique<ConstDataRequier>(_mu));
            iv.emplace_back("H", H, std::make_unique<ConstDataRequier>(_h));
            iv.emplace_back("al", A, std::make_unique<ConstDataRequier>(_alpha));
            iv.emplace_back("bet", B, std::make_unique<ConstDataRequier>(_beta));

            HH::ArgsDef args(std::move(iv), std::move(lv), std::move(I));
            args = HH::concat(args, std::move(orth));
            f = DefaultHyperElasticForce("MooneyOrth2", gen_dir, std::move(args), U, regenerate);
        }

        std::unique_ptr<HyperElasticForceBase> copy() override {
            if (f.force_name == "MooneyOrth2")
                return std::make_unique<MooneyRivlinModel2>(_h, _mu, _alpha, _beta, _Av, f.gen_dir, false);
            else if (f.force_name == "Mooney2")
                return std::make_unique<MooneyRivlinModel2>(_h, _mu, _alpha, f.gen_dir, false);
            else
                throw std::runtime_error("Wrong type");
        }
    };
}
#endif //AORTIC_VALVE_MOONEYRIVLIN_H
