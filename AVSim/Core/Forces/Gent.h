//
// Created by alex on 29.01.2021.
//

#ifndef AORTIC_VALVE_GENT_H
#define AORTIC_VALVE_GENT_H

#include "HyperElastic.h"
namespace World3d{
    class GentModel: public HyperElasticForceBase {
    public:
        DReal _mu, _h, _Jm = 2.3;
        GentModel(double mu, double h, double Jm, std::string gen_dir, bool regenerate = true):
                _mu{mu}, _h{h}, _Jm{Jm}
        {
            type = "GentModel";
            SX Mu = SX::sym("mu");
            SX H  = SX::sym("H" );
            SX jm  = SX::sym("Jm" );
            auto lv = DefaultHyperElasticForce::getDefaultLocalVars();
            SX Ap = lv["Ap"];
            auto I = DefaultHyperElasticForce::getDefaultInvariants();
            SX I1 = I["I1"];
            SX J = I["J"];
            SX U = -Ap * H * Mu * jm/ 2 * SX::log(1 - (I1 + SX::pow(J, -2) - 3) / jm);
            auto iv = DefaultHyperElasticForce::getDefaultAvailableInputVars();
            iv.emplace_back("mu", Mu, std::make_unique<ConstDataRequier>(_mu));
            iv.emplace_back("H", H, std::make_unique<ConstDataRequier>(_h));
            iv.emplace_back("Jm", jm, std::make_unique<ConstDataRequier>(_Jm));
            f = DefaultHyperElasticForce("Gent", gen_dir, std::move(iv), lv, I, U, regenerate);
        }
        std::unique_ptr<HyperElasticForceBase> copy() override{
            return std::make_unique<GentModel>(_mu, _h, _Jm, f.gen_dir, false);
        }
    };
}

#endif //AORTIC_VALVE_GENT_H
