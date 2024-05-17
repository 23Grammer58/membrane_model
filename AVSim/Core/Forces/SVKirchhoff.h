//
// Created by alex on 29.01.2021.
//

#ifndef AORTIC_VALVE_SVKIRCHHOFF_H
#define AORTIC_VALVE_SVKIRCHHOFF_H

#include "HyperElastic.h"
namespace World3d{
    class SVKirchhoffModel: public HyperElasticForceBase{
    public:
        DReal _lambda, _mu, _h;
        SVKirchhoffModel(double h, double lambda, double mu, std::string gen_dir, bool regenerate = true):
                _h{h}, _lambda{lambda}, _mu{mu}
        {
            type = "SVKirchhoffModel";
            SX L = SX::sym("lambda");
            SX Mu = SX::sym("mu");
            SX H  = SX::sym("H" );
            auto lv = DefaultHyperElasticForce::getDefaultLocalVars();
            SX Ap = lv["Ap"];
            auto I = DefaultHyperElasticForce::getDefaultInvariants();
            SX I1 = I["I1"];
//            SX J = I["J"];
            SX I3 = I["I3"];
            SX U = Ap * H * (L / 8 * SX::sq(I1 + 1.0/I3 - 3) + Mu / 4 * (SX::sq(I1 - 1) + SX::sq(1.0/I3 - 1) - 2*I3 + 1));
            auto iv = DefaultHyperElasticForce::getDefaultAvailableInputVars();
            iv.emplace_back("lambda", L, std::make_unique<ConstDataRequier>(_lambda));
            iv.emplace_back("mu", Mu, std::make_unique<ConstDataRequier>(_mu));
            iv.emplace_back("H", H, std::make_unique<ConstDataRequier>(_h));
            f = DefaultHyperElasticForce("SVKirchhoff", gen_dir, std::move(iv), lv, I, U, regenerate);
        }
        std::unique_ptr<HyperElasticForceBase> copy() override{
            return std::make_unique<SVKirchhoffModel>(_h, _lambda, _mu, f.gen_dir, false);
        }
    };

    class SVKirchhoffNoPoissonModel: public HyperElasticForceBase{
    public:
        DReal _lambda, _mu, _h;
        SVKirchhoffNoPoissonModel(double h, double lambda, double mu, std::string gen_dir, bool regenerate = true):
                _h{h}, _lambda{lambda}, _mu{mu}
        {
            type = "SVKirchhoffNoPoissonModel";
            SX L = SX::sym("lambda");
            SX Mu = SX::sym("mu");
            SX H  = SX::sym("H" );
            auto lv = DefaultHyperElasticForce::getDefaultLocalVars();
            SX Ap = lv["Ap"];
            auto I = DefaultHyperElasticForce::getDefaultInvariants();
            SX I1 = I["I1"];
//            SX J = I["J"];
            SX I3 = I["I3"];
            SX U = Ap * H * (Mu * L / (L + 2*Mu) * SX::sq(I1/2 - 1) + Mu/4 * (SX::sq(I1-1) - 2*I3 + 1));
            auto iv = DefaultHyperElasticForce::getDefaultAvailableInputVars();
            iv.emplace_back("lambda", L, std::make_unique<ConstDataRequier>(_lambda));
            iv.emplace_back("mu", Mu, std::make_unique<ConstDataRequier>(_mu));
            iv.emplace_back("H", H, std::make_unique<ConstDataRequier>(_h));
            f = DefaultHyperElasticForce("SVKirchhoffNoPoisson", gen_dir, std::move(iv), lv, I, U, regenerate);
        }
        std::unique_ptr<HyperElasticForceBase> copy() override{
            return std::make_unique<SVKirchhoffNoPoissonModel>(_h, _lambda, _mu, f.gen_dir, false);
        }
    };
}

#endif //AORTIC_VALVE_SVKIRCHHOFF_H
