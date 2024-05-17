//
// Created by alex on 29.01.2021.
//

#ifndef AORTIC_VALVE_HOLZAPFELOGDENGASSER_H
#define AORTIC_VALVE_HOLZAPFELOGDENGASSER_H

#include "HyperElastic.h"
namespace World3d{
    class HGOModel: public HyperElasticForceBase {
    public:
        DReal _h, _c, _k1, _k2, _kappa;
        std::array<DReal, 3> _Av;
        HGOModel(double h, double c, double k1, double k2, double kappa, std::array<double, 3> _Av, std::string gen_dir, bool regenerate = true):
                _h{h}, _c{c}, _k1{k1}, _k2{k2}, _kappa{kappa}, _Av{_Av}
        {
            type = "HGOModel";
            SX H  = SX::sym("H" );
            SX C = SX::sym("C");
            SX K1  = SX::sym("k1" );
            SX K2  = SX::sym("k2" );
            SX kap = SX::sym("kap");
            auto lv = DefaultHyperElasticForce::getDefaultLocalVars();
            SX Ap = lv["Ap"];
            auto I = DefaultHyperElasticForce::getDefaultInvariants();
            std::array<DReal, 3> Av{_Av[0], _Av[1], _Av[2]};
            auto orth = DefaultHyperElasticForce::makeOrthotropicInvariant("I4s", Av);
            SX I1 = I["I1"];
            SX I3 = I["I3"];
            SX I4s = get<2>(orth)["I4s"];
            SX U = Ap * H * (C * (I1 + 1.0/I3 - 3) + K1/(2*K2) * SX::exp(K2*SX::sq(kap*I1 + (1 - 2*kap)*I4s - 1)) - 1);
            auto iv = DefaultHyperElasticForce::getDefaultAvailableInputVars();
            iv.emplace_back("H", H, std::make_unique<ConstDataRequier>(_h));
            iv.emplace_back("C", C, std::make_unique<ConstDataRequier>(_c));
            iv.emplace_back("k1", K1, std::make_unique<ConstDataRequier>(_k1));
            iv.emplace_back("k2", K2, std::make_unique<ConstDataRequier>(_k2));
            iv.emplace_back("kap", kap, std::make_unique<ConstDataRequier>(_kappa));

            HH::ArgsDef args(std::move(iv), std::move(lv), std::move(I));
            args = HH::concat(args, orth);
            f = DefaultHyperElasticForce("HGO0", gen_dir, std::move(args), U, regenerate);
        }
        std::unique_ptr<HyperElasticForceBase> copy() override{
            return std::make_unique<HGOModel>(_h, _c, _k1, _k2, _kappa, _Av, f.gen_dir, false);
        }
    };
}

#endif //AORTIC_VALVE_HOLZAPFELOGDENGASSER_H
