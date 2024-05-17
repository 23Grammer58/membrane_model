//
// Created by alex on 29.01.2021.
//

#ifndef AORTIC_VALVE_MAYYIN_H
#define AORTIC_VALVE_MAYYIN_H

#include "HyperElastic.h"
namespace World3d{
    class MayYinModel: public HyperElasticForceBase {
    public:
        DReal _h, _c0, _c1, _c2;
        std::array<DReal, 3> _Av;
        MayYinModel(double h, double c0, double c1, double c2, std::array<DReal, 3> _Av, std::string gen_dir, bool regenerate = true):
                _h{h}, _c0{c0}, _c1{c1}, _c2{c2}, _Av{_Av}
        {
            type = "MayYinModel";
            SX H  = SX::sym("H" );
            SX C0 = SX::sym("C0");
            SX C1  = SX::sym("C1" );
            SX C2  = SX::sym("C2" );
            auto lv = DefaultHyperElasticForce::getDefaultLocalVars();
            SX Ap = lv["Ap"];
            auto I = DefaultHyperElasticForce::getDefaultInvariants();
            std::array<DReal, 3> Av{_Av[0], _Av[1], _Av[2]};
            auto orth = DefaultHyperElasticForce::makeOrthotropicInvariant("I4s", Av);
            SX I1 = I["I1"];
            SX I3 = I["I3"];
            SX I4s = get<2>(orth)["I4s"];
            SX U = Ap * H * C0 * (SX::exp(C1 * SX::sq(I1 + 1.0/I3 - 3) + C2 * SX::sq(I4s - 1)) - 1);
            auto iv = DefaultHyperElasticForce::getDefaultAvailableInputVars();
            iv.emplace_back("H", H, std::make_unique<ConstDataRequier>(_h));
            iv.emplace_back("C0", C0, std::make_unique<ConstDataRequier>(_c0));
            iv.emplace_back("C1", C1, std::make_unique<ConstDataRequier>(_c1));
            iv.emplace_back("C2", C2, std::make_unique<ConstDataRequier>(_c2));

            HH::ArgsDef args(std::move(iv), std::move(lv), std::move(I));
            args = HH::concat(args, orth);
            f = DefaultHyperElasticForce("MayYin", gen_dir, std::move(args), U, regenerate);
        }
        std::unique_ptr<HyperElasticForceBase> copy() override{
            return std::make_unique<MayYinModel>(_h, _c0, _c1, _c2, _Av, f.gen_dir, false);
        }
    };
}

#endif //AORTIC_VALVE_MAYYIN_H
