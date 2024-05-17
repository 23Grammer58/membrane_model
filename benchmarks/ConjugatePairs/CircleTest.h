//
// Created by Liogky Alexey on 01.10.2022.
//

#ifndef AVSYM_CONJUGATEPAIRS_CIRCLETEST_H
#define AVSYM_CONJUGATEPAIRS_CIRCLETEST_H

#include "../BenchCommon.h"
#include "AVSim/Core/NSWorldWrapper.h"
#include "AVSim/Solvers/NonLinearSolverKinsol.h"
#include "AVSim/Core/NonLinearSolverCustom.h"
#include <Eigen/Dense>
#include "AVSim/Core/MeshGen/Helper.h"

class ZhouModel: public HyperElasticForceBase{
public:
    DReal _mu;
    DReal _h;
    ZhouModel(double mu, double h, std::string gen_dir, bool regenerate = true):
            _mu{mu}, _h{h}
    {
        type = "ZhouModel";
        SX Mu = SX::sym("mu");
        SX H  = SX::sym("H" );
        auto lv = DefaultHyperElasticForce::getDefaultLocalVars();
        SX Ap = lv["Ap"];
        auto I = DefaultHyperElasticForce::getDefaultInvariants();
        SX I1 = I["I1"];
        SX I3 = I["I3"];
        SX U = Ap * H * Mu / 2 * (I1 - SX::log(I3) - 2 + SX::sq(I1 - 2)/2);
        auto iv = DefaultHyperElasticForce::getDefaultAvailableInputVars();
        iv.emplace_back("mu", Mu, std::make_unique<ConstDataRequier>(_mu));
        iv.emplace_back("H", H, std::make_unique<ConstDataRequier>(_h));
        f = DefaultHyperElasticForce("Zhou", gen_dir, std::move(iv), lv, I, U, regenerate);
    }
    std::unique_ptr<HyperElasticForceBase> copy() override{
        return std::make_unique<ZhouModel>(_mu, _h, f.gen_dir, false);
    }
};

class OgdenBPModel: public HyperElasticForceBase{
public:
    DReal _mu;
    DReal _h;
    DReal _alpha;
    OgdenBPModel(double mu, double h, double alpha, std::string gen_dir, bool regenerate = true):
            _mu{mu}, _h{h}, _alpha{alpha}
    {
        type = "OgdenBPModel";
        SX Mu = SX::sym("mu");
        SX H  = SX::sym("H" );
        SX A  = SX::sym("al");
        auto lv = DefaultHyperElasticForce::getDefaultLocalVars();
        SX Ap = lv["Ap"];
        auto I = DefaultHyperElasticForce::getDefaultInvariants();
        SX I1 = I["I1"];
        SX I3 = I["I3"];
        SX D = SX::sq(I1) - 4*I3;
            D = SX::if_else_zero(D > 0, D);
        SX sqrtD = SX::sqrt(D);
        SX sqrl1 = (I1 + sqrtD)/2, sqrl2 = (I1 - sqrtD)/2;
            sqrl2 = SX::if_else_zero(sqrl2 > 0, sqrl2);
        SX W = 2 / SX::sq(A) * (SX::pow(sqrl1, A/2) + SX::pow(sqrl2, A/2) + SX::pow(sqrl1*sqrl2, -A/2) - 3);
        SX U = Ap * H * Mu * W;
        auto iv = DefaultHyperElasticForce::getDefaultAvailableInputVars();
        iv.emplace_back("mu", Mu, std::make_unique<ConstDataRequier>(_mu));
        iv.emplace_back("H", H, std::make_unique<ConstDataRequier>(_h));
        iv.emplace_back("al", A, std::make_unique<ConstDataRequier>(_alpha));
        f = DefaultHyperElasticForce("OgdenBP", gen_dir, std::move(iv), lv, I, U, regenerate);
    }
    std::unique_ptr<HyperElasticForceBase> copy() override{
        return std::make_unique<OgdenBPModel>(_mu, _h, _alpha, f.gen_dir, false);
    }
};

// struct ConjugateValues{
//     std::string x0_tag_name = "v:point",
//                 x_tag_name = "v:x";
//     Object3D& obj;


// };

#endif //AVSYM_CONJUGATEPAIRS_CIRCLETEST_H