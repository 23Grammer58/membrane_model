//
// Created by alex on 23.01.2022.
//

#ifndef AORTIC_VALVE_LINEARELASTIC_H
#define AORTIC_VALVE_LINEARELASTIC_H

#include "../Object3D.h"
#include "HyperElastic.h"

namespace World3d {
    class LinearElasticModel: public ForceBase {
        DReal _mu, _H, _lambda;
        Object3D *_obj = nullptr;
        ConstProperty<F_ind, std::array<Vector, 3>> _Dv;
        ConstProperty<F_ind, DReal> _Ap;
    public:
        void registerObj(Object3D *obj);
        LinearElasticModel(double h, double lambda, double mu): _H{h}, _lambda{lambda}, _mu{mu} {
            type = "LinearElastic";
        }
        std::unique_ptr<ForceBase> clone() override{
            return std::make_unique<LinearElasticModel>(_H, _lambda, _mu);
        }
        int operator()(Object3D& obj) override;
        int fill_matrix(Object3D& obj, SparseMatrix::IndexFrom& m) override;
    };

    class LinearElasticModel1: public HyperElasticForceBase {
        DReal _mu, _H, _lambda;
        Object3D *_obj = nullptr;
        ConstProperty<F_ind, std::array<Vector, 3>> _Dv;
        ConstProperty<F_ind, DReal> _Ap;
    public:
        LinearElasticModel1(double h, double lambda, double mu, std::string gen_dir, bool regenerate = true);
        std::unique_ptr<HyperElasticForceBase> copy() override{
            return std::make_unique<LinearElasticModel1>(_H, _lambda, _mu, f.gen_dir, false);
        }
    };
}

#endif //AORTIC_VALVE_LINEARELASTIC_H
