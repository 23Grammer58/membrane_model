//
// Created by alex on 01.03.2021.
//

#include "Pressure.h"

namespace World3d {
    void SimplePressureLoad::registerObj(Object3D *obj) {
        _obj = obj;
        _S = set_S(*obj);
    }

    int SimplePressureLoad::operator()(Object3D &obj) {
        if (&obj != _obj) registerObj(&obj);
        for (auto f: obj.m_mesh.faces()) {
            auto vv = vert_around(obj.m_mesh, f);
            Vector P_f = _P(obj, f) * _S[f] / 3;
            for (auto v: vv)
                if (obj.is_movable(v))
                    obj.m_F[v] += obj.withBCmask(v, P_f);
        }

        return 0;
    }

    int SimplePressureLoad::element_matrix(Object3D &obj, ForceBase::LocMatrix &matr, F_ind element) {
        matr.resize(3 * 3, 3 * 3);
        auto v = vert_around(obj.m_mesh, element);
        for (int l = 0; l < 3; ++l) {
            auto x = _P(obj, element) / 6 * (obj.m_x[v[(l + 1) % 3]] - obj.m_x[v[(l + 2) % 3]]);
            for (int n = 0; n < 3; ++n)
                matr[n][3 * l + n] = 0;
            matr[0][3 * l + 1] = x[2];
            matr[1][3 * l + 0] = -x[2];
            matr[0][3 * l + 2] = -x[1];
            matr[2][3 * l + 0] = x[1];
            matr[1][3 * l + 2] = x[0];
            matr[2][3 * l + 1] = -x[0];

        }
        for (int n = 1; n < 3; ++n)
            for (int l = 0; l < 3; ++l)
                for (int j = 0; j < 9; ++j)
                    matr[3 * n + l][j] = matr[l][j];
        for (auto i = 0; i < 3; ++i) {
            auto bc = obj.getBC(v[i]);
            for (int l = 0; l < 3; ++l){
                if (!(bc & (1 << l)))
                    matr.ApplyDirichlet(3 * i + l);
            }
        }
        return 0;
    }

    std::unique_ptr <ForceBase> SimplePressureLoad::clone() {
        return std::make_unique<SimplePressureLoad>(_P);
    }


    void SimpleLagrangePressureLoad::registerObj(Object3D *obj) {
        _obj = obj;
        _n0 = set_n0(*obj);
        _Ap = set_Ap(*obj);
    }

    int SimpleLagrangePressureLoad::operator()(Object3D &obj) {
        if (&obj != _obj) registerObj(&obj);
        for (auto f: obj.m_mesh.faces()) {
            auto vv = vert_around(obj.m_mesh, f);
            Vector P_f = _P(obj, f) * _Ap[f] * _n0[f] / 3;
            for (auto v: vv)
                if (obj.is_movable(v))
                    obj.m_F[v] += obj.withBCmask(v, P_f);
        }

        return 0;
    }

    std::unique_ptr<ForceBase> SimpleLagrangePressureLoad::clone() {
        return std::make_unique<SimpleLagrangePressureLoad>(_P);
    }

    int SimpleLagrangePressureLoad::element_matrix(Object3D &obj, ForceBase::LocMatrix &matr, F_ind element) {
        matr.resize(3 * 3, 3 * 3);
        matr.setZero();
        return 0;
    }
};
