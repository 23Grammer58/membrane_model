//
// Created by alex on 29.01.2021.
//

#ifndef AORTIC_VALVE_PRESSURE_INL
#define AORTIC_VALVE_PRESSURE_INL

namespace World3d {
    template<typename RealType>
    void SimpleWrongPressureLoad<RealType>::registerObj(Object3D *obj) {
        _obj = obj;
        _n0 = set_n0(*obj);
        _Ap = set_Ap(*obj);
    }

    template<typename RealType>
    int SimpleWrongPressureLoad<RealType>::operator()(Object3D &obj) {
        if (&obj != _obj) registerObj(&obj);
        for (auto f: obj.m_mesh.faces()) {
            auto vv = vert_around(obj.m_mesh, f);
            Vector P_f = _P * _Ap[f] * _n0[f] / 3;
            for (auto v: vv)
                if (obj.is_movable(v))
                    obj.m_F[v] += obj.withBCmask(v, P_f);
        }

        return 0;
    }

    template<typename RealType>
    std::unique_ptr<ForceBase> SimpleWrongPressureLoad<RealType>::clone() {
        return std::make_unique<SimpleWrongPressureLoad>(_P);
    }

    template<typename RealType>
    int SimpleWrongPressureLoad<RealType>::element_matrix(Object3D &obj, ForceBase::LocMatrix &matr, F_ind element) {
        matr.resize(3 * 3, 3 * 3);
        return 0;
    }
};


#endif //AORTIC_VALVE_PRESSURE_INL
