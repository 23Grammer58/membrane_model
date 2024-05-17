//
// Created by alex on 29.01.2021.
//

#ifndef AORTIC_VALVE_PRESSURE_H
#define AORTIC_VALVE_PRESSURE_H
#include "../Object3D.h"

namespace World3d {
    static double operator""_mmHg(long double P){
        return P * 133.3223684;
    }

    class SimplePressureLoad: public ForceBase {
        using P_func = std::function<DReal(Object3D&, F_ind)>;
        P_func _P;
        UpdatableProperty<F_ind, Vector> _S;
        Object3D *_obj = nullptr;
    public:
        explicit SimplePressureLoad(P_func pf): _P{pf} { type = "SimplePressureLoad"; }
        explicit SimplePressureLoad(DReal pf): _P{[pf](Object3D&, F_ind){ return pf; }} { type = "SimplePressureLoad"; }
        void registerObj(Object3D *obj);
        void setPressure(DReal P){ _P = [P](Object3D&, F_ind){ return P; }; }
        void setPressure(P_func P){ _P = std::move(P); }
        int operator()(Object3D &obj) override;
        std::unique_ptr<ForceBase> clone() override;
        int element_matrix(Object3D &obj, LocMatrix &matr, F_ind element) override;
    };

    template<typename RealType>
    class SimpleWrongPressureLoad : public ForceBase {
        RealType _P;
        ConstProperty<F_ind, Vector> _n0;
        ConstProperty<F_ind, DReal> _Ap;
        Object3D *_obj = nullptr;
    public:
        explicit SimpleWrongPressureLoad(RealType P) : _P{P} { type = "SimpleWrongPressureLoad"; }

        void registerObj(Object3D *obj);
        int operator()(Object3D &obj) override;
        std::unique_ptr<ForceBase> clone() override;
        int element_matrix(Object3D &obj, LocMatrix &matr, F_ind element) override;
    };

    class SimpleLagrangePressureLoad : public ForceBase {
        using P_func = std::function<DReal(Object3D&, F_ind)>;
        P_func _P;
        ConstProperty<F_ind, Vector> _n0;
        ConstProperty<F_ind, DReal> _Ap;
        Object3D *_obj = nullptr;
    public:
        explicit SimpleLagrangePressureLoad(P_func pf): _P{pf} { type = "SimpleLagrangePressureLoad"; }
        explicit SimpleLagrangePressureLoad(DReal pf) : _P{[pf](Object3D&, F_ind){ return pf; }} { type = "SimpleLagrangePressureLoad"; }

        void registerObj(Object3D *obj);
        int operator()(Object3D &obj) override;
        std::unique_ptr<ForceBase> clone() override;
        int element_matrix(Object3D &obj, LocMatrix &matr, F_ind element) override;
    };
}

#include "Pressure.inl"


#endif //AORTIC_VALVE_PRESSURE_H
