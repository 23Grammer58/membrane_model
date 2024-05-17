//
// Created by alex on 29.01.2021.
//

#ifndef AORTIC_VALVE_EDGELOAD_H
#define AORTIC_VALVE_EDGELOAD_H

#include "../Object3D.h"

namespace World3d{
    class BodyEdgeLoad: public ForceBase{
        using P_func = std::function<Vector(Object3D&, E_ind)>;
        P_func _P;
        Object3D *_obj = nullptr;

        std::vector<std::pair<V_ind, Vector>> _moving;
    public:
        explicit BodyEdgeLoad(Vector F): _P([F](Object3D&, E_ind){ return F; }) {type = "BodyEdgeLoad";}
        explicit BodyEdgeLoad(P_func pf): _P{pf} { type = "BodyEdgeLoad"; } 
        
        void registerObj(Object3D *obj);
        int operator()(Object3D &obj) override;
        std::unique_ptr<ForceBase> clone() override { return std::make_unique<BodyEdgeLoad>(_P); }
        int element_matrix(Object3D &obj, LocMatrix &matr, F_ind element) override { matr.setZero(); matr.resize(3 * 3, 3 * 3); return 0; }
        int fill_matrix(Object3D &obj, SparseMatrix::IndexFrom &m) override { return 0; }
    };

    class SimplePressureEdgeLoad : public ForceBase {
        Vector _F;
        Object3D *_obj = nullptr;
        std::vector<std::pair<V_ind, DReal>> _moving;
        std::function<bool(const Object3D&, V_ind)> choose_vertex;
        std::function<bool(const Object3D&, V_ind)> choose_edge_vertex;
        //int b_lbl;
    public:
        explicit SimplePressureEdgeLoad(Vector F, int b_lbl) : _F{F}, choose_vertex{[b_lbl](const Object3D& obj, V_ind v) { return (obj.m_boundary[v] == b_lbl); }},
                                                               choose_edge_vertex{[b_lbl](const Object3D& obj, V_ind v) { return (obj.m_boundary[v] & b_lbl) != 0; }}
         { type = "SimplePressureEdgeLoad"; }
        explicit SimplePressureEdgeLoad(Vector F, std::function<bool(const Object3D&, V_ind)> choose_vertex, std::function<bool(const Object3D&, V_ind)> choose_edge_vertex):
        _F{F}, choose_vertex{choose_vertex}, choose_edge_vertex{choose_edge_vertex} {type = "SimplePressureEdgeLoad";}

        void registerObj(Object3D *obj);

        int operator()(Object3D &obj) override;
        std::unique_ptr<ForceBase> clone() override;
        int element_matrix(Object3D &obj, LocMatrix &matr, F_ind element) override;
        int fill_matrix(Object3D &obj, SparseMatrix::IndexFrom &m) override;
    };

    class SimpleNormalEdgeLoad: public ForceBase{
        DReal _F;
        Object3D* _obj = nullptr;
        std::map<E_ind, DReal> _moving;
        UpdatableProperty<F_ind, Vector> _S;
        int b_lbl;
    public:
        explicit SimpleNormalEdgeLoad(DReal F, int b_lbl): _F{F}, b_lbl{b_lbl} { type = "SimpleNormalEdgeLoad"; }

        void registerObj(Object3D* obj);
        int operator() (Object3D& obj) override;
        std::unique_ptr<ForceBase> clone() override;
        int element_matrix(Object3D& obj, LocMatrix& matr, F_ind element) override;
    };
}


#endif //AORTIC_VALVE_EDGELOAD_H
