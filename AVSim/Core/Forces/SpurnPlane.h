//
// Created by alex on 15.04.2021.
//

#ifndef AORTIC_VALVE_SPURNPLANE_H
#define AORTIC_VALVE_SPURNPLANE_H
#include "../Object3D.h"

namespace World3d {
    class SpurnPlane : public ForceBase {
        UpdatableProperty<F_ind, Vector> _S;
        Object3D *_obj = nullptr;
        std::function<double(double)> m_dist_func;
        std::function<double(double)> m_dist_deriv;
    public:
        using Plane_3 = CGAL::Simple_cartesian<double>::Plane_3;
        double m_P = 0;
        Plane_3 m_pl;

        SpurnPlane(const SpurnPlane&) = default;
        explicit SpurnPlane(double P, double dist_sparse, Plane_3 plane): m_P{ P }, m_pl{plane} {
            set_dist_funcs(dist_sparse);
            type = "SpurnPlaneForce";
        }
        void set_dist_funcs(double dist_sparse);
        void set_dist_funcs(std::function<double(double d)> dist_func, std::function<double(double d)> dist_deriv){
            m_dist_func = std::move(dist_func), m_dist_deriv = std::move(dist_deriv);
        }
        void registerObj(Object3D *obj);
        int operator()(Object3D &obj) override;
        std::unique_ptr<ForceBase> clone() override{ return std::make_unique<SpurnPlane>(*this); }
        int element_matrix(Object3D &obj, LocMatrix &matr, F_ind element) override;
    };

}

#endif //AORTIC_VALVE_SPURNPLANE_H
