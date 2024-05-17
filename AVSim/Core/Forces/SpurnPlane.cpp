//
// Created by alex on 15.04.2021.
//

#include "SpurnPlane.h"

void World3d::SpurnPlane::registerObj(Object3D *obj) {
    _obj = obj;
    _S = set_S(*obj);
}

int World3d::SpurnPlane::operator()(Object3D &obj) {
    if (&obj != _obj) registerObj(&obj);
    for (auto f: obj.m_mesh.faces()) {
        auto vv = vert_around(obj.m_mesh, f);
        auto c = ((obj.m_x[vv[0]] - CGAL::ORIGIN) + (obj.m_x[vv[1]] - CGAL::ORIGIN) + (obj.m_x[vv[2]] - CGAL::ORIGIN))/3;
        double Aq = sqrt(_S[f].squared_length());
        Vector P_f = m_P * Aq * m_pl.orthogonal_vector() / 3;
        for (auto v: vv)
            if (obj.is_movable(v)) {
                auto dist  =  m_pl.orthogonal_vector() * (obj.m_x[v] - m_pl.point());
                double val = m_dist_func(dist);
                obj.m_F[v] += obj.withBCmask(v, val * P_f);
            }
    }

    return 0;
}

int World3d::SpurnPlane::element_matrix(Object3D &obj, ForceBase::LocMatrix &matr, World3d::F_ind element) {
    matr.resize(3 * 3, 3 * 3);
    auto v = vert_around(obj.m_mesh, element);
    auto S = _S[element];
    auto Aq = sqrt(S.squared_length());
    auto s = S / Aq;
    auto orth = m_pl.orthogonal_vector();
    for (int n = 0; n < 3; ++n) {
        auto dist  =  orth * (obj.m_x[v[n]] - m_pl.point());
        double P = m_P * m_dist_func(dist);
        for (int l = 0; l < 3; ++l) {
            auto x = P / 6 * CGAL::cross_product(s, obj.m_x[v[(l + 2) % 3]] - obj.m_x[v[(l + 1) % 3]]);
            for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                matr[3*n + i][3*l+j] = orth[i] * x[j];
        }

        double dPddist = m_P * m_dist_deriv(dist);
        for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            matr[3 * n + i][3 * n + j] += orth[i] * orth[j] *dPddist * Aq / 3;
    }

    for (auto i = 0; i < 3; ++i) {
        auto bc = obj.getBC(v[i]);
        for (int l = 0; l < 3; ++l){
            if (!(bc & (1 << l)))
                matr.ApplyDirichlet(3 * i + l);
        }
    }

    return 0;
}

void World3d::SpurnPlane::set_dist_funcs(double dist_sparse) {
    static const double scl = 0.0;
//    m_dist_func = [width = dist_sparse](double dist) -> double { return (dist > 0) ? exp(-2*dist / width) : (1 - 2*dist / width*(1 - scl*dist / width)); };
//    m_dist_deriv = [width = dist_sparse](double dist) -> double { return (dist > 0) ? -2 * exp(-2*dist / width) / width : -2 / width * (1 - scl*2*dist / width); };
    m_dist_func = [width = dist_sparse](double dist) -> double { return (2*dist < width) ? (1 - 2*dist / width) : 0; };
    m_dist_deriv = [width = dist_sparse](double dist) -> double { return (2*dist < width) ? -2 / width : 0; };
}
