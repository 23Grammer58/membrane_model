//
// Created by alex on 29.01.2021.
//

#include "MSM.h"

void World3d::SimpleMassSpringModel::registerObj(Object3D *obj) {
    _obj = obj;
    _l0 = set_l0(*obj);
    _l = set_l(*obj);
    _k = obj->m_mesh.add_property_map<E_ind, SavedData>("e:SimpleMassSpringModel");
    auto &m = obj->m_mesh;
    for (auto e: m.edges()) {
        std::array<V_ind, 4> v;
        v[0] = m.vertex(e, 0);
        v[1] = m.vertex(e, 1);
        auto vv = vert_opposite(m, e);
        DReal A0 = 0;
        for (int i = 0; i < vv.second; ++i) {
            v[2 + i] = vv.first[i];

            A0 += 0.5 * sqrt(CGAL::cross_product(Vector(obj->m_x0[v[0]], obj->m_x0[v[1]]),
                                                 Vector(obj->m_x0[v[0]],
                                                        obj->m_x0[v[2 + i]])).squared_length());
        }
        _k[e] = _E * (_h * A0) / (_l0[e] * _l0[e]);
    }
}

int World3d::SimpleMassSpringModel::operator()(Object3D &obj) {
    if (&obj != _obj) registerObj(&obj);
    auto &m = obj.m_mesh;
    for (auto e: m.edges()) {
        auto v = vert_around(m, e);
        if (obj.m_boundary[v[0]] && obj.m_boundary[v[1]])
            continue;
        Vector dir = (obj.m_x[v[1]] - obj.m_x[v[0]]) / _l[e];
        DReal dl = _l[e] - _l0[e];
        Vector force = _k[e] * dl * dir;
        if (!obj.m_boundary[v[0]])
            obj.m_F[v[0]] += force;
        if (!obj.m_boundary[v[1]])
            obj.m_F[v[1]] += -force;
    }
    return 0;
}

int World3d::SimpleMassSpringModel::element_matrix(Object3D &obj, ForceBase::LocMatrix &matr, E_ind edge) {
    matr.resize(2 * 3, 2 * 3);
    auto &m = obj.m_mesh;
    auto v = vert_around(m, edge);
    Vector dir = (obj.m_x[v[1]] - obj.m_x[v[0]]) / _l[edge];
    DReal l0 = _l0[edge] / _l[edge];
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            matr[3 + i][3 + j] = matr[i][j] = -_k[edge] * ((i == j) * (1 - l0) + l0 * dir[i] * dir[j]);
            matr[3 + i][j] = matr[i][3 + j] = matr[i][j];
        }
    return 0;
}

std::unique_ptr<World3d::ForceBase> World3d::SimpleMassSpringModel::clone() {
    return std::make_unique<SimpleMassSpringModel>(_E, _h);
}

int World3d::SimpleMassSpringModel::element_matrix(Object3D &obj, ForceBase::LocMatrix &matr, World3d::F_ind element) {
    matr.setZero();
    matr.resize(3 * 3, 3 * 3);
    auto e = edge_around(obj.m_mesh, element);
    LocMatrix locMatrix;
    //TODO: возможно нарушается индексация (как edge согласованы с vertex?)
    //TODO: сделать ApplyDirichlet
    for (int i = 0; i < 3; ++i) {
        double sc = 0.5;
        if (!vert_opposite(obj.m_mesh, e[i], 1).second) sc = 1;
        element_matrix(obj, locMatrix, e[i]);
        int l1 = i, l2 = (i + 1) % 3;
        LocMatrix::DoubleIndex di(3, 3), dl(2, 3);
        for (int j1 = 0; j1 < 3; ++j1)
            for (int j2 = 0; j2 < 3; ++j2) {
                matr[di(l1, j1)][di(l1, j2)] = locMatrix[dl(0, j1)][dl(0, j2)];
                matr[di(l1, j1)][di(l2, j2)] = locMatrix[dl(0, j1)][dl(1, j2)];
                matr[di(l2, j1)][di(l1, j2)] = locMatrix[dl(1, j1)][dl(0, j2)];
                matr[di(l2, j1)][di(l2, j2)] = locMatrix[dl(1, j1)][dl(1, j2)];
            }
    }
    return 0;
}
