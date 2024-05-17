//
// Created by Liogky Alexey on 25.05.2022.
//

#include "VertexLoad.h"

World3d::VertexLoad<void>& World3d::VertexLoad<void>::insertVertices(const std::vector<std::pair<V_ind, double>> &vertices) {
    for (auto& i: vertices) m_vertices.insert(i);
    return *this;
}

World3d::VertexLoad<void>& World3d::VertexLoad<void>::insertVertices(const std::vector<V_ind> &vertices) {
    for (auto& i: vertices) m_vertices.insert({i, 1});
    return *this;
}

World3d::VertexLoad<void>& World3d::VertexLoad<void>::setForceField(std::function<bool(double *, double *)> force,
                                        std::function<bool(double *, double *)> dforce) {
    m_F = force;
    m_F_deriv = dforce;
    return *this;
}

int World3d::VertexLoad<void>::operator()(Object3D &obj) {
    std::array<double, 3> x, F;
    for (auto i: m_vertices) {
        auto v = i.first;
        auto w = i.second;
        for (int k = 0; k < 3; ++k) x[k] = obj.m_x[v][k];
        m_F(x.data(), F.data());
        obj.m_F[v] += obj.withBCmask(v, w * Vector(F[0], F[1], F[2]));
    }
    return 0;
}

int World3d::VertexLoad<void>::fill_matrix(Object3D &obj, SparseMatrix::IndexFrom &m) {
    std::array<double, 3> x;
    std::array<double, 9> J;
    for (auto i: m_vertices){
        auto v = i.first;
        auto w = i.second;
        m_F_deriv(x.data(), J.data());
        std::array<bool, 3> bc = { obj.is_movable(v, 0), obj.is_movable(v, 1), obj.is_movable(v, 2) };
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                if (bc[i] && bc[j]) m(v*3 + i, v*3 + j) += w * J[3*i + j];
    }
    return 0;
}
