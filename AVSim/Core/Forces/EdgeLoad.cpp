//
// Created by alex on 29.01.2021.
//

#include "EdgeLoad.h"

void World3d::BodyEdgeLoad::registerObj(Object3D *obj){
    _obj = obj;
    std::map<V_ind, Vector> moving;
    for (auto e: obj->m_mesh.edges()){
        auto d = _P(*obj, e);
        if (d[0] == 0.0 && d[1] == 0.0 && d[2] == 0.0) continue;
        auto v = vert_around(obj->m_mesh, e);
        auto len = sqrt((obj->m_x0[v[0]] - obj->m_x0[v[1]]).squared_length());
        for (int k = 0; k < 2; ++k){
            auto it = moving.find(v[k]);
            if (it == moving.end())
                moving.insert({v[k], d*len/2});
            else    
                it->second += d*len/2;  
        }
    }
    _moving.resize(0);
    _moving.reserve(moving.size());
    for (auto v: moving)
        _moving.push_back(v);
}
int World3d::BodyEdgeLoad::operator()(Object3D &obj) {
    if (&obj != _obj) registerObj(&obj);
    for (auto i: _moving) 
        obj.m_F[i.first] += obj.withBCmask(i.first, i.second);

    return 0;
}

void World3d::SimplePressureEdgeLoad::registerObj(Object3D *obj) {
    _obj = obj;
    DReal fullLen = 0;
    auto l0 = set_l0(*obj);
    for (auto v: obj->m_mesh.vertices())
        if (choose_vertex(*obj, v)) {
            DReal len = 0;
            for (auto h: CGAL::halfedges_around_target(v, obj->m_mesh)) {
                if (choose_edge_vertex(*obj, target(prev(h, obj->m_mesh), obj->m_mesh)))
                    len += l0[edge(h, obj->m_mesh)] / 2;
            }
            fullLen += len;
            _moving.emplace_back(v, len);
        }

    for (auto &i: _moving)
        i.second /= fullLen;

    if (l0.second) obj->m_mesh.remove_property_map(l0.first);
}

int World3d::SimplePressureEdgeLoad::operator()(Object3D &obj) {
    if (&obj != _obj) registerObj(&obj);
    for (auto i: _moving) {
        obj.m_F[i.first] += obj.withBCmask(i.first, i.second * _F);
    }

    return 0;
}

std::unique_ptr<World3d::ForceBase> World3d::SimplePressureEdgeLoad::clone() {
    return std::make_unique<SimplePressureEdgeLoad>(_F, choose_vertex, choose_edge_vertex);
}

int World3d::SimplePressureEdgeLoad::element_matrix(Object3D &obj, ForceBase::LocMatrix &matr, World3d::F_ind element) {
    matr.setZero();
    matr.resize(3 * 3, 3 * 3);
    return 0;
}

int World3d::SimplePressureEdgeLoad::fill_matrix(Object3D &obj, SparseMatrix::IndexFrom &m) {
    return 0;
}

void World3d::SimpleNormalEdgeLoad::registerObj(Object3D *obj) {
    _obj = obj;
    DReal fullLen = 0;
    auto l0 = set_l0(*obj);
    for (auto v: obj->m_mesh.vertices())
        if (obj->m_boundary[v] == b_lbl) {
            for (auto h: CGAL::halfedges_around_target(v, obj->m_mesh)){
                if (obj->m_boundary[target(prev(h, obj->m_mesh), obj->m_mesh)] & b_lbl)
                    if (_moving.find(edge(h, obj->m_mesh)) == _moving.end()){
                        DReal len = l0[edge(h, obj->m_mesh)];
                        _moving.insert({edge(h, obj->m_mesh), len});
                        fullLen += len;
                    }
            }
        }
    for (auto& i: _moving)
        i.second /= fullLen;

    _S = set_S(*obj);

    if (l0.second) obj->m_mesh.remove_property_map(l0.first);
}

int World3d::SimpleNormalEdgeLoad::operator()(Object3D &obj) {
    if (&obj != _obj) registerObj(&obj);
    for (auto i: _moving) {
        auto f = face_around(obj.m_mesh, i.first);
        Vector force = i.second * _F * obj.m_normal[f[0].first] / 2;
        auto vv = vert_around(obj.m_mesh, i.first);
        for (auto v: vv)
            obj.m_F[v] += force;
    }

    return 0;
}

std::unique_ptr<World3d::ForceBase> World3d::SimpleNormalEdgeLoad::clone() {
    return std::make_unique<SimpleNormalEdgeLoad>(_F, b_lbl);
}

int World3d::SimpleNormalEdgeLoad::element_matrix(Object3D &obj, ForceBase::LocMatrix &matr, World3d::F_ind element) {
    matr.setZero();
    matr.resize(3*3, 3*3);
    auto ee = edge_around(obj.m_mesh, element);
    int ind = -1;
    for (int i = 0; i < 3; ++i)
        if (_moving.find(ee[i]) != _moving.end())
            ind = i;
    if (ind < 0) return 0;
    auto vve = vert_around(obj.m_mesh, ee[ind]);
    auto v = vert_around(obj.m_mesh, element);
    int arr[2] = {0, 1};
    if (vve[0] == v[2] || vve[1] == v[2] || vve[2] == v[2]){
        arr[1] = 2;
        if (vve[0] == v[1] || vve[1] == v[1] || vve[2] == v[1]){
            arr[0] = 1;
        }
    }
    double coef = _moving[ee[ind]] * _F / (2 * 2 * sqrt(_S[element].squared_length()));
    std::array<std::array<std::array<double, 3>, 3>, 3> d2SdQ = {0};
    for (int l = 0; l < 3; ++l){
        auto x = coef*(obj.m_x[v[(l + 1)%3]] - obj.m_x[v[(l + 2)%3]]);
        for (int n = 0; n < 3; ++n)
            d2SdQ[n][l][n] = 0;
        d2SdQ[0][l][1] =  x[2]; d2SdQ[1][l][0] = -x[2];
        d2SdQ[0][l][2] = -x[1]; d2SdQ[2][l][0] =  x[1];
        d2SdQ[1][l][2] =  x[0]; d2SdQ[2][l][1] = -x[0];

    }
    auto& n = obj.m_normal[element];
    for (int i = 0; i < 3; ++i)
        for (int l = 0; l < 3; ++l)
            for (int j = 0; j < 3; ++j) {
                matr[3*arr[0] + i][3 * l + j] = d2SdQ[i][l][j];
                for (int k = 0; k < 3; ++k) {
                    matr[3*arr[0] + i][3 * l + j] += n[i] * n[k] * d2SdQ[k][l][j];
                }
                matr[3*arr[1] + i][3 * l + j] = matr[3*arr[0] + i][3 * l + j];
            }

    return 0;
}
