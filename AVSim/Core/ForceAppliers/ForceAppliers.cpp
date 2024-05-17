//
// Created by alex on 30.06.2020.
//

#include "ForceAppliers.h"

auto World3d::DynamicForceApplier::registerObj(World3d::Object3D *obj) {
    auto t = obj->m_mesh.property_map<V_ind, Vector>("v:velocity");
    if (!t.second)
        throw std::runtime_error("For using DynamicForceApplier Object3D must have \"v:velocity\" property");
    Velocity v = t.first;
    auto t1 = obj->m_mesh.property_map<V_ind, double>("v:mass");
    if (!t1.second)
        throw std::runtime_error("For using DynamicForceApplier Object3D must have \"v:mass\" property");
    Mass m = t1.first;
    return data.insert({obj, {v, m}});
}

int World3d::DynamicForceApplier::operator()(World3d::Object3D &obj) {
    auto it = data.find(&obj);
    if (it == data.end()){
        auto q = registerObj(&obj);
        it = q.first;
    }
    MeshData& md = it->second;
    Velocity v = get<0>(md);
    Mass m = get<1>(md);
    auto& f = obj.m_F;
    for (auto i: obj.m_mesh.vertices()){
        if (obj._is_movable(obj, i)){
            obj.m_next_x[i] = obj.m_x[i] + v[i] * _dt + f[i] * _dt * _dt / (2 * m[i]);
        }
    }
    return 0;
}

int World3d::StaticForceApplier::operator()(World3d::Object3D &obj) {
    for (auto i: obj.m_mesh.vertices()){
        if (obj._is_movable(obj, i)){
            obj.m_next_x[i] = obj.m_x[i] + _delta * obj.withBCmask(i, obj.m_F[i]);
        }
    }
    return 0;
}

int World3d::StaticForceApplierP::operator()(World3d::Object3D &obj) {
    double delta = _delta;
    for (auto& dmp: _dampers) delta = dmp(obj, delta);
    for (auto i: obj.m_mesh.vertices()){
        if (obj._is_movable(obj, i)){
            obj.m_next_x[i] = obj.m_x[i] + delta * obj.withBCmask(i, obj.m_F[i]);
        }
    }
    return 0;
}

double World3d::Damper::operator()(Object3D &obj, double cur_delta) {
    double delta = factor * cur_delta;
    double allow_dx = std::min(max_dx, init_incr * ((init_dx > 0) ? init_dx: max_dx));
    double d2 = delta * delta, r2 = allow_dx * allow_dx;
    bool delta_ok = true;
    bool delta_increase = (factor <= scl * 1);
    if (init_dx > 0) {
        double dl2 = 0;
        for (auto i: obj.m_mesh.vertices()) {
            if (obj._is_movable(obj, i)) {
                auto mv2 = d2 * obj.withBCmask(i, obj.m_F[i]).squared_length();
                if (mv2 > dl2) dl2 = mv2;
                if (mv2 > r2) {
                    d2 /= mv2 / r2;
                    delta_ok = false;
                    delta_increase = false;
                }
                if (delta_increase && mv2 > r2 * scl * scl) delta_increase = false;
            }
        }
        reinit_val = std::max(d2, reinit_val);
    } else {
        for (auto i: obj.m_mesh.vertices()) {
            if (obj._is_movable(obj, i)) {
                auto mv2 = d2 * obj.withBCmask(i, obj.m_F[i]).squared_length();
                init_dx = std::max(mv2, init_dx);
                if (mv2 > r2) {
                    d2 /= mv2 / r2;
                    delta_ok = false;
                    delta_increase = false;
                }
                if (delta_increase && mv2 > r2 * scl * scl) delta_increase = false;
            }
        }
        init_dx = (init_dx > 0) ? sqrt(init_dx): max_dx;
    }
    if (!delta_ok) {
        bad_its++;
        increase_its = 0;
        ok_its = 0;
        delta = sqrt(d2);
    } else {
        ok_its++;
    }
    if (delta_increase){ increase_its++; }
    if (no_damp_it > 0 && increase_its >= no_damp_it){
        factor /= scl;
        std::cout << "Increase factor to " << factor << std::endl;
        delta = factor * cur_delta;
        ok_its = 0; bad_its = 0; increase_its = 0;
    }
    if (bad_its > damp_it){
        if (factor > min_factor) {
            factor *= scl;
            std::cout << "Decrease factor to " << factor << std::endl;
        }
        ok_its = 0; bad_its = 0; increase_its = 0; reinit_val = 2*allow_dx*allow_dx;
    }
    if (re_init_nit > 0 && m_it > 0 && m_it % re_init_nit == 0){
        double max_val = sqrt(reinit_val);
        if (max_val < reinit_min_factor * allow_dx) {
            init_dx = 1.5*max_val;
            std::cout << "Decrease init_dx bound to " << init_dx << std::endl;
        }
        reinit_val = 0;
    }
    m_it++;


    return delta;
}
