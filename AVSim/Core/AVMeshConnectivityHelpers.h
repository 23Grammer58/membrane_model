//
// Created by alex on 30.06.2020.
//

#ifndef AORTIC_VALVE_OBJECT3DINLINEHELPERS_H
#define AORTIC_VALVE_OBJECT3DINLINEHELPERS_H
#include <chrono>
#include "AVMesh.h"

namespace World3d {
    static inline void face_around(const Mesh &mesh, V_ind v, std::vector<F_ind>& res){
        auto it = faces_around_target(mesh.halfedge(v), mesh);
        res.resize(0);
        res.reserve(it.size());
        for (auto f = it.begin(); f != it.end(); ++f)
            if (f->is_valid())
                res.push_back(*f);
    }

    static inline std::vector<F_ind> face_around(const Mesh &mesh, V_ind v){
        auto it = faces_around_target(mesh.halfedge(v), mesh);
        std::vector<F_ind> res;
        res.reserve(it.size());
        for (auto f = it.begin(); f != it.end(); ++f)
            if (f->is_valid())
                res.push_back(*f);
        return res;
    }

    static inline std::array<V_ind, 3> vert_around(const Mesh &mesh, F_ind f) {
        std::array<V_ind, 3> res;
        auto it = vertices_around_face(mesh.halfedge(f), mesh).first;
        for (int i = 0; i < 3; ++i)
            res[i] = *(it++);
        return res;
    }

    static inline void vert_around(const Mesh &mesh, F_ind f, std::vector<V_ind>& res) {
        auto it = vertices_around_face(mesh.halfedge(f), mesh).first;
        for (int i = 0; i < 3; ++i)
            res[i] = *(it++);
    }

    static inline std::array<E_ind, 3> edge_around(const Mesh &mesh, F_ind f) {
        std::array<E_ind, 3> res;
        auto it = edges_around_face(mesh.halfedge(f), mesh).first;
        for (int i = 0; i < 3; ++i)
            res[i] = *(it++);
        return res;
    }

    static inline void edge_around(const Mesh &mesh, V_ind v, std::vector<E_ind>& out){
        auto heds1 = CGAL::halfedges_around_source(v, mesh);
        auto heds2 = CGAL::halfedges_around_target(v, mesh);
        int it = 0;
        out.resize(heds1.size() + heds2.size());
        for (auto h = heds1.begin(); h != heds1.end(); ++h) out[it++] = E_ind(*h);
        std::sort(out.data(), out.data() + it);
        int it1 = it;
        for (auto h = heds2.begin(); h != heds2.end(); ++h){
            E_ind e(*h);
            bool exists = std::binary_search(out.data(), out.data() + it1, e);
            if (!exists) out[it++] = e;
        }
        out.resize(it);
    }

    static inline std::vector<E_ind> edge_around(const Mesh &mesh, V_ind v){
        std::vector<E_ind> res;
        edge_around(mesh, v, res);
        return res;
    }

    static inline std::array<V_ind, 2> vert_around(const Mesh &mesh, E_ind e) {
        std::array<V_ind, 2> res;
        res[0] = mesh.vertex(e, 0);
        res[1] = mesh.vertex(e, 1);
        return res;
    }

    static inline E_ind edge_around(const Mesh &mesh, V_ind v_st, V_ind v_end, F_ind on_face){
        E_ind res(-1);
        auto ee = edge_around(mesh, on_face);
        for (int i = 0; i < ee.size(); ++i){
            E_ind e = ee[i];
            auto ve = vert_around(mesh, e);
            if ((ve[0] == v_st && ve[1] == v_end) || (ve[1] == v_st && ve[0] == v_end)){
                res = e;
                break;
            }
        }
        return res;
    }

    static inline std::array<std::pair<F_ind, bool>, 2> face_around(const Mesh &mesh, E_ind e){
        std::array<std::pair<F_ind, bool>, 2> res;
        auto h = halfedge(e, mesh);
        if (h.is_valid()) {
            auto f = face(h, mesh);
            if (f.is_valid())
                res[0] = std::pair<F_ind, bool>{f, true};
        }
        auto oh = opposite(h, mesh);
        if (oh.is_valid()){
            auto f = face(oh, mesh);
            if (f.is_valid())
                res[1] = std::pair<F_ind, bool>{f, true};
        }
        if (res[1].second && !res[0].second) swap(res[0], res[1]);
        return res;
    }

    static inline std::pair<V_ind, bool> vert_opposite(const Mesh &mesh, E_ind e, int i) {
        assert(i == 0 || i == 1);
        V_ind v[2];
        std::pair<V_ind, bool> res;
        res.second = false;
        v[0] = mesh.vertex(e, 0);
        v[1] = mesh.vertex(e, 1);
        {
            auto h = mesh.halfedge(e);
            auto f = CGAL::face(h, mesh);
            auto opf = CGAL::face(CGAL::opposite(h, mesh), mesh);
            std::array<F_ind, 2> ff = {f, opf};
            if (!ff[i].is_valid()){
                return res;
            } else {
                auto vv = vert_around(mesh, ff[i]);
                for (int i = 0; i < 3; ++i){
                    if (vv[i] != v[0] && vv[i] != v[1]){
                        res.first = vv[i];
                        res.second = true;
                        break;
                    }
                }
            }
        }
//        auto it = vertices_around_face(mesh.halfedge(e, i), mesh);
//        if (it.size() >= 3) {
//            auto i = it.first;
//            for (unsigned j = 0; j < 3; ++j)
//                if (*i != v[0] && *i != v[1]) {
//                    res.first = *i;
//                    res.second = true;
//                    break;
//                } else ++i;
//        }

        return res;
    }

    static inline std::pair<std::array<V_ind, 2>, int> vert_opposite(const Mesh &mesh, E_ind e) {
        V_ind v[2];
        std::pair<std::array<V_ind, 2>, int> res;
        res.second = 0;
        v[0] = mesh.vertex(e, 0);
        v[1] = mesh.vertex(e, 1);
        {
            auto h = mesh.halfedge(e);
            auto f = CGAL::face(h, mesh);
            auto opf = CGAL::face(CGAL::opposite(h, mesh), mesh);
            std::array<F_ind, 2> ff = {f, opf};
            res.second = 2;
            if (!f.is_valid()){
                res.second--;
                std::swap(ff[0], ff[1]);
            }
            if (!opf.is_valid()){
                res.second--;
            }
            for (int l = 0; l < res.second; ++l){
                auto vv = vert_around(mesh, ff[l]);
                for (int i = 0; i < 3; ++i){
                    if (vv[i] != v[0] && vv[i] != v[1]){
                        res.first[l] = vv[i];
                        break;
                    }
                }
            }
        }

//        for (int l = 0; l < 2; ++l) {
//            auto it = vertices_around_face(mesh.halfedge(e, l), mesh);
//            if (it.size() >= 3) {
//                auto i = it.first;
//                for (unsigned j = 0; j < 3; ++j)
//                    if (*i != v[0] && *i != v[1]) {
//                        res.first[res.second++] = *i;
//                        break;
//                    } else ++i;
//            }
//        }

        return res;
    }

    static inline std::array<std::pair<V_ind, bool>, 3> vert_opposite(const Mesh &mesh, F_ind f){
        std::array<std::pair<V_ind, bool>, 3>  res;
        auto p = halfedges_around_face(mesh.halfedge(f), mesh);
        auto it = p.begin();
        V_ind existing;
        bool bexist = false;
        constexpr int shft = 1;
        for (int i = 0; i < 3; ++i){
            res[(i+shft)%3].second = false;
            V_ind v0 = CGAL::target(*(it), mesh);
            auto op = CGAL::opposite(*it, mesh);
            V_ind v1 = CGAL::target(op, mesh);
            F_ind opf = CGAL::face(op, mesh);
            if (!opf.is_valid()) {
                ++it;
                continue;
            }
            auto opit = vertices_around_face(mesh.halfedge(opf), mesh).first;
            for (int l = 0; l < 3; ++l){
                if (v0 != *opit && v1 != *opit) {
                    res[(i+shft)%3] = {*opit, true};
                    if (!bexist){
                        existing = *opit;
                        bexist = true;
                    }
                    break;
                }
                ++opit;
            }
            ++it;
        }
        if (bexist)
            for (int i = 0; i < 3; ++i)
                if (!res[i].second) res[i].first = existing;

        return res;
    }

    class Timer
    {
    private:
        using clock_t = std::chrono::high_resolution_clock;
        using second_t = std::chrono::duration<double, std::ratio<1> >;

        std::chrono::time_point<clock_t> m_beg;

    public:
        Timer() : m_beg(clock_t::now())
        {
        }

        void reset()
        {
            m_beg = clock_t::now();
        }

        double elapsed() const
        {
            return std::chrono::duration_cast<second_t>(clock_t::now() - m_beg).count();
        }
    };
}



#endif //AORTIC_VALVE_OBJECT3DINLINEHELPERS_H
