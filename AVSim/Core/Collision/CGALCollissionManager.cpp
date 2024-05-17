//
// Created by alex on 17.11.2021.
//

#include "CGALCollissionManager.h"

int World3d::CGALColManVFfull::default_solve_collision_algo(std::map<ObjectID, ObjColInfo> &shape) {
    for (auto& ob: shape){
        Object3D* to = ob.second.obj;
        for (auto col: ob.second.cd){
            Object3D* from = col.from->obj;
            if (col.from->status & 1) { //soft-soft contact
                Vector nr = col.nr;
                auto fvv = vert_around(from->m_mesh, col.f);
                auto w = col.w;
                double ma = 0.5, mb = 0.5;
                const Point p = CGAL::ORIGIN + BaryEval(from->m_next_x[fvv[0]] - CGAL::ORIGIN,from->m_next_x[fvv[1]] - CGAL::ORIGIN,from->m_next_x[fvv[2]] - CGAL::ORIGIN,w);
                const Point q = CGAL::ORIGIN + BaryEval(from->m_x[fvv[0]] - CGAL::ORIGIN, from->m_x[fvv[1]] - CGAL::ORIGIN, from->m_x[fvv[2]] - CGAL::ORIGIN, w);
                const Vector vr = (to->m_next_x[col.v] - to->m_x[col.v]) - (p - q);
                Vector corr{CGAL::NULL_VECTOR};
                double dot = vr * nr;
                if (dot < 0) {
                    const double j = col.margin - nr * (to->m_next_x[col.v] - p);
                    corr += j * nr;

                    if (to->is_movable(col.v)
                        || !from->is_movable(fvv[0]) || !from->is_movable(fvv[1]) || !from->is_movable(fvv[2])) {
                        auto dif = ((nr * from->m_normal[col.f] >= 0) ? 1.0 : -1.0) * ma * corr;
                        to->m_next_x[col.v] += to->withBCmask(col.v, dif);
                    }
                    for (unsigned j = 0; j < 3; ++j)
                        if (from->is_movable(fvv[j]) || !to->is_movable(col.v)) {
                            auto dif = -mb * static_cast<double>(w[j]) * corr;
                            from->m_next_x[fvv[j]] += from->withBCmask(fvv[j], dif);
                        }
                }
            } else { //soft-rigid contact
                if (!to->is_movable(col.v)) continue;
                Point &next = to->m_next_x[col.v];
                const Vector vr = to->m_next_x[col.v] - to->m_x[col.v];
                Vector corr{CGAL::NULL_VECTOR};
                double dot = vr * col.nr, dott = col.nr * from->m_normal[col.f];
                if (dot < 0 && dott >= 0) {
                    const double j = col.margin - col.nr * (next - col.proj);
                    corr += j * col.nr;
                } else if (dott < 0) {
                    const double j = col.margin + col.nr * (next - col.proj);
                    corr -= j * col.nr;
                }
                next += corr;
            }
        }
    }
    return 0;
}

bool World3d::CGALColManVFfull::set_margin(ObjectID id, double margin) {
    auto it = m_shape.find(id);
    if (it == m_shape.end()) return false;
    it->second.margin = 2*margin;
    return true;
}

void World3d::CGALColManVFfull::findCollisions() {
    //std::vector<VBBox> vb;
    //std::vector<FBBox> fb;
//    static int it = 0;
//    static double time = 0;

    for (auto& i: m_shape){
        auto& obj = *i.second.obj;
        std::function<Point(V_ind)> getP = [&obj](V_ind v) { return obj.m_next_x[v]; };
        if ((i.second.status & 1) == 0) getP = [&obj](V_ind v) { return obj.m_x[v]; };
        auto mrg = i.second.margin;
        auto& fb = i.second.fb;
        if (fb.empty()) {
            fb.resize(obj.m_mesh.num_faces());
            int ind = 0;
            for (auto f: obj.m_mesh.faces()) {
                auto v = vert_around(obj.m_mesh, f);
                fb[ind++] = FBBox(make_bbox(getP(v[0]), getP(v[1]), getP(v[2]), mrg), f);
            }
        } else {
            for (auto& fd: fb) {
                auto f = F_ind(fd.handle());
                auto v = vert_around(obj.m_mesh, f);
                fd = FBBox(make_bbox(getP(v[0]), getP(v[1]), getP(v[2]), mrg), f);
            }
        }
        if ((i.second.status & 1) != 0){
            auto& vb = i.second.vb;
            if (vb.empty()){
                vb.resize(obj.m_mesh.num_vertices());
                int ind = 0;
                for (auto v: obj.m_mesh.vertices()) {
                    vb[ind++] = VBBox(make_bbox(getP(v), mrg), v);
                }
            } else {
                for (auto& vd: vb) {
                    auto v = V_ind(vd.handle());
                    vd = VBBox(make_bbox(getP(v), mrg), v);
                }
            }
        }
    }
    for (auto i1 = m_shape.begin(); i1 != m_shape.end(); ++i1) { //from
        auto &fb = i1->second.fb;
        for (auto i2 = m_shape.begin(); i2 != m_shape.end(); ++i2) { //to
            if (i1 == i2) continue;
            if ((i2->second.status & 1) == 0) continue;
            auto &vb = i1->second.vb;
            double mrg = (i1->second.margin + i2->second.margin) / 2;
            auto check_collision = [ob1 = i1->second.obj, ob2 = i2->second.obj, mrg, &contacts = i2->second.cd, &it = i1](
                    const VBBox &vb, const FBBox &fb) {
                auto fv = World3d::vert_around(ob1->m_mesh, F_ind(fb.handle()));
                std::array<Point, 3> tri = {ob1->m_next_x[fv[0]], ob1->m_next_x[fv[1]], ob1->m_next_x[fv[2]]};
                Point p = CGAL::ORIGIN;
                double sqrd = DBL_MAX;
                auto vv = V_ind(vb.handle());
                ProjectOnTriangle(tri[0], tri[1], tri[2], ob2->m_next_x[vv], p, sqrd);
                double m = mrg + 2 * sqrt((ob2->m_next_x[vv] - ob2->m_x[vv]).squared_length());
                if (sqrd < m * m) {
                    ContactData cd;
                    cd.v = vv;
                    cd.f = F_ind(fb.handle());
                    cd.dist2 = sqrd;
                    cd.proj = p;
                    cd.from = &it->second;
                    cd.margin = m;
                    cd.w = BaryCoord(tri[0], tri[1], tri[2], p) - CGAL::ORIGIN;
                    cd.nr = (ob2->m_next_x[vv] - p) / sqrt(sqrd);
                    contacts.push_back(cd);
                }
            };
//            World3d::Timer tm;
//            tm.reset();
            CGAL::box_intersection_d(vb.begin(), vb.end(), fb.begin(), fb.end(), check_collision);
//            time += tm.elapsed();
        }
    }

//    it++;
//    if (it >= 1000){
//        std::cout << it << ": intersection time = " << time << std::endl;
//        exit(-1);
//    }
}
