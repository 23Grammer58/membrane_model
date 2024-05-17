//
// Created by Liogky Alexey on 17.06.2022.
//

#include "ContactForce.h"
#include "../Collision/BulletVolumeTree.h"

World3d::ContactForce::ContactForce(){
    type = "ContactForce";
    m_vf = [](){ return std::make_shared<BulletVolumeTree>(); };
}

World3d::ContactForce::ContactForce(int contact_type){
    type = "ContactForce";
    m_vf = [](){ return std::make_shared<BulletVolumeTree>(); };
    m_ContactType = contact_type;
};
World3d::ContactForce::ContactForce(int contact_type, DistFunc f, DistFunc df){
    type = "ContactForce";
    m_vf = [](){ return std::make_shared<BulletVolumeTree>(); };
    m_ContactType = contact_type;
    m_F = std::move(f);
    m_dF = std::move(df);
}

void World3d::ContactForce::ObjectWrapper::updateTree(bool enforce) {
    if (m_obj == nullptr) return;
    if (!m_fd || m_obj->m_mesh.num_faces() != m_fd->size()){
        if (!m_fd) m_fd = std::make_shared<std::vector<FaceData>>();
        m_fd->reserve(m_obj->m_mesh.num_faces());
        for (auto f: m_obj->m_mesh.faces()) {
            m_fd->emplace_back(f);
            auto v = vert_around(m_obj->m_mesh, f);
            double loc_thickness = shape_thickness(f.idx(), CollisionThicknessBase::FACE);
            m_fd->back().m_leaf = m_fbvt->insert(BoxVol::FromPoints({m_obj->m_x[v[0]], m_obj->m_x[v[1]], m_obj->m_x[v[2]]}).Expand(loc_thickness), &m_fd->back());
        }
    } else if (is_dynamic || enforce) {
        for (auto& fd: *m_fd) {
            auto v = vert_around(m_obj->m_mesh, fd.id);
            double loc_thickness = shape_thickness(fd.id.idx(), CollisionThicknessBase::FACE);
            auto newvol = BoxVol::FromPoints({m_obj->m_x[v[0]], m_obj->m_x[v[1]], m_obj->m_x[v[2]]}).Expand(loc_thickness);
            m_fbvt->update(fd.m_leaf, newvol);
        }
    }
}

bool World3d::ContactForce::registerWorld(World3d::World &w) {
    if (m_w == &w) return true;
    m_obws.clear();
    m_w = &w;
    unsigned N = 0;
    for (auto& o: w.objs()){
        ObjectWrapper ow;
        ow.m_obj = &o.second.first;
        ow.is_dynamic = o.second.second > 0;
        ow.shape_thickness = CollisionThicknessConst(m_DefThickParam * compute_mid_edge_len(o.second.first));
        if (o.second.second) {
            ow.setVertIndOff(N);
            N += 3 * o.second.first.m_mesh.num_vertices();
        }
        m_obws.insert({o.first, std::move(ow)}).first->second.setVT(m_vf());
    }
    return true;
}

bool World3d::ContactForce::addObj(Object3D &obj, const World3d::ObjectID &id, bool is_dynamic) {
    ObjectWrapper ow;
    ow.m_obj = &obj;
    ow.is_dynamic = is_dynamic;
    ow.shape_thickness = CollisionThicknessConst(m_DefThickParam * compute_mid_edge_len(obj));
    auto it = m_obws.insert({id, std::move(ow)});
    if (it.second) {
        it.first->second.setVT(m_vf());
        auto& objs = m_w->objs();
        unsigned N = 0;
        for (auto it = objs.begin(); it != objs.end(); ++it){
            if (!it->second.second) continue;
            auto& o = it->second.first;
            auto jt = m_obws.find(it->first);
            if (jt != m_obws.end()) jt->second.setVertIndOff(N);
            N += 3 * o.m_mesh.num_vertices();
        }
    }
    return it.second;
}

bool World3d::ContactForce::setShapeThickness(const World3d::ObjectID &id, CollisionThickness thickness) {
    auto it = m_obws.find(id);
    if (it != m_obws.end()) return it->second.shape_thickness = thickness, true;
    return false;
}

void World3d::ContactForce::removeObj(const World3d::ObjectID &id) { m_obws.erase(id); }

int World3d::ContactForce::operator()() {
    assert(m_F && "Force function was not set");
    updateProximityTree();
    for (auto& i: m_contacts) i.resize(0);
    performCulling();
    contactRemoveDuplicates();
    setAllForces();

    return 0;
}

int World3d::ContactForce::fill_matrix(SparseMatrix *sm) {
    assert(m_dF && "Derivative force function was not set");
    setAlldForces(sm);
    return 0;
}

std::unique_ptr<World3d::WorldForceBase> World3d::ContactForce::clone() {
    auto other = std::make_unique<ContactForce>();
    other->m_vf = m_vf;
    other->m_w = m_w;
    other->m_obws = m_obws;
    other->m_DefThickParam = m_DefThickParam;
    other->m_ContactType = m_ContactType;
    return other;
}

void World3d::ContactForce::performCulling() {
    if (m_ContactType & SelfMask) {
        for (auto& obwt: m_obws) {
            auto& obw = obwt.second;
            if (!obw.is_dynamic) continue;
            VolumeIntersect(obw.getVT(), obw.getVT(),
                            [this, &obw](const VolumeTreeBase::LeafIndex l0, void *data0,
                                         const VolumeTreeBase::LeafIndex l1, void *data1) -> void {
                                auto *fd0 = static_cast<ObjectWrapper::FaceData *>(data0), *fd1 = static_cast<ObjectWrapper::FaceData *>(data1);
                                if (fd0 == fd1) return;
                                CollideElem ce0{&obw, fd0->id}, ce1{&obw, fd1->id};
                                this->checkAndCollectContact({ce0, ce1});
                            });
        }
    }
    if ((m_ContactType & OtherMask) && (m_obws.size() >= 2)) {
        auto it = m_obws.begin();
        auto st = it++;
        for (; it != m_obws.end(); ++it)
        for (auto jt = st; jt != it; ++jt){
                auto &obw0 = jt->second, &obw1 = it->second;
                if (!obw0.is_dynamic && !obw1.is_dynamic) continue;
                VolumeIntersect(obw0.getVT(), obw1.getVT(),
                                [this, &obw0, &obw1](const VolumeTreeBase::LeafIndex l0, void *data0,
                                                     const VolumeTreeBase::LeafIndex l1, void *data1) -> void {
                                    auto *fd0 = static_cast<ObjectWrapper::FaceData *>(data0), *fd1 = static_cast<ObjectWrapper::FaceData *>(data1);
                                    CollideElem ce0{&obw0, fd0->id}, ce1{&obw1, fd1->id};
                                    this->checkAndCollectContact({ce0, ce1});
                                });
            }
    }
}

bool World3d::ContactForce::computeVV(World3d::ContactForce::Contact &c, World3d::DReal h,
                                      World3d::ContactForce::ObjectWrapper *o0, World3d::V_ind v0,
                                      World3d::ContactForce::ObjectWrapper *o1, World3d::V_ind v1) {
    c.diff = o1->m_obj->m_x[v1] - o0->m_obj->m_x[v0];
    c.sqd = c.diff.squared_length();
    if (c.sqd > h * h) return false;
    c.w[0] = 1, c.w[1] = 1;
    c.elems[0] = Contact::Elem{o0, static_cast<int>(v0.idx()), o0->m_obj->m_x[v0]};
    c.elems[1] = Contact::Elem{o1, static_cast<int>(v1.idx()), o1->m_obj->m_x[v1]};
    if (c.sqd > 2 * (abs(c.diff[0]) + abs(c.diff[1]) + abs(c.diff[2])) * std::numeric_limits<DReal>::epsilon()){
        c.d = sqrt(c.sqd);
        c.n = c.diff / c.d;
    }
    c.loc_thickness = h;
    return true;
}

static inline void setForceOnObject(World3d::Object3D* obj, bool is_moving, World3d::V_ind v, const World3d::Vector F){
    if (!is_moving) return;
    obj->m_F[v] += obj->withBCmask(v, F);
}

void World3d::ContactForce::setVVForce(World3d::ContactForce::Contact &c, bool is_self_inter) {
    double sF = m_F(is_self_inter ? VertVertSelf : VertVert, c);
    setForceOnObject(c.elems[0].obw->m_obj, c.elems[0].obw->is_dynamic, V_ind(c.elems[0].geom_elem_id), -sF * c.n);
    setForceOnObject(c.elems[1].obw->m_obj, c.elems[1].obw->is_dynamic, V_ind(c.elems[1].geom_elem_id), +sF * c.n);
}

void World3d::ContactForce::setVEForce(World3d::ContactForce::Contact &c, bool is_self_inter) {
    double sF = m_F(is_self_inter ? VertEdgeSelf : VertEdge, c);
    setForceOnObject(c.elems[0].obw->m_obj, c.elems[0].obw->is_dynamic, V_ind(c.elems[0].geom_elem_id), -sF * c.n);
    auto vv = vert_around(c.elems[1].obw->m_obj->m_mesh, E_ind(c.elems[1].geom_elem_id));
    setForceOnObject(c.elems[1].obw->m_obj, c.elems[1].obw->is_dynamic, vv[0], +c.w[1 + 0] * sF * c.n);
    setForceOnObject(c.elems[1].obw->m_obj, c.elems[1].obw->is_dynamic, vv[1], +c.w[1 + 1] * sF * c.n);
}

void World3d::ContactForce::setVFForce(World3d::ContactForce::Contact &c, bool is_self_inter) {
    double sF = m_F(is_self_inter ? VertFaceSelf : VertFace, c);
    setForceOnObject(c.elems[0].obw->m_obj, c.elems[0].obw->is_dynamic, V_ind(c.elems[0].geom_elem_id), -sF * c.n);
    auto vv = vert_around(c.elems[1].obw->m_obj->m_mesh, F_ind(c.elems[1].geom_elem_id));
    setForceOnObject(c.elems[1].obw->m_obj, c.elems[1].obw->is_dynamic, vv[0], +c.w[1 + 0] * sF * c.n);
    setForceOnObject(c.elems[1].obw->m_obj, c.elems[1].obw->is_dynamic, vv[1], +c.w[1 + 1] * sF * c.n);
    setForceOnObject(c.elems[1].obw->m_obj, c.elems[1].obw->is_dynamic, vv[2], +c.w[1 + 2] * sF * c.n);
}

void World3d::ContactForce::setEEForce(World3d::ContactForce::Contact &c, bool is_self_inter) {
    double sF = m_F(is_self_inter ? EdgeEdgeSelf : EdgeEdge, c);
    auto vv0 = vert_around(c.elems[0].obw->m_obj->m_mesh, E_ind(c.elems[0].geom_elem_id)),
         vv1 = vert_around(c.elems[1].obw->m_obj->m_mesh, E_ind(c.elems[1].geom_elem_id));
    setForceOnObject(c.elems[0].obw->m_obj, c.elems[0].obw->is_dynamic, vv0[0], -c.w[0 + 0] * sF * c.n);
    setForceOnObject(c.elems[0].obw->m_obj, c.elems[0].obw->is_dynamic, vv0[1], -c.w[0 + 1] * sF * c.n);
    setForceOnObject(c.elems[1].obw->m_obj, c.elems[1].obw->is_dynamic, vv1[0], +c.w[2 + 0] * sF * c.n);
    setForceOnObject(c.elems[1].obw->m_obj, c.elems[1].obw->is_dynamic, vv1[1], +c.w[2 + 1] * sF * c.n);
}

void World3d::ContactForce::setXXdForce(SparseMatrix *sm, World3d::ContactForce::Contact &c, int contact_type, int nel1,
                                        int nel2, const std::array<V_ind, 4>& vv) {
    if (c.d == 0) return;
    std::array<unsigned char, 4> o;
    std::fill(o.data(), o.data() + nel1, 0);
    std::fill(o.data() + nel1, o.data() + nel1 + nel2, 1);
    std::array<char, 4> r;
    std::fill(r.data(), r.data() + nel1, -1);
    std::fill(r.data() + nel1, r.data() + nel1 + nel2, 1);
    //dF_{qi} / dQ_{lj} = s_q * (df/d(d) * d(d)/dQ_{lj} * n_i + f(d) * dn_i / dQ_{lj}); q,l = 0->3, i,j = 0->2
    //s_q = r_q * w_q
    //d(d)/dQ_{lj} = r_l * w_l * n_j
    //dn_i / dQ_{lj} = r_l * w_l * (δ_{i,j} - n_i*n_j) / d
    double F = m_F(contact_type, c);
    double dF = m_dF(contact_type, c);
    int nobj = nel1 + nel2;
    for (int q = 0; q < nobj; ++q) if (c.elems[o[q]].obw->is_dynamic)
    for (int i = 0; i < 3; ++i) if (c.elems[o[q]].obw->m_obj->is_movable(vv[q], i))
    for (int l = 0; l < nobj; ++l) if (c.elems[o[l]].obw->is_dynamic)
    for (int j = 0; j < 3; ++j) if (c.elems[o[q]].obw->m_obj->is_movable(vv[l], j)){
        auto val = r[q]*c.w[q]*r[l]*c.w[l] * (dF * c.n[i] * c.n[j] + F * (((i==j) ? 1 : 0) - c.n[i] * c.n[j]) / c.d);
        if (val != 0)
            (*sm)(c.elems[o[q]].obw->getVIndOff() + 3 * vv[q].idx() + i, c.elems[o[l]].obw->getVIndOff() + 3 * vv[l].idx() + j) += val;
    }

}

void World3d::ContactForce::setAllForces() {
    for (auto& c: m_contacts[0]) setVVForce(c, false);
    for (auto& c: m_contacts[1]) setVEForce(c, false);
    for (auto& c: m_contacts[2]) setVFForce(c, false);
    for (auto& c: m_contacts[3]) setEEForce(c, false);
    for (auto& c: m_contacts[4]) setVVForce(c, true);
    for (auto& c: m_contacts[5]) setVEForce(c, true);
    for (auto& c: m_contacts[6]) setVFForce(c, true);
    for (auto& c: m_contacts[7]) setEEForce(c, true);
}

void World3d::ContactForce::setAlldForces(SparseMatrix *sm) {
    for (auto& c: m_contacts[0]) setVVdForce(sm, c, false);
    for (auto& c: m_contacts[1]) setVEdForce(sm, c, false);
    for (auto& c: m_contacts[2]) setVFdForce(sm, c, false);
    for (auto& c: m_contacts[3]) setEEdForce(sm, c, false);
    for (auto& c: m_contacts[4]) setVVdForce(sm, c, true);
    for (auto& c: m_contacts[5]) setVEdForce(sm, c, true);
    for (auto& c: m_contacts[6]) setVFdForce(sm, c, true);
    for (auto& c: m_contacts[7]) setEEdForce(sm, c, true);
}

void World3d::ContactForce::setForce(World3d::ContactForce::Contact &c, int contact_type) {
    switch (contact_type) {
        case VertVert: setVVForce(c, false); break;
        case VertEdge: setVEForce(c, false); break;
        case VertFace: setVFForce(c, false); break;
        case EdgeEdge: setEEForce(c, false); break;
        case VertVertSelf: setVVForce(c, true); break;
        case VertEdgeSelf: setVEForce(c, true); break;
        case VertFaceSelf: setVFForce(c, true); break;
        case EdgeEdgeSelf: setEEForce(c, true); break;
        default: throw std::runtime_error("Faced unknown contact type = " + std::to_string(contact_type));
    }
}

void World3d::ContactForce::setdForce(SparseMatrix *sm, World3d::ContactForce::Contact &c, int contact_type) {
    switch (contact_type) {
        case VertVert: setVVdForce(sm, c, false); break;
        case VertEdge: setVEdForce(sm, c, false); break;
        case VertFace: setVFdForce(sm, c, false); break;
        case EdgeEdge: setEEdForce(sm, c, false); break;
        case VertVertSelf: setVVdForce(sm, c, true); break;
        case VertEdgeSelf: setVEdForce(sm, c, true); break;
        case VertFaceSelf: setVFdForce(sm, c, true); break;
        case EdgeEdgeSelf: setEEdForce(sm, c, true); break;
        default: throw std::runtime_error("Faced unknown contact type = " + std::to_string(contact_type));
    }
}

bool World3d::ContactForce::computeVE(World3d::ContactForce::Contact &c, World3d::DReal h,
                                      World3d::ContactForce::ObjectWrapper *o0, World3d::V_ind v0,
                                      World3d::ContactForce::ObjectWrapper *o1, World3d::E_ind e1) {
    DReal max_sqd = h * h;
    auto v1 = vert_around(o1->m_obj->m_mesh, e1);
    if (!GeomProjs::ProjectOnSegment(o1->m_obj->m_x[v1[0]], o1->m_obj->m_x[v1[1]], o0->m_obj->m_x[v0], c.elems[1].proj, c.w[1], c.sqd, &max_sqd))
        return false;
    c.diff = c.elems[1].proj - o0->m_obj->m_x[v0];
    c.w[0] = 1; c.w[2] = 1 - c.w[1];
    c.elems[0].proj = o0->m_obj->m_x[v0], c.elems[0].obw = o0, c.elems[0].geom_elem_id = v0.idx();
    c.elems[1].obw = o1, c.elems[1].geom_elem_id = e1.idx();
    if (c.sqd > 2 * (abs(c.diff[0]) + abs(c.diff[1]) + abs(c.diff[2])) * std::numeric_limits<DReal>::epsilon()){
        c.d = sqrt(c.sqd);
        c.n = c.diff / c.d;
    }
    c.loc_thickness = h;
    return true;
}

bool World3d::ContactForce::computeVF(World3d::ContactForce::Contact &c, World3d::DReal h,
                                      World3d::ContactForce::ObjectWrapper *o0, World3d::V_ind v0,
                                      World3d::ContactForce::ObjectWrapper *o1, World3d::F_ind f1) {
    DReal max_sqd = h * h;
    auto v1 = vert_around(o1->m_obj->m_mesh, f1);
    if (!GeomProjs::ProjectOnTriangle(o1->m_obj->m_x[v1[0]], o1->m_obj->m_x[v1[1]], o1->m_obj->m_x[v1[2]], o0->m_obj->m_x[v0], c.elems[1].proj, c.sqd, &max_sqd))
        return false;
    c.diff = c.elems[1].proj - o0->m_obj->m_x[v0];
    c.w[0] = 1;
    auto ww = GeomProjs::BaryCoord(o1->m_obj->m_x[v1[0]], o1->m_obj->m_x[v1[1]], o1->m_obj->m_x[v1[2]], c.elems[1].proj);
    for (int k = 0; k < 3; ++k) c.w[k+1] = ww[k];
    c.elems[0].proj = o0->m_obj->m_x[v0], c.elems[0].obw = o0, c.elems[0].geom_elem_id = v0.idx();
    c.elems[1].obw = o1, c.elems[1].geom_elem_id = f1.idx();
    if (c.sqd > 2 * (abs(c.diff[0]) + abs(c.diff[1]) + abs(c.diff[2])) * std::numeric_limits<DReal>::epsilon()){
        c.d = sqrt(c.sqd);
        c.n = c.diff / c.d;
    }
    c.loc_thickness = h;
    return true;
}

bool World3d::ContactForce::computeEE(World3d::ContactForce::Contact &c, World3d::DReal h,
                                      World3d::ContactForce::ObjectWrapper *o0, World3d::E_ind e0,
                                      World3d::ContactForce::ObjectWrapper *o1, World3d::E_ind e1) {
    DReal max_sqd = h * h;
    auto v0 = vert_around(o0->m_obj->m_mesh, e0), v1 = vert_around(o1->m_obj->m_mesh, e1);
    Vector vprj0, vprj1;
    Point o = CGAL::ORIGIN;
    if (!GeomProjs::SegmentsProjects(o1->m_obj->m_x[v1[0]]-o, o1->m_obj->m_x[v1[1]]-o, o0->m_obj->m_x[v0[0]]-o, o0->m_obj->m_x[v0[1]]-o, vprj1, vprj0, c.w[3], c.w[1], c.sqd, &max_sqd))
        return false;
    c.elems[0].proj = o + vprj0, c.elems[1].proj = o + vprj1;
    c.diff = c.elems[1].proj - c.elems[0].proj;
    c.w[0] = 1 - c.w[1], c.w[2] = 1 - c.w[3];
    c.elems[0].obw = o0, c.elems[0].geom_elem_id = e0.idx();
    c.elems[1].obw = o1, c.elems[1].geom_elem_id = e1.idx();
    if (c.sqd > 2 * (abs(c.diff[0]) + abs(c.diff[1]) + abs(c.diff[2])) * std::numeric_limits<DReal>::epsilon()){
        c.d = sqrt(c.sqd);
        c.n = c.diff / c.d;
    }
    c.loc_thickness = h;
    return true;
}

void World3d::ContactForce::collectVV(World3d::ContactForce::ObjectWrapper *o0, World3d::V_ind v0,
                                      World3d::ContactForce::ObjectWrapper *o1, World3d::V_ind v1) {
    DReal h = eval_thickness(o0, o1, v0.idx(), v1.idx(), CollisionThicknessBase::VERTEX, CollisionThicknessBase::VERTEX);
    Contact c;
    if (o0 > o1) { std::swap(o0, o1), std::swap(v0, v1); }
    if (computeVV(c, h, o0, v0, o1, v1)) m_contacts[0].push_back(c);
}

void World3d::ContactForce::collectVE(World3d::ContactForce::ObjectWrapper *o0, World3d::V_ind v0,
                                      World3d::ContactForce::ObjectWrapper *o1, World3d::E_ind e1) {
    DReal h = eval_thickness(o0, o1, v0.idx(), e1.idx(), CollisionThicknessBase::VERTEX, CollisionThicknessBase::EDGE);
    Contact c;
    if (computeVE(c, h, o0, v0, o1, e1)) m_contacts[1].push_back(c);
}

void World3d::ContactForce::collectVF(World3d::ContactForce::ObjectWrapper *o0, World3d::V_ind v0,
                                      World3d::ContactForce::ObjectWrapper *o1, World3d::F_ind f1) {
    DReal h = eval_thickness(o0, o1, v0.idx(), f1.idx(), CollisionThicknessBase::VERTEX, CollisionThicknessBase::FACE);
    Contact c;
    if (computeVF(c, h, o0, v0, o1, f1)) m_contacts[2].push_back(c);
}

void World3d::ContactForce::collectEE(World3d::ContactForce::ObjectWrapper *o0, World3d::E_ind e0,
                                      World3d::ContactForce::ObjectWrapper *o1, World3d::E_ind e1) {
    DReal h = eval_thickness(o0, o1, e0.idx(), e1.idx(), CollisionThicknessBase::EDGE, CollisionThicknessBase::EDGE);
    Contact c;
    if (o0 > o1) { std::swap(o0, o1), std::swap(e0, e1); }
    if (computeEE(c, h, o0, e0, o1, e1)) m_contacts[3].push_back(c);
}

void World3d::ContactForce::collectVVS(World3d::ContactForce::ObjectWrapper *o, World3d::V_ind v0, World3d::V_ind v1) {
    if (v1 < v0) std::swap(v1, v0);
    DReal h = eval_thickness(o, o, v0.idx(), v1.idx(), CollisionThicknessBase::VERTEX, CollisionThicknessBase::VERTEX);
    Contact c;
    if (computeVV(c, h, o, v0, o, v1)) m_contacts[4].push_back(c);
}

void World3d::ContactForce::collectVES(World3d::ContactForce::ObjectWrapper *o, World3d::V_ind v0, World3d::E_ind e1) {
    DReal h = eval_thickness(o, o, v0.idx(), e1.idx(), CollisionThicknessBase::VERTEX, CollisionThicknessBase::EDGE);
    Contact c;
    if (computeVE(c, h, o, v0, o, e1)) m_contacts[5].push_back(c);
}

void World3d::ContactForce::collectVFS(World3d::ContactForce::ObjectWrapper *o, World3d::V_ind v0, World3d::F_ind f1) {
    DReal h = eval_thickness(o, o, v0.idx(), f1.idx(), CollisionThicknessBase::VERTEX, CollisionThicknessBase::FACE);
    Contact c;
    if (computeVF(c, h, o, v0, o, f1)) m_contacts[6].push_back(c);
}

void World3d::ContactForce::collectEES(World3d::ContactForce::ObjectWrapper *o, World3d::E_ind e0, World3d::E_ind e1) {
    if (e1 < e0) std::swap(e1, e0);
    DReal h = eval_thickness(o, o, e0.idx(), e1.idx(), CollisionThicknessBase::EDGE, CollisionThicknessBase::EDGE);
    Contact c;
    if (computeEE(c, h, o, e0, o, e1)) m_contacts[7].push_back(c);
}

void World3d::ContactForce::checkAndCollectContact(const std::array<CollideElem, 2> &cp) {
    if (cp[0].obw == cp[1].obw){
        if (!(m_ContactType & SelfMask)) return;
        if (cp[0].id == cp[1].id) return;
        std::array<V_ind, 3> vv[2] = {vert_around(cp[0].obw->m_obj->m_mesh, cp[0].id), vert_around(cp[1].obw->m_obj->m_mesh, cp[1].id)};
        std::array<E_ind, 3> ee[2];
        std::array<V_ind, 6> ve[2];
        if (m_ContactType & (VertEdgeSelf | EdgeEdgeSelf)){
            for (int i = 0; i < 2; ++i){
                ee[i] = edge_around(cp[i].obw->m_obj->m_mesh, cp[i].id);
                for (int k = 0; k < 3; ++k){
                    std::array<V_ind, 2> vve = vert_around(cp[i].obw->m_obj->m_mesh, ee[i][k]);
                    ve[i][2*k + 0] = vve[0], ve[i][2*k + 1] = vve[1];
                }
            }
        }
        std::array<char, 3> link[2];
        int nocommon = 3;
        {
            for (char i = 0; i < 3; ++i)
                link[0][i] = link[1][i] = -1;
            for (char i = 0; i < 3; ++i) {
                for (char j = 0; j < 3; ++j)
                    if (vv[0][i] == vv[1][j]) {
                        link[0][i] = j;
                        link[1][j] = i;
                        --nocommon;
                        break;
                    }
            }
        }
        if (m_ContactType & VertVertSelf) {
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j)
                    if (vv[0][i] != vv[1][j]) collectVVS(cp[0].obw, vv[0][i], vv[1][j]);
        }
        if (m_ContactType & VertEdgeSelf){
            for (int k = 0; k < 2; ++k)
                for (int i = 0; i < 3; ++i)
                    for (int j = 0; j < 3; ++j)
                        if (vv[k][i] != ve[(k+1)%2][2*j + 0] && vv[k][i] != ve[(k+1)%2][2*j + 1])
                            collectVES(cp[k].obw, vv[k][i], ee[(k+1)%2][j]);
        }
        if (m_ContactType & VertFaceSelf){
            for (int k = 0; k < 2; ++k)
                for (int i = 0; i < 3; ++i)
                    if (link[k][i] < 0) collectVFS(cp[k].obw, vv[k][i], cp[(k+1)%2].id);
        }
        if (m_ContactType & EdgeEdgeSelf){
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j) {
                    if (ve[0][2*i + 0] != ve[1][2*i + 0] && ve[0][2*i + 0] != ve[1][2*i + 1] &&
                        ve[0][2*i + 1] != ve[1][2*i + 0] && ve[0][2*i + 1] != ve[1][2*i + 1])
                        collectEES(cp[0].obw, ee[0][i], ee[1][j]);
                }
        }
    } else {
        if (!(m_ContactType & OtherMask)) return;
        std::array<V_ind, 3> vv[2] = {vert_around(cp[0].obw->m_obj->m_mesh, cp[0].id), vert_around(cp[1].obw->m_obj->m_mesh, cp[1].id)};
        std::array<E_ind, 3> ee[2];
        std::array<V_ind, 6> ve[2];
        if (m_ContactType & (VertEdgeSelf | EdgeEdgeSelf)){
            for (int i = 0; i < 2; ++i){
                ee[i] = edge_around(cp[i].obw->m_obj->m_mesh, cp[i].id);
                for (int k = 0; k < 3; ++k){
                    std::array<V_ind, 2> vve = vert_around(cp[i].obw->m_obj->m_mesh, ee[i][k]);
                    ve[i][2*k + 0] = vve[0], ve[i][2*k + 1] = vve[1];
                }
            }
        }
        if (m_ContactType & VertVert)
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j)
                    collectVV(cp[0].obw, vv[0][i], cp[1].obw, vv[1][j]);
        if (m_ContactType & VertEdge)
            for (int k = 0; k < 2; ++k)
                for (int i = 0; i < 3; ++i)
                    for (int j = 0; j < 3; ++j)
                        collectVE(cp[k].obw, vv[k][i], cp[(k+1)%2].obw, ee[(k+1)%2][j]);
        if (m_ContactType & VertFace){
            for (int k = 0; k < 2; ++k)
                for (int i = 0; i < 3; ++i)
                    collectVF(cp[k].obw, vv[k][i], cp[(k+1)%2].obw, cp[(k+1)%2].id);
        }
        if (m_ContactType & EdgeEdge){
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j)
                    collectEE(cp[0].obw, ee[0][i], cp[0].obw, ee[1][j]);
        }
    }
}

void World3d::ContactForce::contactRemoveDuplicates() {
    for (auto& i: m_contacts){ //remove duplicates
        if (i.empty()) continue;
        std::sort(i.begin(), i.end(), [](const Contact& a, const Contact& b) {
            if (a.elems[0].obw != b.elems[0].obw) return a.elems[0].obw < b.elems[0].obw;
            if (a.elems[1].obw != b.elems[1].obw) return a.elems[1].obw < b.elems[1].obw;
            if (a.elems[0].geom_elem_id != b.elems[0].geom_elem_id) return a.elems[0].geom_elem_id < b.elems[0].geom_elem_id;
            return a.elems[1].geom_elem_id < b.elems[1].geom_elem_id;
        });
        auto is_eq = [](const Contact& a, const Contact& b){
            return a.elems[0].obw == b.elems[0].obw && a.elems[1].obw == b.elems[1].obw &&
                   a.elems[0].geom_elem_id != b.elems[0].geom_elem_id && a.elems[1].geom_elem_id == b.elems[1].geom_elem_id;
        };
        auto first = i.begin(), last = i.end();
        auto result = first;
        while (++first != last) {
            if (!is_eq(*result, *first) && ++result != first) {
                *result = std::move(*first);
            }
        }
        i.resize(std::distance(++result, first));
    }
}

World3d::DReal World3d::ContactForce::compute_mid_edge_len(const Object3D &o) {
    DReal l = 0;
    for (auto e: o.m_mesh.edges()) {
        auto v = vert_around(o.m_mesh, e);
        l += sqrt((o.m_x[v[1]] - o.m_x[v[0]]).squared_length());
    }
    if (o.m_mesh.num_edges() > 0) l /= o.m_mesh.num_edges();
    return l;
}

void World3d::ContactForce::setVVdForce(SparseMatrix *sm, World3d::ContactForce::Contact &c, bool is_self_inter) {
    std::array<V_ind, 4> vv{V_ind(c.elems[0].geom_elem_id), V_ind(c.elems[1].geom_elem_id)};
    setXXdForce(sm, c, is_self_inter ? VertVertSelf : VertVert, 1, 1, vv);
}

void World3d::ContactForce::setVEdForce(SparseMatrix *sm, World3d::ContactForce::Contact &c, bool is_self_inter) {
    auto ve = vert_around(c.elems[1].obw->m_obj->m_mesh, E_ind(c.elems[1].geom_elem_id));
    std::array<V_ind, 4> vv{V_ind(c.elems[0].geom_elem_id), ve[0], ve[1]};
    setXXdForce(sm, c, is_self_inter ? VertEdgeSelf : VertEdge, 1, 2, vv);
}

void World3d::ContactForce::setVFdForce(SparseMatrix *sm, World3d::ContactForce::Contact &c, bool is_self_inter) {
    auto vf = vert_around(c.elems[1].obw->m_obj->m_mesh, F_ind(c.elems[1].geom_elem_id));
    std::array<V_ind, 4> vv{V_ind(c.elems[0].geom_elem_id), vf[0], vf[1], vf[2]};
    setXXdForce(sm, c, is_self_inter ? VertFaceSelf : VertFace, 1, 3, vv);
}

void World3d::ContactForce::setEEdForce(SparseMatrix *sm, World3d::ContactForce::Contact &c, bool is_self_inter) {
    auto ve0 = vert_around(c.elems[0].obw->m_obj->m_mesh, E_ind(c.elems[0].geom_elem_id)),
            ve1 = vert_around(c.elems[1].obw->m_obj->m_mesh, E_ind(c.elems[1].geom_elem_id));
    std::array<V_ind, 4> vv{ve0[0], ve0[1], ve1[0], ve1[1]};
    setXXdForce(sm, c, is_self_inter ? EdgeEdgeSelf : EdgeEdge, 2, 2, vv);
}

World3d::DReal World3d::ContactForce::eval_thickness(World3d::ContactForce::ObjectWrapper* o0, World3d::ContactForce::ObjectWrapper* o1,
                                    int id0, int id1, 
                                    int primitive_type0, int primitive_type1){
    if (o0->m_obj == o1->m_obj) return o0->shape_thickness(id0, id1, primitive_type0, primitive_type1);
    else return (o0->shape_thickness(id0, primitive_type0) + o0->shape_thickness(id1, primitive_type1))/2;
}
