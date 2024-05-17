//
// Created by Liogky Alexey on 06.06.2022.
//

#include "BridsonCollisionManager.h"
#include "BulletVolumeTree.h"

using namespace World3d;

void BridsonCollisionManager::SpringAndFrictionParams::_init_data_from_BridsonCollisionManager(BridsonCollisionManager* self, ObjectID id0, ObjectID id1){
    k_spr = self->m_Kspr(id0, id1);
    mu = self->m_friction_mu(id0, id1);
    eps = self->m_RepulsSprCoef(id0, id1);
}

BridsonCollisionManager::BridsonCollisionManager(){
    m_vtfactory = [](){ return std::make_shared<BulletVolumeTree>(); };
}

void BridsonCollisionManager::addObj(Object3D& obj, ObjectID id, bool is_dynamic, double thickness){
    if (is_dynamic){
        auto pinfo = std::make_shared<DynamicObject3DCI>();
        DynamicObject3DCI& info = *pinfo;
        info.m_obj = &obj;
        info.m_id = id;
        info.m_thickness = CollisionThicknessConst(thickness);
        info.m_thickness.registerObj(&obj);
        auto m_xit = obj.m_mesh.add_property_map<V_ind, Point>("v:collision_prev_x");
        info.m_x = m_xit.first;
        if (m_xit.second) for (auto v: obj.m_mesh.vertices()) info.m_x[v] = obj.m_x[v];
        info.m_next_x = obj.m_x;
        info.m_dx = obj.m_mesh.add_property_map<V_ind, Vector>("v:dx").first;
        for (auto v: obj.m_mesh.vertices()) info.m_dx[v] = info.m_next_x[v] - info.m_x[v];
        info.m_v = obj.m_mesh.add_property_map<V_ind, Vector>("v:velocity").first;
        info.m_dv = obj.m_mesh.add_property_map<V_ind, Vector>("v:diff_velocity").first;
        auto it = obj.m_mesh.add_property_map<V_ind, DReal>("v:mass");
        info.m_m = it.first;
        if (it.second) for (auto v: obj.m_mesh.vertices()) info.m_m[v] = 1;
        info.fdat.reserve(obj.m_mesh.num_faces());
        info.m_f0bvt = m_vtfactory();
        info.m_febvt = m_vtfactory();
        for (auto f: obj.m_mesh.faces()){
            auto v = vert_around(obj.m_mesh, f);
            DynamicObject3DCI::FaceData fd;
            fd.id = f;
            fd.backptr = pinfo.get();
            info.fdat.push_back(fd);
            auto& rfd = info.fdat.back();
            rfd.m_f0leaf = info.m_f0bvt->insert(BoxVol::FromPoints({info.m_x[v[0]], info.m_x[v[1]], info.m_x[v[2]]}), &rfd);
            rfd.m_feleaf = info.m_febvt->insert(BoxVol::FromPoints({info.m_x[v[0]], info.m_x[v[1]], info.m_x[v[2]],
                                                                    /*info.m_x[v[0]] , info.m_x[v[1]] , info.m_x[v[2]]*/ }), &rfd);
        }
        m_dobjs.push_back(std::move(pinfo));
    } else
        throw std::runtime_error("Is not implemented");
}

void BridsonCollisionManager::removeObj(ObjectID id){
    int i = 0, j = 0;
    for (; i < m_dobjs.size(); ++i){
        if (m_dobjs[i]->m_id != id) {
            if (j != i) m_dobjs[j] = std::move(m_dobjs[i]);
            ++j;
        }
    }
    if (!m_dobjs.empty() && j != i) m_dobjs.resize(m_dobjs.size() - 1);
}

void BridsonCollisionManager::findCollisions(){ 
    update(); 
}

void BridsonCollisionManager::solveCollisions(){
    applyRelaxRepulsion(m_MaxRelaxIts);
    applyRepulsionImpulses();
    applyRelaxRepulsion(m_MaxRelaxIts);
    if (!applyCCDRepulsions(m_maxCCDIts)) {
        std::cout << "Warrning: Perform impact zone!\n"; //TODO: remove this
        applyRelaxRepulsion(m_MaxRelaxIts);
        applyFailSafe(m_maxFailSafeIts);
    }
    applyDx();
}

void BridsonCollisionManager::updateCCD(){
    for (auto& o: m_dobjs) {
        for (auto& fd: o->fdat) {
            auto v = vert_around(o->m_obj->m_mesh, fd.id);
            auto newvol = BoxVol::FromPoints({o->m_x[v[0]], o->m_x[v[1]], o->m_x[v[2]], o->m_x[v[0]] + o->m_dx[v[0]], o->m_x[v[1]] + o->m_dx[v[1]], o->m_x[v[2]] + o->m_dx[v[2]]}).ScaleExpand(1 + 1e-6);
            o->m_febvt->update(fd.m_feleaf, newvol);
        }
        o->m_febvt->optimize();
    }
}

void BridsonCollisionManager::update(){
    for (auto& o: m_dobjs) {
        for (auto v: o->m_obj->m_mesh.vertices())
            o->m_dx[v] = o->m_next_x[v] - o->m_x[v];
        for (auto& fd: o->fdat) {
            auto v = vert_around(o->m_obj->m_mesh, fd.id);
            auto newvol = BoxVol::FromPoints({o->m_x[v[0]], o->m_x[v[1]], o->m_x[v[2]]}).Expand(o->m_thickness(fd.id, CollisionThicknessBase::FACE)/2);
            o->m_f0bvt->update(fd.m_f0leaf, newvol);
        }
        o->m_f0bvt->optimize();
        o->m_thickness.setupIndex(CTB::EDGE, CTB::EDGE);
        o->m_thickness.setupIndex(CTB::VERTEX, CTB::FACE);
    }
    //updateCCD();
}

void BridsonCollisionManager::setAsInitialNonCollidingState(){
    for (auto& ob: m_dobjs)
        for (auto v: ob->m_obj->m_mesh.vertices())
            ob->m_x[v] = ob->m_obj->m_x[v];   
}

void BridsonCollisionManager::applyDx(){
    for (auto& ob: m_dobjs)
        for (auto v: ob->m_obj->m_mesh.vertices()){
            ob->m_x[v] = ob->m_next_x[v] = ob->m_x[v] + ob->m_obj->withBCmask(v, ob->m_dx[v]);
        }
}

void BridsonCollisionManager::redistributeImpulseV(Vector impulse, V_ind v, DynamicObject3DCI *obj){
    obj->m_dx[v] += impulse / obj->m_m[v];
}
void BridsonCollisionManager::redistributeImpulseE(Vector impulse, V_ind* v/*[2]*/, double* w/*[2]*/, DynamicObject3DCI *obj){
    for (int k = 0; k < 2; ++k)
        obj->m_dx[v[k]] += impulse * w[k] / obj->m_m[v[k]];
}
void BridsonCollisionManager::redistributeImpulseT(Vector impulse, V_ind* v/*[3]*/, double* w/*[3]*/, DynamicObject3DCI *obj){
    for (int k = 0; k < 3; ++k)
        obj->m_dx[v[k]] += impulse * w[k] / obj->m_m[v[k]];
}
double BridsonCollisionManager::computedx_nVF(CollisionVF& col){
    return (col.pobj->m_dx[col.vp] - GeomProjs::BaryEval(col.tobj->m_dx[col.vt[0]], col.tobj->m_dx[col.vt[1]], col.tobj->m_dx[col.vt[2]], col.wt)) * col.n;
}
double BridsonCollisionManager::computedx_nEE(CollisionEE& col){
    return -( (1-col.w[0])*col.obj[0]->m_dx[col.v[0][0]] + col.w[0]*col.obj[0]->m_dx[col.v[0][1]] -
                (1-col.w[1])*col.obj[1]->m_dx[col.v[1][0]] - col.w[1]*col.obj[1]->m_dx[col.v[1][1]]) * col.n;
}
double BridsonCollisionManager::applyInelasticRepulsionVF(CollisionVF& col){
    auto& n = col.n;
    double dxn = computedx_nVF(col);
    if (dxn <= 0) return 0;
    Vector I = dxn / (1 / col.pobj->m_m[col.vp] + col.wt[0]*col.wt[0]/col.tobj->m_m[col.vt[0]] + col.wt[1]*col.wt[1]/col.tobj->m_m[col.vt[1]] + col.wt[2]*col.wt[2]/col.tobj->m_m[col.vt[2]]) * n;
    redistributeImpulseT(I, col.vt.data(), col.wt.data(), col.tobj);
    redistributeImpulseV(-I, col.vp, col.pobj);
    return dxn;
}
double BridsonCollisionManager::applyInelasticRepulsionEE(CollisionEE& col){
    auto& n = col.n;
    double dxn = computedx_nEE(col);
    if (dxn <= 0) return 0;
    auto sq = [](auto t) {return t*t; };
    Vector I = dxn / (sq(1-col.w[0])/col.obj[0]->m_m[col.v[0][0]] + sq(col.w[0])/col.obj[0]->m_m[col.v[0][1]] + sq(1-col.w[1])/col.obj[1]->m_m[col.v[1][0]] + sq(col.w[1])/col.obj[1]->m_m[col.v[1][1]]) * n;
    std::array<double, 2> ww{1-col.w[0], col.w[0]};
    redistributeImpulseE(I, col.v[0].data(), ww.data(), col.obj[0]);
    ww[0] = 1-col.w[1], ww[1] = col.w[1];
    redistributeImpulseE(-I, col.v[1].data(), ww.data(), col.obj[1]);
    return dxn;
}
double BridsonCollisionManager::applySpringRepulsionVF(CollisionVF& col){
    double eps = col.sfp.eps;
    double k_spr = col.sfp.k_spr; //TODO: k = P*h^2/H
    double d = col.sfp.h - (col.prj - col.pobj->m_x[col.vp]) * col.n;
    double dxn = computedx_nVF(col);
    if (dxn > -eps*d) {
        double j = eps * d + dxn;
        double spr = k_spr * d * m_dt * m_dt;
        double Ir = std::min(j, spr);
        redistributeImpulseT(Ir * col.n, col.vt.data(), col.wt.data(), col.tobj);
        redistributeImpulseV(-Ir * col.n, col.vp, col.pobj);
        return  dxn - computedx_nVF(col);
    }
    return 0;
}
void BridsonCollisionManager::applyFrictionRepulsionVF(CollisionVF& col, double ddxn){
    if (col.sfp.mu <= 0) return;
    Vector dxr = col.pobj->m_dx[col.vp] - GeomProjs::BaryEval(col.tobj->m_dx[col.vt[0]], col.tobj->m_dx[col.vt[1]], col.tobj->m_dx[col.vt[2]], col.wt);
    Vector dxt = dxr - (dxr * col.n) * col.n;
    double dxt_sqr_len = dxt.squared_length();
    if (dxt_sqr_len < 1e-20) return;
    double vpre = sqrt(dxt.squared_length());
    Vector t = dxt / vpre;
    double I = (vpre - std::max(vpre - abs(ddxn)*col.sfp.mu, 0.0)) / (1 / col.pobj->m_m[col.vp] + col.wt[0]*col.wt[0]/col.tobj->m_m[col.vt[0]] + col.wt[1]*col.wt[1]/col.tobj->m_m[col.vt[1]] + col.wt[2]*col.wt[2]/col.tobj->m_m[col.vt[2]]);
    redistributeImpulseT(I*t, col.vt.data(), col.wt.data(), col.tobj);
    redistributeImpulseV(-I*t, col.vp, col.pobj);
}
void BridsonCollisionManager::applyFrictionRepulsionEE(CollisionEE& col, double ddxn){
    if (col.sfp.mu <= 0 || ddxn == 0.0) return;
    Vector dxr = -( (1-col.w[0])*col.obj[0]->m_dx[col.v[0][0]] + col.w[0]*col.obj[0]->m_dx[col.v[0][1]] -
                    (1-col.w[1])*col.obj[1]->m_dx[col.v[1][0]] - col.w[1]*col.obj[1]->m_dx[col.v[1][1]]);
    Vector dxt = dxr - (dxr * col.n) * col.n;
    double dxt_sqr_len = dxt.squared_length();
    if (dxt_sqr_len < 1e-20) return;
    double vpre = sqrt(dxt.squared_length());
    Vector t = dxt / vpre;
    auto sq = [](auto t) {return t*t; };
    double I = (vpre - std::max(vpre - abs(ddxn)*col.sfp.mu, 0.0)) / (sq(1-col.w[0])/col.obj[0]->m_m[col.v[0][0]] + sq(col.w[0])/col.obj[0]->m_m[col.v[0][1]] + sq(1-col.w[1])/col.obj[1]->m_m[col.v[1][0]] + sq(col.w[1])/col.obj[1]->m_m[col.v[1][1]]);
    std::array<double, 2> ww{1-col.w[0], col.w[0]};
    redistributeImpulseE(I*t, col.v[0].data(), ww.data(), col.obj[0]);
    ww[0] = 1-col.w[1], ww[1] = col.w[1];
    redistributeImpulseE(-I*t, col.v[1].data(), ww.data(), col.obj[1]);
}
double BridsonCollisionManager::applySpringRepulsionEE(CollisionEE& col){
    double eps = col.sfp.eps;
    double k_spr = col.sfp.k_spr; //TODO: k = P*h^2/H
    double d = col.sfp.h - (col.prj[0] - col.prj[1]) * col.n;
    double dxn = computedx_nEE(col);
    if (dxn > -eps*d) {
        double j = eps * d + dxn;
        double spr = k_spr * d * m_dt * m_dt;
        double Ir = std::min(j, spr);
        std::array<double, 2> ww{1 - col.w[0], col.w[0]};
        redistributeImpulseE(Ir * col.n, col.v[0].data(), ww.data(), col.obj[0]);
        ww[0] = 1 - col.w[1], ww[1] = col.w[1];
        redistributeImpulseE(-Ir * col.n, col.v[1].data(), ww.data(), col.obj[1]);
        return dxn - computedx_nEE(col);
    }
    return 0;
}
void BridsonCollisionManager::applyRepulsionVF(CollisionVF& col){
    double ddxn = applyInelasticRepulsionVF(col);
    ddxn += applySpringRepulsionVF(col);
    applyFrictionRepulsionVF(col, ddxn);
}
void BridsonCollisionManager::applyRepulsionEE(CollisionEE& col){
    double ddxn = applyInelasticRepulsionEE(col);
    ddxn += applySpringRepulsionEE(col);
    applyFrictionRepulsionEE(col, ddxn);
}

void BridsonCollisionManager::applyRepulsionImpulses(){
    struct NarrowPhase{
        Point prj0, prj1;
        double sqd = 0;
        int nPrim0, nPrim1;
        DynamicObject3DCI::FaceData* fd[2];
        std::array<V_ind, 3> vv[2];
        double w0, w1;
        double thickness;

        bool InitVF(DynamicObject3DCI::FaceData** _fd, std::array<V_ind, 3>* _vv, int nTri, int nPonTri){
            fd[0] = _fd[0], fd[1] = _fd[1];
            vv[0] = _vv[0], vv[1] = _vv[1];
            nPrim0 = nTri, nPrim1 = nPonTri;
            auto vp = vv[nTri][nPonTri];
            auto ft = fd[(nTri+1)%2]->id;
            auto& vt = vv[(nTri+1)%2];
            auto& t_x = fd[(nTri+1)%2]->backptr->m_x, &p_x = fd[nTri]->backptr->m_x;
            prj1 = p_x[vp];
            thickness = 0;
            if (fd[0]->backptr->m_obj != fd[1]->backptr->m_obj)
                thickness = (fd[nTri]->backptr->m_thickness(vp.idx(), CTB::VERTEX) + fd[(nTri+1)%2]->backptr->m_thickness(ft.idx(), CTB::FACE)) / 2;
            else 
                thickness = fd[0]->backptr->m_thickness(vp.idx(), ft.idx(), CTB::VERTEX, CTB::FACE);
            double max_sqr_dist = thickness * thickness;
            bool is_close = GeomProjs::ProjectOnTriangle(t_x[vt[0]], t_x[vt[1]], t_x[vt[2]], prj1, prj0, sqd, &max_sqr_dist);
            return is_close;
        }
        bool InitEE(DynamicObject3DCI::FaceData** _fd, std::array<V_ind, 3>* _vv, int nE0, int nE1){
            fd[0] = _fd[0], fd[1] = _fd[1];
            vv[0] = _vv[0], vv[1] = _vv[1];
            nPrim0 = nE0, nPrim1 = nE1;
            auto& e0_x = fd[0]->backptr->m_x, &e1_x = fd[1]->backptr->m_x;
            E_ind e0 = edge_around(fd[0]->backptr->m_obj->m_mesh, vv[0][nE0], vv[0][(nE0+1)%3], fd[0]->id),
                  e1 = edge_around(fd[1]->backptr->m_obj->m_mesh, vv[1][nE1], vv[1][(nE1+1)%3], fd[1]->id);  
            
            if (fd[0]->backptr->m_obj != fd[1]->backptr->m_obj)
                thickness = (fd[0]->backptr->m_thickness(e0.idx(), CTB::EDGE) + fd[1]->backptr->m_thickness(e1.idx(), CTB::EDGE)) / 2;
            else
                thickness = fd[0]->backptr->m_thickness(e0.idx(), e1.idx(), CTB::EDGE, CTB::EDGE);
            double max_sqr_dist = thickness * thickness;
            w0 = -1, w1 = -1;
            Vector vprj0, vprj1;
            bool is_close = GeomProjs::SegmentsProjects(e0_x[vv[0][nE0]] - CGAL::ORIGIN, e0_x[vv[0][(nE0+1)%3]] - CGAL::ORIGIN,
                                                        e1_x[vv[1][nE1]] - CGAL::ORIGIN, e1_x[vv[1][(nE1+1)%3]] - CGAL::ORIGIN,
                                                        vprj0, vprj1, w0, w1, sqd, &max_sqr_dist);
            prj0 = CGAL::ORIGIN + vprj0, prj1 = CGAL::ORIGIN + vprj1;
            return is_close;
        }
    };
    struct Handler{
        NarrowPhase* np = nullptr;
        BridsonCollisionManager* self = nullptr;
        void InitVF(NarrowPhase& _np, BridsonCollisionManager* _ptr){
            np = &_np;
            self = _ptr;
        }
        void InitEE(NarrowPhase& _np, BridsonCollisionManager* _ptr){
            np = &_np;
            self = _ptr;
        }
        bool applyVF(){
            CollisionVF col;
            col.tobj = np->fd[(np->nPrim0+1)%2]->backptr; col.pobj = np->fd[np->nPrim0]->backptr;
            col.ft = np->fd[(np->nPrim0+1)%2]->id; col.vt = np->vv[(np->nPrim0+1)%2]; col.vp = np->vv[np->nPrim0][np->nPrim1];
            col.prj = np->prj0; col.sqd = np->sqd;
            // auto ww = GeomProjs::BaryCoord(col.tobj->m_x[col.vt[0]], col.tobj->m_x[col.vt[1]], col.tobj->m_x[col.vt[2]], col.pobj->m_x[col.vp]);
            auto ww = GeomProjs::BaryCoord(col.tobj->m_x[col.vt[0]], col.tobj->m_x[col.vt[1]], col.tobj->m_x[col.vt[2]], col.prj);
            col.wt = std::array<double, 3>{ww[0], ww[1], ww[2]};
            col.n = (col.prj - col.pobj->m_x[col.vp]) / sqrt(col.sqd);
            col.sfp.h = np->thickness; 
            col.sfp._init_data_from_BridsonCollisionManager(self, np->fd[0]->backptr->m_id, np->fd[1]->backptr->m_id);
            ///handle point vs triangle self-intersection
            self->applyRepulsionVF(col);
            return true;
        }
        bool applyEE(){
            CollisionEE col;
            col.obj[0] = np->fd[0]->backptr; col.obj[1] = np->fd[1]->backptr;
            col.v[0] = std::array<V_ind, 2>{np->vv[0][np->nPrim0], np->vv[0][(np->nPrim0+1)%3]}, col.v[1] = std::array<V_ind, 2>{np->vv[1][np->nPrim1], np->vv[1][(np->nPrim1+1)%3]};
            col.prj[0] = np->prj0, col.prj[1] = np->prj1;
            col.w[0] = np->w0, col.w[1] = np->w1;
            col.sqd = np->sqd; col.n = (col.prj[0] - col.prj[1]) / sqrt(np->sqd);
            col.sfp.h = np->thickness;
            col.sfp._init_data_from_BridsonCollisionManager(self, np->fd[0]->backptr->m_id, np->fd[1]->backptr->m_id);
            ///handle edge vs edge self-intersection
            self->applyRepulsionEE(col);
            return true;
        }
    };
    //self-collision
    if (m_collision_status & SelfCollision){
        for (int i = 0; i < m_dobjs.size(); ++i){
            auto &obw = m_dobjs[i];
            _performBroadPhaseSelf<NarrowPhase, Handler>(*(obw->m_f0bvt));
        }
    }
    //different objects
    if (m_collision_status & DynamicVsDynamicCollission){
        for (int i = 1; i < m_dobjs.size(); ++i)
        for (int j = 0; j < i; ++j){
            auto &obw0 = m_dobjs[j], &obw1 = m_dobjs[i];
            _performBroadPhaseDiff<NarrowPhase, Handler>(*(obw0->m_f0bvt), *(obw1->m_f0bvt));
        }
    }
}
int BridsonCollisionManager::applyCCDRepulsion(){
    CHECKCOLON();
    struct NarrowPhaseCCD{
        int nPrim0, nPrim1;
        DynamicObject3DCI::FaceData* fd[2];
        std::array<V_ind, 3> vv[2];
        double m_t;
        bool InitVF(DynamicObject3DCI::FaceData** _fd, std::array<V_ind, 3>* _vv, int nTri, int nPonTri){
            fd[0] = _fd[0], fd[1] = _fd[1];
            vv[0] = _vv[0], vv[1] = _vv[1];
            nPrim0 = nTri, nPrim1 = nPonTri;
            auto ft = fd[(nTri+1)%2]->id;
            auto vp = vv[nTri][nPonTri];
            auto& vt = vv[(nTri+1)%2];
            auto& t_x = fd[(nTri+1)%2]->backptr->m_x, &p_x = fd[nTri]->backptr->m_x;
            auto& t_dx = fd[(nTri+1)%2]->backptr->m_dx, &p_dx = fd[nTri]->backptr->m_dx;
            Point o = t_x[vt[0]];
            bool has_intersection = GeomProjs::PerformVFCCD(Vector(CGAL::NULL_VECTOR), t_x[vt[1]] - o, t_x[vt[2]] - o, p_x[vp] - o,
                                                            t_dx[vt[0]], t_x[vt[1]] - o + t_dx[vt[1]], t_x[vt[2]] - o + t_dx[vt[2]], p_x[vp] - o + p_dx[vp], &m_t);
            return has_intersection;
        }
        bool InitEE(DynamicObject3DCI::FaceData** _fd, std::array<V_ind, 3>* _vv, int nE0, int nE1){
            fd[0] = _fd[0], fd[1] = _fd[1];
            vv[0] = _vv[0], vv[1] = _vv[1];
            nPrim0 = nE0, nPrim1 = nE1;
            auto& e0_x = fd[0]->backptr->m_x, &e1_x = fd[1]->backptr->m_x;
            auto& e0_dx = fd[0]->backptr->m_dx, &e1_dx = fd[1]->backptr->m_dx;
            Point o = e0_x[vv[0][nE0]];
            bool has_intersection = GeomProjs::PerformEECCD(Vector(CGAL::NULL_VECTOR), e0_x[vv[0][(nE0+1)%3]] - o, e1_x[vv[1][nE1]] - o, e1_x[vv[1][(nE1+1)%3]] - o,
                                                            e0_dx[vv[0][nE0]], e0_x[vv[0][(nE0+1)%3]] - o + e0_dx[vv[0][(nE0+1)%3]], e1_x[vv[1][nE1]] - o + e1_dx[vv[1][nE1]], e1_x[vv[1][(nE1+1)%3]] - o + e1_dx[vv[1][(nE1+1)%3]], &m_t);
            return has_intersection;
        }
    };
    struct HandlerCCD {
        NarrowPhaseCCD *np = nullptr;
        BridsonCollisionManager *self = nullptr;

        void InitVF(NarrowPhaseCCD &_np, BridsonCollisionManager *_ptr) {
            np = &_np;
            self = _ptr;
        }

        void InitEE(NarrowPhaseCCD &_np, BridsonCollisionManager *_ptr) {
            np = &_np;
            self = _ptr;
        }
        bool applyVF(){
            auto vp = np->vv[np->nPrim0][np->nPrim1];
            auto& vt = np->vv[(np->nPrim0+1)%2];
            auto& t_x = np->fd[(np->nPrim0+1)%2]->backptr->m_x, &p_x = np->fd[np->nPrim0]->backptr->m_x;
            auto& t_dx = np->fd[(np->nPrim0+1)%2]->backptr->m_dx, &p_dx = np->fd[np->nPrim0]->backptr->m_dx;

            Point pc = p_x[vp] + np->m_t * p_dx[vp];
            Point ac = t_x[vt[0]] + np->m_t * t_dx[vt[0]], bc = t_x[vt[1]] + np->m_t * t_dx[vt[1]], cc = t_x[vt[2]] + np->m_t * t_dx[vt[2]];
            auto ww = GeomProjs::BaryCoord(ac, bc, cc, pc);
            {//TODO: remove this
                auto sqdt = (GeomProjs::BaryEval(ac, bc, cc, ww) - pc).squared_length();
                if (sqdt > 1e-6) {
                    std::cout << "Warning: Wrong VF intersect decision sqd = "<< sqdt << std::endl;
                    std::cout << "\tT: \n" << std::hexfloat  << 
                                 "\t\tv0: " << t_x[vt[0]] << " -> " << t_dx[vt[0]] << "\n"
                                 "\t\tv1: " << t_x[vt[1]] << " -> " << t_dx[vt[1]] << "\n"
                                 "\t\tv2: " << t_x[vt[2]] << " -> " << t_dx[vt[2]] << "\n"
                                 "\tP: \n\t\tv3: " << p_x[vp] << " -> " << p_dx[vp] << "\n";
                    auto pc_proj = GeomProjs::BaryEval(ac, bc, cc, ww);
                    std::cout << "\tt_eval = " << std::defaultfloat << std::setprecision(16) << np->m_t << " where sqd = " << sqdt << "\n"
                                << "\tT (at t_eval): \n"
                                << "\t\tv0: " << ac << "\n"
                                << "\t\tv1: " << bc << "\n"
                                << "\t\tv2: " << cc << "\n"
                                << "\tP (at t_eval): \n\t\tv3: " << pc << "\n"
                                << "\tP projection (at t_eval): \n\t\t    " << pc_proj << "\n";
                    std::cout << "\tIn plane test: " << CGAL::cross_product(bc - ac, cc - ac) * (pc - ac) << "\n";
                    std::cout << "\tIn triangle test: " 
                        << CGAL::cross_product(ac - pc, bc - pc) * CGAL::cross_product(ac - cc, bc - cc) << " " 
                        << CGAL::cross_product(bc - pc, cc - pc) * CGAL::cross_product(bc - ac, cc - ac) << " " 
                        << CGAL::cross_product(cc - pc, ac - pc) * CGAL::cross_product(cc - bc, ac - bc) << "\n\n"; 
                    bool tmp = np->InitVF(np->fd, np->vv, np->nPrim0, np->nPrim1);
                    std::cout << "tmp = " << tmp << std::endl;
                    return false;
                }
            }
            auto w = std::array<double, 3>{ww[0], ww[1], ww[2]};
            auto sqd0 = (GeomProjs::BaryEval(t_x[vt[0]], t_x[vt[1]], t_x[vt[2]], w) - p_x[vp]).squared_length();
            auto d0 = sqrt(sqd0);
            double h = 0;
            if (np->fd[0]->backptr->m_obj != np->fd[1]->backptr->m_obj) 
                h = (np->fd[np->nPrim0]->backptr->m_thickness(vp, CTB::VERTEX) + np->fd[(np->nPrim0+1)%2]->backptr->m_thickness(np->fd[(np->nPrim0+1)%2]->id, CTB::FACE)) / 2;
            else
                h = np->fd[0]->backptr->m_thickness(vp, np->fd[(np->nPrim0+1)%2]->id, CTB::VERTEX, CTB::FACE);
            ObjectID oid[2] = {np->fd[0]->backptr->m_id, np->fd[1]->backptr->m_id}; 
            double dmin = h * (1 - self->m_ccdPenetrateDepthCoef(oid[0], oid[1]));
            if (d0 > dmin){
                //make proxy positions after some free movement
                double alpha = (d0 - dmin) / d0;
                std::array<Point, 4> tmp_x = {t_x[vt[0]], t_x[vt[1]], t_x[vt[2]], p_x[vp]};
                std::array<Vector, 4> tmp_dx = {alpha * np->m_t*t_dx[vt[0]], alpha * np->m_t*t_dx[vt[1]], alpha * np->m_t*t_dx[vt[2]], alpha * np->m_t*p_dx[vp]};
                p_x[vp] += alpha * np->m_t * p_dx[vp]; p_dx[vp] *= 1 - alpha * np->m_t;
                for (int i = 0; i < 3; ++i) {
                    t_x[vt[i]] += alpha * np->m_t * t_dx[vt[i]];
                    t_dx[vt[i]] *= 1 - alpha * np->m_t;
                }

                CollisionVF col;
                col.tobj = np->fd[(np->nPrim0+1)%2]->backptr; col.pobj = np->fd[np->nPrim0]->backptr;
                col.ft = np->fd[(np->nPrim0+1)%2]->id; col.vt = np->vv[(np->nPrim0+1)%2]; col.vp = np->vv[np->nPrim0][np->nPrim1];
                col.prj = GeomProjs::BaryEval(t_x[vt[0]] , t_x[vt[1]], t_x[vt[2]], w); col.sqd = dmin*dmin;
                col.wt = w;
                col.n = (col.prj - col.pobj->m_x[col.vp]) / dmin;
                col.sfp.h = h;
                col.sfp._init_data_from_BridsonCollisionManager(self, oid[0], oid[1]);
                self->applyInelasticRepulsionVF(col);
                //return to initial positions before free movement
                p_x[vp] = tmp_x[3]; p_dx[vp] += tmp_dx[3];
                for (int i = 0; i < 3; ++i) {
                    t_x[vt[i]] = tmp_x[i];
                    t_dx[vt[i]] += tmp_dx[i];
                }
            } else {
                CollisionVF col;
                col.tobj = np->fd[(np->nPrim0+1)%2]->backptr; col.pobj = np->fd[np->nPrim0]->backptr;
                col.ft = np->fd[(np->nPrim0+1)%2]->id; col.vt = np->vv[(np->nPrim0+1)%2]; col.vp = np->vv[np->nPrim0][np->nPrim1];
                col.prj = GeomProjs::BaryEval(t_x[vt[0]] , t_x[vt[1]], t_x[vt[2]], w); col.sqd = sqd0;
                col.wt = w;
                col.n = (col.prj - col.pobj->m_x[col.vp]) / d0;
                col.sfp.h = h;
                col.sfp._init_data_from_BridsonCollisionManager(self, np->fd[0]->backptr->m_id, np->fd[1]->backptr->m_id);
                self->applyInelasticRepulsionVF(col);
                self->applySpringRepulsionVF(col);
            }
            return true;
        }
        bool applyEE(){
            auto& e0_x = np->fd[0]->backptr->m_x, &e1_x = np->fd[1]->backptr->m_x;
            auto& e0_dx = np->fd[0]->backptr->m_dx, &e1_dx = np->fd[1]->backptr->m_dx;

            std::array<std::array<V_ind, 2>, 2> lvv{std::array<V_ind, 2>{np->vv[0][np->nPrim0], np->vv[0][(np->nPrim0+1)%3]},
                                                    std::array<V_ind, 2>{np->vv[1][np->nPrim1], np->vv[1][(np->nPrim1+1)%3]}};
            E_ind e0 = edge_around(np->fd[0]->backptr->m_obj->m_mesh, lvv[0][0], lvv[0][1], np->fd[0]->id),
                  e1 = edge_around(np->fd[1]->backptr->m_obj->m_mesh, lvv[1][0], lvv[1][1], np->fd[1]->id);
            double h = 0;
            if (np->fd[0]->backptr->m_obj != np->fd[1]->backptr->m_obj)
                h = (np->fd[0]->backptr->m_thickness(e0.idx(), CTB::EDGE) + np->fd[1]->backptr->m_thickness(e1.idx(), CTB::EDGE)) / 2;
            else
                h = np->fd[0]->backptr->m_thickness(e0.idx(), e1.idx(), CTB::EDGE, CTB::EDGE);
            ObjectID oid[2] = {np->fd[0]->backptr->m_id, np->fd[1]->backptr->m_id}; 
            double dmin = h * (1 - self->m_ccdPenetrateDepthCoef(oid[0], oid[1]));

            Point a0t = e0_x[lvv[0][0]] + np->m_t * e0_dx[lvv[0][0]], a1t = e0_x[lvv[0][1]] + np->m_t * e0_dx[lvv[0][1]];
            Point b0t = e1_x[lvv[1][0]] + np->m_t * e1_dx[lvv[1][0]], b1t = e1_x[lvv[1][1]] + np->m_t * e1_dx[lvv[1][1]];
            Vector prj0, prj1; double w[2]; double sqd; Point o = CGAL::ORIGIN;
            GeomProjs::SegmentsProjects(a0t-o, a1t-o, b0t-o, b1t-o, prj0, prj1, w[0], w[1], sqd);
            
            {//TODO: remove this
                if (sqd > 1e-6) {
                    std::cout << "Warning: Wrong EE intersect decision sqd = "<< sqd << std::endl;
                    std::cout << "\tE0: \n" << std::hexfloat <<
                                 "\t\tv0: " << e0_x[lvv[0][0]] << " -> " << e0_dx[lvv[0][0]] << "\n"
                                 "\t\tv1: " << e0_x[lvv[0][1]] << " -> " << e0_dx[lvv[0][1]] << "\n"
                                 "\tE1: \n" << 
                                 "\t\tv0: " << e1_x[lvv[1][0]] << " -> " << e1_dx[lvv[1][0]] << "\n"
                                 "\t\tv1: " << e1_x[lvv[1][1]] << " -> " << e1_dx[lvv[1][1]] << "\n";
                    std::cout << "\tt_eval = " << std::defaultfloat << std::setprecision(16) <<  np->m_t << " where sqd = " << sqd << "\n"
                                << "\tE0 (at t_eval): \n"
                                << "\t\tv0:  " << a0t << "\n"
                                << "\t\tv1:  " << a1t << "\n"
                                << "\t\tprj: " << prj0 << "\n"
                                << "\tE1 (at t_eval): \n"
                                << "\t\tv0: " << b0t << "\n"
                                << "\t\tv1: " << b1t << "\n"
                                << "\t\tprj: " << prj1 << "\n";
                    std::cout << "\tIn plane test: " << CGAL::cross_product(a0t - b0t, a1t - b0t) * (b1t - b0t) << "\n";
                    std::cout << "\tSign test: " 
                        << -CGAL::cross_product(a0t - b1t, a1t - b1t) * CGAL::cross_product(a0t - b0t, a1t - b0t) << " " 
                        << CGAL::cross_product(a1t - b1t, b0t - b1t) * CGAL::cross_product(a1t - a0t, b0t - a0t) << " " 
                        << CGAL::cross_product(b0t - b1t, a0t - b1t) * CGAL::cross_product(b0t - a1t, a0t - a1t) << "\n\n"; 
                    bool tmp = np->InitEE(np->fd, np->vv, np->nPrim0, np->nPrim1);
                    std::cout << "tmp = " << tmp << std::endl;
                    return false;
                }
            }
            auto sqd0 = (GeomProjs::Lerp(e0_x[lvv[0][0]], e0_x[lvv[0][1]], w[0]) - GeomProjs::Lerp(e1_x[lvv[1][0]], e1_x[lvv[1][1]], w[1])).squared_length();
            auto d0 = sqrt(sqd0);
            if (d0 > dmin){
                //make proxy positions after some free movement
                double alpha = (d0 - dmin) / d0;
                std::array<Point, 4> tmp_x = {e0_x[lvv[0][0]], e0_x[lvv[0][1]], e1_x[lvv[1][0]], e1_x[lvv[1][1]]};
                std::array<Vector, 4> tmp_dx = {e0_dx[lvv[0][0]], e0_dx[lvv[0][1]], e1_dx[lvv[1][0]], e1_dx[lvv[1][1]]};
                for (auto& dx: tmp_dx) dx *= alpha * np->m_t;
                for (int i = 0; i < 2; ++i){
                    e0_x[lvv[0][i]] += tmp_dx[0 + i]; e0_dx[lvv[0][i]] *= 1 - alpha * np->m_t;
                    e1_x[lvv[1][i]] += tmp_dx[2 + i]; e1_dx[lvv[1][i]] *= 1 - alpha * np->m_t;
                }
                CollisionEE col;
                col.obj[0] = np->fd[0]->backptr; col.obj[1] = np->fd[1]->backptr;
                col.v[0] = lvv[0], col.v[1] = lvv[1];
                col.prj[0] = GeomProjs::Lerp(e0_x[lvv[0][0]], e0_x[lvv[0][1]], w[0]), col.prj[1] = GeomProjs::Lerp(e1_x[lvv[1][0]], e1_x[lvv[1][1]], w[1]);
                col.w[0] = w[0], col.w[1] = w[1];
                col.sqd = dmin*dmin; col.n = (col.prj[0] - col.prj[1]) / dmin;
                col.sfp.h = h;
                col.sfp._init_data_from_BridsonCollisionManager(self, oid[0], oid[1]);
                self->applyInelasticRepulsionEE(col);
                //return to initial positions before free movement
                for (int i = 0; i < 2; ++i){
                    e0_x[lvv[0][i]] = tmp_x[0 + i]; e0_dx[lvv[0][i]] += tmp_dx[0 + i];
                    e1_x[lvv[1][i]] = tmp_x[2 + i]; e1_dx[lvv[1][i]] += tmp_dx[2 + i];
                }
            } else {
                CollisionEE col;
                col.obj[0] = np->fd[0]->backptr; col.obj[1] = np->fd[1]->backptr;
                col.v[0] = lvv[0], col.v[1] = lvv[1];
                col.prj[0] = GeomProjs::Lerp(e0_x[lvv[0][0]], e0_x[lvv[0][1]], w[0]), col.prj[1] = GeomProjs::Lerp(e1_x[lvv[1][0]], e1_x[lvv[1][1]], w[1]);
                col.w[0] = w[0], col.w[1] = w[1];
                col.sqd = sqd0; col.n = (col.prj[0] - col.prj[1]) / d0;
                col.sfp.h = h;
                col.sfp._init_data_from_BridsonCollisionManager(self, oid[0], oid[1]);
                self->applyInelasticRepulsionEE(col);
                self->applySpringRepulsionEE(col);
            }
            return true;
        }
    };
    int ncols = 0;
    //self-collision
    struct TempOpts{
        static DReal GetColTime(NarrowPhaseCCD& h){ return h.m_t; }
    };
    if (m_collision_status & SelfCollision){
        for (int i = 0; i < m_dobjs.size(); ++i){
            auto &obw = m_dobjs[i];
            ncols += _performBroadPhaseSelf<NarrowPhaseCCD, HandlerCCD, TempOpts>(*(obw->m_febvt));
        }
    }
    //different objects
    if (m_collision_status & DynamicVsDynamicCollission){
        for (int i = 1; i < m_dobjs.size(); ++i)
            for (int j = 0; j < i; ++j){
                auto &obw0 = m_dobjs[j], &obw1 = m_dobjs[i];
                ncols += _performBroadPhaseDiff<NarrowPhaseCCD, HandlerCCD, TempOpts>(*(obw0->m_febvt), *(obw1->m_febvt));
            }
    }
    CHECKCOLOFF();
    return ncols;
}

bool BridsonCollisionManager::applyRelaxRepulsion(int maxRelaxIts){
    for (int it = 0; it < maxRelaxIts; ++it) {
        int nrelax = 0;
        for (auto& ob: m_dobjs) {
            auto obj = ob->m_obj;
            auto reps = m_RelaxEps(ob->m_id);
            for (auto e: obj->m_mesh.edges()) {
                auto v = vert_around(obj->m_mesh, e);
                auto l0 = ob->m_x[v[1]] - ob->m_x[v[0]];
                auto dl = ob->m_dx[v[1]] - ob->m_dx[v[0]];
                auto l1 = l0 + dl;
                DReal l0_2 = l0.squared_length(), l1_2 = l1.squared_length();
                bool small = (l1_2 - l0_2 < (reps - 2 - 1e-6) * reps * l0_2);
                bool big = (l1_2 - l0_2 > (reps + 2 + 1e-6) * reps * l0_2);
                if (!small && !big) continue;
                DReal dl_2 = dl.squared_length(), l0dl = dl * l0;
                Vector dxc = (ob->m_dx[v[1]] * ob->m_m[v[1]] + ob->m_dx[v[0]] * ob->m_m[v[0]]) /
                                (ob->m_m[v[1]] + ob->m_m[v[0]]);
                DReal k = l0dl / dl_2, q =
                        -(small ? (reps - 2) * reps : (reps + 2) * reps) * l0_2 / dl_2;
                DReal D = sqrt(k * k - q);
                DReal alpha = -k + D;
                if (alpha > 1) alpha = -k - D;
                ob->m_dx[v[0]] = dxc + (ob->m_dx[v[0]] - dxc) * alpha;
                ob->m_dx[v[1]] = dxc + (ob->m_dx[v[1]] - dxc) * alpha;
                nrelax++;
            }
        }
        if (nrelax <= 0) return true;
    }
    return false;
}

bool BridsonCollisionManager::applyCCDRepulsions(int maxCCDIts){
    for (int it = 0; it < maxCCDIts; ++it) {
        updateCCD();
        int ncol = applyCCDRepulsion();
        if (!ncol) return true;
    }
    return false;
}

void BridsonCollisionManager::applyFailSafe(int maxFailSafeIts){
    std::vector<ImpactZone> izs;
    int it = 0;
    for (; it < maxFailSafeIts; ++it){
        // std::cout << "FailSafe: it = " << it << "\n"; //TODO: remove this
        updateCCD();
        // std::cout << "\t updateCCD done\n"; //TODO: remove this
        int ncol = collectImpactZones(izs);
        // std::cout << "\t collectImpactZones done, ncol = " << ncol << "\n"; //TODO: remove this
        if (ncol <= 0) return;
        for (auto& iz: izs) {
            handleImpactZone(iz);
        }
        // std::cout << "\t handleImpactZone done" << std::endl; //TODO: remove this
    }
    if (it == maxFailSafeIts)
        std::cerr << "[COLLISION_ALGO WARNING] BridsonCollisionManager::applyFailSafe\n"
                    "\tFailSafe algorithm reached maximum fail-safe algo number of iterations("<<maxFailSafeIts << ")."
                    "It may lead to nonfree-intersection state. Try to set a larger value of maxFailSafeIts using set_MaxFailSafeIts(...) method\n";
}

static std::ostream& print_vec(std::vector<unsigned> v, std::string name = "") {
    static int n=0;
    if (name.empty()) name = std::string("vec_") + std::to_string(n);
    std::cout << name << "[" << v.size() << "]{ ";
    for (int i = 0; i < v.size() - 1; ++i)
        std::cout << v[i] << ", ";
    if (!v.empty()) std::cout << v.back();
    std::cout << " }";
    ++n;
    return std::cout;
}

static std::ostream& print_vec(std::vector<std::pair<unsigned, unsigned>> v, std::string name = "") {
    static int n=0;
    if (name.empty()) name = std::string("vec_") + std::to_string(n);
    std::cout << name << "[" << v.size() << "]{ ";
    for (int i = 0; i < v.size() - 1; ++i)
        std::cout << "(" << v[i].first << ", " << v[i].second << ")" << ", ";
    if (!v.empty()) std::cout << "(" << v.back().first << ", " << v.back().second << ")";
    std::cout << " }";
    ++n;
    return std::cout;
}

int BridsonCollisionManager::collectImpactZones(std::vector<ImpactZone>& izs){
    struct FullVInd{
        std::map<DynamicObject3DCI*, unsigned> obj_num;
        std::map<unsigned, DynamicObject3DCI*> obj_rev;
        std::vector<unsigned int> obj_id_shift;
        FullVInd(std::vector<std::shared_ptr<DynamicObject3DCI>>& objs){
            obj_id_shift.resize(objs.size()+1, 0);
            for (int i = 1; i < objs.size()+1; ++i) obj_id_shift[i] = obj_id_shift[i-1] + objs[i-1]->m_obj->m_mesh.num_vertices();
            for (int i = 0; i < objs.size(); ++i) {
                obj_num.insert({objs[i].get(), i});
                obj_rev.insert({i, objs[i].get()});
            }
        }
        unsigned Ind(unsigned obj_num, V_ind vid) { return obj_id_shift[obj_num] + vid.idx(); }
        unsigned Ind(DynamicObject3DCI* obj, V_ind vid) {
            auto it = obj_num.find(obj);
            if (it != obj_num.end()){
                return Ind(it->second, vid);
            } else {
                throw std::runtime_error("Can't find object");
            }
        }
        std::tuple<DynamicObject3DCI*, unsigned, V_ind> LocInd(unsigned gind){
            auto it = std::upper_bound(obj_id_shift.begin(), obj_id_shift.end(), gind);
            if (it == obj_id_shift.end()) throw std::runtime_error("wrong global index");
            int onum = std::distance(obj_id_shift.begin(), it) - 1;
            auto jt = obj_rev.find(onum);
            if (jt != obj_rev.end()) return {jt->second, jt->first, V_ind(gind - obj_id_shift[onum])};
            else throw std::runtime_error("wrong global index");
        }
    };
    FullVInd fvid(m_dobjs);
    std::vector<unsigned> impact_nodes; impact_nodes.reserve(fvid.obj_id_shift.back() / 20);
    std::vector<std::pair<unsigned, unsigned>> impact_links; impact_links.reserve(fvid.obj_id_shift.back() / 20 * 3);
    using CollsionCont = std::vector<std::pair<std::array<unsigned, 4>, int>>;
    CollsionCont collisions;
    struct SaveZonesStructure{
        FullVInd* pfvid;
        std::vector<std::pair<unsigned, unsigned>>* pimpact_links;
        std::vector<unsigned>* pimpact_nodes;
        CollsionCont* collisions;
        BridsonCollisionManager *self;

        bool InitVF(DynamicObject3DCI::FaceData** fd, std::array<V_ind, 3>* vv, int nTri, int nPonTri){
            auto ft = fd[(nTri+1)%2]->id;
            auto vp = vv[nTri][nPonTri];
            auto& vt = vv[(nTri+1)%2];
            auto& t_x = fd[(nTri+1)%2]->backptr->m_x, &p_x = fd[nTri]->backptr->m_x;
            auto& t_dx = fd[(nTri+1)%2]->backptr->m_dx, &p_dx = fd[nTri]->backptr->m_dx;
            CGAL::Origin o;
            bool has_intersection = GeomProjs::PerformVFCCD(t_x[vt[0]] - o, t_x[vt[1]] - o, t_x[vt[2]] - o, p_x[vp] - o,
                                                            t_x[vt[0]] + t_dx[vt[0]] - o, t_x[vt[1]] + t_dx[vt[1]] - o, t_x[vt[2]] + t_dx[vt[2]] - o, p_x[vp] + p_dx[vp] - o);
            if (has_intersection){
                std::array<unsigned, 4> ss;
                for (int i = 0; i < 3; ++i) ss[i] = pfvid->Ind(fd[(nTri+1)%2]->backptr, vt[i]);
                ss[3] = pfvid->Ind(fd[nTri]->backptr, vp);
                for (int i = 0; i < 4; ++i) pimpact_nodes->push_back(ss[i]);
                for (int i = 0; i < 4; ++i)
                for (int j = 0; j < 4; ++j)
                    if (i != j) pimpact_links->push_back({ss[i], ss[j]});
                collisions->emplace_back(std::move(ss), 1);
            }
            return has_intersection;
        }
        bool InitEE(DynamicObject3DCI::FaceData** fd, std::array<V_ind, 3>* vv, int nE0, int nE1){
            auto& e0_x = fd[0]->backptr->m_x, &e1_x = fd[1]->backptr->m_x;
            auto& e0_dx = fd[0]->backptr->m_dx, &e1_dx = fd[1]->backptr->m_dx;
            CGAL::Origin o;
            bool has_intersection = GeomProjs::PerformEECCD(e0_x[vv[0][nE0]] - o, e0_x[vv[0][(nE0+1)%3]] - o, e1_x[vv[1][nE1]] - o, e1_x[vv[1][(nE1+1)%3]] - o,
                                                            e0_x[vv[0][nE0]] + e0_dx[vv[0][nE0]] - o, e0_x[vv[0][(nE0+1)%3]] + e0_dx[vv[0][(nE0+1)%3]] - o, e1_x[vv[1][nE1]] + e1_dx[vv[1][nE1]] - o, e1_x[vv[1][(nE1+1)%3]] + e1_dx[vv[1][(nE1+1)%3]] - o);
            if (has_intersection){
                std::array<unsigned, 4> ss;
                ss[0] = pfvid->Ind(fd[0]->backptr, vv[0][nE0]);
                ss[1] = pfvid->Ind(fd[0]->backptr, vv[0][(nE0+1)%3]);
                ss[2] = pfvid->Ind(fd[1]->backptr, vv[1][nE1]);
                ss[3] = pfvid->Ind(fd[1]->backptr, vv[1][(nE1+1)%3]);
                for (int i = 0; i < 4; ++i) pimpact_nodes->push_back(ss[i]);
                for (int i = 0; i < 4; ++i)
                for (int j = 0; j < 4; ++j)
                    if (i != j) pimpact_links->push_back({ss[i], ss[j]});
                collisions->emplace_back(std::move(ss), 2);
            }
            return has_intersection;
        }
    };
    struct NoOp{
        void InitVF(SaveZonesStructure &_np, BridsonCollisionManager *_ptr){}
        void InitEE(SaveZonesStructure &_np, BridsonCollisionManager *_ptr){}
        bool applyVF(){ return true; }
        bool applyEE(){ return true; }
    };
    SaveZonesStructure szs;
    szs.pfvid = &fvid;
    szs.pimpact_nodes = &impact_nodes;
    szs.pimpact_links = &impact_links;
    szs.collisions = &collisions;
    szs.self = this;
    NoOp nope;
    int ncols = 0;
    //self-collision
    if (m_collision_status & SelfCollision){
        for (int i = 0; i < m_dobjs.size(); ++i){
            auto &obw = m_dobjs[i];
            ncols += _performBroadPhaseSelfExt(*(obw->m_febvt), szs, nope);
        }
    }
    //different objects
    if (m_collision_status & DynamicVsDynamicCollission){
        for (int i = 1; i < m_dobjs.size(); ++i)
            for (int j = 0; j < i; ++j){
                auto &obw0 = m_dobjs[j], &obw1 = m_dobjs[i];
                ncols += _performBroadPhaseDiffExt(*(obw0->m_febvt), *(obw1->m_febvt), szs, nope);
            }
    }
    if (ncols <= 0) return ncols;
    for (auto& iz: izs) iz.reset_dx();
    auto lunique = [](auto& v){ //remove duplicates
        auto first = v.begin();
        auto last = v.end();
        if (first == last)
            return;

        auto result = first;
        while (++first != last) {
            if (!(*result == *first) && ++result != first) {
                *result = std::move(*first);
            }
        }
        v.resize(std::distance(v.begin(), ++result));
    };
    {   //insert impact zones from previous computation
        unsigned n = impact_nodes.size();
        unsigned add_n = std::accumulate(izs.begin(), izs.end(), 0, [](const unsigned init, const auto& iz) { return init + iz.iz.size(); });
        impact_nodes.resize(n + add_n);
        for (auto& iz: izs) for (int i = 0; i < iz.iz.size(); ++i) impact_nodes[n++] = fvid.Ind(iz.iz[i].first, iz.iz[i].second);
    }
    {   //insert impact links and collisions from previous computation
        int n = impact_links.size();
        int add_n = std::accumulate(izs.begin(), izs.end(), 0, [](const unsigned init, const auto& iz) { return init + iz.collisions.size(); });
        impact_links.resize(n + 12*add_n);
        int nc = collisions.size();
        collisions.resize(nc + add_n);
        for (auto& iz: izs)
            for (int c = 0; c < iz.collisions.size(); ++c) {
                auto& v = iz.collisions[c].v;
                std::array<DynamicObject3DCI*, 4> objs{iz.collisions[c].obj[0], iz.collisions[c].obj[0], iz.collisions[c].obj[0], iz.collisions[c].obj[1]};
                if (iz.collisions[c].type == ImpactZone::ColType::EdgeEdge)
                    objs[2] = objs[3];
                std::array<unsigned, 4> ss{fvid.Ind(objs[0], v[0]), fvid.Ind(objs[1], v[1]), fvid.Ind(objs[2], v[2]), fvid.Ind(objs[3], v[3])};
                collisions[nc++] = std::pair<std::array<unsigned, 4>, int>{ss, iz.collisions[c].type == ImpactZone::ColType::EdgeEdge ? 2 : 1};
                for (int i = 0; i < 4; ++i)
                for (int j = 0; j < 4; ++j)
                    if (i != j) impact_links[n++] = std::pair<unsigned, unsigned>{ss[i], ss[j]};
            }
    }
    std::sort(impact_nodes.begin(), impact_nodes.end());
    // print_vec(impact_nodes, "impact_nodes_nouniq") << std::endl;
    lunique(impact_nodes);
    // print_vec(impact_nodes, "impact_nodes") << std::endl;
    std::sort(impact_links.begin(), impact_links.end());
    // print_vec(impact_links, "impact_links_init_nouniq") << std::endl;
    lunique(impact_links);
    // print_vec(impact_links, "impact_links_init") << std::endl;
    std::vector<unsigned> st_links(impact_nodes.size() + 1, 0);
    st_links[impact_nodes.size()] = impact_links.size();
    for (int i = 1; i < impact_nodes.size(); ++i)
    for (unsigned j = st_links[i-1]+1; st_links[i] <= st_links[i-1]; ++j)
        if (impact_nodes[i-1] != impact_links[j].first) st_links[i] = j;
    //make remap on sequential index
    for (int i = 0; i < impact_nodes.size(); ++i) {
        auto it = impact_nodes.begin();
        for (int j = st_links[i]; j < st_links[i + 1]; ++j) {
            auto &l = impact_links[j];
            l.first = i;
            it = std::lower_bound(it, impact_nodes.end(), l.second);
            l.second = std::distance(impact_nodes.begin(), it);
        }
    }
    std::vector<unsigned> iz_mark(impact_nodes.size(), 0);
    int niz = 0;
    std::vector<unsigned> niz_sizes;
    // print_vec(st_links, "st_links") << std::endl;
    // print_vec(impact_links, "impact_links") << std::endl;
    // std::cout << "st_links.size() = " << st_links.size() << std::endl;
    // std::cout << "iz_mark.size() = " << iz_mark.size() << std::endl;
    // std::cout << "impact_links.size() = " << impact_links.size() << std::endl;
    // std::cout << "Before DFS" << std::endl;
    {//dfs
        std::stack<unsigned> q;
        for (unsigned v = 0; v < impact_nodes.size(); ++v){
            if (iz_mark[v] > 0) continue;
            iz_mark[v] = ++niz;
            q.push(v);
            // std::cout << "Push v = " << v << " niz = " << niz << std::endl;
            niz_sizes.push_back(1);
            while (!q.empty()){
                unsigned vq = q.top();
                q.pop();
                // std::cout << "\t vq=" << vq << std::endl;
                for (auto i = st_links[vq]; i < st_links[vq+1]; ++i)
                {
                    // std::cout << "\t\t i = " << i;
                    auto n = impact_links[i].second;
                    if (iz_mark[n] == 0){
                        // std::cout << " n=" << n << std::endl;
                        iz_mark[n] = niz;
                        ++niz_sizes.back();
                        q.push(n);
                    } 
                    // else 
                        // std::cout << " iz_mark=" << iz_mark[n] << " n=" << n << std::endl;
                }
            }
            // std::cout << "q.size() = " << q.size() << std::endl;
        }
    }
    // print_vec(niz_sizes, "niz_sizes") << std::endl;
    // std::cout << "After DFS" << std::endl;
    izs.resize(niz);
    for (int i = 0; i < niz; ++i) izs[i].iz.reserve(niz_sizes[i]);
    for (int k = 0; k < impact_nodes.size(); ++k) {
        auto rev = fvid.LocInd(impact_nodes[k]);
        izs[iz_mark[k] - 1].iz.emplace_back(get<0>(rev), get<2>(rev));
    }
    // std::cout << "After IZS computation" << std::endl;
    for (auto & collision : collisions){
        auto n = std::lower_bound(impact_nodes.begin(), impact_nodes.end(), collision.first[0]) - impact_nodes.begin();
        ImpactZone::ImpactColInfo ci;
        auto vi0 = fvid.LocInd(collision.first[0]), vi3 = fvid.LocInd(collision.first[3]);
        ci.obj[0] = get<0>(vi0), ci.obj[1] = get<0>(vi3);
        ci.type = (collision.second == 1) ? ImpactZone::ColType::VertexFace : ImpactZone::ColType::EdgeEdge;
        ci.v[0] = get<2>(vi0), ci.v[1] = get<2>(fvid.LocInd(collision.first[1])), ci.v[2] = get<2>(fvid.LocInd(collision.first[2])), ci.v[3] = get<2>(vi3);
        izs[iz_mark[n] - 1].collisions.push_back(ci);
    }
    // std::cout << "Before SetDX computation" << std::endl;
    for (auto& iz: izs) iz.set_dx();
    // std::cout << "After SetDX computation" << std::endl;

    return ncols;
}

void BridsonCollisionManager::handleImpactZone(ImpactZone& iz){
    Vector xc = CGAL::NULL_VECTOR;
    Vector dxc = CGAL::NULL_VECTOR;
    DReal mc = 0;
    for (auto & i : iz.iz){
        auto& obj = i.first;
        auto v = i.second;
        xc += (obj->m_x[v] - CGAL::ORIGIN) * obj->m_m[v];
        dxc += obj->m_dx[v] * obj->m_m[v];
        mc += obj->m_m[v];
    }
    xc /= mc, dxc /= mc;
    Vector L = CGAL::NULL_VECTOR;
    std::array<DReal, 6> Idat{0};
    for (auto & r : iz.iz) {
        auto &obj = r.first;
        auto v = r.second;
        auto x = obj->m_x[v] - CGAL::ORIGIN - xc;
        L += obj->m_m[v] * cross_product(x, obj->m_dx[v] - dxc);
        double xsqr = x.squared_length();
        for (int i = 0, k = 0; i < 3; ++i)
        for (int j = 0; j <= i; ++j)
            Idat[k++] += obj->m_m[v] * ((i == j) ? xsqr : 0 - x[i] * x[j]);
    }
    typedef CGAL::Simple_cartesian<DReal>::Aff_transformation_3 Aff;
    Aff I(Idat[0], Idat[1], Idat[3],
            Idat[1], Idat[2], Idat[4],
            Idat[3], Idat[4], Idat[5]);
    Vector wn = I.inverse().transform(L);
    DReal w = sqrt(wn.squared_length());
    Vector n = CGAL::NULL_VECTOR;
    if (w > 1e-15) n = wn / w;
    for (auto & r : iz.iz) {
        auto &obj = r.first;
        auto v = r.second;
        auto x = obj->m_x[v] - CGAL::ORIGIN - xc;
        Vector xf = (x * n) * n;
        Vector xr = x - xf;
        obj->m_dx[v] = xc + dxc + xf + cos(w) * xr + sin(w) * cross_product(n, xr) - (obj->m_x[v] - CGAL::ORIGIN);
    }
}
