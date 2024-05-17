//
// Created by alex on 06.07.2020.
//

#include "BulletCollisionManager.h"

namespace World3d {
    using namespace std;

    btVector3 p2Bt(const Vector &p) { return btVector3(p.x(), p.y(), p.z()); }

    btVector3 p2Bt(const Point &p) { return btVector3(p.x(), p.y(), p.z()); }

    Vector Bt2p(const btVector3 &p) { return Vector(p.x(), p.y(), p.z()); }

    Point Bt2pp(const btVector3 &p) { return Point(p.x(), p.y(), p.z()); }

    static inline btDbvtVolume VolumeOf(const Object3D &m, const F_ind f,
                                        double margin) {
        auto v = vert_around(m.m_mesh, f);
        const btVector3 tmp[] = {p2Bt(m.m_x[v[0]]),
                                 p2Bt(m.m_x[v[1]]),
                                 p2Bt(m.m_x[v[2]])};
        const btVector3 *pts[] = {&tmp[0], &tmp[1], &tmp[2]};
        btDbvtVolume vol = btDbvtVolume::FromPoints(pts, 3);
        vol.Expand(btVector3(margin, margin, margin));
        return (vol);
    }

    BulletObject::BulletCollisionObject(Object3D& net, Bullet_Wrapper *body): m_body {body}, m_net{net}{
        m_nleaf.resize(net.m_mesh.num_vertices());
        _m_data1.reserve(net.m_mesh.num_vertices());
        for (auto v: net.m_mesh.vertices()){
            _m_data1.push_back(pair<V_ind, Bullet_Wrapper *>{v, body});
            m_nleaf[v] = m_ndbvt.insert(btDbvtVolume::FromCR(p2Bt(net.m_x[v]), m_body->getMargin()), &_m_data1.back());
        }

        m_fleaf.resize(net.m_mesh.num_faces());
        _m_data2.reserve(net.m_mesh.num_faces());
        for (auto f: net.m_mesh.faces()){
            _m_data2.push_back(pair<F_ind, Bullet_Wrapper *>{f, body});
            m_fleaf[f] = m_fdbvt.insert(VolumeOf(net, f, m_body->getMargin()), &_m_data2.back());
        }
        m_worldTransform.setIdentity();

        m_collisionShape = new BulletCollisionShape(this);
        m_collisionShape->setMargin(m_body->getMargin());
        m_collisionFlags &= ~CF_STATIC_OBJECT;
        updateBounds();
    }

void BulletObject::updateBounds() {
    if (m_fdbvt.m_root) {
        const btVector3 &mins = m_fdbvt.m_root->volume.Mins();
        const btVector3 &maxs = m_fdbvt.m_root->volume.Maxs();
        const btScalar csm = getCollisionShape()->getMargin();
        const btVector3 mrg = btVector3(csm, csm, csm) * 1;
        m_bounds[0] = mins - mrg;
        m_bounds[1] = maxs + mrg;
    } else if (m_ndbvt.m_root){
        const btVector3 &mins = m_ndbvt.m_root->volume.Mins();
        const btVector3 &maxs = m_ndbvt.m_root->volume.Maxs();
        const btScalar csm = getCollisionShape()->getMargin();
        const btVector3 mrg = btVector3(csm, csm, csm) * 1;
        m_bounds[0] = mins - mrg;
        m_bounds[1] = maxs + mrg;
    } else{
        m_bounds[0] = m_bounds[1] = btVector3(0, 0, 0);
    }

}

void BulletObject::updateCollisionInfo() {
    ATTRIBUTE_ALIGNED16(btDbvtVolume) vol;
    const static double extend_multiplier = 3.0; //recommend be more 1.0
    const static btScalar margin_scale = 0.001;//0.25 should be > 0
    for (auto v: m_net.m_mesh.vertices()) {
        vol = btDbvtVolume::FromCR(p2Bt(m_net.m_x[v]), getMargin());
        btDbvtNode *leaf = m_nleaf[v];
        btVector3 velocity = p2Bt(m_net.m_next_x[v] - m_net.m_x[v]) * extend_multiplier;
        btScalar margin = getMargin() * static_cast<double>(margin_scale);
        m_ndbvt.update(leaf, vol, velocity, margin);
    }

    for (auto f: m_net.m_mesh.faces()) {
        auto v = vert_around(m_net.m_mesh, f);
        const btVector3 tmp[] = {p2Bt(m_net.m_x[v[0]]),
                                 p2Bt(m_net.m_x[v[1]]),
                                 p2Bt(m_net.m_x[v[2]]),
                                 p2Bt(m_net.m_next_x[v[0]]),
                                 p2Bt(m_net.m_next_x[v[1]]),
                                 p2Bt(m_net.m_next_x[v[2]])};
        vol = btDbvtVolume::FromPoints(tmp, 6);
        btScalar mrg = getMargin();
        vol.Expand(btVector3(mrg, mrg, mrg));
        m_fdbvt.update(m_fleaf[f], vol, margin_scale*getMargin());
    }

    m_ndbvt.optimizeIncremental(1);
    m_fdbvt.optimizeIncremental(1);

    updateBounds();
}

void BulletObject::CollisionHandler(BulletObject *psb) {
    if (this->m_body == psb->m_body) return;
    BulletObject::Collider docol;
    docol.mrg = getMargin() + psb->getMargin();
    /* psb0 nodes vs psb1 faces	*/
    if (state & 1) {
        docol.psb[0] = this;
        docol.psb[1] = psb;
        docol.psb[0]->m_ndbvt.collideTT(docol.psb[0]->m_ndbvt.m_root,
                                        docol.psb[1]->m_fdbvt.m_root,
                                        docol);
    }
    /* psb1 nodes vs psb0 faces	*/
    if (psb->state & 1) {
        docol.psb[0] = psb;
        docol.psb[1] = this;
        docol.psb[0]->m_ndbvt.collideTT(docol.psb[0]->m_ndbvt.m_root,
                                        docol.psb[1]->m_fdbvt.m_root,
                                        docol);
    }
}

void BulletObject::CollisionHandler(const btCollisionObject *pco, btVector3 *triangle) {

    BulletObject::ColliderStatic docollide;

    ATTRIBUTE_ALIGNED16(btDbvtVolume) volume;
    volume = btDbvtVolume::FromPoints(triangle, 3);
    const btScalar margin = pco->getCollisionShape()->getMargin();
    volume.Expand(btVector3(margin, margin, margin));
    docollide.psb = this;
    for (int i = 0; i < 3; ++i) docollide.m_triangle[i] = triangle[i];
    docollide.mrg = pco->getCollisionShape()->getMargin();

    m_ndbvt.collideTV(m_ndbvt.m_root, volume, docollide);
}

void BulletObject::ColliderStatic::Process(const btDbvtNode *leaf) {
    using NodeData = std::pair<V_ind, Bullet_Wrapper*>;
    NodeData *data =
            static_cast<NodeData *>(leaf->data);
    auto pnt = data->second->m_net.m_x[data->first];
    auto btpnt = p2Bt(pnt);
    btVector3 prj(0, 0, 0);
    btScalar sqd{SIMD_INFINITY};
    ProjectOnTriangle(m_triangle[0], m_triangle[1], m_triangle[2], btpnt, prj, sqd);

    const double m = mrg + sqrt((data->second->m_net.m_next_x[data->first] - pnt).squared_length());
    if (static_cast<double>(sqd) < (m * m)) {
        double len = sqrt(sqd);

        Bullet_Wrapper::RContact c;
        c.m_normal = (btpnt - prj) / len;
        c.m_torthogonal = btCross(m_triangle[1] - m_triangle[0], m_triangle[2] - m_triangle[0]);
        c.m_proj = prj;
        c.m_margin = m;
        c.m_node = data->first;
        c.m_obstr = data->second;
        psb->m_body->m_rcontacts.push_back(c);
    }
}

void BulletObject::Collider::Process(const btDbvtNode *lnode,
                                     const btDbvtNode *lface) {
    using NodeData = std::pair<V_ind, Bullet_Wrapper*>;
    using FaceData = std::pair<F_ind, Bullet_Wrapper*>;
    NodeData *data = static_cast<NodeData *>(lnode->data);
    FaceData *edata = static_cast<FaceData *>(lface->data);
    auto v = data->second->m_net.m_x[data->first];
    auto vv = vert_around(edata->second->m_net.m_mesh, edata->first);
    auto &p = edata->second->m_net.m_x;
    btVector3 btpnt = p2Bt(v);
    btVector3 prj(0, 0, 0);
    btScalar sqd{SIMD_INFINITY};
    btVector3 face[] = {
            p2Bt(p[vv[0]]),
            p2Bt(p[vv[1]]),
            p2Bt(p[vv[2]])
    };
    ProjectOnTriangle(face[0], face[1], face[2], btpnt, prj, sqd);
    const double m = mrg + 2 * sqrt((data->second->m_net.m_next_x[data->first] - v).squared_length());

    if (static_cast<double>(sqd) < m * m) {
        const btVector3 w = BaryCoord(face[0], face[1], face[2], prj);
        const double ma = 1;
        double mb = 1;
        if (!(psb[1]->state & 1))
            mb = 0;
        const double ms = ma + mb;
        if (ms > 0) {
            Bullet_Wrapper::SContact c;
            c.m_normal = (btpnt - prj) / sqrt(sqd);
            c.m_margin = m;
            c.m_node = data->first;
            c.m_face = edata->first;
            c.m_weights = w;
            c.m_cfm[0] = ma / ms * 1;
            c.m_cfm[1] = mb / ms * 1;
            c.m_nobstr = edata->second;
            c.m_fobstr = data->second;
            psb[0]->m_body->m_scontacts.push_back(c);
        }
    }
}

Bullet_Wrapper::Bullet_Wrapper(Object3D &net, double margin) : m_net{net}, m_mrg{margin} {
    m_objs.push_back(new BulletObject(m_net, this));
}

void Bullet_Wrapper::addObjs2World(btCollisionWorld *wrld) {
    for (int i = 0, bc = m_objs.size(); i < bc; ++i)
        wrld->addCollisionObject(m_objs[i], btBroadphaseProxy::DefaultFilter, btBroadphaseProxy::AllFilter);
}
void Bullet_Wrapper::removeObjsFromWorld(btCollisionWorld* wrld){
    for (int i = 0, bc = m_objs.size(); i < bc; ++i)
        wrld->removeCollisionObject(m_objs[i]);
}

void Bullet_Wrapper::updateCollisionInfo() {
    for (int i = 0, bc = m_objs.size(); i < bc; ++i)
        m_objs[i]->updateCollisionInfo();

    m_scontacts.resize(0);
    m_rcontacts.resize(0);
}

void Bullet_Wrapper::solveRContacts(Bullet_Wrapper *psb) {
    for (int i = 0, ni = psb->m_rcontacts.size(); i < ni; ++i) {
        const RContact &c = psb->m_rcontacts[i];
        const Vector &nr = Bt2p(c.m_normal);
        const Vector &nt = Bt2p(c.m_torthogonal);
        const Point &p = Bt2pp(c.m_proj);
        auto &n = c.m_node;
        auto &no = c.m_obstr->m_net;
        if (!no._is_movable(no, n)) continue;
        Point &next = no.m_next_x[n];
        const Vector vr = no.m_next_x[n] - no.m_x[n];
        Vector corr{CGAL::NULL_VECTOR};
        double dot = vr * nr, dott = nr * nt;
        if (dot < 0 && dott >= 0) {
            const double j = c.m_margin - nr * (next - p);
            corr += j * nr;
        } else if (dott < 0) {
            const double j = c.m_margin + nr * (next - p);
            corr -= j * nr;
        }
        next += corr;
    }
}

void Bullet_Wrapper::solveSContacts(Bullet_Wrapper *psb) {
    for (int i = 0, ni = psb->m_scontacts.size(); i < ni; ++i) {
        const SContact &c = psb->m_scontacts[i];
        Vector nr = Bt2p(c.m_normal);
        V_ind n = c.m_node;
        F_ind f = c.m_face;
        Object3D &fo = c.m_nobstr->m_net;
        Object3D &no = c.m_fobstr->m_net;
        auto fvv = vert_around(fo.m_mesh, f);
        const Point p = CGAL::ORIGIN + BaryEval(fo.m_next_x[fvv[0]] - CGAL::ORIGIN,
                                                fo.m_next_x[fvv[1]] - CGAL::ORIGIN,
                                                fo.m_next_x[fvv[2]] - CGAL::ORIGIN,
                                                c.m_weights);
        const Point q = CGAL::ORIGIN + BaryEval(fo.m_x[fvv[0]] - CGAL::ORIGIN,
                                                fo.m_x[fvv[1]] - CGAL::ORIGIN,
                                                fo.m_x[fvv[2]] - CGAL::ORIGIN,
                                                c.m_weights);
        bool penetrate_node = false;
//            {
//                Vector vnormal{CGAL::NULL_VECTOR};
//                auto& mesh = no.m_mesh;
//                auto it = faces_around_target(mesh.halfedge(n), mesh);
//                int cnt = 0;
//                for (auto fit = it.begin(); fit != it.end(); ++fit)
//                    if (fit->is_valid()) {
//                        auto _f = *fit;
//                        auto vv = vert_around(mesh, _f);
//                        vnormal += CGAL::cross_product(Vector(no.m_x[vv[0]], no.m_x[vv[1]]),
//                                                       Vector(no.m_x[vv[0]], no.m_x[vv[2]]));
//                        cnt++;
//                    }
//                if (cnt) vnormal /= sqrt(vnormal.squared_length());
//                else vnormal = -nr;
//                if (vnormal * nr > 0) penetrate_node = true;
//                if (penetrate_node) nr *= -1;
//            }
        const Vector vr = (no.m_next_x[n] - no.m_x[n]) - (p - q);
        Vector corr{CGAL::NULL_VECTOR};
        double dot = vr * nr;
        if (dot < 0 || penetrate_node) {
            const double j = c.m_margin - nr * (no.m_next_x[n] - p);
            corr += j * nr;

            if (no._is_movable(no, n)
                || !fo._is_movable(fo, fvv[0]) || !fo._is_movable(fo, fvv[1]) || !fo._is_movable(fo, fvv[2])) {
                auto dif = c.m_cfm[0] * corr;
                no.m_next_x[n] += no.withBCmask(n, dif);
            }
            for (unsigned j = 0; j < 3; ++j)
                if (fo._is_movable(fo, fvv[j]) || !no._is_movable(no, n)) {
                    auto dif = -c.m_cfm[1] * static_cast<double>(c.m_weights[j]) * corr;
                    fo.m_next_x[fvv[j]] += fo.withBCmask(fvv[j], dif);
                }
        }

//            const Vector vr = (no.m_next_x[n] - no.m_x[n]) - (p - q);
//            Vector corr{CGAL::NULL_VECTOR};
//            double dot = vr * nr;
//            if (dot < 0) {
//                const double j = c.m_margin - nr * (no.m_next_x[n] - p);
//                corr += j * nr;
//
//                if (no._is_movable(no, n)
//                    || !fo._is_movable(fo, fvv[0]) || !fo._is_movable(fo, fvv[1]) || !fo._is_movable(fo, fvv[2])) {
//                    auto dif = (!(nr * fo.m_normal[f] < 0) ? 1.0 : -1.0) * c.m_cfm[0] * corr;
//                    no.m_next_x[n] += no.withBCmask(n, dif);
//                }
//                for (unsigned j = 0; j < 3; ++j)
//                    if (fo._is_movable(fo, fvv[j]) || !no._is_movable(no, n)) {
//                        auto dif = -c.m_cfm[1] * static_cast<double>(c.m_weights[j]) * corr;
//                        fo.m_next_x[fvv[j]] += fo.withBCmask(fvv[j], dif);
//                    }
//            }
    }
}

Bullet_Wrapper::~Bullet_Wrapper() {
    for (int i = 0; i < m_objs.size(); ++i)
        delete m_objs[i];
}

void BulletCollisionManager::addObj(Object3D& newobj, ObjectID id, bool is_dynamic, double margin) {
    if (is_dynamic){
        Bullet_Wrapper* bw = new Bullet_Wrapper(newobj, margin);
        bw->addObjs2World(_collisionWorld.get());
        _collisions.insert({id, {bw, true}});
    }
    else{
        auto obj = convert_Object3D_to_btTriangleMesh(newobj, id, margin);
        _collisionWorld->addCollisionObject(obj, btBroadphaseProxy::StaticFilter, btBroadphaseProxy::AllFilter ^ btBroadphaseProxy::StaticFilter);
        _collisions.insert({id, {obj, false}});
    }
}

void BulletCollisionManager::removeObj(ObjectID id) {
    auto it = _collisions.find(id);
    if (it == _collisions.end()) return;
    auto& obj = it->second;
    if (obj.second){
        Bullet_Wrapper* bw = static_cast<Bullet_Wrapper*>(obj.first);
        bw->removeObjsFromWorld(_collisionWorld.get());
        delete bw;
    }else{
        btCollisionObject* bo = static_cast<btCollisionObject*>(obj.first);
        _collisionWorld->removeCollisionObject(bo);
        _internal_data.erase(it->first);
    }
    _collisions.erase(it);
}

void BulletCollisionManager::findCollisions(){
    for (auto& i: _collisions){
        if (i.second.second){
            Bullet_Wrapper* body = static_cast<Bullet_Wrapper*>(i.second.first);
            body->updateCollisionInfo();
        }
    }

    _collisionWorld->performDiscreteCollisionDetection();
}

void BulletCollisionManager::set_margin(ObjectID id, double margin) {
    auto it = _collisions.find(id);
    if (it == _collisions.end()) return;
    if (it->second.second){
        Bullet_Wrapper* bw = static_cast<Bullet_Wrapper*>(it->second.first);
        bw->setMargin(margin);
    } else{
        get<3>(_internal_data[id]).get()->setMargin(margin);
    }
}

void BulletCollisionManager::set_margin(double margin) {
    for (auto it: _collisions){
        if (it.second.second){
            Bullet_Wrapper* bw = static_cast<Bullet_Wrapper*>(it.second.first);
            bw->setMargin(margin);
        } else{
            get<3>(_internal_data[it.first]).get()->setMargin(margin);
        }
    }
}

double BulletCollisionManager::get_margin(ObjectID id) {
    auto it = _collisions.find(id);
    if (it == _collisions.end()) return -1.0;
    if (it->second.second){
        Bullet_Wrapper* bw = static_cast<Bullet_Wrapper*>(it->second.first);
        return bw->getMargin();
    } else{
        return get<3>(_internal_data[id]).get()->getMargin();
    }
}

void BulletCollisionManager::registerCollisionWorld(std::pair<Vector, Vector> min_max_Aabb, int maxProxies) {
    btVector3 worldAabbMin(min_max_Aabb.first.x(), min_max_Aabb.first.y(), min_max_Aabb.first.z());
    btVector3 worldAabbMax(min_max_Aabb.second.x(), min_max_Aabb.second.y(), min_max_Aabb.second.z());

    _broadphase = std::make_unique<btAxisSweep3>(worldAabbMin, worldAabbMax, maxProxies);
    _collisionConfig = std::make_unique<MyCollisionConfiguration>(btDefaultCollisionConstructionInfo());
    _dispatcher = std::make_unique<btCollisionDispatcher>(_collisionConfig.get());
    _collisionWorld = std::make_unique<btCollisionWorld>(_dispatcher.get(), _broadphase.get(), _collisionConfig.get());
}

btCollisionObject *
BulletCollisionManager::convert_Object3D_to_btTriangleMesh(const Object3D& st, const ObjectID id, double mrg) {
    using namespace std;
    vector<btVector3> vertices(st.m_mesh.num_vertices());
    for (auto v: st.m_mesh.vertices())
        vertices[v] = p2Bt(st.m_x[v]);
    vector<int> triag_indeces(st.m_mesh.num_faces()*3);
    for(auto f: st.m_mesh.faces()){
        auto v = vert_around(st.m_mesh, f);
        for (auto j = 0; j < 3; ++j)
            triag_indeces[3*f + j] = v[j];
    }
    const int vertStride = sizeof(btVector3);
    const int indexStride = 3 * sizeof(int);
    auto indexVertexArrays = make_unique<btTriangleIndexVertexArray>( st.m_mesh.num_faces(),
                                                                      triag_indeces.data(),
                                                                      indexStride,
                                                                      st.m_mesh.num_vertices(),
                                                                      reinterpret_cast<btScalar*> ((vertices.data())),
                                                                      vertStride);
    bool useQuantizedAabbCompression = true;
    auto meshShape = make_unique<btBvhTriangleMeshShape>(indexVertexArrays.get(), useQuantizedAabbCompression);
    meshShape->setMargin(mrg);//16);
    meshShape->buildOptimizedBvh();
    auto newOb = make_unique<btCollisionObject>();
    btTransform tr;
    tr.setIdentity();
    newOb->setWorldTransform(tr);
    newOb->setInterpolationWorldTransform(tr);
    newOb->setCollisionShape(meshShape.get());

    auto it = _internal_data.insert(pair<ObjectID, StaticInfo>{id, StaticInfo{move(vertices), move(triag_indeces), move(indexVertexArrays), move(meshShape), move(newOb)}});
    return get<4>(it.first->second).get();
}

void BulletCollisionManager::solveCollisions(){ m_solve_collission(_collisions); }
}

#include "BulletCollisionAlgorithms.h"

namespace World3d{
    MyCollisionConfiguration::MyCollisionConfiguration(const btDefaultCollisionConstructionInfo& constructionInfo)
            : btDefaultCollisionConfiguration(constructionInfo)
    {
        void* mem;
        mem = btAlignedAlloc(sizeof(NetNetCollisionAlgorithm::CreateFunc), 16);
        m_softSoftCreateFunc = new (mem) NetNetCollisionAlgorithm::CreateFunc;

        mem = btAlignedAlloc(sizeof(NetConcaveCollisionAlgorithm::CreateFunc), 16);
        m_softBodyConcaveCreateFunc = new (mem) NetConcaveCollisionAlgorithm::CreateFunc;

        mem = btAlignedAlloc(sizeof(NetConcaveCollisionAlgorithm::CreateFunc), 16);
        m_swappedSoftBodyConcaveCreateFunc = new (mem) NetConcaveCollisionAlgorithm::SwappedCreateFunc;
        m_swappedSoftBodyConcaveCreateFunc->m_swapped = true;
    }

    MyCollisionConfiguration::~MyCollisionConfiguration()
    {
        m_softSoftCreateFunc->~btCollisionAlgorithmCreateFunc();
        btAlignedFree(m_softSoftCreateFunc);
        m_softBodyConcaveCreateFunc->~btCollisionAlgorithmCreateFunc();
        btAlignedFree(m_softBodyConcaveCreateFunc);
        m_swappedSoftBodyConcaveCreateFunc->~btCollisionAlgorithmCreateFunc();
        btAlignedFree(m_swappedSoftBodyConcaveCreateFunc);
    }

    btCollisionAlgorithmCreateFunc* MyCollisionConfiguration::getCollisionAlgorithmCreateFunc(int proxyType0, int proxyType1)
    {
        ///softbody versus softbody
        if ((proxyType0 == SOFTBODY_SHAPE_PROXYTYPE) && (proxyType1 == SOFTBODY_SHAPE_PROXYTYPE))
        {
            return m_softSoftCreateFunc;
        }

        ///softbody versus convex
        if (proxyType0 == SOFTBODY_SHAPE_PROXYTYPE && btBroadphaseProxy::isConcave(proxyType1))
        {
            return m_softBodyConcaveCreateFunc;
        }

        ///convex versus soft body
        if (btBroadphaseProxy::isConcave(proxyType0) && proxyType1 == SOFTBODY_SHAPE_PROXYTYPE)
        {
            return m_swappedSoftBodyConcaveCreateFunc;
        }
        return btDefaultCollisionConfiguration::getCollisionAlgorithmCreateFunc(proxyType0, proxyType1);
    }

    int PenaltyCollisionAlgo::default_penalty_force(const PenaltyCollisionAlgo::PenaltyContact &c) {
        auto v = c.m_node;
        Object3D* obj = c.m_obj;
        Vector vnormal{CGAL::NULL_VECTOR};
        {
            auto& mesh = obj->m_mesh;
            auto it = faces_around_target(mesh.halfedge(v), mesh);
            int cnt = 0;
            for (auto fit = it.begin(); fit != it.end(); ++fit)
                if (fit->is_valid()) {
                    auto f = *fit;
                    auto vv = vert_around(obj->m_mesh, f);
                    vnormal += CGAL::cross_product(Vector(obj->m_x[vv[0]], obj->m_x[vv[1]]),
                                                   Vector(obj->m_x[vv[0]], obj->m_x[vv[2]]));
                    cnt++;
                }
            if (cnt) vnormal /= sqrt(vnormal.squared_length());
            else vnormal = -c.m_proj_dir;
        }
        Vector tnormal = c.m_normal;
        if (vnormal*tnormal > 0) tnormal *= -1;

        double fn = obj->m_F[v] * tnormal;
        double de = (obj->m_x[v] - c.m_proj)*tnormal;
        double d0 = c.m_margin/2;
        if (de > 2*d0 || (de > 0 && fn > 0)) return 0;
        auto sq = [](auto x) { return x*x; };
        double s = ( (de >= d0) ? sq(2 - de/d0) : 3-2*de/d0 );
        if (de < 0 && fn > 0) fn *= -1;
        Vector Fn = -fn * tnormal;
        obj->m_F[v] += Fn*s;

        return 0;
    }

    int PenaltyCollisionAlgo::operator()(map<ObjectID, std::pair<void *, bool>> &collisions) {
//        static int it = 0;
//        static double t[10] = {0};
//        World3d::Timer tm, _tm[10];
        for (auto i: collisions){
            if (i.second.second){
                Bullet_Wrapper* body = static_cast<Bullet_Wrapper*>(i.second.first);
                {
                    for (int j = 0, ni = body->m_scontacts.size(); j < ni; ++j) {
//                        tm.reset();
//                        _tm[0].reset();
                        const auto &s = body->m_scontacts[j];
                        V_ind n = s.m_node;
                        F_ind f = s.m_face;
                        Object3D &fo = s.m_nobstr->m_net;
                        Object3D &no = s.m_fobstr->m_net;
                        if (!no.is_movable(n)) continue;
                        auto fvv = vert_around(fo.m_mesh, f);
//                        t[4] += _tm[0].elapsed();
//
//                        _tm[1].reset();
                        std::array<double, 3> proj = {0};
                        for (int l = 0; l < 3; ++l) {
                            auto& p = fo.m_x[fvv[l]];
                            for (int k = 0; k < 3; ++k) proj[k] += s.m_weights[l] * p[k];
                        }
                        PenaltyContact c{n,
                                         &no,
                                         Point(proj[0], proj[1], proj[2]),
                                         Bt2p(s.m_normal),
                                         fo.m_normal[f],
                                         s.m_margin};
//                        t[5] += _tm[1].elapsed();
//                        t[0] += tm.elapsed(); tm.reset();

                        set_force(c);
//                        t[1] += tm.elapsed();
                    }
                    for (int j = 0, ni = body->m_rcontacts.size(); j < ni; ++j) {
//                        tm.reset();
                        const auto &s = body->m_rcontacts[j];
                        auto &n = s.m_node;
                        auto &no = s.m_obstr->m_net;
                        if (!no.is_movable(n)) continue;

                        PenaltyContact c;
                        c.m_node = n;
                        c.m_obj = &no;
                        c.m_proj = CGAL::ORIGIN + Bt2p(s.m_proj);
                        c.m_proj_dir = Bt2p(s.m_normal);
                        c.m_normal = Bt2p(s.m_torthogonal.normalized());
                        c.m_margin = s.m_margin;
//                        t[2] += tm.elapsed(); tm.reset();

                        set_force(c);
//                        t[3] += tm.elapsed();
                    }
                }
            }
        }
//        it++;
//        if (it % 1000 == 0) {
//            std::cout << "PenaltyCollisionAlgo perform (" << it << ") :";
//            for (int i = 0; i < 10; ++i) std::cout << " " << t[i];
//            std::cout << std::endl;
//        }

        for (auto i: collisions){
            if (!i.second.second) continue;
            Bullet_Wrapper* body = static_cast<Bullet_Wrapper*>(i.second.first);
            if (body->m_scontacts.size() == 0 && body->m_rcontacts.size() == 0) continue;
            body->m_net.update_next_x();
        }
        return 0;
    }

}



