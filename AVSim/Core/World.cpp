//
// Created by alex on 30.06.2020.
//

#include "World.h"

void setPrevWatch(World3d::Object3D& obj){
    auto it = obj.m_mesh.add_property_map<World3d::V_ind, World3d::Point>("v:prev_x_1");
    for (auto v: obj.m_mesh.vertices())
        it.first[v] = obj.m_x[v];
    auto prev_x_updater = [](World3d::Object3D &obj) -> int {
        auto it = obj.m_mesh.property_map<World3d::V_ind, World3d::Point>("v:prev_x_1");
        if (!it.second) it = obj.m_mesh.add_property_map<World3d::V_ind, World3d::Point>("v:prev_x_1");
        std::copy(obj.m_x.begin(), obj.m_x.end(), it.first.begin());
        return 0;
    };
    obj.add_updater(prev_x_updater);
}
//
//void updateShiftWatch(World3d::Object3D& obj){
//    auto it = obj.m_mesh.property_map<World3d::V_ind, World3d::Vector>("v:shift");
//    for (auto v: obj.m_mesh.vertices())
//        it.first[v] =
//}

World3d::ObjectID World3d::World::addObject3D(Object3D &&obj, int status) {
    static int objectID = 1;
    auto it = _objects.insert({objectID, {obj, status}});
    setPrevWatch(it.first->second.first);
    _renderer->set_renderer(objectID, it.first->second.first, it.first->second.second != 0);
    it.first->second.first.__id = objectID;
    _collisionManager.get()->addObj(it.first->second.first, it.first->first, status > 0);
    for (auto& f: m_forces) f.second->addObj(it.first->second.first, it.first->first, status > 0);
    return objectID++;
}

void World3d::World::setCollider(std::unique_ptr<CollisionManagerBase>&& collider){
    _collisionManager = move(collider);
    for (auto& o: _objects)
        _collisionManager.get()->addObj(o.second.first, o.first, o.second.second > 0);
}

World3d::Object3D World3d::World::releaseObj(ObjectID oid){
    auto it = _objects.find(oid);
    if (it == _objects.end()) throw std::runtime_error("Object not found");
    Object3D obj = std::move(it->second.first);
    removeObj(oid);
    return obj;
}

void World3d::World::removeObj(ObjectID oid){
    auto it = _objects.find(oid);
    if (it == _objects.end()) return;
    _collisionManager.get()->removeObj(oid);
    _renderer.get()->delete_renderer(oid, it->second.first);
    for (auto& f: m_forces) f.second->removeObj(oid);
    _objects.erase(oid);
}

// update data in updatable properties of objects
// F = F(x) - compute forces in current position
int World3d::World::UpdateForces() {
    int status = 0;
    for (auto& i: _objects){
        if (i.second.second) {
            i.second.first.apply_updaters();
            if ((status = i.second.first.apply_forces()) < 0){
                std::cout << "Error in apply_forces(...) on object \"" << i.second.first.name << "\" with id = " << i.first << std::endl;
                return status;
            }
        }
    }
    for (auto& f: m_forces)
        if ((status = f.second.operator()()) < 0){
            std::cout << "Error in WorldForce  with id = " << f.first << std::endl;
            return status;
        }

    return status;
}

// x_next = x_cur + Δx(F)
int World3d::World::ApplyForces() {
    for (auto& i: _objects){
        if (i.second.second)
            i.second.first.update_next_x();
    }
    return 0;
}

// x_next = x_next + collision_shift(x_next, x_cur)
int World3d::World::ApplyCollision(){
    _collisionManager.get()->findCollisions();
    _collisionManager.get()->solveCollisions();
    return 0;
}

int World3d::World::StepSimulation() {
    static Timer t;

    int status = _stepAlgo(this);

    if (t.elapsed() > 1.0/60) {
        Render();
        t.reset();
    }

    return status;
}

int World3d::World::Simulation(std::function<bool(StepSimInfo&, World*)> stopCondition){
    StepSimInfo ssi;
    ssi.status = 0;
    for (ssi.it = 0; !stopCondition(ssi, this); ++ssi.it){
        ssi.status = StepSimulation();
    }
    return ssi.status;
}

int World3d::World::Render() {
    for (auto& i: _objects)
        i.second.first.apply_render();
    _renderer->render();
    return 0;
}

void World3d::World::setRenderer(Renderer&& renderer){
    _renderer = move(renderer);
    for (auto& o: _objects)
        _renderer->set_renderer(o.first, o.second.first, o.second.second != 0);
}

World3d::World::World() {
//    _renderer = std::make_unique<World3d::DefaultRenderer>();
    _renderer = std::make_unique<World3d::RendererBase>();
    _collisionManager = std::make_unique<CollisionManagerBase>();
    _stepAlgo = [](World* w)->int{
        int status = 0;
        if ((status = w->UpdateForces()) < 0) return status;
        if ((status = w->ApplyForces()) < 0) return status;
        if ((status = w->ApplyCollision()) < 0) return status;
        if ((status = w->ApplyNext()) < 0) return status;
        return status;
    };
}

// x_cur = next(x_next, x_cur), typically next(x_next, x_cur) = x_next
int World3d::World::ApplyNext() {
    for (auto& i: _objects)
        if(i.second.second)
            i.second.first.apply_next_x();
    return 0;
}

World3d::Force_ID World3d::World::addForce(World3d::Force f, World3d::ObjectID oid) {
    return _objects[oid].first.add_force(f);
}

void World3d::World::setForceApplier(World3d::ForceApplier f, World3d::ObjectID oid) {
    _objects[oid].first.set_NextXUpdater(f);
}

std::map<World3d::ObjectID, World3d::Force_ID> World3d::World::addForce(World3d::Force f) {
    std::map<ObjectID, Force_ID> res;
    for (auto& i: _objects){
        if (i.second.second)
            res.insert({i.first, i.second.first.add_force(f)});
    }
    return res;
}

void World3d::World::setForceApplier(World3d::ForceApplier f) {
    for (auto& i: _objects){
        if (i.second.second)
            i.second.first.set_NextXUpdater(f);
    }
}

World3d::Object3D& World3d::World::obj(ObjectID id) {
    return _objects[id].first;
}

World3d::WorldForce& World3d::World::force(WorldForce_ID id) {
    return m_forces[id.id];
}

double World3d::World::getShift(ObjectID id, unsigned normOut, unsigned normIn){
    ConstProperty<V_ind, Point> it = obj(id).m_mesh.property_map<V_ind, Point>("v:prev_x_1");

    struct ShiftRange{
        ConstProperty<V_ind, Point>* prop;
        Object3D* obj;
        class iterator: public std::iterator<std::output_iterator_tag, Vector>{
            Mesh::Vertex_iterator vit;
            ConstProperty<V_ind, Point>* prop;
            Object3D* obj;
            mutable Vector res;
        public:
            explicit iterator(Mesh::Vertex_iterator vit, ConstProperty<V_ind, Point>* prop, Object3D* obj) : vit(vit), prop(prop), obj(obj) {}
            iterator& operator++() {vit++; return *this;}
            iterator operator++(int) {iterator retval = *this; ++(*this); return retval;}
            bool operator==(iterator other) const {return vit == other.vit;}
            bool operator!=(iterator other) const {return !(*this == other);}
            reference operator*() const {
                res = Vector(prop->first[*vit], obj->m_x[*vit]);
                return res;
            }
        };
        iterator begin() const {return iterator(obj->m_mesh.vertices().begin(), prop, obj);}
        iterator end() const {return iterator(obj->m_mesh.vertices().end(), prop, obj);}
        ShiftRange(ConstProperty<V_ind, Point>* prop, Object3D* obj): prop(prop), obj(obj) {}
    };

    return obj(id).residual<Vector>(ShiftRange(&it, &obj(id)), normOut, normIn);
}

double World3d::World::getNorm(ObjectID id, std::string on_tag, unsigned normOut, unsigned normIn){
    return obj(id).residual(normOut, normIn, on_tag);
}

double World3d::World::getWorldNorm(std::string on_tag, unsigned normObj, unsigned normOut, unsigned normIn){
    if (_objects.size() == 0) return 0;
    if (_objects.size() == 1) return _objects.begin()->second.first.residual(normOut, normIn, on_tag);
    double res = 0;
    if (normObj == 0){
        for (auto& i: _objects)
            res = std::max(res, i.second.first.residual(normOut, normIn, on_tag));
    } else if (normObj == 1){
        for (auto& i: _objects)
            res += i.second.first.residual(normOut, normIn, on_tag);
    } else if (normObj == 2){
        for (auto& i: _objects) {
            double q = i.second.first.residual(normOut, normIn, on_tag);
            res += q * q;
        }
        res = sqrt(res);
    } else {
        for (auto& i: _objects) {
            double q = i.second.first.residual(normOut, normIn, on_tag);
            res += pow(q, normObj);
        }
        res = pow(res, 1.0/normObj);
    }
    return res;
}

double World3d::World::getWorldShift(unsigned normObj, unsigned normOut, unsigned normIn){
    if (_objects.size() == 0) return 0;
    if (_objects.size() == 1) return getShift(_objects.begin()->first, normOut, normIn);
    double res = 0;
    if (normObj == 0){
        for (auto i: _objects)
            res = std::max(res, getShift(i.first, normOut, normIn));
    } else if (normObj == 1){
        for (auto i: _objects)
            res += getShift(i.first, normOut, normIn);
    } else if (normObj == 2){
        for (auto i: _objects) {
            double q = getShift(i.first, normOut, normIn);
            res += q * q;
        }
        res = sqrt(res);
    } else {
        for (auto i: _objects) {
            double q = getShift(i.first, normOut, normIn);
            res += pow(q, normObj);
        }
        res = pow(res, 1.0/normObj);
    }
    return res;
}

World3d::World::~World() {
    for (auto& o: _objects){
        _renderer->delete_renderer(o.first, o.second.first);
    }
    _renderer->render();
}

void World3d::World::clearObjects() {
    for (auto& o: _objects){
        _renderer->delete_renderer(o.first, o.second.first);
        _collisionManager.get()->removeObj(o.first);
    }
    _renderer->render();
    _objects.clear();
}

void World3d::World::clearRenderer() {
    for (auto& o: _objects){
        _renderer->delete_renderer(o.first, o.second.first);
    }
    _renderer->render();
    _renderer = std::make_unique<World3d::RendererBase>();
}

void World3d::World::clearStepAlgo() {
    _stepAlgo = [](World* w)->int{
        int status = 0;
        if ((status = w->UpdateForces()) < 0) return status;
        if ((status = w->ApplyForces()) < 0) return status;
        if ((status = w->ApplyCollision()) < 0) return status;
        if ((status = w->ApplyNext()) < 0) return status;
        return status;
    };
}

void World3d::World::clearCollisionManager() { setCollider(std::make_unique<CollisionManagerBase>()); }

void World3d::World::clear() {
    clearObjects();
    clearRenderer();
    clearCollisionManager();
    clearStepAlgo();
}

World3d::WorldForce_ID World3d::World::addForce(World3d::WorldForce f) {
    static int force_id = 0;
    auto it = m_forces.insert({force_id, std::move(f)});
    it.first->second->registerWorld(*this);
    return WorldForce_ID(force_id++);
}

int World3d::World::compStateSize() const {
    int N = 0;
    for (const auto &i: _objects)
        if (i.second.second)
            N += 3 * i.second.first.m_mesh.num_vertices();
    return N;
}

int World3d::World::compResidual(double *R) {
    auto &objs = _objects;
    if (!full_residual_updated){
        UpdateForces();
        full_residual_updated = true;
    }

    for (auto &i: objs)
        if (i.second.second) {
            Object3D &obj = i.second.first;
            for (auto v: obj.m_mesh.vertices()) {
                auto bc = obj.getBC(v);
                for (int j = 0; j < 3; ++j)
                    *(R++) = (bc & (1 << j)) ? obj.m_F[v][j] : 0.0;
            }
        }
    return 0;
}

int World3d::World::compJacobian(SparseMatrix *sm) {
    auto &objs = _objects;
    int N = 0;
    int stat = 0;
    for (auto &i: objs)
        if (i.second.second) {
            Object3D &obj = i.second.first;
            SparseMatrix::IndexFrom locm(*sm, N, N);
            for (auto &F: obj.m_forces) {
                if ((stat = F.second->fill_matrix(obj, locm)) < 0) return stat;
            }
            N += 3 * i.second.first.m_mesh.num_vertices();
        }
    for (auto& f: m_forces)
        if ((stat = f.second->fill_matrix(sm)) < 0) return stat;
    return stat;
}

int World3d::World::compStateVector(double *X) const {
    for (auto &i: _objects)
        if (i.second.second)
            for (auto v: i.second.first.m_mesh.vertices()) {
                for (int k = 0; k < 3; ++k)
                    *(X++) = i.second.first.m_x[v][k];
            }
    return 0;
}

int World3d::World::setStateVector(const double *X, bool update_all_forces) {
    for (auto &i: _objects)
        if (i.second.second)
            for (auto v: i.second.first.m_mesh.vertices()) {
                auto bc = i.second.first.getBC(v);
                auto &p = i.second.first.m_x[v];
                p = Point((bc & 1) ? *(X + 0) : p[0], (bc & 2) ? *(X + 1) : p[1], (bc & 4) ? *(X + 2) : p[2]);
                if (std::isnan(X[0]) || std::isnan(X[1]) || std::isnan(X[2]))
                    throw std::runtime_error("ERROR: setStateVector(): set NAN Point at v = " + std::to_string(v.idx()) + " of obj = \"" + i.second.first.name + "\"");
                X += 3;
            }
    int stat = 0;
    if (update_all_forces){        
        stat = UpdateForces();
        full_residual_updated = true;
    } else 
        for (auto& i: _objects)
            i.second.first.apply_updaters();    
    return stat;
}

int World3d::World::compFResidual(double* R, ObjectID oid, Force_ID fid){
    full_residual_updated = false;
    Object3D* pobj = nullptr;
    auto &objs = _objects;

    for (auto it = objs.begin(); it != objs.end(); ++it){
        auto& i = *it;
        if (i.second.second) {
            Object3D &obj = i.second.first;
            if (i.first != oid) {
                R += 3*obj.m_mesh.num_vertices();
                continue;
            }
            pobj = &obj;
            break;
        }
    }
    pobj->set_zero_residual();
    int status = 0;
    if ((status = pobj->apply_force(fid)) < 0) return status;
    for (auto v: pobj->m_mesh.vertices()) {
        auto bc = pobj->getBC(v);
        for (int j = 0; j < 3; ++j)
            *(R++) = (bc & (1 << j)) ? pobj->m_F[v][j] : 0.0;
    }
    return status;
}

int World3d::World::compWFResidual(double* R, WorldForce_ID wfid){
    full_residual_updated = false;
    auto& f = force(wfid);

    for (auto& obj: objs()) obj.second.first.set_zero_residual();
    
    int status = f.operator()();
    if (status < 0){
        std::cout << "Error in WorldForce  with id = " << wfid.id << std::endl;
    }
    for (auto& obj: objs()){
        if (!obj.second.second) continue;
        Object3D* pobj = &obj.second.first;
        for (auto v: pobj->m_mesh.vertices()) {
            auto bc = pobj->getBC(v);
            for (int j = 0; j < 3; ++j)
                *(R++) = (bc & (1 << j)) ? pobj->m_F[v][j] : 0.0;
        }
    }
    return status;
}

int World3d::World::compFJacobian(SparseMatrix *sm, ObjectID oid, Force_ID fid){
    int N = 0;
    int stat = 0;
    for (auto& i: objs())
        if (i.second.second){
            if (i.first == oid){
                Object3D &obj = i.second.first;
                SparseMatrix::IndexFrom locm(*sm, N, N);
                auto Fit = obj.m_forces.find(fid);
                if (Fit == obj.m_forces.end())
                    throw std::runtime_error("Force id = " + std::to_string(fid) + " on object " + std::to_string(oid.id) + " not found");
                auto& F = Fit->second;
                if ((stat = F->fill_matrix(obj, locm)) < 0) return stat;
                break;
            }
            N += 3 * i.second.first.m_mesh.num_vertices();
        }
    return stat;    
}

int World3d::World::compWFJacobian(SparseMatrix *sm, WorldForce_ID wfid){
    return force(wfid)->fill_matrix(sm);
}

//void World3d::NewtonSolverStepAlgo::operator()(World3d::World *w) {
//    //removed timers
////    static Timer tt;
////    auto print_T = [](std::string s = ""){ std::cout << "\ttime = " << tt.elapsed() << " " << s << std::endl; };
////    print_T("Before compute forces");
//    w->UpdateForces();
////    print_T("After compute forces");
//
//    auto& objs = w->objs();
//    int N = 0;
//    for (auto& i: objs)
//        if (i.second.second)
//            N += 3 * i.second.first.m_mesh.num_vertices();
//    auto& solver = _solver;
//
//    solver.sm.clear();
//    solver.sm.resize(N);
//    solver.rhs.resize(N);
//    std::fill(solver.x.begin(), solver.x.end(), 0.0);
//    solver.x.resize(N);
////    print_T("After resize internal solver data");
//
//    N = 0;
//    for (auto& i: objs)
//        if (i.second.second){
//            Object3D& obj = i.second.first;
//            SparseMatrix::IndexFrom locm(solver.sm, N, N);
////            print_T("Before assemble matrix");
//            for (auto& F: obj.m_forces){
////                print_T("Before assemble matrix from " + F.second->type);
//                F.second->fill_matrix(obj, locm);
//            }
////            print_T("Before copy RHS");
//            for (auto v: obj.m_mesh.vertices()) {
//                auto bc = obj.getBC(v);
//                for (int j = 0; j < 3; ++j)
//                    solver.rhs[N + 3 * v.idx() + j] = (bc & (1 << j)) ? obj.m_F[v][j] : 0.0;
//            }
//
//            N += 3 * i.second.first.m_mesh.num_vertices();
//        }
//    N = 0;
//    double tau = tau_algo(m_it, w);
//    m_it++;
////    print_T("After assembling data");
//
//    solver.solver.SetMatrix(solver.sm);
////    print_T("After set matrix");
//    solver.solver.Solve(solver.rhs, solver.x);
////    print_T("After solve matrix");
//    for (auto& i: objs)
//        if (i.second.second) {
//            Object3D &obj = i.second.first;
//            for (auto i: obj.m_mesh.vertices()){
//                auto bc = obj.getBC(i);
//                if (obj._is_movable(obj, i)){
//                    Vector dx((bc & 1) ? solver.x[N+3*i] : 0.0, (bc & 2) ? solver.x[N+3*i + 1] : 0.0, (bc & 4) ? solver.x[N+3*i + 2] : 0.0);
//                    obj.m_next_x[i] = obj.m_x[i] - tau * dx;
//                }
//            }
//            N += 3 * i.second.first.m_mesh.num_vertices();
//        }
//
//    w->ApplyCollision();
//    w->ApplyNext();
////    print_T("After apply result");
//}
bool World3d::StaticStopCondition::operator()(World3d::StepSimInfo &info, World3d::World *w) {
    auto printStepSimulationError = [](std::ostream& out, StepSimInfo& info, World3d::Timer* time){
        out << "During computation of rhs error occured with status = " << info.status << " at it = " << info.it;
        if (time) out<<" time = " << time -> elapsed();
        out << std::endl;
    };
    auto printCurrentState = [](std::ostream& out, int it, double eps, double resid, World3d::Timer* time){
        out << "it " << it << ": rel_resid = " << eps << " abs_resid = " << resid;
        if (time) out<<" time = " << time -> elapsed();
        out << std::endl;
    };
    auto printConverge = [](std::ostream& out, int it, double eps, double resid, World3d::Timer* time){
        out << "Algorithm is converged: it = " << it << " rel_resid = " << eps << " abs_resid = " << resid;
        if (time) out<<" time = " << time -> elapsed();
        out << std::endl;
    };
    auto printDiverge = [](std::ostream& out, int it, double eps, double resid, World3d::Timer* time){
        out << "Algorithm is diverged: it = " << it << " rel_resid = " << eps << " abs_resid = " << resid;
        if (time) out<<" time = " << time -> elapsed();
        out << std::endl;
    };
    auto printMaxIts = [](std::ostream& out, int it, int maxits, double eps, double resid, World3d::Timer* time){
        out << "Algorithm is reached maximum of iterations (" << maxits << "): it = " << it << " rel_resid = " << eps << " abs_resid = " << resid;
        if (time) out<<" time = " << time -> elapsed();
        out << std::endl;
    };
    if (m_external_call) m_external_call(info, w);
    if (info.status < 0){
        if (use_cout) printStepSimulationError(std::cout, info, p_time);
        if (p_lgfile) printStepSimulationError(*p_lgfile, info, p_time);
        last_ssi = info;
        return true;
    }
    if (info.it == 1) {
        m_resid_init = w->getWorldNorm("v:force");
        if (m_resid_init == 0) {
            if (use_cout) printConverge(std::cout, info.it, 0, m_resid_init, p_time);
            if (p_lgfile) printConverge(*p_lgfile, info.it, 0, m_resid_init, p_time);
            info.status = 0;
            last_ssi = info;
            return true;
        }
    }
    if (info.it > 0 && (info.it % (p_freq ? *p_freq : 1) == 0 || info.it >= (p_maxits ? *p_maxits : INT_MAX))) {
        double resid = w->getWorldNorm("v:force");
        double eps = resid / m_resid_init;

        if (use_cout) printCurrentState(std::cout, info.it, eps, resid, p_time);
        if (p_lgfile) printCurrentState(*p_lgfile, info.it, eps, resid, p_time);
        if ((p_err && eps < *p_err) || (p_abs_err && resid < *p_abs_err)) {
            if (use_cout) printConverge(std::cout, info.it, eps, resid, p_time);
            if (p_lgfile) printConverge(*p_lgfile, info.it, eps, resid, p_time);
            info.status = 0;
            last_ssi = info;
            return true;
        }
        if (std::isnan(eps) || (p_diverge_lim && eps >= *p_diverge_lim)) {
            if (use_cout) printDiverge(std::cout, info.it, eps, resid, p_time);
            if (p_lgfile) printDiverge(*p_lgfile, info.it, eps, resid, p_time);
            info.status = -10;
            last_ssi = info;
            return true;
        }
        if (info.it >= (p_maxits ? *p_maxits : INT_MAX)) {
            if (use_cout) printMaxIts(std::cout, info.it, (p_maxits ? *p_maxits : INT_MAX), eps, resid, p_time);
            if (p_lgfile) printMaxIts(*p_lgfile, info.it, (p_maxits ? *p_maxits : INT_MAX), eps, resid, p_time);
            info.status = 2;
            last_ssi = info;
            return true;
        }
    }
    last_ssi = info;
    return false;
}
