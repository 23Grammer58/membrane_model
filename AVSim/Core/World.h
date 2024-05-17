//
// Created by alex on 30.06.2020.
//

#ifndef AORTIC_VALVE_WORLD_H
#define AORTIC_VALVE_WORLD_H
#include <utility>
#include <map>
#include <vector>
#include <algorithm>
#include <memory>

#include "Object3D.h"
#include "AVSim/Solvers/LinearSolver.h"
#include "RendererInterface.h"
#include "CollisionManagerInterface.h"
#include "ForceAppliersInterface.h"

namespace World3d {
    class World;
    typedef int ObjectLbl;
    typedef std::map<ObjectID, std::pair<Object3D, ObjectLbl>> ObjectContainer;
    typedef std::function<int(World*)> StepSimulationAlgo;

    struct SolverInternalData{
    public:
        LinearSolver solver;
        SparseMatrix sm;
        std::vector<double> rhs, x;
    };

    typedef std::unique_ptr<CollisionManagerBase> Collider;
    typedef std::unique_ptr<RendererBase> Renderer;

    struct StepSimInfo{
        int it = 0;
        int status = 0;
    };
    struct WorldForce_ID{
        int id = -1;
        WorldForce_ID(int id = -1): id{id} {}
    };

    class World {
        ObjectContainer _objects;
        Renderer _renderer;
        Collider _collisionManager;
        StepSimulationAlgo _stepAlgo;
        std::map<int, WorldForce> m_forces;

    public:
        World();
        virtual ~World();
        ObjectID addObject3D(Object3D&& obj, int status);
        Object3D& obj(ObjectID id);
        WorldForce& force(WorldForce_ID id);
        ObjectContainer& objs() { return _objects; }
        Object3D releaseObj(ObjectID oid);
        void removeObj(ObjectID oid);
        Force_ID addForce(Force f, ObjectID oid);
        void removeForce(ObjectID oid, Force_ID fid) { obj(oid).remove_force(fid); }
        void setForceApplier(ForceApplier f, ObjectID oid);
        void setForceApplier(ForceApplier f);
        std::map<ObjectID, Force_ID> addForce(Force f);
        WorldForce_ID addForce(WorldForce f);
        void removeForce(WorldForce_ID id) { m_forces.erase(id.id); }
        void setCollider(Collider&& collider);
        CollisionManagerBase* getCollider() { return _collisionManager.get(); }
        void setRenderer(Renderer&& renderer);
        void setStepSimulationAlgo(StepSimulationAlgo&& ssa) { _stepAlgo = std::move(ssa); }
        void setStepSimulationAlgo(const StepSimulationAlgo& ssa) { _stepAlgo = ssa; }
        void clearObjects();
        void clearRenderer();
        void clearCollisionManager();
        void clearStepAlgo();
        void clear();

        int UpdateForces();
        int ApplyForces();
        int ApplyCollision();
        int ApplyNext();
        int Render();

        int StepSimulation();

        int Simulation(std::function<bool(StepSimInfo&, World*)> stopCondition = [](StepSimInfo&, World*){ return false; });

        double getNorm(ObjectID id, std::string on_tag = "v:force", unsigned normOut = 2, unsigned normIn = 2);
        double getWorldNorm(std::string on_tag = "v:force", unsigned normObj = 2, unsigned normOut = 2, unsigned normIn = 2);
        double getShift(ObjectID id, unsigned normOut = 2, unsigned normIn = 2);
        double getWorldShift(unsigned normObj = 2, unsigned normOut = 2, unsigned normIn = 2);

        int compStateSize() const;
        int compStateVector(double *X) const;
        int setStateVector(const double* X, bool update_all_forces = true);
        int compResidual(double *R);
        int compJacobian(SparseMatrix *sm); 
        int compFResidual(double* R, ObjectID oid, Force_ID fid);
        int compWFResidual(double* R, WorldForce_ID wfid);
        int compFJacobian(SparseMatrix *sm, ObjectID oid, Force_ID fid);
        int compWFJacobian(SparseMatrix *sm, WorldForce_ID wfid);
        private: 
        bool full_residual_updated = false;
    };



    struct StaticStopCondition{
        double* p_err, *p_abs_err = nullptr;
        double* p_diverge_lim;
        int* p_maxits;
        int* p_freq;
        World3d::Timer* p_time;
        std::ostream* p_lgfile;
        bool use_cout;
        StepSimInfo last_ssi;
        std::function<bool(StepSimInfo& info, World *w)> m_external_call = nullptr;

        StaticStopCondition() = default;
        StaticStopCondition(const StaticStopCondition& ) = default;
        StaticStopCondition(StaticStopCondition&& ) noexcept = default;
        StaticStopCondition& operator=(const StaticStopCondition&) = default;
        StaticStopCondition& operator=(StaticStopCondition&&) = default;
        StaticStopCondition(double* err, double* diverge_lim, int* maxits, int* freq = nullptr, World3d::Timer* time = nullptr, std::ostream* lgfile = nullptr, bool use_cout = true):
                p_err(err), p_diverge_lim(diverge_lim), p_maxits(maxits), p_freq(freq), p_time(time), p_lgfile(lgfile), use_cout(use_cout) { }
        StaticStopCondition(double& err, double& diverge_lim, int& maxits, int* freq = nullptr, World3d::Timer* time = nullptr, std::ostream* lgfile = nullptr, bool use_cout = true):
                p_err(&err), p_diverge_lim(&diverge_lim), p_maxits(&maxits), p_freq(freq), p_time(time), p_lgfile(lgfile), use_cout(use_cout) { }
        StaticStopCondition& setErr(double& err){ return p_err = &err, *this; }
        StaticStopCondition& setAbsErr(double& abs_err){ return p_abs_err = &abs_err, *this; }
        StaticStopCondition& setDivergeLim(double& diverge_lim){ return p_diverge_lim = &diverge_lim, *this; }
        StaticStopCondition& setMaxits(int& maxits){ return p_maxits = &maxits, *this; }
        StaticStopCondition& setFreq(int& freq){ return p_freq = &freq, *this; }
        StaticStopCondition& setTime(World3d::Timer& time){ return p_time = &time, *this; }
        StaticStopCondition& setLogOut(std::ostream& lgfile){ return p_lgfile = &lgfile, *this; }
        StaticStopCondition& isUseCout(bool use){ return use_cout = use, *this; }
        StaticStopCondition& addExternalCall(std::function<bool(StepSimInfo& info, World *w)> external){ return m_external_call = external, *this; }
        bool operator()(StepSimInfo& info, World *w);
    private:
        double m_resid_init = -1;
    };

    struct NewtonSolverStepAlgo{
        int m_it = 0;
        std::function<int(NewtonSolverStepAlgo& ns, World* w)> newton_algo = [](NewtonSolverStepAlgo& ns, World* w){
                int stat = 0;
                if ((stat = ns.updateMeshData(w)) < 0) return stat;
                int N = ns.getNodfs(w);
                if ((stat = ns.assembleSystem(w, N)) < 0) return stat;
                if (!ns.Solve()) return (stat = -4);
                double tau = 1.0;
                if ((stat = ns.applySolution(w, tau)) < 0) return stat;
                ns.m_it++;
                return stat;
        };
        SolverInternalData _solver;


        int getNodfs(World3d::World *w){
            return w->compStateSize();
        }
        int updateMeshData(World3d::World *w){
            return w->UpdateForces();
        }
        int assembleSystem(World3d::World *w, int Nodfs){
            auto& objs = w->objs();
            auto& solver = _solver;
            solver.sm.clear();
            solver.sm.resize(Nodfs);
            solver.rhs.resize(Nodfs);

            int stat = 0;
            if ((stat = w->compResidual(solver.rhs.data())) < 0) return stat;
            if ((stat = w->compJacobian(&solver.sm)) < 0) return stat;

            return stat;
        }
        void assembleMatrix(World3d::World *w, int Nodfs){
            auto& objs = w->objs();
            auto& solver = _solver;
            solver.sm.clear();
            solver.sm.resize(Nodfs);

            w->compJacobian(&solver.sm);
        }
        void assembleRHS(World3d::World *w, int Nodfs){
            auto& objs = w->objs();
            auto& solver = _solver;
            solver.rhs.resize(Nodfs);

            w->compResidual(solver.rhs.data());
        }
        double computeRhsNorm(unsigned int N) { return computeRhsNorm(N, N); }
        double computeRhsNorm(unsigned int N1, unsigned int N2){
            auto& rhs = _solver.rhs;
            auto get_loc_norm = [N = N2](auto begin, auto end){
                if (N == 0)
                    return *std::max_element(begin, end);
                else if (N == 1){
                    double res = 0;
                    for (auto i = begin; i != end; ++i) res += fabs(*i);
                    return res;
                } else {
                    double res = 0;
                    for (auto i = begin; i != end; ++i) res += pow(fabs(*i), N);
                    return pow(res, 1.0/N);
                }
            };
            if (N1 == N2){
                return get_loc_norm(rhs.begin(), rhs.end());
            } else {
                auto get_glob_norm = [N = N1, get_loc_norm](std::vector<double>& v){
                    int K = v.size()/3;
                    double res = 0;
                    if (N == 0){
                        for (int i = 0; i < K; ++i) res = std::max(res, get_loc_norm(v.data()+3*i, v.data() + 3*(i+1)));
                    } else if (N == 1){
                        for (int i = 0; i < K; ++i) res += get_loc_norm(v.data()+3*i, v.data() + 3*(i+1));
                    } else {
                        for (int i = 0; i < K; ++i) res += pow(get_loc_norm(v.data()+3*i, v.data() + 3*(i+1)), N);
                        res = pow(res, 1.0/N);
                    }
                    return res;
                };
                return get_glob_norm(rhs);
            }
        }
        bool Solve(){
            auto& solver = _solver;
            if (!solver.solver.SetMatrix(solver.sm)) return false;
            std::fill(solver.x.begin(), solver.x.end(), 0.0);
            solver.x.resize(solver.rhs.size(), 0.0);
            return solver.solver.Solve(solver.rhs, solver.x);
        }
        int applySolution(World3d::World *w, double tau){
            int N = 0;
            auto& solver = _solver;
            auto& objs = w->objs();
            for (auto& i: objs)
                if (i.second.second) {
                    Object3D &obj = i.second.first;
                    for (auto i: obj.m_mesh.vertices()){
                        auto bc = obj.getBC(i);
                        if (obj._is_movable(obj, i)){
                            Vector dx((bc & 1) ? solver.x[N+3*i] : 0.0, (bc & 2) ? solver.x[N+3*i + 1] : 0.0, (bc & 4) ? solver.x[N+3*i + 2] : 0.0);
                            obj.m_next_x[i] = obj.m_x[i] - tau * dx;
                        }
                    }
                    N += 3 * i.second.first.m_mesh.num_vertices();
                }
            int stat = 0;
            if ((stat = w->ApplyCollision()) < 0) return stat;
            if ((stat = w->ApplyNext()) < 0) return stat;
            return stat;
        }
        int operator()(World* w){
            return newton_algo(*this, w);
        }
    };
    struct ChangableNewtonStepAlgo{
        using TauAlgo = std::function<double(int it, double rel_resid, double abs_resid, NewtonSolverStepAlgo& ns)>;
        using Self = ChangableNewtonStepAlgo;
        double resid_init = -1;
        double resid_now = -1;
        unsigned norm = 1;
        bool use_cout = true;
        std::ostream* logout = nullptr;
        TauAlgo tau_algo = [](int, double, double, NewtonSolverStepAlgo&) -> double{ return 1; };

        Self& setTauAlgo(TauAlgo ta) { return tau_algo = std::move(ta), *this; }
        Self& setLogOut(std::ostream& log) { return logout = &log, *this; }
        int operator()(NewtonSolverStepAlgo& ns, World* w){
            int stat = 0;
            if ((stat = ns.updateMeshData(w)) < 0) return stat;
            int N = ns.getNodfs(w);
            if ((stat = ns.assembleSystem(w, N)) < 0) return stat;
            resid_now = ns.computeRhsNorm(norm);
            if (ns.m_it == 0) resid_init = resid_now;
            double tau = tau_algo(ns.m_it, resid_now/resid_init, resid_now, ns);
            if (!ns.Solve()) {
                if (use_cout) std::cout << "Linear solver doesn't solve the system." << std::endl;
                if (logout) (*logout) << "Linear solver doesn't solve the system." << std::endl;
                return (stat = -4);
            }
            if ((stat = ns.applySolution(w, tau)) < 0) return stat;
            ns.m_it++;
            return stat;
        }
    };
}



#endif //AORTIC_VALVE_WORLD_H
