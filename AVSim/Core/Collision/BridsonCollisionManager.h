//
// Created by Liogky Alexey on 06.06.2022.
//

#ifndef AVSIM_BRIDSONCOLLISIONMANAGER_H
#define AVSIM_BRIDSONCOLLISIONMANAGER_H

#include "../CollisionManagerInterface.h"
#include "CollisionThickness.h"
#include <stack>

#define CHECKCOLLISIONREMOVE

#ifdef CHECKCOLLISIONREMOVE
    #define CHECKCOLON() m_make_check = true
#else
    #define CHECKCOLON
#endif

#ifdef CHECKCOLLISIONREMOVE
    #define CHECKCOLOFF() m_make_check = false
#else
    #define CHECKCOLOFF
#endif

#ifdef CHECKCOLLISIONREMOVE
    #define CHECKCOLDEF() bool m_make_check = false
#else
    #define CHECKCOLDEF
#endif 

#ifdef CHECKCOLLISIONREMOVE
    #define CHECKCOLIF(X) if (m_make_check) { X; }
    #define CHECKSTATE(X) X
#else
    #define CHECKCOLIF(X)
    #define CHECKSTATE(X)
#endif

#ifdef CHECKCOLLISIONREMOVE
    #define CHECKCOLIFBEGIN() if (m_make_check) {
#else
    #define CHECKCOLIFBEGIN() if (false) {
#endif

#ifdef CHECKCOLLISIONREMOVE
    #define CHECKCOLIFEND() }
#else
    #define CHECKCOLIFEND() }
#endif

namespace World3d {
    class BridsonCollisionManager: public CollisionManagerBase {
        using CTB = CollisionThicknessBase;
        struct DynamicObject3DCI{
            Object3D* m_obj;
            ObjectID m_id;
            Mesh::Property_map<V_ind, Point>  m_x;  ///<actual current position
            Mesh::Property_map<V_ind, Point>  m_next_x;
            Mesh::Property_map<V_ind, Vector> m_dx; ///<shift to next position
            Mesh::Property_map<V_ind, Vector> m_v;  ///<cur velocity
            Mesh::Property_map<V_ind, Vector> m_dv; ///<shift to next velocity
            Mesh::Property_map<V_ind, DReal>  m_m;  ///<node mass
            CollisionThickness m_thickness = CollisionThicknessConst(0.5);
            std::shared_ptr<VolumeTreeBase> m_f0bvt; ///<face (Tprev) volume tree
            std::shared_ptr<VolumeTreeBase> m_febvt; ///<face extended (Tprev + Tnext) volume tree
            struct FaceData{
                F_ind id;
                VolumeTreeBase::LeafIndex m_f0leaf;
                VolumeTreeBase::LeafIndex m_feleaf;
                DynamicObject3DCI* backptr;
            };
            std::vector<FaceData> fdat;
        };
//      to handle fixed objects
        struct StaticObject3DCI {
            Object3D* m_obj;
            ObjectID m_id;
            Mesh::Property_map<V_ind, Point> m_x; ///<actual position
            CollisionThickness m_thickness = CollisionThicknessConst(0.5);
        };
        ///To store values dependings on object
        template<typename ValT>
        struct ObjectValues{
            ValT def; 
            std::map<ObjectID, ValT> val;
            ObjectValues(ValT def = ValT()): def{def} {}
            
            ValT operator()() const { return def; }
            ValT operator()(ObjectID id0) const{
                auto it = val.find(id0);
                if (it != val.end()) return it->second;
                else return def;
            }
            ObjectValues& set_default(const ValT& v){ return def = v, *this; }
            ObjectValues& insert(const ValT& v, ObjectID id0){ return val[id0] = v, *this; }
        };
        ///To store values dependings on pair of objects
        template<typename ValT>
        struct ObjectInteractionValues{
            ValT def; 
            std::map<std::pair<ObjectID, ObjectID>, ValT> val;
            ObjectInteractionValues(ValT def = ValT()): def{def} {}
            
            ValT operator()() const { return def; }
            ValT operator()(ObjectID id0, ObjectID id1) const{
                if (id0 > id1) std::swap(id0, id1);
                auto it = val.find(std::pair<ObjectID, ObjectID>(id0, id1));
                if (it != val.end()) return it->second;
                else return def;
            }
            ObjectInteractionValues& set_default(const ValT& v){ return def = v, *this; }
            ObjectInteractionValues& insert(const ValT& v, ObjectID id0, ObjectID id1){
                if (id0 > id1) std::swap(id0, id1);
                val[std::pair<ObjectID, ObjectID>(id0, id1)] = v;
                return *this;
            }
        };
        std::vector<std::shared_ptr<DynamicObject3DCI>> m_dobjs;
        DReal m_dt; ///< time step
        ObjectInteractionValues<DReal> m_Kspr; ///< spring repulsion param
        ObjectInteractionValues<DReal> m_friction_mu;   ///< friction coefficient
        ObjectInteractionValues<DReal> m_RepulsSprCoef = ObjectInteractionValues<DReal>(0.1); ///< max repulsion shift relative to thickness per iteration
        ObjectValues<DReal> m_RelaxEps = ObjectValues<DReal>(0.15); ///<in [0, 1] - max changing of edge length per collision iteration
        int m_MaxRelaxIts = 3;
        int m_maxCCDIts = 5;
        int m_maxFailSafeIts = 2;
        ObjectInteractionValues<DReal>  m_ccdPenetrateDepthCoef = ObjectInteractionValues<DReal>(0.4); //in [0, 1]
        VTFactry m_vtfactory;
        int m_collision_status = SelfCollision | DynamicVsDynamicCollission | DynamicVsStaticCollission;
    public:
        static const int NoSolveCollision = 0x0;
        static const int SelfCollision = 0x1;
        static const int DynamicVsDynamicCollission = 0x2;
        static const int DynamicVsStaticCollission = 0x4;

        CHECKCOLDEF();
        BridsonCollisionManager();
        void addObj(Object3D& obj, ObjectID id, bool is_dynamic, double thickness) override;
        void set_thickness(ObjectID id, CollisionThickness thickness){ for (auto& o: m_dobjs) if (o->m_id == id) { (o->m_thickness = thickness).registerObj(o->m_obj); return; }}
        void set_thickness(ObjectID id, DReal thickness){ return set_thickness(id, CollisionThickness(CollisionThicknessConst(thickness))); }
        void set_Kspring(DReal k) { m_Kspr.set_default(k); }
        void set_Kspring(DReal k, ObjectID id0, ObjectID id1) { m_Kspr.insert(k, id0, id1); }
        void set_FrictionCoef(DReal mu = 0.0) { m_friction_mu.set_default(mu); }
        void set_FrictionCoef(DReal mu, ObjectID id0, ObjectID id1) { m_friction_mu.insert(mu, id0, id1); }
        void set_RepulsSprCoef(DReal eps = 0.1) { m_RepulsSprCoef.set_default(eps); }
        void set_RepulsSprCoef(DReal eps, ObjectID id0, ObjectID id1) { m_RepulsSprCoef.insert(eps, id0, id1); }
        void set_Dt(DReal dt) { m_dt = dt; }
        void setVTFactory(VTFactry vf) { m_vtfactory = std::move(vf); }
        void set_RelaxEps(DReal eps = 0.15) { m_RelaxEps.set_default(eps); }
        void set_RelaxEps(DReal eps, ObjectID id) { m_RelaxEps.insert(eps, id); }
        void set_PenetrateDepthCoef(DReal coef = 0.4) { m_ccdPenetrateDepthCoef.set_default(coef); }
        void set_PenetrateDepthCoef(DReal coef, ObjectID id0, ObjectID id1) { m_ccdPenetrateDepthCoef.insert(coef, id0, id1); }
        void set_MaxRelaxIts(int its = 3){ m_MaxRelaxIts = its; }
        void set_MaxCCDIts(int its = 5){ m_maxCCDIts = its; }
        void set_MaxFailSafeIts(int its = 2) { m_maxFailSafeIts = its; }
        void setCollisionStatus(int status = SelfCollision | DynamicVsDynamicCollission | DynamicVsStaticCollission){ m_collision_status = status; }
        void removeObj(ObjectID id) override;
        void findCollisions() override;
        void solveCollisions() override;
        void updateCCD();
        void update();
        void setAsInitialNonCollidingState();
        void applyDx();
        struct SpringAndFrictionParams{
            double h;
            double k_spr;
            double eps = 0.1;
            double mu = 0;
            void _init_data_from_BridsonCollisionManager(BridsonCollisionManager* self, ObjectID id0, ObjectID id1);
        };
        friend class SpringAndFrictionParams;
        struct CollisionVF{
            DynamicObject3DCI *tobj, *pobj;
            F_ind ft;
            std::array<V_ind, 3> vt;
            V_ind vp;
            Point prj;
            double sqd;
            std::array<double, 3> wt;
            Vector n;
            SpringAndFrictionParams sfp;
        };
        struct CollisionEE{
            std::array<DynamicObject3DCI*, 2> obj;
            std::array<std::array<V_ind, 2>, 2> v;
            std::array<Point, 2> prj;
            std::array<double, 2> w;
            double sqd;
            Vector n;
            SpringAndFrictionParams sfp;
        };
        void redistributeImpulseV(Vector impulse, V_ind v, DynamicObject3DCI *obj);
        void redistributeImpulseE(Vector impulse, V_ind* v/*[2]*/, double* w/*[2]*/, DynamicObject3DCI *obj);
        void redistributeImpulseT(Vector impulse, V_ind* v/*[3]*/, double* w/*[3]*/, DynamicObject3DCI *obj);
        double computedx_nVF(CollisionVF& col);
        double computedx_nEE(CollisionEE& col);
        double applyInelasticRepulsionVF(CollisionVF& col);
        double applyInelasticRepulsionEE(CollisionEE& col);
        double applySpringRepulsionVF(CollisionVF& col);
        double applySpringRepulsionEE(CollisionEE& col);
        void applyFrictionRepulsionVF(CollisionVF& col, double ddxn);
        void applyFrictionRepulsionEE(CollisionEE& col, double ddxn);
        void applyRepulsionVF(CollisionVF& col);
        void applyRepulsionEE(CollisionEE& col);

        void applyRepulsionImpulses();
        int applyCCDRepulsion();
        bool applyRelaxRepulsion(int maxRelaxIts = 1);
        bool applyCCDRepulsions(int maxCCDIts = 1);
        struct ImpactZone{
            enum ColType{
                VertexFace = 1,
                EdgeEdge = 2
            };
            struct ImpactColInfo{
                std::array<V_ind, 4> v;
                std::array<DynamicObject3DCI*, 2> obj;
                ColType type;
            };
            std::vector<std::pair<DynamicObject3DCI*, V_ind>> iz;
            std::vector<Vector> init_dx;
            std::vector<ImpactColInfo> collisions;
            void reset_dx() { for (int i = 0; i < iz.size(); ++i) iz[i].first->m_dx[iz[i].second] = init_dx[i]; }
            void set_dx() { init_dx.resize(iz.size()); for (int i = 0; i < iz.size(); ++i) init_dx[i] = iz[i].first->m_dx[iz[i].second];}
        };
        void applyFailSafe(int maxFailSafeIts = 1);
        int collectImpactZones(std::vector<ImpactZone>& izs);
        void handleImpactZone(ImpactZone& iz);

    private:
        struct CollisionStatistics{
            int nVFBroadCases = 0;
            int nEEBroadCases = 0;
            int nVFNarrowCases = 0;
            int nEENarrowCases = 0;
        };

        template<typename NarrowPhase>
        struct DefaultOpts { 
            static DReal GetColTime(NarrowPhase& h){ return 1; }
        };
        template<typename NarrowPhase, typename Handler, typename Opt = DefaultOpts<NarrowPhase>>
        int _performBroadPhaseSelf(VolumeTreeBase& a);

        template<typename NarrowPhase, typename Handler, typename Opt = DefaultOpts<NarrowPhase>>
        int _performBroadPhaseSelfExt(VolumeTreeBase& a, NarrowPhase& np, Handler& handler);

        template<typename NarrowPhase, typename Handler, typename Opt = DefaultOpts<NarrowPhase>>
        int _performBroadPhaseDiff(VolumeTreeBase& a, VolumeTreeBase& b);

        template<typename NarrowPhase, typename Handler, typename Opt = DefaultOpts<NarrowPhase>>
        int _performBroadPhaseDiffExt(VolumeTreeBase& a, VolumeTreeBase& b, NarrowPhase& np, Handler& handler);

        template<typename NarrowPhase, typename Handler, typename Opt = DefaultOpts<NarrowPhase>>
        int _performBroadPhase(VolumeTreeBase& a, VolumeTreeBase& b);

        template<typename NarrowPhase, typename Handler, typename Opt = DefaultOpts<NarrowPhase>>
        int _performBroadPhaseExt(VolumeTreeBase& a, VolumeTreeBase& b, NarrowPhase& np, Handler& handler);
    };
}
#include "BridsonCollisionManager.inl"


#endif //AVSIM_BRIDSONCOLLISIONMANAGER_H
