//
// Created by alex on 06.07.2020.
//

#ifndef AORTIC_VALVE_COLLISIONMANAGER_H
#define AORTIC_VALVE_COLLISIONMANAGER_H

#include "../CollisionManagerInterface.h"

#include "BulletCollision/CollisionShapes/btCollisionShape.h"
#include "BulletCollision/CollisionShapes/btTriangleIndexVertexArray.h"
#include "BulletCollision/CollisionShapes/btBvhTriangleMeshShape.h"
#include "BulletCollision/CollisionDispatch/btCollisionObject.h"
#include "BulletCollision/CollisionDispatch/btDefaultCollisionConfiguration.h"
#include "BulletCollision/CollisionDispatch/btCollisionWorld.h"
#include "BulletCollision/BroadphaseCollision/btAxisSweep3.h"

#include "BulletCollision/CollisionShapes/btCollisionShape.h"
#include "BulletCollision/CollisionShapes/btConcaveShape.h"
#include "LinearMath/btQuickprof.h"
#include "LinearMath/btPolarDecomposition.h"
#include "BulletCollision/BroadphaseCollision/btBroadphaseInterface.h"
#include "BulletCollision/CollisionDispatch/btCollisionDispatcher.h"
#include "BulletCollision/CollisionShapes/btConvexInternalShape.h"
#include "BulletCollision/NarrowPhaseCollision/btGjkEpa2.h"

namespace World3d {
    btVector3 p2Bt (const Vector & p);
    Vector Bt2p (const btVector3& p);
    btVector3 p2Bt(const Point &p);
    Point Bt2pp(const btVector3 &p);

    class BulletCollisionManager;

    class Bullet_Wrapper
    {

    public:
        class BulletCollisionObject: public btCollisionObject
        {
        public:
            struct Collider : public btDbvt::ICollide
            {
                BulletCollisionObject* psb[2];
                double mrg;

                virtual void Process (const btDbvtNode* lnode,
                                      const btDbvtNode* lface) override;
            };


            struct ColliderStatic : public btDbvt::ICollide
            {
                BulletCollisionObject* psb;
                btVector3 m_triangle[3];
                double mrg;

                virtual void Process (const btDbvtNode* leaf) override;
            };
            typedef btAlignedObjectArray<btDbvtNode*> tDbvtArray;


        public:
            Bullet_Wrapper* m_body;
            Object3D& m_net;
            tDbvtArray m_nleaf;
            tDbvtArray m_fleaf;
            btDbvt m_ndbvt;                    // Nodes tree
            btDbvt m_fdbvt;                    // Faces tree
            btVector3 m_bounds[2];
            char state = 1;

            BulletCollisionObject(Object3D& net, Bullet_Wrapper* backPtr);
            virtual ~BulletCollisionObject() {
                delete m_collisionShape;
            }
            void CollisionHandler(BulletCollisionObject* psb);
            void CollisionHandler(const btCollisionObject* pco, btVector3* triangle);
            double getMargin() { return m_body->getMargin(); }
            void updateCollisionInfo();
            virtual void getAabb(btVector3& aabbMin, btVector3& aabbMax) const
            {
                aabbMin = m_bounds[0];
                aabbMax = m_bounds[1];
            }
            void updateBounds();
        private:
            std::vector<std::pair<V_ind, Bullet_Wrapper*>> _m_data1;
            std::vector<std::pair<F_ind, Bullet_Wrapper*>> _m_data2;
        };

        struct SContact
        {
            V_ind m_node;         // Node
            F_ind m_face;         // Face
            btVector3 m_weights;      // Weigths
            btVector3 m_normal;       // Normal
            double m_margin;        // Margin
            double m_cfm[2];        // Constraint force mixing
            Bullet_Wrapper* m_nobstr;   // Node body
            Bullet_Wrapper* m_fobstr;   // Face body
        };

        struct RContact
        {
            using Node = V_ind;
            Node m_node;         // Node
            btVector3 m_proj;         // Projection
            btVector3 m_normal;       // Normal
            btVector3 m_torthogonal;      // Orthogonal vector of triangle codirected with normal
            double m_margin;        // Margin
            Bullet_Wrapper* m_obstr;    // Obstruction
        };

        typedef btAlignedObjectArray<SContact> tSContactArray;
        typedef btAlignedObjectArray<RContact> tRContactArray;
        typedef btAlignedObjectArray<BulletCollisionObject*> ObjArray;

    public:
        Bullet_Wrapper(Object3D& net, double margin = 0.1);
        ~Bullet_Wrapper();

        double getMargin() { return m_mrg; }
        void setMargin(double margin) { m_mrg = margin; }
        void addObjs2World(btCollisionWorld* wrld);
        void removeObjsFromWorld(btCollisionWorld* wrld);
        void updateCollisionInfo();
        static void solveSContacts(Bullet_Wrapper* body);
        static void solveRContacts(Bullet_Wrapper* body);

        tSContactArray m_scontacts;        // Soft contacts
        tRContactArray m_rcontacts;
        Object3D& m_net;
    private:
        ObjArray m_objs;
        double m_mrg = 0.1;

        friend class BulletCollisionManager;
    };

    using BulletObject = Bullet_Wrapper::BulletCollisionObject;

    class BulletCollisionShape : public btConcaveShape
    {
    public:
        BulletObject* m_body;

        BulletCollisionShape(BulletObject* backptr)
        {
            m_shapeType = SOFTBODY_SHAPE_PROXYTYPE;
            m_body = backptr;
        }
        virtual ~BulletCollisionShape() {}

        void processAllTriangles(btTriangleCallback* /*callback*/, const btVector3& /*aabbMin*/, const btVector3& /*aabbMax*/) const{/*not yet*/assert(0 && "processAllTriangles");}
        virtual void setLocalScaling(const btVector3& /*scaling*/) { /*na*/ }
        virtual const btVector3& getLocalScaling() const { static const btVector3 dummy(1, 1, 1); return dummy; }
        virtual void calculateLocalInertia(btScalar /*mass*/, btVector3& /*inertia*/) const{/*not yet*/ assert(0 && "calculateLocalInertia"); }
        virtual const char* getName() const { return "Bullet_Wrapper::BulletCollisionObject"; }
        virtual void getAabb(const btTransform& t, btVector3& aabbMin, btVector3& aabbMax) const
        {
            //m_body->getAabb(aabbMin, aabbMax);
            const btVector3 mins= m_body->m_bounds[0];
            const btVector3 maxs= m_body->m_bounds[1];
            const btVector3 crns[]={t*btVector3(mins.x(),mins.y(),mins.z()),
                                    t*btVector3(maxs.x(),mins.y(),mins.z()),
                                    t*btVector3(maxs.x(),maxs.y(),mins.z()),
                                    t*btVector3(mins.x(),maxs.y(),mins.z()),
                                    t*btVector3(mins.x(),mins.y(),maxs.z()),
                                    t*btVector3(maxs.x(),mins.y(),maxs.z()),
                                    t*btVector3(maxs.x(),maxs.y(),maxs.z()),
                                    t*btVector3(mins.x(),maxs.y(),maxs.z())};
            aabbMin=aabbMax=crns[0];
            for(int i=1;i<8;++i)
            {
                aabbMin.setMin(crns[i]);
                aabbMax.setMax(crns[i]);
            }
        }
    };

    template <typename T>
    static inline T Clamp(const T& x, const T& l, const T& h)
    {
        return (x < l ? l : x > h ? h : x);
    }

    static inline btVector3 ProjectOnAxis(const btVector3& v, //v supposed normalized
                                          const btVector3& a)
    {
        return (a * btDot(v, a));
    }
//
    static inline btVector3 ProjectOnPlane(const btVector3& v, //v supposed normalized
                                           const btVector3& a)
    {
        return (v - ProjectOnAxis(v, a));
    }

//
    static inline void ProjectOrigin(const btVector3& a,
                                     const btVector3& b,
                                     btVector3& prj,
                                     btScalar& sqd)
    {
        const btVector3 d = b - a;
        const btScalar m2 = d.length2();
        if (m2 > SIMD_EPSILON)
        {
            const btScalar t = Clamp<btScalar>(-btDot(a, d) / m2, 0, 1);
            const btVector3 p = a + d * t;
            const btScalar l2 = p.length2();
            if (l2 < sqd)
            {
                prj = p;
                sqd = l2;
            }
        }
        else {
            const btScalar l2 = a.length2();
            if (l2 < sqd)
            {
                prj = a;
                sqd = l2;
            }
        }
    }
//
//initial value of sqd set maximal distance of searching, if search is infinite then should be set as SIMD_INFINITY
    static inline void ProjectOrigin(const btVector3& a,
                                     const btVector3& b,
                                     const btVector3& c,
                                     btVector3& prj,
                                     btScalar& sqd)
    {
        const btVector3& q = btCross(b - a, c - a);
        const btScalar m2 = q.length2();
        if (m2 > SIMD_EPSILON)
        {
            const btVector3 n = q / btSqrt(m2);
            const btScalar k = btDot(a, n);
            const btScalar k2 = k * k;
            if (k2 < sqd)
            {
                const btVector3 p = n * k;
                if ((btDot(btCross(a - p, b - p), q) > 0) &&
                    (btDot(btCross(b - p, c - p), q) > 0) &&
                    (btDot(btCross(c - p, a - p), q) > 0))
                {
                    prj = p;
                    sqd = k2;
                }
                else
                {
                    ProjectOrigin(a, b, prj, sqd);
                    ProjectOrigin(b, c, prj, sqd);
                    ProjectOrigin(c, a, prj, sqd);
                }
            }
        }
        else
            ProjectOrigin(a, b, prj, sqd);
    }

    static inline void ProjectOnTriangle(const btVector3& a,
                                     const btVector3& b,
                                     const btVector3& c,
                                     const btVector3& o,
                                     btVector3& prj,
                                     btScalar& sqd)
    {
        sqd = SIMD_INFINITY;
        ProjectOrigin(a - o, b - o, c - o, prj, sqd);
        prj += o;
    }

//
    template <typename T, typename R>
    static inline T BaryEval(const T& a,
                             const T& b,
                             const T& c,
                             const R& coord)
    {
        return (a * coord.x() + b * coord.y() + c * coord.z());
    }
//
    static inline btVector3 BaryCoord(const btVector3& a,
                                      const btVector3& b,
                                      const btVector3& c,
                                      const btVector3& p)
    {
        const btScalar w[] = {btCross(a - p, b - p).length(),
                              btCross(b - p, c - p).length(),
                              btCross(c - p, a - p).length()};
        const btScalar isum = 1 / (w[0] + w[1] + w[2]);
        return (btVector3(w[1] * isum, w[2] * isum, w[0] * isum));
    }

    class MyCollisionConfiguration : public btDefaultCollisionConfiguration
    {
        btCollisionAlgorithmCreateFunc* m_softSoftCreateFunc;
        btCollisionAlgorithmCreateFunc* m_softBodyConcaveCreateFunc;
        btCollisionAlgorithmCreateFunc* m_swappedSoftBodyConcaveCreateFunc;
    public:
        MyCollisionConfiguration(const btDefaultCollisionConstructionInfo& constructionInfo = btDefaultCollisionConstructionInfo());
        virtual ~MyCollisionConfiguration();
        ///creation of soft-soft and soft-rigid, and otherwise fallback to base class implementation
        virtual btCollisionAlgorithmCreateFunc* getCollisionAlgorithmCreateFunc(int proxyType0, int proxyType1) override;
    };

    class BulletCollisionManager: public CollisionManagerBase {
    public:
//        struct CollisionTimers{
//            double* m_broad_phase = nullptr;
//            double* m_narrow_phase = nullptr;
//            double* m_process = nullptr;
//            CollisionTimers() {}
//            CollisionTimers(double* broad_phase, double* narrow_phase, double* process): m_broad_phase{broad_phase}, m_narrow_phase{narrow_phase}, m_process{process} {}
//            World3d::Timer m_timer;
//        };
        using SolveColAlgo = std::function<int(std::map<ObjectID, std::pair<void*, bool>>&)>;
        BulletCollisionManager(std::pair<Vector, Vector> min_max_Aabb = {{-100, -100, -100}, {100, 100, 100}},
                            int maxProxies = 1000){ registerCollisionWorld(min_max_Aabb, maxProxies); }
        virtual void addObj(Object3D& obj, ObjectID id, bool is_dynamic = true, double margin = 0.1) override ;
        virtual void removeObj(ObjectID id) override;
        virtual void findCollisions() override;
        virtual void solveCollisions() override;
        void setCustomSolveColAlgo(SolveColAlgo algo) { m_solve_collission = std::move(algo); }
        void set_margin(ObjectID id, double margin);
        void set_margin(double margin);
        double get_margin(ObjectID id);
    private:
        void registerCollisionWorld(std::pair<Vector, Vector> min_max_Aabb, int maxProxies);
        btCollisionObject* convert_Object3D_to_btTriangleMesh(const Object3D& st, const ObjectID id, double mrg);

        std::unique_ptr<btCollisionWorld> _collisionWorld;
        std::map<ObjectID, std::pair<void*, bool>> _collisions;
        //if true - Bullet_Wrapper* else btCollisionObject*
        std::unique_ptr<btBroadphaseInterface> _broadphase;
        std::unique_ptr<btCollisionConfiguration> _collisionConfig;
        std::unique_ptr<btDispatcher> _dispatcher;
        using StaticInfo = std::tuple<
                std::vector<btVector3>,
                std::vector<int>,
                std::unique_ptr<btTriangleIndexVertexArray>,
                std::unique_ptr<btBvhTriangleMeshShape>,
                std::unique_ptr<btCollisionObject>>;
        std::map<ObjectID, StaticInfo> _internal_data;
        SolveColAlgo m_solve_collission = [](std::map<ObjectID, std::pair<void*, bool>>& collisions){
            for (auto i: collisions){
                if (i.second.second){
                    Bullet_Wrapper* body = static_cast<Bullet_Wrapper*>(i.second.first);
                    Bullet_Wrapper::solveSContacts(body);
                    Bullet_Wrapper::solveRContacts(body);
                }
            }
            return 0;
        };
//        CollisionTimers m_ctm;
    };

    struct PenaltyCollisionAlgo{
        struct PenaltyContact{
            V_ind m_node;
            Object3D* m_obj;
            Point m_proj;
            Vector m_proj_dir;
            Vector m_normal;
            double m_margin;
        };
        static int default_penalty_force(const PenaltyContact& c);
        std::function<int(const PenaltyContact& )> set_force = [](const PenaltyContact& c){ return default_penalty_force(c); };
        int operator()(std::map<ObjectID, std::pair<void*, bool>>& collisions);
    };

}


#endif //AORTIC_VALVE_COLLISIONMANAGER_H
