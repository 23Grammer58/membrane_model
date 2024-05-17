//
// Created by alex on 07.07.2020.
//

#ifndef AORTIC_VALVE_BULLETCOLLISIONALGORITHMS_H
#define AORTIC_VALVE_BULLETCOLLISIONALGORITHMS_H

#include "BulletCollision/BroadphaseCollision/btCollisionAlgorithm.h"
#include "BulletCollision/BroadphaseCollision/btDispatcher.h"
#include "BulletCollision/CollisionDispatch/btCollisionObjectWrapper.h"
#include "BulletCollision/CollisionShapes/btCollisionShape.h"
#include "BulletCollision/CollisionShapes/btConcaveShape.h"
#include "BulletCollision/CollisionDispatch/btCollisionObject.h"
#include "BulletCollision/BroadphaseCollision/btBroadphaseInterface.h"
#include "BulletCollision/CollisionShapes/btTriangleCallback.h"
#include "BulletCollision/NarrowPhaseCollision/btPersistentManifold.h"
class btDispatcher;
#include "BulletCollision/BroadphaseCollision/btBroadphaseProxy.h"
#include "BulletCollision/CollisionDispatch/btCollisionCreateFunc.h"
#include "BulletCollisionManager.h"
class btCollisionShape;

#include "LinearMath/btHashMap.h"

#include "BulletCollision/BroadphaseCollision/btQuantizedBvh.h"
using namespace World3d;

struct btTriIndex
{
    int m_PartIdTriangleIndex;
    class btCollisionShape* m_childShape;

    btTriIndex(int partId, int triangleIndex, btCollisionShape* shape)
    {
        m_PartIdTriangleIndex = (partId << (31 - MAX_NUM_PARTS_IN_BITS)) | triangleIndex;
        m_childShape = shape;
    }

    int getTriangleIndex() const
    {
        // Get only the lower bits where the triangle index is stored
        unsigned int x = 0;
        unsigned int y = (~(x & 0)) << (31 - MAX_NUM_PARTS_IN_BITS);
        return (m_PartIdTriangleIndex & ~(y));
    }
    int getPartId() const
    {
        // Get only the highest bits where the part index is stored
        return (m_PartIdTriangleIndex >> (31 - MAX_NUM_PARTS_IN_BITS));
    }
    int getUid() const
    {
        return m_PartIdTriangleIndex;
    }
};

///For each triangle in the concave mesh that overlaps with the AABB of a soft body (m_softBody), processTriangle is called.
class NetTriangleCallback : public btTriangleCallback
{
    BulletObject* m_softBody;
    const btCollisionObject* m_triBody;

    btVector3 m_aabbMin;
    btVector3 m_aabbMax;

    btManifoldResult* m_resultOut;

    btDispatcher* m_dispatcher;
    const btDispatcherInfo* m_dispatchInfoPtr;
    btScalar m_collisionMarginTriangle;

    btHashMap<btHashKey<btTriIndex>, btTriIndex> m_shapeCache;

public:
    int m_triangleCount;

    NetTriangleCallback(btDispatcher* dispatcher, const btCollisionObjectWrapper* body0Wrap, const btCollisionObjectWrapper* body1Wrap, bool isSwapped): m_dispatcher(dispatcher),
                                                                                                                                                         m_dispatchInfoPtr(nullptr)
    {
        m_softBody = static_cast<BulletObject*>(const_cast<btCollisionObject*>(isSwapped ? body1Wrap->getCollisionObject() : body0Wrap->getCollisionObject()));
        m_triBody = isSwapped ? body0Wrap->getCollisionObject() : body1Wrap->getCollisionObject();

        clearCache();
    }

    void setTimeStepAndCounters(btScalar collisionMarginTriangle, const btCollisionObjectWrapper* triBodyWrap, const btDispatcherInfo& dispatchInfo, btManifoldResult* resultOut){
        m_dispatchInfoPtr = &dispatchInfo;
        m_collisionMarginTriangle = collisionMarginTriangle;
        m_resultOut = resultOut;

        btVector3 aabbWorldSpaceMin, aabbWorldSpaceMax;
        m_softBody->getAabb(aabbWorldSpaceMin, aabbWorldSpaceMax);
        btVector3 halfExtents = (aabbWorldSpaceMax - aabbWorldSpaceMin) * btScalar(0.5);
        btVector3 softBodyCenter = (aabbWorldSpaceMax + aabbWorldSpaceMin) * btScalar(0.5);

        btTransform softTransform;
        softTransform.setIdentity();
        softTransform.setOrigin(softBodyCenter);

        btTransform convexInTriangleSpace;
        convexInTriangleSpace = triBodyWrap->getWorldTransform().inverse() * softTransform;
        btTransformAabb(halfExtents, m_collisionMarginTriangle, convexInTriangleSpace, m_aabbMin, m_aabbMax);
    }

    virtual ~NetTriangleCallback() { clearCache(); }

    virtual void processTriangle(btVector3* triangle, int partId, int triangleIndex) {
        (void) partId, (void) triangleIndex;
        m_softBody->CollisionHandler(m_triBody, triangle);
    }


        void clearCache()
    {
        for (int i = 0; i < m_shapeCache.size(); i++)
        {
            btTriIndex* tmp = m_shapeCache.getAtIndex(i);
            btAssert(tmp);
            btAssert(tmp->m_childShape);
            //m_softBody->m_sparsesdf.RemoveReferences(tmp->m_childShape);  //necessary?
            delete tmp->m_childShape;
        }
        m_shapeCache.clear();
    }

    SIMD_FORCE_INLINE const btVector3& getAabbMin() const
    {
        return m_aabbMin;
    }
    SIMD_FORCE_INLINE const btVector3& getAabbMax() const
    {
        return m_aabbMax;
    }
};



class NetConcaveCollisionAlgorithm : public btCollisionAlgorithm
{
    bool m_isSwapped;

    NetTriangleCallback m_NetTriangleCallback;

public:
    NetConcaveCollisionAlgorithm(const btCollisionAlgorithmConstructionInfo& ci, const btCollisionObjectWrapper* body0Wrap, const btCollisionObjectWrapper* body1Wrap, bool isSwapped)
            : btCollisionAlgorithm(ci),
              m_isSwapped(isSwapped),
              m_NetTriangleCallback(ci.m_dispatcher1, body0Wrap, body1Wrap, isSwapped)
    {}

    virtual ~NetConcaveCollisionAlgorithm() {};

    virtual void processCollision(const btCollisionObjectWrapper* body0Wrap, const btCollisionObjectWrapper* body1Wrap, const btDispatcherInfo& dispatchInfo, btManifoldResult* resultOut){
        //btCollisionObject* convexBody = m_isSwapped ? body1 : body0;
        const btCollisionObjectWrapper* triBody = m_isSwapped ? body0Wrap : body1Wrap;

        if (triBody->getCollisionShape()->isConcave())
        {
            const btCollisionObject* triOb = triBody->getCollisionObject();
            const btConcaveShape* concaveShape = static_cast<const btConcaveShape*>(triOb->getCollisionShape());

            //	if (convexBody->getCollisionShape()->isConvex())
            {
                btScalar collisionMarginTriangle = concaveShape->getMargin();

                //			resultOut->setPersistentManifold(m_btSoftBodyTriangleCallback.m_manifoldPtr);
                m_NetTriangleCallback.setTimeStepAndCounters(collisionMarginTriangle, triBody, dispatchInfo, resultOut);

                concaveShape->processAllTriangles(&m_NetTriangleCallback, m_NetTriangleCallback.getAabbMin(), m_NetTriangleCallback.getAabbMax());

                //	resultOut->refreshContactPoints();
            }
        }
    }

    btScalar calculateTimeOfImpact(btCollisionObject* body0, btCollisionObject* body1, const btDispatcherInfo& dispatchInfo, btManifoldResult* resultOut)
    {
        (void) body0, (void) body1, (void) dispatchInfo, (void) resultOut;
        return 1.0f; //do Nothing
    }

    virtual void getAllContactManifolds(btManifoldArray& manifoldArray)
    {
        (void) manifoldArray;
        //we don't add any manifolds
    }

    void clearCache()
    {
        m_NetTriangleCallback.clearCache();
    }

    struct CreateFunc : public btCollisionAlgorithmCreateFunc
    {
        virtual btCollisionAlgorithm* CreateCollisionAlgorithm(btCollisionAlgorithmConstructionInfo& ci, const btCollisionObjectWrapper* body0Wrap, const btCollisionObjectWrapper* body1Wrap)
        {
            void* mem = btAlignedAlloc(static_cast<size_t>(sizeof(NetConcaveCollisionAlgorithm)), 16);//ci.m_dispatcher1->allocateCollisionAlgorithm(sizeof(NetConcaveCollisionAlgorithm));
            return new (mem) NetConcaveCollisionAlgorithm(ci, body0Wrap, body1Wrap, false);
        }
    };

    struct SwappedCreateFunc : public btCollisionAlgorithmCreateFunc
    {
        virtual btCollisionAlgorithm* CreateCollisionAlgorithm(btCollisionAlgorithmConstructionInfo& ci, const btCollisionObjectWrapper* body0Wrap, const btCollisionObjectWrapper* body1Wrap)
        {
            void* mem = btAlignedAlloc(static_cast<size_t>(sizeof(NetConcaveCollisionAlgorithm)), 16);
            //ci.m_dispatcher1->allocateCollisionAlgorithm(sizeof(NetConcaveCollisionAlgorithm));
            return new (mem) NetConcaveCollisionAlgorithm(ci, body0Wrap, body1Wrap, true);
        }
    };
};

#include "BulletCollision/BroadphaseCollision/btCollisionAlgorithm.h"
#include "BulletCollision/CollisionDispatch/btCollisionCreateFunc.h"
#include "BulletCollision/CollisionDispatch/btCollisionObjectWrapper.h"
#include "BulletCollision/BroadphaseCollision/btDispatcher.h"

class btPersistentManifold;

class NetNetCollisionAlgorithm : public btCollisionAlgorithm
{
    bool m_ownManifold;
    btPersistentManifold* m_manifoldPtr;

public:
    NetNetCollisionAlgorithm(const btCollisionAlgorithmConstructionInfo& ci)
            : btCollisionAlgorithm(ci) {}

    virtual void processCollision(const btCollisionObjectWrapper* body0Wrap, const btCollisionObjectWrapper* body1Wrap, const btDispatcherInfo& /*dispatchInfo*/, btManifoldResult* /*resultOut*/)
    {
        BulletObject* soft0 = static_cast<BulletObject*>(const_cast<btCollisionObject*>(body0Wrap->getCollisionObject()));
        BulletObject* soft1 = static_cast<BulletObject*>(const_cast<btCollisionObject*>(body1Wrap->getCollisionObject()));
        soft0->CollisionHandler(soft1);
    }

    virtual btScalar calculateTimeOfImpact(btCollisionObject* body0, btCollisionObject* body1, const btDispatcherInfo& dispatchInfo, btManifoldResult* resultOut)
    {
        (void) body0, (void) body1, (void) dispatchInfo, (void) resultOut;
        return 1.f; //do Nothing
    }

    virtual void getAllContactManifolds(btManifoldArray& manifoldArray)
    {
        if (m_manifoldPtr && m_ownManifold)
            manifoldArray.push_back(m_manifoldPtr);
    }

    NetNetCollisionAlgorithm(btPersistentManifold* mf, const btCollisionAlgorithmConstructionInfo& ci, const btCollisionObjectWrapper* body0Wrap, const btCollisionObjectWrapper* body1Wrap)
            :NetNetCollisionAlgorithm(ci) { (void) mf, (void) body0Wrap, (void)body1Wrap; }

    virtual ~NetNetCollisionAlgorithm() {}

    struct CreateFunc : public btCollisionAlgorithmCreateFunc
    {
        virtual btCollisionAlgorithm* CreateCollisionAlgorithm(btCollisionAlgorithmConstructionInfo& ci, const btCollisionObjectWrapper* body0Wrap, const btCollisionObjectWrapper* body1Wrap)
        {
            int bbsize = sizeof(NetNetCollisionAlgorithm);
            void* ptr = ci.m_dispatcher1->allocateCollisionAlgorithm(bbsize);
            return new (ptr) NetNetCollisionAlgorithm(nullptr, ci, body0Wrap, body1Wrap);
        }
    };
};

#endif //AORTIC_VALVE_BULLETCOLLISIONALGORITHMS_H
