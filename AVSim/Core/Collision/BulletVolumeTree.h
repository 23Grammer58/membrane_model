//
// Created by Liogky Alexey on 31.05.2022.
//

#ifndef AVSIM_BULLETVOLUMETREE_H
#define AVSIM_BULLETVOLUMETREE_H

#include "../CollisionManagerInterface.h"
#include "BulletCollision/BroadphaseCollision/btDbvt.h"

namespace World3d {
    class BulletVolumeTree : public VolumeTreeBase {
    public:
        struct Collider : public btDbvt::ICollide {
            const VolumeTreeBase::CollideProcessor *cp;

            void Process(const btDbvtNode *lfirst,
                         const btDbvtNode *lsecond) override {
                LeafIndex l1, l2;
                l1.ptr = const_cast<btDbvtNode *>(lfirst), l2.ptr = const_cast<btDbvtNode *>(lsecond);

                cp->operator()(l1, lfirst->data, l2, lsecond->data);
            }
        };

        BulletVolumeTree() { drive_id = 1; }
        void clear() override { tree.clear(); }
        bool empty() const override { return tree.empty(); };
        LeafIndex insert(const BoxVol &box, void *data) override {
            auto c2p = [](auto p) { return btVector3(p[0], p[1], p[2]); };
            LeafIndex li;
            li.ptr = tree.insert(btDbvtVolume::FromMM(c2p(box.Mins()), c2p(box.Maxs())), data);
            return li;
        };
        void update(LeafIndex id, BoxVol &new_volume) override {
            auto c2p = [](auto p) { return btVector3(p[0], p[1], p[2]); };
            auto vol = btDbvtVolume::FromMM(c2p(new_volume.Mins()), c2p(new_volume.Maxs()));
            tree.update(static_cast<btDbvtNode *>(id.ptr), vol);
        };
        void remove(LeafIndex id) override { tree.remove(static_cast<btDbvtNode *>(id.ptr)); };
        void optimize() override { tree.optimizeIncremental(1); };
        [[nodiscard]] std::shared_ptr<VolumeTreeBase> copy_empty() const override { return std::make_shared<BulletVolumeTree>(); }

        btDbvt tree;
    private:
        virtual void collide(VolumeTreeBase &other, const VolumeTreeBase::CollideProcessor &cp) {
            Collider cl;
            cl.cp = &cp;
            tree.collideTT(tree.m_root, static_cast<BulletVolumeTree &>(other).tree.m_root, cl);
        };
    };
}

#endif //AVSIM_BULLETVOLUMETREE_H
