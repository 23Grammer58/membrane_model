//
// Created by Liogky Alexey on 17.06.2022.
//

#ifndef AVSIM_CONTACTFORCE_H
#define AVSIM_CONTACTFORCE_H

#include "../CollisionManagerInterface.h"
#include "../World.h"
#include "../Collision/CollisionThickness.h"
namespace World3d {
    class ContactForce : public WorldForceBase{
        public:
        struct ObjectWrapper {
            Object3D *m_obj = nullptr;
            bool is_dynamic = true;
            CollisionThickness shape_thickness;
            void setVT(std::shared_ptr<VolumeTreeBase> vt) { m_fbvt = std::move(vt); }
            void setVertIndOff(unsigned off) { m_vert_ind_offset = off; }
            unsigned getVIndOff() const { return m_vert_ind_offset; }
            VolumeTreeBase& getVT() { return *m_fbvt; }
            void updateTree(bool enforce = false);
            struct FaceData{
                F_ind id;
                VolumeTreeBase::LeafIndex m_leaf;
                FaceData(F_ind id): id{id} {}
                FaceData(F_ind id, VolumeTreeBase::LeafIndex leaf): id{id}, m_leaf{leaf} {}
            };
        private:
            std::shared_ptr<VolumeTreeBase> m_fbvt;
            std::shared_ptr<std::vector<FaceData>> m_fd;
            unsigned m_vert_ind_offset = 0;
        };
        static const int NoneContacts = 0x0;
        static const int VertVert = 0x1;
        static const int VertEdge = 0x2;
        static const int VertFace = 0x4;
        static const int EdgeEdge = 0x8;
        static const int VertVertSelf = 0x10;
        static const int VertEdgeSelf = 0x20;
        static const int VertFaceSelf = 0x40;
        static const int EdgeEdgeSelf = 0x80;

        static const int OtherMask = 0xF;
        static const int SelfMask = 0xF0;
        static const int NumContactTypes = 8;

        struct Contact{
            struct Elem{
                ObjectWrapper* obw;
                int geom_elem_id;
                Point proj;
            };
            std::array<Elem, 2> elems;
            std::array<double, 4> w;
            Vector diff = CGAL::NULL_VECTOR, n = CGAL::NULL_VECTOR;
            DReal sqd = 0, d = 0;
            DReal loc_thickness = 0;
        };
        using DistFunc = std::function<double(int contact_type, const Contact& c)>;

        std::array<std::vector<Contact>, NumContactTypes> m_contacts;
        std::map<ObjectID, ObjectWrapper> m_obws;
        DistFunc m_F; ///<distance force function
        DistFunc m_dF; ///<derivative of distance force function
        World* m_w = nullptr;
        VTFactry m_vf;
        DReal m_DefThickParam = 0.1;
        int m_ContactType = VertFace | EdgeEdge | VertFaceSelf | EdgeEdgeSelf;

        ContactForce();
        ContactForce(int contact_type);
        ContactForce(int contact_type, DistFunc f, DistFunc df);
        bool registerWorld(World& w) override;
        bool addObj(Object3D& obj, const ObjectID& id, bool is_dynamic) override;
        bool setShapeThickness(const ObjectID& id, CollisionThickness thickness);
        void setVTFactory(VTFactry vf) { m_vf = std::move(vf); }
        void setContactType(int contact_type) { m_ContactType = contact_type; }
        void setDistanceFunc(DistFunc f, DistFunc df) { m_F = std::move(f), m_dF = std::move(df); }
        void removeObj(const ObjectID &id) override;
        int operator()() override;
        int fill_matrix(SparseMatrix *sm) override;
        std::unique_ptr<WorldForceBase> clone() override;
    private:
        struct CollideElem{
            ObjectWrapper* obw;
            F_ind id;
        };
        void updateProximityTree(){ for (auto& o: m_obws) o.second.updateTree(); }
        void performCulling();

        void setVVForce(Contact& c, bool is_self_inter = false);
        void setVEForce(Contact& c, bool is_self_inter = false);
        void setVFForce(Contact& c, bool is_self_inter = false);
        void setEEForce(Contact& c, bool is_self_inter = false);
        void setForce(Contact& c, int contact_type);
        void setAllForces();
        void setXXdForce(SparseMatrix *sm, Contact& c, int contact_type, int nel1, int nel2, const std::array<V_ind, 4>& vv);
        void setVVdForce(SparseMatrix *sm, Contact& c, bool is_self_inter = false);
        void setVEdForce(SparseMatrix *sm, Contact& c, bool is_self_inter = false);
        void setVFdForce(SparseMatrix *sm, Contact& c, bool is_self_inter = false);
        void setEEdForce(SparseMatrix *sm, Contact& c, bool is_self_inter = false);
        void setdForce(SparseMatrix *sm, Contact& c, int contact_type);
        void setAlldForces(SparseMatrix *sm);
        static bool computeVV(Contact& c, DReal h, ObjectWrapper* o0, V_ind v0, ObjectWrapper* o1, V_ind v1);
        static bool computeVE(Contact& c, DReal h, ObjectWrapper* o0, V_ind v0, ObjectWrapper* o1, E_ind e1);
        static bool computeVF(Contact& c, DReal h, ObjectWrapper* o0, V_ind v0, ObjectWrapper* o1, F_ind f1);
        static bool computeEE(Contact& c, DReal h, ObjectWrapper* o0, E_ind e0, ObjectWrapper* o1, E_ind e1);
        void collectVV(ObjectWrapper* o0, V_ind v0, ObjectWrapper* o1, V_ind v1);
        void collectVE(ObjectWrapper* o0, V_ind v0, ObjectWrapper* o1, E_ind e1);
        void collectVF(ObjectWrapper* o0, V_ind v0, ObjectWrapper* o1, F_ind f1);
        void collectEE(ObjectWrapper* o0, E_ind e0, ObjectWrapper* o1, E_ind e1);
        void collectVVS(ObjectWrapper* o, V_ind v0, V_ind v1);
        void collectVES(ObjectWrapper* o, V_ind v0, E_ind e1);
        void collectVFS(ObjectWrapper* o, V_ind v0, F_ind f1);
        void collectEES(ObjectWrapper* o, E_ind e0, E_ind e1);
        void checkAndCollectContact(const std::array<CollideElem, 2>& cp);
        void contactRemoveDuplicates();
        static DReal compute_mid_edge_len(const Object3D& o);
        static DReal eval_thickness(World3d::ContactForce::ObjectWrapper* o0, World3d::ContactForce::ObjectWrapper* o1,
                                    int id0, int id1, 
                                    int primitive_type0, int primitive_type1);
    };
}


#endif //AVSIM_CONTACTFORCE_H
