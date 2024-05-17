#ifndef AORTIC_VALVE_COLLISSIONTHICKESS_H
#define AORTIC_VALVE_COLLISSIONTHICKESS_H

#include "../CollisionManagerInterface.h"

namespace World3d {
    ///When calculating self-collisions of the body, it is required to estimate the distance between the parts of the body 
    ///at which these parts begin to come into contact. When using a mesh with a step less than the thickness of the body, 
    ///as such a distance, it is necessary to use not the thickness of the body, but some smaller value.
    ///This class provides methods to evaluate effective collision distance between initially close parts of the mesh
    struct CollisionThicknessBase{
        ///Types of mesh primitives
        static const int VERTEX = 0;
        static const int EDGE = 1;
        static const int FACE = 2;

        Object3D* m_obj;
        ///@param id identification number of the primitive 
        ///@param mesh_primitive_type one of values VERTEX, EDGE, FACE
        ///@return thickness on the mesh primitive for collision with different objects
        virtual DReal operator()(int id, int mesh_primitive_type) = 0;
        ///@param id0,id1 identification numbers of primitives of object 
        ///@param mesh_primitive_type0,mesh_primitive_type1 one of values VERTEX, EDGE, FACE
        ///@return effective self-collision thickness between parts of the object
        virtual DReal operator()(int id0, int id1, int mesh_primitive_type0, int mesh_primitive_type1) = 0;
        ///Set some internal data to speed-up querying for specified types of mesh parts. Should be called after registerObj(...)
        virtual void setupIndex(int type0, int type1) = 0;
        virtual void registerObj(Object3D* obj){ m_obj = obj; }
        virtual void removeObj(){ m_obj = nullptr; }
        virtual std::shared_ptr<CollisionThicknessBase> copyEmpty() = 0;
        virtual ~CollisionThicknessBase(){ removeObj(); }
    };

    struct CollisionThicknessConst: public CollisionThicknessBase{
        DReal thickness = 0.5;

        CollisionThicknessConst() = default;
        CollisionThicknessConst(DReal thickness): thickness{thickness} {}
        CollisionThicknessConst& setThickness(DReal thickness){ CollisionThicknessConst::thickness = thickness; return *this; }
        DReal operator()(int id, int mesh_primitive_type) override { return thickness; }
        DReal operator()(int id0, int id1, int mesh_primitive_type0, int mesh_primitive_type1) override { return thickness; }
        void setupIndex(int type0, int type1) override {}
        std::shared_ptr<CollisionThicknessBase> copyEmpty() override{ return std::make_shared<CollisionThicknessConst>(); }
    };
    ///Class to give effective self-collision thickness between parts of the object
    template<class ThicknessT, typename std::enable_if<std::is_base_of<CollisionThicknessBase, ThicknessT>::value>::type* = nullptr>
    struct CollisionThicknessCompressed: public ThicknessT{
    protected:
        using CTB = CollisionThicknessBase;
        std::array<const void*, 3*(3+1)/2> m_distance_map_matrix = {nullptr};
        DReal max_compression = 3;
        std::string m_init_x_tag_name = "v:point";
    public:
        CollisionThicknessCompressed() = default;
        CollisionThicknessCompressed(ThicknessT t): ThicknessT(t) {}
        CollisionThicknessCompressed& setMaxCompression(DReal compress) { max_compression = compress > 1 ? compress : 1; return *this; }
        CollisionThicknessCompressed& setInitCoordTagName(std::string name){ return m_init_x_tag_name = name, *this; }
        DReal operator()(int id0, int mesh_primitive_type0) override{ return ThicknessT::operator()(id0, mesh_primitive_type0); }
        DReal operator()(int id0, int id1, int mesh_primitive_type0, int mesh_primitive_type1) override;
        void setupIndex(int type0, int type1) override;
        void registerObj(Object3D* obj) override;
        void removeObj() override;
        std::shared_ptr<CollisionThicknessBase> copyEmpty() override{ return std::make_shared<CollisionThicknessCompressed<ThicknessT>>(); }
    };

    struct CollisionThickness{
        std::shared_ptr<CollisionThicknessBase> m_ptr;

        CollisionThickness() {}
        template<class ThicknessT, typename std::enable_if<std::is_base_of<CollisionThicknessBase, ThicknessT>::value>::type* = nullptr>
        CollisionThickness(ThicknessT f): m_ptr{std::make_shared<ThicknessT>(f)} {}
        template<class ThicknessT, typename std::enable_if<std::is_base_of<CollisionThicknessBase, ThicknessT>::value>::type* = nullptr>
        CollisionThickness& operator=(ThicknessT f){ m_ptr = std::make_shared<ThicknessT>(f); return *this; }
        
        DReal operator()(int id, int mesh_primitive_type){ return m_ptr->operator()(id, mesh_primitive_type); }
        DReal operator()(int id0, int id1, int mesh_primitive_type0, int mesh_primitive_type1){
            return m_ptr->operator()(id0, id1, mesh_primitive_type0, mesh_primitive_type1);
        }
        void setupIndex(int type0, int type1) { return m_ptr->setupIndex(type0, type1); }
        void registerObj(Object3D* obj){ return m_ptr->registerObj(obj); }
        void removeObj(){ return m_ptr->removeObj(); }
    };
}

#include "CollisionThickness.inl"

#endif