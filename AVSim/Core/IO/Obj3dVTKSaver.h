//
// Created by alex on 29.04.2021.
//

#ifndef AORTIC_VALVE_OBJ3DVTKSAVER_H
#define AORTIC_VALVE_OBJ3DVTKSAVER_H

#include "../Object3D.h"

namespace World3d{
struct VTKSaver{
    enum Sparsity{
        NODE = 1,
        EDGE = 2,
        FACE = 3,
        HALFEDGE = 4,
    };
    enum DType{
        REAL_TP,
        INT_TP
    };
    struct TagWrapper{
    private:
        std::vector<char> data;
        std::function<void (int, void*)> compute;
    public:
        std::string name;
        Sparsity sp;
        DType tp;
        int dim;

        void Init(DType _tp, int _dim, Sparsity _sp, std::function<void (int, void*)> _compute, std::string _name){
            tp = _tp, dim = _dim, sp = _sp, compute = std::move(_compute), name = std::move(_name);
            switch (tp) {
                case REAL_TP: data.resize(_dim * sizeof(double)); break;
                case INT_TP: data.resize(_dim * sizeof(int)); break;
            }
        }
        template<typename TYPE>
        TYPE* evaluate(int id) { compute(id, data.data()); return reinterpret_cast<TYPE*>(data.data()); }
    };
    Object3D* m_obj = nullptr;
    std::function<std::array<double, 3>(V_ind)> m_x;
    std::vector<TagWrapper> m_tags;

    VTKSaver() = default;
    VTKSaver& setObj(const Object3D& obj){ return m_obj = const_cast<Object3D*>(&obj), *this; }
    VTKSaver& setX(const std::function<std::array<double, 3>(V_ind)>& x) { return m_x = x, *this; }
    VTKSaver& setX(const Mesh::Property_map<V_ind, Point>& x) { return m_x = [&x](V_ind v){ return std::array<double, 3>{x[v][0], x[v][1], x[v][2]}; }, *this; }

#define SETTAG(TYPE, DTYPE, SPARSITY, SP_TYPE) \
    template<int n> \
    VTKSaver& setTag(const std::function<std::array<TYPE, n>(SP_TYPE)>& x, std::string name) { \
        TagWrapper tw; \
        tw.Init(DTYPE, n, SPARSITY, [x](int v, void* d) { auto arr = x(SP_TYPE(v)); std::copy(arr.begin(), arr.end(), static_cast<TYPE*>(d)); }, name); \
        m_tags.push_back(tw);\
        return *this; \
    }

    SETTAG(double, REAL_TP, NODE, V_ind)
    SETTAG(double, REAL_TP, EDGE, E_ind)
    SETTAG(double, REAL_TP, FACE, F_ind)
    SETTAG(double, REAL_TP, HALFEDGE, HE_ind)
    SETTAG(int, INT_TP, NODE, V_ind)
    SETTAG(int, INT_TP, EDGE, E_ind)
    SETTAG(int, INT_TP, FACE, F_ind)
    SETTAG(int, INT_TP, HALFEDGE, HE_ind)

#undef SETTAG
#define SETSPECTAG(SP_TYPE)    \
    VTKSaver& setTag(const Mesh::Property_map<SP_TYPE, Point>& x, std::string name) { \
        return setTag<3>([x](SP_TYPE v){ return std::array<double, 3>{x[v][0], x[v][1], x[v][2]}; }, name);\
    }\
    VTKSaver& setTag(const Mesh::Property_map<SP_TYPE, Vector>& x, std::string name) {\
        return setTag<3>([x](SP_TYPE v){ return std::array<double, 3>{x[v][0], x[v][1], x[v][2]}; }, name);\
    }\
    VTKSaver& setTag(const Mesh::Property_map<SP_TYPE, double>& x, std::string name) {\
        return setTag<1>([x](SP_TYPE v){ return std::array<double, 1>{x[v]}; }, name);\
    }\
    VTKSaver& setTag(const Mesh::Property_map<SP_TYPE, int>& x, std::string name) {\
        return setTag<1>([x](SP_TYPE v){ return std::array<int, 1>{x[v]}; }, name);\
    }\
    template<int n>\
    VTKSaver& setTag(const Mesh::Property_map<SP_TYPE, std::array<double, n>>& x, std::string name) {\
        return setTag<n>([x](SP_TYPE v){ return x[v]; }, name);\
    }                          \
    template<int n>\
    VTKSaver& setTag(const Mesh::Property_map<SP_TYPE, std::array<int, n>>& x, std::string name) {\
        return setTag<n>([x](SP_TYPE v){ return x[v]; }, name);\
    }
    SETSPECTAG(V_ind)
    SETSPECTAG(E_ind)
    SETSPECTAG(F_ind)
    SETSPECTAG(HE_ind)
#undef SETSPECTAG

    bool setTagByName(std::string tag_name, std::string save_name = "") {
        if (!m_obj) return false;
        if (save_name.empty()) save_name = tag_name;
#define CHECKANDSETTAG(SP_TYPE, TYPE)\
        {   \
            auto it = m_obj->m_mesh.property_map<SP_TYPE, TYPE>(tag_name);\
            if (it.second) {\
                setTag(it.first, save_name);\
                return true;\
            }\
        }
#define CHECKANDSETTAGARR(SP_TYPE, TYPE, n)\
        {   \
            auto it = m_obj->m_mesh.property_map<SP_TYPE, std::array<TYPE, n>>(tag_name);\
            if (it.second) {\
                setTag<n>(it.first, save_name);\
                return true;\
            }\
        }
#define CHECKTAGONSPARSITY(SP_TYPE) \
        if (m_obj->m_mesh.property_type<SP_TYPE>(tag_name) != typeid(void)){\
        CHECKANDSETTAG(SP_TYPE, Point);\
        CHECKANDSETTAG(SP_TYPE, Vector);\
        CHECKANDSETTAG(SP_TYPE, double);\
        CHECKANDSETTAG(SP_TYPE, int);\
        CHECKANDSETTAGARR(SP_TYPE, double, 2);\
        CHECKANDSETTAGARR(SP_TYPE, double, 3);\
        CHECKANDSETTAGARR(SP_TYPE, double, 6);\
        CHECKANDSETTAGARR(SP_TYPE, double, 9);\
        CHECKANDSETTAGARR(SP_TYPE, int, 2);\
        CHECKANDSETTAGARR(SP_TYPE, int, 3);\
        CHECKANDSETTAGARR(SP_TYPE, int, 6);\
        CHECKANDSETTAGARR(SP_TYPE, int, 9);\
        }

        CHECKTAGONSPARSITY(V_ind);
        CHECKTAGONSPARSITY(E_ind);
        CHECKTAGONSPARSITY(F_ind);
        CHECKTAGONSPARSITY(HE_ind);
#undef CHECKTAGONSPARSITY
#undef CHECKANDSETTAGARR
#undef CHECKANDSETTAG

        return false;
    }

    int save(std::string file);
};

}


#endif //AORTIC_VALVE_OBJ3DVTKSAVER_H
