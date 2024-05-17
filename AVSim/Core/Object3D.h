//
// Created by alex on 13.06.2020.
//

#ifndef AORTIC_VALVE_OBJECT3D_H
#define AORTIC_VALVE_OBJECT3D_H
#include <functional>
#include <vector>
#include <stdexcept>
#include <fstream>
#include <algorithm>
#include "AVSim/Solvers/LinearSolverInterface.h"
#include "AVMesh.h"
#include "ForceInterface.h"

namespace World3d {
    class Object3D;
    typedef std::function<int(Object3D&)> Processor;
    typedef std::function<int(Object3D&)> Updater;
    template <typename CheckObj>
        using Checker = std::function<bool(Object3D&, const CheckObj&)>;
    using DirichletCond = std::function<unsigned char(Object3D&, const V_ind&)>;
    typedef int Updater_ID;
    typedef int Force_ID;
    //typedef int ObjectID;
    struct ObjectID{
        int id;
        ObjectID(int id = -1): id(id) {};
        operator int() const { return id; }
    };
    typedef Processor NextXUpdater;
    template<class Ind, class Type>
    struct ConstProperty: public std::pair<Mesh::Property_map<Ind, Type>, bool>{
        Type& operator[](const Ind& ind) { return this->first[ind]; }
        ConstProperty() {}
        ConstProperty(std::pair<Mesh::Property_map<Ind, Type>, bool>&& o):
        std::pair<Mesh::Property_map<Ind, Type>, bool>(o)
        {}
    };
    template<class Ind, class Type>
    struct UpdatableProperty: public std::pair<Mesh::Property_map<Ind, Type>, Updater_ID>{
        Type& operator[](Ind& ind) { return this->first[ind]; }
        UpdatableProperty() {}
        UpdatableProperty(std::pair<Mesh::Property_map<Ind, Type>, Updater_ID>&& o):
        std::pair<Mesh::Property_map<Ind, Type>, Updater_ID>(o)
        {}
    };

    class Object3D {
    public:
        std::string name = "";
        ObjectID __id;
        Mesh m_mesh;
        std::map<int, Force> m_forces;
        std::map<int, Processor> m_updaters;    //при копировании объекта UpdatableProperty не копируются!!!
        NextXUpdater _next_x_updater = [](Object3D& obj) -> int {
            throw std::runtime_error("For object " + obj.name + " was not setted NextXUpdater, but was called");
            return 0;
        };
        Processor m_apply_next_x = [](Object3D& obj) -> int{
            for (auto v: obj.m_mesh.vertices()) {
                obj.m_x[v] = obj.m_next_x[v];
            }
            return 0;
        };
        UpdatableProperty<F_ind, Vector> m_normal;
        Mesh::Property_map<V_ind, Vector> m_F;
        Mesh::Property_map<V_ind, Point> m_x0;
        Mesh::Property_map<V_ind, Point> m_x;
        Mesh::Property_map<V_ind, Point> m_next_x;
        Mesh::Property_map<V_ind, int> m_boundary;
        Processor* _renderer = nullptr;
        DirichletCond _bc = [](Object3D& obj, const V_ind& v) -> unsigned char{
            if (obj.m_boundary[v] == 0) return 7;
            else return 0;
        };
        bool _is_movable(Object3D& obj, const V_ind& v) { return _bc(obj, v); }

        Object3D() { init(); };
        Object3D(Mesh&& mesh): m_mesh(mesh) { init(); };
        Object3D(const Mesh& mesh): m_mesh(mesh) { init(); };
        Object3D(std::string filename) { init(); read(filename);  }
        Object3D(Object3D&& obj);
        Object3D(const Object3D& obj);
        Object3D& operator=(Object3D&& obj);
        Object3D& operator=(const Object3D& obj);

        void init();
        void reset_mesh();

        void setDirichletCond(DirichletCond bc) { _bc = bc; }
        unsigned char getBC(V_ind vi) { return _bc(*this, vi); }
        template<typename Vtype>
        Vtype withBCmask(V_ind vi, const Vtype& v) { auto bc = _bc(*this, vi); return Vtype((bc & 1) ? v[0] : 0, (bc & 2) ? v[1] : 0, (bc & 4) ? v[2] : 0); }
        void set_is_movable_checker(Checker<V_ind> is_movable) { _bc = [is_movable](Object3D& obj, const V_ind& v) -> unsigned char{ return (is_movable(obj, v)) ? 7 : 0; }; }
        bool is_movable(V_ind v) { return _is_movable(*this, v); }
        bool is_movable(const V_ind& v, int odf) { return _bc(*this, v) & (1 << odf); }
        void set_renderer(Processor& r) { _renderer = &r; }
        void save(std::string file, const Mesh::Property_map<V_ind, Point>& view_x, const Mesh::Property_map<F_ind, Vector>& view_n) const ;
        void save(std::string file, const Mesh::Property_map<V_ind, Point>& view_x) const;
        void save(std::string file) const;
        bool read(std::string file);
        Force_ID add_force(Force f);
        Updater_ID add_updater(Updater f);
        void set_NextXUpdater(NextXUpdater next_x_updater) {_next_x_updater = next_x_updater; }
        void remove_force(Force_ID force_id);
        void remove_updater(Updater_ID updater_id);
        int apply_forces();
        int apply_force(int fid);
        void apply_updaters();
        void apply_render();
        void apply_next_x();
        int update_next_x() { return _next_x_updater(*this); }
        double residual(unsigned normOut = 2, unsigned normIn = 2, std::string on_tag = "v:force");
        template<typename V, typename Iterable>
        double residual(const Iterable& v, unsigned normOut = 2, unsigned normIn = 2);
        ObjectID get_id() { return __id; }
        Object3D& setName(const std::string& name){ return this->name = name, *this;}
        void set_zero_residual();
    private:
        void repair_connectivity();
    };

    void save_mesh(std::string filename, Mesh& m);

    bool check_flat_initial_template(Object3D& obj, double tolerance = 1e-4);

    ConstProperty<E_ind, DReal> set_l0(Object3D& obj, bool forced = false);
    UpdatableProperty<E_ind, DReal> set_l(Object3D& obj);
    ConstProperty<F_ind, DReal> set_Ap(Object3D& obj, bool forced = false);
    ConstProperty<F_ind, Vector> set_n0(Object3D& obj, bool forced = false);
    UpdatableProperty<F_ind, Vector> set_S(Object3D& obj);
    ConstProperty<F_ind, std::array<Vector, 3>> set_D_vecs(Object3D& obj, bool forced = false);
    ConstProperty<F_ind, std::array<DReal, 6>> set_DD_matrix(Object3D& obj, bool forced = false);
}

#include "AVMeshConnectivityHelpers.h"

#define HAS_METHOD(POSTFIX, MNAME, CHECK) \
template<typename T, typename = void>\
struct __HasMethod_##POSTFIX: std::false_type\
{};\
template<typename T>\
struct __HasMethod_##POSTFIX<T, std::enable_if_t<CHECK<decltype(std::declval<T>().MNAME())>::value>>: std::true_type\
{};\
template<typename T>\
constexpr bool __HasMethod_v_##POSTFIX = __HasMethod_##POSTFIX<T>::value;

HAS_METHOD(dim, dimension, std::is_integral)
HAS_METHOD(size, size, std::is_integral)

template<typename T, typename = void>
struct __HasMethod_brack: std::false_type
{};
template<typename T>
struct __HasMethod_brack<T, std::enable_if_t<std::is_convertible<decltype(std::declval<T>().operator[](0)), double>::value>>: std::true_type
{};
template<typename T>
constexpr bool __HasMethod_v_brack = __HasMethod_brack<T>::value;
//HAS_METHOD(brack, operator[], std::is_floating_point)

template<typename V, typename Iterable>
static std::function<double(const Iterable& it)>
        get_norm_function(unsigned norm1, unsigned norm2){

    static auto NV = [](const V& v) -> size_t {
        if constexpr (__HasMethod_v_dim<V>) { return v.dimension(); }
        else if constexpr (__HasMethod_v_size<V>) { return v.size(); }
        else { return 1; }
    };
    static auto get_v = [](const V& v, int i) -> double{
        if constexpr (__HasMethod_v_brack<V>)
            return static_cast<double>(v[i]);
        else
            return static_cast<double>(v);
    };
    std::function<double(const Iterable&)> nrm;
    if (norm1 != norm2) {
        std::function<double(const Iterable&)>& nrm1 = nrm;
        std::function<double(const V &)> nrm2;
        switch (norm2) {
            case 0:
                nrm2 = [](const V &v) {
                    double res = 0;
                    for (int i = 0; i < NV(v); ++i)
                            res = std::max(res, fabs(get_v(v, i)));
                    return res;
                };
                break;
            case 1:
                nrm2 = [](const V &v) {
                    double res = 0;
                    for (int i = 0; i < NV(v); ++i)
                        res += fabs(get_v(v, i));
                    return res;
                };
                break;
            case 2:
                nrm2 = [](const V &v) {
                    double res = 0;
                    for (int i = 0; i < NV(v); ++i) {
                        res += get_v(v, i) * get_v(v, i);
                    }
                    return std::sqrt(res);
                };
                break;
            default:
                nrm2 = [norm2](const V &v) {
                    double res = 0;
                    for (int i = 0; i < NV(v); ++i)
                        res += pow(fabs(get_v(v, i)), norm2);
                    return pow(res, 1.0 / norm2);
                };
        }
        switch (norm1) {
            case 0:
                nrm1 = [&nrm2](const Iterable& v) {
                    double res = 0;
                    for (auto &i: v)
                        res = std::max(res, fabs(nrm2(i)));
                    return res;
                };
                break;
            case 1:
                nrm1 = [&nrm2](const Iterable& v) {
                    double res = 0;
                    for (auto &i: v)
                        res += fabs(nrm2(i));
                    return res;
                };
                break;
            case 2:
                nrm1 = [&nrm2](const Iterable& v) {
                    double res = 0;
                    for (auto &i: v) {
                        double q = fabs(nrm2(i));
                        res += q * q;
                    }
                    return sqrt(res);
                };
                break;
            default:
                nrm1 = [norm1, &nrm2](const Iterable& v) {
                    double res = 0;
                    for (auto &i: v) {
                        double q = nrm2(i);
                        res += pow(q, norm1);
                    }
                    return pow(res, 1.0 / norm1);
                };
                break;
        }
    } else {
        switch (norm1) {
            case 0:
                nrm = [](const Iterable& v) {
                    double res = 0;
                    for (auto &i: v)
                        for (int j = 0; j < NV(i); ++j)
                            res = std::max(res, (fabs(get_v(i, j))));
                    return res;
                };
                break;
            case 1:
                nrm = [](const Iterable& v) {
                    double res = 0;
                    for (auto &i: v)
                        for (int j = 0; j < NV(i); ++j)
                            res += fabs(get_v(i, j));
                    return res;
                };
                break;
            case 2:
                nrm = [](const Iterable& v) {
                    double res = 0;
                    for (auto &i: v) {
                        for (int j = 0; j < NV(i); ++j)
                            res += get_v(i, j) * get_v(i, j);
                    }
                    return sqrt(res);
                };
                break;
            default:
                nrm = [norm1](const Iterable& v) {
                    double res = 0;
                    for (auto &i: v) {
                        for (int j = 0; j < NV(i); ++j)
                            res += pow(fabs(get_v(i, j)), norm1);
                    }
                    return pow(res, 1.0 / norm1);
                };
                break;
        }
    }
    return nrm;
}
#undef HAS_METHOD

template<typename V, typename Iterable>
double World3d::Object3D::residual(const Iterable& v, unsigned normOut, unsigned normIn){
    auto nrm = get_norm_function<V, Iterable>(normOut, normIn);
    return nrm(v);
}

/*
 * "v:point" - exists by default
 * "v:connectivity", "h:connectivity", and "f:connectivity"
 * "v:removed", "e:removed", and "f:removed"
 * "f:normal"
 * "v:boundary_lbl"
 * "v:velocity"
 * "v:mass"
 * "v:force"
 */


#endif //AORTIC_VALVE_OBJECT3D_H
