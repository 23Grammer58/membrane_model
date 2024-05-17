#ifndef MINENERGYDEFORMATOR_H
#define MINENERGYDEFORMATOR_H

#include <cstring>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <sstream>
#include <vector>
#include <set>
#include <algorithm>
#include <gsl/gsl_multimin.h>
#include <chrono>
#include <memory>
#include <Eigen/Dense>
#include "AVSim/Core/Object3D.h"
#include "AVSim/Core/World.h"

struct EnergyDeformatorMesh{
    typedef double DReal;
    typedef unsigned VIND;
    typedef unsigned EIND;
    typedef unsigned FIND;
    typedef Eigen::Matrix<DReal, 3, 1> Vector;
    typedef std::array<VIND, 2> Edge;
    typedef std::array<VIND, 3> Face;
    typedef std::pair<std::array<VIND, 4>, bool> Digedral;

    virtual Vector x(VIND v) const = 0;
    virtual void set_x(VIND v, Vector x) = 0;
    virtual Vector x0(VIND v) const = 0;
    virtual bool is_fix(VIND v) const = 0;
    virtual Edge edge(EIND e) const  = 0;
    virtual Face face(FIND f) const = 0;
    virtual Digedral digedral(EIND e) const = 0;
    virtual unsigned num_vertices() const = 0;
    virtual unsigned num_edges() const = 0;
    virtual unsigned num_faces() const = 0;
    virtual void updateToNextIteration(){}

    virtual DReal l(EIND e) const {
        Edge v = edge(e);
        return (x(v[1]) - x(v[0])).stableNorm();
    }
    virtual DReal l0(EIND e) const {
        Edge v = edge(e);
        return (x0(v[1]) - x0(v[0])).stableNorm();
    }
    virtual Vector Sor(FIND f) const {
        Face v = face(f);
        return (x(v[1]) - x(v[0])).cross(x(v[2]) - x(v[0]))/2;
    }
    virtual Vector n(FIND f) const {
        Face v = face(f);
        return (x(v[1]) - x(v[0])).cross(x(v[2]) - x(v[0])).normalized().eval();
    }
    virtual DReal S(FIND f) const {
        return Sor(f).stableNorm();
    }
    virtual DReal S0(FIND f) const {
        Face v = face(f);
        return (x0(v[1]) - x0(v[0])).cross(x0(v[2]) - x0(v[0])).stableNorm()/2;
    }
};

struct EnergyDeformatorOb3D: public EnergyDeformatorMesh{
    using Renderer = std::shared_ptr<World3d::RendererBase>;
private:
    using Object3D = World3d::Object3D;
    using VIND = unsigned;
    using EIND = unsigned;
    using FIND = unsigned;
    using V_ind = World3d::V_ind;
    using E_ind = World3d::E_ind;
    using F_ind = World3d::F_ind;
    using Point = World3d::Point;

    using Storage = World3d::Mesh::Property_map<V_ind, Point>;
    Object3D& obj;
    int loc_id;
    Storage sx, sx0;

    Renderer _renderer;
    World3d::Timer rtt;

    void render(){
        obj.apply_render();
        _renderer->render();
    }

public:
    EnergyDeformatorOb3D(Object3D& obj, std::string x_tag = "v:x", std::string x0_tag = "v:point"): obj{obj}{
        static int id = 0;
        loc_id = id++;
        auto sx_s = obj.m_mesh.property_map<V_ind, Point>(x_tag);
        auto sx0_s = obj.m_mesh.property_map<V_ind, Point>(x0_tag);
        if (!sx_s.second || !sx0_s.second)
            throw std::runtime_error("In object is not setted \"" + x_tag + "\" or \"" + x0_tag + "\" tag");
        sx = sx_s.first;
        sx0 = sx0_s.first;

        _renderer = std::make_shared<World3d::RendererBase>();
    }

    void setRenderer(Renderer&& renderer){
        _renderer = move(renderer);
        _renderer->set_renderer(loc_id, obj);
    }
    ~EnergyDeformatorOb3D(){
        _renderer->delete_renderer(loc_id, obj);
    }

    Vector x(VIND vu) const override{ V_ind v(vu); return Vector{sx [v][0], sx [v][1], sx [v][2]}; }
    Vector x0(VIND vu) const override{ V_ind v(vu); return Vector{sx0[v][0], sx0[v][1], sx0[v][2]}; }
    void set_x(VIND vu, Vector x) override { sx[V_ind(vu)] = Point(x[0], x[1], x[2]); }
    bool is_fix(VIND vu) const override { V_ind v(vu); return !obj.is_movable(v); }
    Edge edge(EIND eu) const override {
        E_ind e(eu);
        auto v = World3d::vert_around(obj.m_mesh, e);
        return Edge{v[0], v[1]};
    }
    Face face(FIND fu) const override {
        F_ind f(fu);
        auto v = World3d::vert_around(obj.m_mesh, f);
        return Face{v[0], v[1], v[2]};
    }
    Digedral digedral(EIND eu) const override {
        E_ind e(eu);
        Digedral res;
        auto f = World3d::face_around(obj.m_mesh, e);
        if (!f[1].second){
            res.second = false;
            auto v = World3d::vert_around(obj.m_mesh, f[0].first);
            res.first = {v[0], v[1], v[2], v[0]};
        } else {
            auto ve =  World3d::vert_around(obj.m_mesh, e);
            auto v = std::array<std::array<V_ind, 3>, 2>{World3d::vert_around(obj.m_mesh, f[0].first), World3d::vert_around(obj.m_mesh, f[1].first)};
            int id[2] = {0, 0};
            bool reversed = false;
            for (int fi = 0; fi < 2; ++fi) {
                auto& fv = v[fi];
                for (int i = 0; i < 3; ++i)
                    if ((fv[i] == ve[0] && fv[(i + 1)%3] == ve[1])) {
                        id[fi] = i;
                        break;
                    } else if (fv[i] == ve[1] && fv[(i + 1)%3] == ve[0]){
                        id[fi] = i;
                        if (fi == 0) reversed = true;
                        break;
                    }
            }
            if (reversed) {
                std::swap(id[0], id[1]);
                std::swap(v[0], v[1]);
            }
            res.second = true;
            res.first = {v[0][id[0]], v[0][(id[0] + 1)%3], v[0][(id[0] + 2)%3], v[1][(id[1] + 2)%3]};
        }
        return res;
    }
    unsigned num_vertices() const override { return obj.m_mesh.num_vertices(); }
    unsigned num_edges() const override { return obj.m_mesh.num_edges(); }
    unsigned num_faces() const override { return obj.m_mesh.num_faces(); }
    void updateToNextIteration() override{
        if (rtt.elapsed() > 1.0/60){
            render();
            rtt.reset();
        }
    }
};


class MinEnergyDeformator{
    using EDM = EnergyDeformatorOb3D;
    using EFDM = EnergyDeformatorMesh;
private:
    EDM m_mesh;

public:
    class EnergyFunctor{
        static inline double zero(const EFDM& mesh, int node_id, double* df) { return 0; }
    private:
        using Energy = double(*)(const EFDM& mesh, int node_id, double* df);
        Energy func;
    public:
        using Vector = EFDM::Vector;
        EnergyFunctor( double(*func)(const EFDM& mesh, int node_id, double* df) = zero): func{func} {};
        virtual double operator()(const EFDM& mesh, int node_id, double* df)
        { return func(mesh, node_id, df); }
    };

    // this functions should save gradient into df, if df != nullptr
    using NodeEnergy =  std::shared_ptr<EnergyFunctor>; //df[3]
    using EdgeEnergy = std::shared_ptr<EnergyFunctor>; //df[6]
    using DigedralEnergy = std::shared_ptr<EnergyFunctor>; //df[12]

private:
    std::vector<NodeEnergy> m_node_energy;
    std::vector<EdgeEnergy> m_edge_energy;
    std::vector<DigedralEnergy> m_digedral_energy;

    std::vector<int> m_node_map;
    int m_N_vars = 0;

public:
    MinEnergyDeformator(EDM mesh): m_mesh{mesh}
    {
        set_node_map();
    }

    void addNodeEnergyComponent(NodeEnergy f) { m_node_energy.push_back(f); }
    void addEdgeEnergyComponent(EdgeEnergy f) { m_edge_energy.push_back(f); }
    void addDigedralEnergyComponent(DigedralEnergy f) { m_digedral_energy.push_back(f); }
    double getEnergy() { return compute(); }
private:

    inline void set_coords_to_mesh(const gsl_vector* X){
        for (int j = 0, i = 0; j < m_mesh.num_vertices(); ++j){
            if (m_mesh.is_fix(j)){
                continue;
            }
            EFDM::Vector v;
            for (int k = 0; k < 3; ++k)
                v[k] = gsl_vector_get(X, i++);
            m_mesh.set_x(j, v);
        }
    }

    inline void set_node_map(){
        int i = 0;
        m_node_map.resize(m_mesh.num_vertices());
        for (int j = 0; j < m_mesh.num_vertices(); ++j){
            if (m_mesh.is_fix(j)){
                m_node_map[j] = -1;
                continue;
            }
            m_node_map[j] = i;
            i++;
        }
        m_N_vars = 3 * i;
    }

    double compute(){
        double f = 0;
        for (int i = 0; i < m_mesh.num_edges(); ++i){
            auto [ve, is_degedral] = m_mesh.digedral(i);
            if (m_mesh.is_fix(ve[0]) && m_mesh.is_fix(ve[1])) continue;
            for (auto& en: m_edge_energy)
                f += en->operator()(m_mesh, i, nullptr);
            if (is_degedral)
                for (auto& en: m_digedral_energy)
                    f += en->operator()(m_mesh, i, nullptr);
        }

        for (int i = 0; i < m_mesh.num_vertices(); ++i)
            for (auto& en: m_node_energy)
                f += en->operator()(m_mesh, i, nullptr);
        return f;
    }

    inline void set_templated(int* id, int n, gsl_vector *df, double* ddf){
        for (int i = 0; i < n; ++i)
            if (id[i] >= 0)
                for (int j = 0; j < 3; ++j) {
                    assert(!std::isnan(ddf[i * 3 + j]));
                    gsl_vector_set(df, id[i] * 3 + j, gsl_vector_get(df, id[i] * 3 + j) + ddf[i * 3 + j]);
                }
    }

    void compute_df(gsl_vector *df){
        gsl_vector_set_all(df, 0);

        for (int i = 0; i < m_mesh.num_edges(); ++i) {
//            if (i == 8)
//                std::cout << "Now we get nan!" << std::endl;
            auto [ve, is_degedral] = m_mesh.digedral(i);
            {
                double loc_df[6] = {};
                for (auto &en: m_edge_energy)
                    en->operator()(m_mesh, i, loc_df);
                int id[2] = {m_node_map[ve[0]], m_node_map[ve[1]]};
                set_templated(id, 2, df, loc_df);
            }

            if (is_degedral) {
                double loc_df[12] = {};
                for (auto& en: m_digedral_energy)
                    en->operator()(m_mesh, i, loc_df);
                int id[4] = {m_node_map[ve[2]], m_node_map[ve[3]],
                             m_node_map[ve[0]], m_node_map[ve[1]]};
                set_templated(id, 4, df, loc_df);
            }
        }

        for (int i = 0; i < m_mesh.num_vertices(); ++i){
            double loc_df[3] = { 0 };
            for (auto& en: m_node_energy)
                en->operator()(m_mesh, i, loc_df);
            set_templated(&m_node_map[i], 1, df, loc_df);
        }
    }

    static void _fdf4gsl(const gsl_vector *x, void *params, double *f, gsl_vector *df)
    {
        MinEnergyDeformator* m = (MinEnergyDeformator*)params;
        m->set_coords_to_mesh(x);
        *f = m->compute();
        m->compute_df(df);
    }

    static double _f4gsl (const gsl_vector *x, void *params)
    {
        MinEnergyDeformator* m = (MinEnergyDeformator*)params;
        m->set_coords_to_mesh(x);
        return m->compute();
    }

    static void _df4gsl (const gsl_vector *x, void *params, gsl_vector *df)
    {
        MinEnergyDeformator* m = (MinEnergyDeformator*)params;
        m->set_coords_to_mesh(x);
        m->compute_df(df);
    }

public:
    int find_minimum_energy_df(int freq = 1, double step_sz = 1.e-4, double tol = 1.e-4, double epsabs = 1.e-2, int maxits = 50000, double time/*secs*/ = 150){
        size_t iter = 0;
        int status;

        const gsl_multimin_fdfminimizer_type *T;
        gsl_multimin_fdfminimizer *s;

        void* par = this;

        int g_N = m_N_vars;
        gsl_vector *x = gsl_vector_alloc (g_N);
        for (int j = 0, i = 0; j < m_mesh.num_vertices(); ++j){
            if (m_mesh.is_fix(j)){
                continue;
            }
            EFDM::Vector dat = m_mesh.x(j);
            for (int k = 0; k < 3; ++k)
                gsl_vector_set(x, i++, dat[k]);// + (k == 2) ? 1.0 : 0.0);
        }

        gsl_multimin_function_fdf my_func;

        my_func.n = g_N;
        my_func.f = _f4gsl;
        my_func.df = _df4gsl;
        my_func.fdf = _fdf4gsl;
        my_func.params = par;

        T = gsl_multimin_fdfminimizer_conjugate_fr;
        s = gsl_multimin_fdfminimizer_alloc (T, g_N);

        gsl_multimin_fdfminimizer_set (s, &my_func, x, step_sz, tol);

        class _Timer
        {
        private:
            using clock_t = std::chrono::high_resolution_clock;
            using second_t = std::chrono::duration<double, std::ratio<1> >;
            std::chrono::time_point<clock_t> m_beg;
        public:
            _Timer() : m_beg(clock_t::now()){}
            void reset() { m_beg = clock_t::now(); }

            double elapsed() const
            {
                return std::chrono::duration_cast<second_t>(clock_t::now() - m_beg).count();
            }
        };

        _Timer t;
        do
        {
            m_mesh.updateToNextIteration();

            iter++;
            status = gsl_multimin_fdfminimizer_iterate (s);
            if (status)
                break;

            status = gsl_multimin_test_gradient (s->gradient, epsabs);
            if (status == GSL_SUCCESS)
                std::cout << "Minimum found at: " << iter << ": f() = " << s->f << std::endl;

            if ((freq > 0) && !(iter % freq))
                std::cout << iter << ": f() = " << s->f << std::endl;
        }
        while (status == GSL_CONTINUE && iter < maxits && t.elapsed() < time);
        if (status != GSL_SUCCESS)
            std::cout << "Success didn't reached status = " << status << ": " <<  iter << ": f() = " << s->f << std::endl;
        m_mesh.updateToNextIteration();


        gsl_multimin_fdfminimizer_free (s);
        gsl_vector_free (x);

        return status;
    }

};

void set_default_length_constr(MinEnergyDeformator& m, double weight);
void set_default_digedral_angle_constr(MinEnergyDeformator& m, double weight, double scale = 1.0, bool convexity = true);
void set_plane_constr(MinEnergyDeformator& m, double weight, const Eigen::Vector3d& normal, double b);
void set_isotrop_force(MinEnergyDeformator& m, double weight, const Eigen::Vector3d& force);

#endif