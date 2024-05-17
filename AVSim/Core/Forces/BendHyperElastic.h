//
// Created by alex on 29.01.2021.
//

#ifndef AORTIC_VALVE_BENDHYPERELASTIC_H
#define AORTIC_VALVE_BENDHYPERELASTIC_H

#include "HyperElastic.h"
#include "Integrator1D.h"

//#define USE_BST
#define USE_DEBUG_WATCHES

namespace World3d {
    using SX_Case = std::pair<DReal, casadi::SX>;
    template<typename T, typename ...Ts>
    casadi::SX __SX_switch(casadi::SX& expr, casadi::SX& def, T t, Ts ... ts){
        if constexpr (sizeof...(ts) > 0) {
            return casadi::SX::if_else(expr == t.first, t.second, __SX_switch(expr, def, ts...));
        }
        else return casadi::SX::if_else(expr == t.first, t.second, def);
    }

    template<typename ...Ts>
    casadi::SX SX_switch(casadi::SX& expr, casadi::SX& def, Ts... ts){
        if constexpr (sizeof...(ts) > 0) {
            return __SX_switch(expr, def, ts...);
        } else return def;
    }

    template<typename T, typename ...Ts>
    casadi::SX __SX_switch(double tol, casadi::SX& expr, casadi::SX& def, T t, Ts ... ts){
        if constexpr (sizeof...(ts) > 0) {
            return casadi::SX::if_else(casadi::SX::abs(expr - t.first) < tol, t.second, __SX_switch(expr, def, ts...));
        }
        else return casadi::SX::if_else(casadi::SX::abs(expr - t.first) < tol, t.second, def);
    }

    template<typename ...Ts>
    casadi::SX SX_switch(double tol, casadi::SX& expr, casadi::SX& def, Ts... ts){
        if constexpr (sizeof...(ts) > 0) {
            return __SX_switch(tol, expr, def, ts...);
        } else return def;
    }

    class BendingForce: public ForceBase{
    public:
        using SX = casadi::SX;
        using SXVector = casadi::SXVector;
        using string = std::string;
        using Function = casadi::Function;
        template <typename T> using vector = std::vector<T>;

        using DataRequier = HH::DataRequier;
        using ConstDataRequier = HH::ConstDataRequier;
        using ConstDataRequierExt = HH::ConstDataRequierExt;
        using DefaultGenInputVar = HH::DefaultGenInputVar;
        using Requier = std::unique_ptr<DataRequier>;
        using LocalVar = DefaultHyperElasticForce::LocalVar;
        template<class T> using VarContainer = HH::VarContainer<T>;
        using GenFunction = HH::GenFunction;
        using Invariant = HH::Invariant;

        struct _BoundaryData{
            enum BTypes{
                FREE = 1,
                CLAMPED = 2
            };
            struct Entry{
                DReal btype = 0;
                DReal normal[3] = {0, 0, 0};
                Entry() = default;
                Entry(DReal btype): btype{btype} {}
                Entry(DReal btype, DReal nx, DReal ny, DReal nz): btype{btype} { normal[0] = nx, normal[1] = ny, normal[2] = nz;}
                Entry(DReal btype, DReal* norm): btype{btype} { std::copy(norm, norm + 3, normal); }
                Entry(DReal btype, std::array<DReal, 3> norm): btype{btype} { std::copy(norm.begin(), norm.end(), normal); }
            };
            struct InternalEntry{
                Entry entry;
                DReal loc_normal[2];
                DReal l0;
            };
            using FbndData = std::array<InternalEntry, 3>;

            FbndData data;
            bool prepared = false;
            std::function<InternalEntry(E_ind)> boundary;
            Object3D* _obj = nullptr;

            void set_boundary(Object3D& obj, const std::map<E_ind, Entry>& bdata);
            void generate_data(Object3D& obj);
            int check_boundary(Object3D& obj);
            void repair_boundary(Object3D& obj);
            void repair_if_necessary(Object3D& obj);
            void generate_data(Object3D& obj, const std::function<int(V_ind)>& conv);
            const DReal* get_data(F_ind f, std::vector<V_ind>& v);
        };
        using BoundaryData = _BoundaryData;

        struct BoundaryRequier: public DataRequier{
            BoundaryData* p;
            explicit BoundaryRequier(BoundaryData* bp): p{bp} {}
            void registerObj(Object3D* obj) override;
            void saveData(const DReal** target, F_ind f, std::vector<V_ind>& around) override
                {*target = p->get_data(f, around);}
            [[nodiscard]] std::unique_ptr<DataRequier> copy() const override { return std::make_unique<BoundaryRequier>(*this); }
        };
        struct k0Requier: public DataRequier{
            string gen_dir = "../generated";
            bool regenerate = true;
            BoundaryData* m_p = nullptr;
            ConstProperty<World3d::F_ind, std::array<DReal, 3>> k0;
            void registerObj(Object3D* obj) override;
            void saveData(const DReal** target, F_ind f, std::vector<V_ind>& around) override { *target = &k0[f][0]; }
            [[nodiscard]] std::unique_ptr<DataRequier> copy() const override { return std::make_unique<k0Requier>(*this); }
        };
        string force_name;
        string gen_dir;
        SX U;

        Object3D *_obj = nullptr, *_jobj = nullptr;
        GenFunction force_int, jacobian_int;
        GenFunction force_bnd, jacobian_bnd;
        int force_verb = 4, jacob_verb = 4;
#ifdef USE_DEBUG_WATCHES
        std::array<std::vector<GenFunction>, 2> watches; //{watches_int, watches_bnd}
        vector<DefaultGenInputVar> m_input_vars;
        VarContainer<LocalVar> m_loc_vars;
        VarContainer<Invariant> m_invariants;
#endif
        BoundaryData bdata;

        BendingForce(string f_name, string save_dir, vector<DefaultGenInputVar> input_vars, VarContainer<LocalVar>& loc_vars, VarContainer<Invariant>& invariants, SX U, bool regenerate = true);

        void setVerbosity(int force_verbosity, int jacobian_verbosity){ force_verb = force_verbosity, jacob_verb = jacobian_verbosity; }
        void setVerbosity(int verbosity) { setVerbosity(verbosity, verbosity); }
        void compute_forces(vector<DefaultGenInputVar>& generatedInput, VarContainer<LocalVar>& loc_vars, VarContainer<Invariant>& invariants);
        SX compute_force_expression(VarContainer<LocalVar>& loc_vars, VarContainer<Invariant>& invariants, bool is_bnd = true);
        void actualize_internal_data();
        BendingForce(const BendingForce& b){ *this = b; }
        std::unique_ptr<ForceBase> clone() override{ return std::make_unique<BendingForce>(*this); }
        BendingForce(BendingForce&& b) noexcept{ *this = std::move(b); }
        BendingForce& operator=(const BendingForce& a);
        BendingForce& operator=(BendingForce&& a) noexcept;

    private:
        static vector<DefaultGenInputVar> getDefaultAvailableInputVarsSimple();
        void Init(vector<DefaultGenInputVar> input_vars, VarContainer<LocalVar>& loc_vars, VarContainer<Invariant>& invariants, bool regenerate);
        std::vector<V_ind> _vdat;
    public:
        void prepareJacobianFunction(bool regenerate = true);
        explicit BendingForce(const DefaultHyperElasticForce& f, bool regenerate = true);
        void set_boundary_data(Object3D& obj, const std::map<E_ind, BoundaryData::Entry>& data)
            { bdata.set_boundary(obj, data); }
        void generate();
        void link_function();
        void registerObj(Object3D* obj);
        void registerjObj(Object3D* obj);
        int operator() (Object3D& obj) override;
        int element_matrix(Object3D& obj, ForceBase::LocMatrix& matr, F_ind f) override;

        BendingForce& add_watch(const std::function<SX(VarContainer<LocalVar>&, VarContainer<Invariant>&, casadi::SX U, bool is_bnd)>& watch_expr);
#ifdef USE_DEBUG_WATCHES
        Eigen::MatrixXd eval_watch(Object3D& obj, F_ind f, GenFunction* watch_int, GenFunction* watch_bnd, bool print = false);
        std::vector<Eigen::MatrixXd> eval_watches(Object3D& obj, F_ind f, bool print = false);
#endif
        static SX compute_current_curvature_expr(VarContainer<LocalVar> &loc_vars, VarContainer<Invariant> &invariants, SX u, bool is_bnd);

#ifndef USE_BST
        static VarContainer<LocalVar> getDefaultLocalVars();
#else
        static VarContainer<LocalVar> getDefaultLocalVarsBST();
#endif
        static VarContainer<Invariant> getDefaultInvariants(){ return DefaultHyperElasticForce::getDefaultInvariants();}
        static vector<DefaultGenInputVar> getDefaultAvailableInputVars();
    };
}


#endif //AORTIC_VALVE_BENDHYPERELASTIC_H
