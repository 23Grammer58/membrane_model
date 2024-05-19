//
// Created by alex on 29.01.2021.
//

#ifndef AORTIC_VALVE_HYPERELASTIC_H
#define AORTIC_VALVE_HYPERELASTIC_H


#include "../Object3D.h"
#include "casadi/casadi.hpp"
#include <string>
#include "ForcesCommon.h"

#if __cplusplus >= 201703L
#include <filesystem>
#else
#include <experimental/filesystem>
#endif

#define USE_DEBUG_WATCHES

namespace World3d{
    namespace HyperElasticHelpers{
        using SX = casadi::SX;
        using SXVector = casadi::SXVector;
        using string = std::string;
        using Function = casadi::Function;
        template <typename T>
        using vector = std::vector<T>;

        struct DataRequier{
            virtual void registerObj(Object3D* obj) {}
            virtual void saveData(const DReal**  target, F_ind f, std::vector<V_ind>& around) = 0;// {};
            [[nodiscard]] virtual std::unique_ptr<DataRequier> copy() const = 0;// { return std::make_unique<DataRequier>(*this); };
        };

        struct ConstDataRequier: public DataRequier{
            const DReal data;
            void registerObj(Object3D* obj) override {}
            void saveData(const DReal**  target, F_ind f, std::vector<V_ind>& around) override { *target = &data; }
            [[nodiscard]] std::unique_ptr<DataRequier> copy() const override { return std::make_unique<ConstDataRequier>(*this); }
            explicit ConstDataRequier(DReal data): data{data} {};
            ConstDataRequier& operator=(const ConstDataRequier&) = default;
        };
        template<typename T>
        struct ConstDataTagRequier: public DataRequier{
            ConstProperty<F_ind, T> A;
            DReal data;
            string tagName;
            void registerObj(Object3D* obj) override {
                A = obj->m_mesh.property_map<F_ind, T>(tagName);
                if (!A.second) throw std::runtime_error("Tag \"" +  tagName + "\" is not exists on the object \"" + obj->name + "\"");
            }
            void saveData(const DReal**  target, F_ind f, std::vector<V_ind>& around) override { data = A[f]; *target = &data; }
            [[nodiscard]] std::unique_ptr<DataRequier> copy() const override { return std::make_unique<ConstDataTagRequier>(*this); }
            explicit ConstDataTagRequier(std::string tag): tagName{tag} {};
            ConstDataTagRequier& operator=(const ConstDataTagRequier&) = default;
        };
        template<typename T = DReal>
        struct ConstDataTagOrValRequier: public DataRequier{
            ConstProperty<F_ind, T> A;
            DReal data;
            T val;
            std::string tagName;
            bool data_from_tag = false;

            void registerObj(Object3D* obj) override {
                A = obj->m_mesh.property_map<F_ind, T>(tagName);
                data_from_tag = A.second;
            }
            void saveData(const DReal**  target, F_ind f, std::vector<V_ind>& around) override { if (data_from_tag) data = A[f]; else data = val; *target = &data; }
            [[nodiscard]] std::unique_ptr<DataRequier> copy() const override { return std::make_unique<ConstDataTagOrValRequier>(*this); }
            ConstDataTagOrValRequier(std::string tag, DReal val): tagName{tag}, val{val} {};
            ConstDataTagOrValRequier& operator=(const ConstDataTagOrValRequier&) = default;
        };
        struct ConstDataRequierExt: DataRequier{
            const DReal* data;
            void registerObj(Object3D* obj) override {}
            void saveData(const DReal**  target, F_ind f, std::vector<V_ind>& around) override { *target = data; }
            [[nodiscard]] std::unique_ptr<DataRequier> copy() const override { return std::make_unique<ConstDataRequierExt>(*this); }
            explicit ConstDataRequierExt(const DReal* data): data{data} {};
        };

        struct PntRequier: public DataRequier{
            const int _N;
            Object3D* _obj = nullptr;
            Vector vdat;
            void registerObj(Object3D* obj) override { _obj = obj; }
            void saveData(const DReal** target, F_ind f, std::vector<V_ind>& around) override {
                vdat = _obj->m_x[around[_N]] - _obj->m_x[around[0]];
                *target = &vdat[0];
            }
            [[nodiscard]] std::unique_ptr<DataRequier> copy() const override { return std::make_unique<PntRequier>(*this); }
            explicit PntRequier(int n): _N{n} {}
        };
        struct PntInitRequier: public DataRequier{
            const int _N;
            Object3D* _obj = nullptr;
            Vector vdat;
            void registerObj(Object3D* obj) override { _obj = obj; }
            void saveData(const DReal** target, F_ind f, std::vector<V_ind>& around) override {
                vdat = _obj->m_x0[around[_N]]  - _obj->m_x0[around[0]];
                *target = &vdat[0];
            }
            [[nodiscard]] std::unique_ptr<DataRequier> copy() const override { return std::make_unique<PntInitRequier>(*this); }
            explicit PntInitRequier(int n): _N{n} {}
        };

        template<int N>
        struct VectorRequier: public DataRequier{
            ConstProperty<F_ind, std::array<DReal, N>> A;
            string tagName;
            void registerObj(Object3D* obj) override {
                A = obj->m_mesh.property_map<F_ind, std::array<DReal, N>>(tagName);
                if (!A.second) throw std::runtime_error("Tag \"" +  tagName + "\" is not exists on the object \"" + obj->name + "\"");
            }
            void saveData(const DReal** target, F_ind f, std::vector<V_ind>& around) override { *target = &A[f][0]; }
            VectorRequier(string tagName): tagName{tagName} {}
            [[nodiscard]] std::unique_ptr<DataRequier> copy() const override { return std::make_unique<VectorRequier>(*this); }
        };

        template<int N>
        struct ConstVectorRequier: public DataRequier{
            std::array<DReal, N> dat;
            void registerObj(Object3D* obj) override { }
            void saveData(const DReal** target, F_ind f, std::vector<V_ind>& around) override { *target = dat.data(); }
            ConstVectorRequier(std::array<DReal, N> arr): dat{arr} {}
            [[nodiscard]] std::unique_ptr<DataRequier> copy() const override { return std::make_unique<ConstVectorRequier>(*this); }
        };

        struct DefaultGenInputVar{
            using Requier = std::unique_ptr<DataRequier>;
            string name;
            SX var;
            Requier require;
            void setDataRequier(Requier&& req) { require = move(req); }
            explicit DefaultGenInputVar(const string& name, int nrow = 1, int ncol = 1): name{name}, var{SX::sym(name, nrow, ncol)} {}
            DefaultGenInputVar(const string& name, int nrow, int ncol, Requier&& require): name{name}, var{SX::sym(name, nrow, ncol)}, require{std::move(require)} {}
            DefaultGenInputVar(string name, const SX& sym): name{std::move(name)}, var{sym} {}
            DefaultGenInputVar(string name, const SX& sym, Requier&& require): name{std::move(name)}, var{sym}, require{std::move(require)} {}
            DefaultGenInputVar(const DefaultGenInputVar& v): name{v.name}, var{v.var}, require{v.require->copy()} {}
            DefaultGenInputVar() = default;
            DefaultGenInputVar& operator=(const DefaultGenInputVar& v);
            SX get_sym() const { return var; }
            string getName() const { return name; }
            DefaultGenInputVar& operator=(DefaultGenInputVar&& v) noexcept = default;
            DefaultGenInputVar(DefaultGenInputVar&& v) noexcept = default;
        };
        struct GenFunction{
#define CHECK_GENERATE_NANS
            vector<DefaultGenInputVar> generatedInput;
            vector<const DReal*> funcRealInput;
            SX func_expr;
            std::shared_ptr<Function> func;
            vector<int> dense_remap;

            GenFunction() = default;
            GenFunction(GenFunction&&) noexcept = default;
            GenFunction(const GenFunction&) = default;
            GenFunction& operator=(const GenFunction& ) = default;
            GenFunction& operator=(GenFunction&& ) noexcept = default;
            GenFunction(const SX& expr, vector<DefaultGenInputVar>& input);
            GenFunction(vector<DefaultGenInputVar>& input, const string& fname, const string& gen_dir = "");

            GenFunction(const SX& expr, vector<DefaultGenInputVar>& input, const string& fname, const string& gen_dir = "", bool regenerate = true);

            void set_generatedInput(vector<DefaultGenInputVar>& input){ generatedInput = input; }
            void set_func_expression(const SX& expr){ func_expr = expr; }
            void update_input();
            void generate(const string& fname, string gen_dir = "");
            void link_function(const string& fname, const string& gen_dir = "");
            void registerObj(Object3D* obj);
            int operator()(DReal** out, F_ind f, std::vector<V_ind>& v, int debugVar = 0);
            void print_input();

#ifdef CHECK_GENERATE_NANS
        private:
            static std::ostream& printMatrix(std::ostream& out, DReal* arr, int lda, int ldb);
            static bool arrayHasNaN(DReal* arr, int lda, int ldb);
            bool funcRealInputHasNaN();
            std::ostream& print_funcRealInput(std::ostream& out, bool oneline = true);
#endif
        };
        GenFunction makeJacobianFunction(string f_name, string save_dir, GenFunction func, SX var, bool regenerate = true, bool print_gen = false);

        struct LocalVar{
            string name;
            SX sym;
            SX val;
            LocalVar() = default;
            LocalVar(const LocalVar& v) = default;
            LocalVar(string name, int sz1 = 1, int sz2 = 1): sym(SX::sym(name, sz1, sz2)), val(SX::sym(name, sz1, sz2)), name{name} {}
            LocalVar(string name, SX val): sym(SX::sym(name, val.size1(), val.size2())), val{val}, name{name} { }
            LocalVar(string name, SX sym, SX val): name{name}, sym{sym}, val{val} {}
            string getName() const { return name; }
            SX get_sym() const { return sym; }
            void substituteIn(SX& ex) const { ex = SX::substitute(ex, sym, val); }
        };
        struct Invariant{
            string name;
            SX sym;
            SX val;
            SX Qderiv;
            SX C2dval;
            SX C2dderiv;
            Invariant(string name, SX val, SX C2dval, SX Qderiv, SX C2dderiv);
            Invariant(string name, SX sym, SX val, SX C2dval, SX Qderiv, SX C2dderiv);
            string getName() const { return name; }
            SX get_sym() const { return sym; }
        };
        template<class Nameable>
        struct VarContainer{
            std::vector<Nameable> seq;
            std::map<string, int> streq;
            SX& operator[](string name) { return seq[streq[name]].sym; }
            SX& operator[](int id) { return seq[id].sym; }
            const SX& operator[](string name) const { auto it = streq.find(name); return seq[it->second].sym; }
            const SX& operator[](int id) const { return seq[id].sym; }

            auto find(string name);
            const auto find(string name) const;
            Nameable& at(string name);
            Nameable& at(int id) { return seq.at(id); }
            const Nameable& at(string name) const;
            const Nameable& at(int id) const { return seq.at(id); }
            void push_back(Nameable var);
            void push_back(std::initializer_list<Nameable> vars);
            void remove(string name);
            void remove(int id);
            int size() const { return seq.size(); }
            auto begin() { return seq.begin(); }
            auto end() { return seq.end(); }
            auto begin() const { return seq.begin(); }
            auto end() const { return seq.end(); }
            auto cbegin() const { return seq.cbegin(); }
            auto cend() const { return seq.cend(); }
            //updates string iterator
            void update();
            void test();
        };


        using ArgsDef = std::tuple<vector<DefaultGenInputVar>, VarContainer<LocalVar>, VarContainer<Invariant>>;
        void printVars(const vector<DefaultGenInputVar>& vars, string prefix = "");
        void printVars(const VarContainer<LocalVar>& vars, string prefix = "");
        void printVars(const VarContainer<Invariant>& vars, string prefix = "");
        void printVars(const vector<DefaultGenInputVar>& vg, const VarContainer<LocalVar>& vl, const VarContainer<Invariant>& vi, string prefix = "");
        void printVars(const ArgsDef& v, string prefix = "");

        ArgsDef concat(ArgsDef& a, ArgsDef b);

        void correlateInputAndVariables(const vector<DefaultGenInputVar>& generatedInput, VarContainer<LocalVar>& localVars, VarContainer<Invariant>& invariants);
        void correlateInputAndVariables(const vector<DefaultGenInputVar>& generatedInput, VarContainer<LocalVar>& localVars);
        SX applyVariablesToExpr(SX expr, VarContainer<LocalVar>& loc_vars);

        VarContainer<LocalVar> getDefaultLocalVars();
        VarContainer<Invariant> getDefaultInvariants();
        vector<DefaultGenInputVar> getDefaultAvailableInputVars();

        vector<DefaultGenInputVar> makeOrthotropicInvariantInput(string vecName, string vec2dName, string vecdtAName, string vectorTag);
        vector<DefaultGenInputVar> makeOrthotropicInvariantInput(string vecName, string vec2dName, string vecdtAName, std::array<DReal, 3> Av);
        std::pair<VarContainer<LocalVar>, VarContainer<Invariant>> makeOrthotropicInvariantVars(string InvName, string vecName, string vec2dName, string vecdtAName);
        std::pair<VarContainer<LocalVar>, VarContainer<Invariant>> makeAnisotropicInvariantVars(string InvName, std::array<string, 2> vecName, std::array<string, 2> vec2dName, std::array<string, 2> vecdtAName);

        ArgsDef makeOrthotropicInvariant(string InvName, string vectorTag);
        ArgsDef makeOrthotropicInvariant(string InvName, std::array<DReal, 3> A);
        ArgsDef makeAnisotropicInvariant(string InvName, string vectorTag1, string vectorTag2);
    };

    namespace HH = HyperElasticHelpers;

    class DefaultHyperElasticForce{
    public:
        using SX = casadi::SX;
        using SXVector = casadi::SXVector;
        using string = std::string;
        using Function = casadi::Function;
        template <typename T> using vector = std::vector<T>;
        using Requier = std::unique_ptr<HH::DataRequier>;
        using ConstDataRequier = HH::ConstDataRequier;
        using ConstDataRequierExt = HH::ConstDataRequierExt;
        using PntRequier = HH::PntRequier;
        using PntInitRequier = HH::PntInitRequier;
        template<int N> using VectorRequier = HH::VectorRequier<N>;
        template<int N> using ConstVectorRequier = HH::ConstVectorRequier<N>;
        using DefaultGenInputVar = HH::DefaultGenInputVar;
        using GenFunction = HH::GenFunction;
        using LocalVar = HH::LocalVar;
        using Invariant = HH::Invariant;
        template<typename T> using VarContainer = HH::VarContainer<T>;
        using ArgsDef = HH::ArgsDef;

        vector<DefaultGenInputVar> inputVars;
        VarContainer<LocalVar> localVars;
        VarContainer<Invariant> invariants;
        SX U;
        GenFunction genForce, genJacobian;
#ifdef USE_DEBUG_WATCHES
        std::vector<GenFunction> m_watches;
#endif

        string gen_dir = "", force_name = "";
        Object3D *_obj = nullptr, *_jobj = nullptr;
        int force_verb = 4, jacob_verb = 4;

        std::vector<V_ind> _vdat;

        DefaultHyperElasticForce() = default;
        DefaultHyperElasticForce(DefaultHyperElasticForce&& f) = default;
        DefaultHyperElasticForce& operator=(DefaultHyperElasticForce&& f) = default;
        DefaultHyperElasticForce(const DefaultHyperElasticForce& f);
        DefaultHyperElasticForce& operator=(const DefaultHyperElasticForce& f);

        DefaultHyperElasticForce(string f_name, string save_dir, ArgsDef args, SX U, bool regenerateRHS, bool regenerateJ):
                DefaultHyperElasticForce(f_name, save_dir, args, U, regenerateRHS){ prepareJacobianFunction(regenerateJ); }
        DefaultHyperElasticForce(string f_name, string save_dir, ArgsDef args, SX U, bool regenerate = true):
                DefaultHyperElasticForce(f_name, save_dir, std::move(std::get<0>(args)), std::move(std::get<1>(args)), std::move(std::get<2>(args)), U, regenerate) {}
        DefaultHyperElasticForce(string f_name, string save_dir, vector<DefaultGenInputVar> input_vars, VarContainer<LocalVar> loc_vars, VarContainer<Invariant> invariants, SX U, bool regenerateRHS, bool regenerateJ):
                DefaultHyperElasticForce(f_name, save_dir, input_vars, loc_vars, invariants, U, regenerateRHS) { prepareJacobianFunction(regenerateJ); }
        DefaultHyperElasticForce(string f_name, string save_dir, vector<DefaultGenInputVar> input_vars, VarContainer<LocalVar> loc_vars, VarContainer<Invariant> invariants, SX U, bool regenerate = true);
        void prepareJacobianFunction(bool regenerate = true);

        void setVerbosity(int force_verbosity, int jacobian_verbosity){ force_verb = force_verbosity, jacob_verb = jacobian_verbosity; }
        void setVerbosity(int verbosity) { setVerbosity(verbosity, verbosity); }
        void correlateInputAndVariables() { HH::correlateInputAndVariables(genForce.generatedInput, localVars, invariants); }
        SX compute_force();
        SX compute_jacobian();
        void generate();
        void generate_jacobian();
        void link_function() { genForce.link_function(force_name, gen_dir); }
        void link_jacobian(){ genJacobian.link_function(force_name + "_Jacob", gen_dir); }
        void registerObj(Object3D* obj) { genForce.registerObj(_obj = obj); }
        void registerjObj(Object3D* obj){ genJacobian.registerObj(_jobj = obj); }
        int operator() (Object3D& obj);
        int element_matrix(Object3D& obj, ForceBase::LocMatrix& matr, F_ind element);

        static VarContainer<LocalVar> getDefaultLocalVars() { return HH::getDefaultLocalVars(); }
        static VarContainer<Invariant> getDefaultInvariants() { return HH::getDefaultInvariants(); }
        static vector<DefaultGenInputVar> getDefaultAvailableInputVars() { return HH::getDefaultAvailableInputVars(); }
        static ArgsDef makeOrthotropicInvariant(string InvName, string vectorTag){
            return HH::makeOrthotropicInvariant(std::move(InvName), std::move(vectorTag));
        }
        static ArgsDef makeOrthotropicInvariant(string InvName, std::array<DReal, 3> A) {
            return HH::makeOrthotropicInvariant(std::move(InvName), A);
        }

        int add_watch(const std::function<SX(VarContainer<LocalVar>&, VarContainer<Invariant>&, casadi::SX U)>& watch_expr, bool quiet = false, bool use_random_prefix = false, const std::string& prefix = "he_watch");
#ifdef USE_DEBUG_WATCHES
        Eigen::MatrixXd eval_watch(Object3D& obj, F_ind f, int watch_num, bool print = false);
        Eigen::MatrixXd eval_watch(Object3D& obj, F_ind f, GenFunction* watch, bool print = false);
        std::vector<Eigen::MatrixXd> eval_watches(Object3D& obj, F_ind f, bool print = false);
#endif
        static SX compute_HS_tensor_expr(VarContainer<LocalVar>& lVars, VarContainer<Invariant>& lInvs, casadi::SX U);
        static SX compute_S_tensor_expr(VarContainer<LocalVar>& lVars, VarContainer<Invariant>& lInvs, casadi::SX U);
    };

    struct HyperElasticForceBase: public ForceBase{
        using SX = casadi::SX;
        using Requier = DefaultHyperElasticForce::Requier;
        using ConstDataRequier = HH::ConstDataRequier;
        template<typename T> using ConstDataTagRequier = HH::ConstDataTagRequier<T>;
        using IV = HH::DefaultGenInputVar;

        DefaultHyperElasticForce f;
        bool _with_jacob = false;

        void prepareJacobianFunction(bool regenerate = true);
        int operator() (Object3D& obj) override { return f(obj); }
        int element_matrix(Object3D& obj, ForceBase::LocMatrix& matr, F_ind element) override{
            return f.element_matrix(obj, matr, element);
        }
        virtual std::unique_ptr<HyperElasticForceBase> copy() = 0;
        std::unique_ptr<ForceBase> clone() final;
    };

    namespace HyperElasticHelpers{
        ArgsDef makeDefaultInits(const DefaultHyperElasticForce& f, const VarContainer<LocalVar>& local_vars);
    }
}

#include "HyperElastic.inl"


#endif //AORTIC_VALVE_HYPERELASTIC_H
