//
// Created by alex on 29.01.2021.
//

#include "HyperElastic.h"

using namespace World3d;
using namespace HH;
using namespace std;

DefaultGenInputVar & DefaultGenInputVar::operator=(const DefaultGenInputVar &v) {
    if (this == &v) return *this;
    name = v.name;
    var = v.var;
    require = v.require->copy();
    return *this;
}

GenFunction::GenFunction(const SX &expr, vector<DefaultGenInputVar> &input): func_expr{expr}, generatedInput{input}{
    update_input();
}

GenFunction::GenFunction(vector<DefaultGenInputVar> &input, const string &fname, const string &gen_dir) : generatedInput{input}{
    link_function(fname, gen_dir);
}

GenFunction::GenFunction(const SX &expr, vector<DefaultGenInputVar> &input,
                         const string &fname, const string &gen_dir,bool regenerate)
                         : func_expr{expr}, generatedInput{input} {
    update_input();
    if (regenerate)
        generate(fname, gen_dir);
    link_function(fname, gen_dir);
}

void GenFunction::update_input() {
    vector<DefaultGenInputVar> updated;
    for (auto& i: generatedInput)
        if (SX::depends_on(func_expr, i.var))
            updated.emplace_back(std::move(i));
    generatedInput = std::move(updated);
}

void GenFunction::registerObj(Object3D *obj) {
    funcRealInput.resize(generatedInput.size());
    for (auto& i: generatedInput) {
        i.require.get()->registerObj(obj);
    }
}


int GenFunction::operator()(World3d::DReal **out, World3d::F_ind f, std::vector<V_ind> &v,
                                                      int verb) {
    for (int i = 0; i < generatedInput.size(); ++i)
        generatedInput[i].require.get()->saveData(&funcRealInput[i], f, v);
#ifdef CHECK_GENERATE_NANS
    if (funcRealInputHasNaN()){
        if (verb >= 1) print_funcRealInput(std::cout, false);
        if (verb >= 1) std::cout << "Input is broken at f = " << f <<  std::endl;
        return -2;
    }
#endif
    func->operator()(reinterpret_cast<const double**>(funcRealInput.data()), reinterpret_cast<double**>(out), nullptr, nullptr);

    if (dense_remap.size()){
        for (int i = dense_remap.size() - 1; i >= 0 && dense_remap[i] != i; --i){
            out[0][dense_remap[i]] = out[0][i];
            out[0][i] = 0;
        }
    }

#ifdef CHECK_GENERATE_NANS
    if (arrayHasNaN(*out, func->size1_out(0), func->size2_out(0))){
        if (verb >= 3) print_funcRealInput(std::cout, false);
        if (verb >= 3) printMatrix(std::cout, *out, func->size1_out(0), func->size2_out(0));
        if (verb >= 2) std::cout << "Result is broken at f = " << f << std::endl;
        return -1;
    }

    if (verb == 10){
        print_funcRealInput(std::cout, false);
        printMatrix(std::cout, *out, func->size1_out(0), func->size2_out(0));
        std::cout << "Result is printed for f = " << f << std::endl;
        std::cout << std::endl;
    }
#endif
    return 0;
}

void GenFunction::print_input() {
    std::cout << "\t" << "input:\n";
    for (auto& i: generatedInput)
        std::cout << "\t\t" << i.name << ": " << i.var.size1() << "x" << i.var.size2() << "\n";
}

void GenFunction::generate(const string &fname, string gen_dir) {
    SXVector input;
    input.reserve(generatedInput.size());
    for (auto& i: generatedInput)
        input.push_back(i.var);
    casadi::Dict dict;
#if __cplusplus >= 201703L
#define CONSTEXPR constexpr
#else
#define CONSTEXPR
#endif
    if CONSTEXPR (std::is_same<DReal, long double>())
        dict.insert({"casadi_real", "long double"});
    else if CONSTEXPR (std::is_same<DReal, double>())
        dict.insert({"casadi_real", "double"});
    else if CONSTEXPR (std::is_same<DReal, float>())
        dict.insert({"casadi_real", "float"});
    else throw std::runtime_error("Generated type for this DReal type is not defined");
#undef CONSTEXPR

#if __cplusplus >= 201703L
    using namespace std::filesystem;
#else
    using namespace std::experimental::filesystem;
#endif
    path path = gen_dir;
    auto s = status(path);
    if (!status_known(s) || s.type() != file_type::directory)
        create_directory(gen_dir);

    Function(fname, input, {func_expr}).generate(fname, dict);
    std::string mv = "";
    if (gen_dir != "") {
        mv = "mv " + fname + ".c " + gen_dir + ";";
        if (gen_dir[gen_dir.size()-1] != '/')
            gen_dir += "/";
    }
    std::string compile_command = " gcc -fPIC -shared -O3 "+gen_dir+fname+".c -o "+gen_dir+fname+".so";
    int flag = system((mv + compile_command).c_str());
    casadi_assert(flag==0, "Compilation failed");
}

void GenFunction::link_function(const string &fname, const string &gen_dir) {
    func = std::make_shared<Function>(casadi::external(fname, gen_dir+fname + ".so"));

    casadi::Sparsity s = func->sparsity_out(0);

    bool full = (s.nnz() == s.size1() * s.size2());
    if (!full) {
        for (int j = 0; j < s.size2(); ++j)
            for (int i = 0; i < s.size1(); ++i) {
                int ind = s.get_nz(i, j);
                if (ind != -1)
                    dense_remap.push_back(i + j * s.size1());
            }
    }
    else dense_remap.resize(0);
}

std::ostream & GenFunction::printMatrix(std::ostream &out, World3d::DReal *arr, int lda, int ldb) {
    using namespace std;
    vector<string> sm(lda * ldb);
    size_t maxlen = 0;
    for (int i = 0; i < lda; ++i)
        for (int j = 0; j < ldb; ++j) {
            std::ostringstream streamObj;
            streamObj << std::setprecision(std::numeric_limits<DReal>::digits10);
            streamObj << arr[j * lda + i];
            maxlen = max(maxlen, (sm[i * ldb + j] = streamObj.str()).length());
        }
    maxlen += 1;
    out << std::setprecision(std::numeric_limits<DReal>::digits10);
    for (int i = 0; i < lda; ++i) {
        for (int j = 0; j < ldb; ++j) {
            for (int k = 0, dif = maxlen - sm[i * ldb + j].length(); k < dif; ++k)
                out << " ";
            out << sm[i * ldb + j];
        }
        out << endl;
    }
    return out;
}

bool GenFunction::arrayHasNaN(World3d::DReal *arr, int lda, int ldb) {
    for (int i = 0; i < lda; ++i)
        for (int j = 0; j < ldb; ++j)
            if (std::isnan(arr[i * ldb + j]))
                return true;
    return false;
}

std::ostream &GenFunction::print_funcRealInput(std::ostream &out, bool oneline) {
    out << "RealInput[" << funcRealInput.size() << "]: ";
    if (!oneline) out << "\n";
    if (funcRealInput.size() != generatedInput.size()) {
        out << "Error: funcRealInput.size() != generatedInput.size()";
        if (!oneline) out << "\n";
        return out;
    }
    int ncnt = std::min(funcRealInput.size(), generatedInput.size());
    for (int i = 0; i < ncnt; ++i) {
        out << generatedInput[i].name << "[" << generatedInput[i].var.size1() << "]["
            << generatedInput[i].var.size2() << "] = {" << std::flush;
        for (int j = 0; j < generatedInput[i].var.size1(); ++j) {
            out << "{";
            if (generatedInput[i].var.size2())
                out << funcRealInput[i][0 * generatedInput[i].var.size1() + j];
            for (int k = 1; k < generatedInput[i].var.size2(); ++k)
                out << ", " << funcRealInput[i][k * generatedInput[i].var.size1() + j];
            out << "}";
            if (j != generatedInput[i].var.size1() - 1)
                out << ", ";
        }
        out << "}";
        if (i != ncnt - 1) out << (oneline ? ", " : "\n");
    }
    if (!oneline) out << "\n";
    return out;
}

bool GenFunction::funcRealInputHasNaN() {
    for (int i = 0, ncnt = funcRealInput.size(); i < ncnt; ++i) {
        for (int j = 0; j < generatedInput[i].var.size1(); ++j)
            for (int k = 0; k < generatedInput[i].var.size2(); ++k)
                if (std::isnan(funcRealInput[i][j * generatedInput[i].var.size2() + k]))
                    return true;
    }
    return false;
}

Invariant::Invariant(string name, SX val, SX C2dval, SX Qderiv, SX C2dderiv) :
        name{name},
        sym(SX::sym(name, val.size1(), val.size2())),
        val{val}, C2dval{C2dval},
        Qderiv{Qderiv}, C2dderiv{C2dderiv} {}

Invariant::Invariant(string name, SX sym, SX val, SX C2dval, SX Qderiv, SX C2dderiv) :
        name{name}, sym{sym}, val{val}, C2dval{C2dval}, Qderiv{Qderiv}, C2dderiv{C2dderiv} {}

void HH::printVars(const vector<DefaultGenInputVar> &vars, string prefix) {
    for (auto& i: vars) {
        SX sym = i.get_sym();
        std::cout << prefix << i.getName() << ": " << sym.size1() << " x " << sym.size2() << std::endl;
    }
}

void HH::printVars(const VarContainer<LocalVar> &vars, string prefix) {
    for (auto& i: vars) {
        SX sym = i.get_sym();
        std::cout << prefix << i.getName() << ": " << sym.size1() << " x " << sym.size2() << std::endl;
    }
}

void HH::printVars(const VarContainer<Invariant> &vars, string prefix) {
    for (auto& i: vars) {
        SX sym = i.get_sym();
        std::cout << prefix << i.getName() << ": " << sym.size1() << " x " << sym.size2() << std::endl;
    }
}

void HH::printVars(const vector<DefaultGenInputVar> &vg, const VarContainer<LocalVar> &vl, const VarContainer<Invariant> &vi, string prefix) {
    std::cout << prefix << "InputVars:\n";
    printVars(vg, prefix+"\t");
    std::cout << prefix << "LocalVars:\n";
    printVars(vl, prefix+"\t");
    std::cout << prefix << "Invariants:\n";
    printVars(vi, prefix+"\t");
}

void HH::printVars(const ArgsDef &v, string prefix) {
    printVars(get<0>(v), get<1>(v), get<2>(v), prefix);
}

ArgsDef HH::concat(ArgsDef &a, ArgsDef b) {
    ArgsDef res(std::move(a));
    for (auto& i: get<0>(b))
        get<0>(res).push_back(std::move(i));
    for (auto& i: get<1>(b))
        get<1>(res).push_back(std::move(i));
    for (auto& i: get<2>(b))
        get<2>(res).push_back(std::move(i));
    return res;
}

VarContainer<LocalVar> HH::getDefaultLocalVars() {
    static VarContainer<LocalVar> res;
    static bool req = false;
    if (!req){
        LocalVar Ap("Ap");
        LocalVar DD_in("DD_in", 6);
        LocalVar Pnt0("Pnt0", 3);
        LocalVar Pnt1("Pnt1", 3);
        LocalVar Pnt2("Pnt2", 3);
        LocalVar Pnt0_init("Pnt0_init", 3);
        LocalVar Pnt1_init("Pnt1_init", 3);
        LocalVar Pnt2_init("Pnt2_init", 3);
        LocalVar D("D", 3, 3); //D1, D2, D3
        LocalVar L_in("L_in", 6);

        SX S2 = SX::cross(Pnt1.sym - Pnt0.sym, Pnt2.sym - Pnt0.sym);
        SX S2n = SX::norm_2(S2);
        LocalVar Aq("Aq");//, S2n/2);
        LocalVar n("n", 3);//, S2/S2n);
        LocalVar np("np", 3);

        LocalVar Q("Q", SX::horzcat({Pnt0.sym, Pnt1.sym, Pnt2.sym}));
        LocalVar P("P", SX::horzcat({Pnt0_init.sym, Pnt1_init.sym, Pnt2_init.sym}));
        SX DD_val(3, 3);
        for (auto i = 0 ; i < 3; ++i)
            for (auto j = i; j < 3; ++j)
                DD_val(j, i) = DD_val(i, j) = DD_in.sym(j + (i + 1) * (i > 0), 0);
        LocalVar DD("DD", DD_val);
        LocalVar F("F", SX::mtimes(Q.sym, D.sym.T()));
        LocalVar C("C", SX::mtimes(F.sym.T(), F.sym));
        SX L_v(2, 3);
        for (int i = 0; i < 6; ++i)
            L_v(i%2, i/2) = L_in.sym(i, 0);
        LocalVar S2d("S2d", 2, 3);
        LocalVar L("L", L_v);
        LocalVar F2d("F2d", SX::mtimes({S2d.sym, Q.sym,  L.sym.T()})); //TODO: wrong value!!! S2d should be changable Euler shift matrix
        SX F2dC = SX::mtimes({Q.sym,  L.sym.T()});

        LocalVar C2d("C2d", SX::mtimes({F2dC.T(), F2dC})); //F2d.sym.T(),  F2d.sym

        res.push_back({Pnt0_init, Pnt1_init, Pnt2_init, Ap, DD_in, Pnt0, Pnt1, Pnt2, Aq, n, np,
                       D, L_in, Q, P, DD, F, S2d, C, L, F2d, C2d});

        req = true;
    }
    return res;
}

VarContainer<Invariant> HH::getDefaultInvariants() {
    static VarContainer<Invariant> res;
    static bool req = false;
    if (!req){
        VarContainer<LocalVar> v = getDefaultLocalVars();
        SX I1_v = SX::trace(SX::mtimes({v["DD"], v["Q"].T(), v["Q"]}));
        SX I1_C = SX::trace(v["C2d"]);
        SX I1_dQ = 2 * SX::mtimes({v["Q"], v["DD"].T()});
        SX I1_dC = SX::eye(2);
        res.push_back(Invariant{"I1", I1_v, I1_C, I1_dQ, I1_dC});
        SX J_v = v["Aq"] / v["Ap"];
        SX J_C = SX::sqrt(SX::det(v["C2d"]));
        SX Qnm0 = SX::cross(v["n"], v["Pnt2"] - v["Pnt1"]);
        SX Qnm1 = SX::cross(v["n"], v["Pnt0"] - v["Pnt2"]);
        SX Qnm2 = -(Qnm1 + Qnm0);
        SX Qnm = SX::horzcat({Qnm0, Qnm1, Qnm2});
//                SX Qnm = SX::horzcat({SX::cross(v["n"], v["Pnt2"] - v["Pnt1"]),
//                                      SX::cross(v["n"], v["Pnt0"] - v["Pnt2"]),
//                                      SX::cross(v["n"], v["Pnt1"] - v["Pnt0"])});
        SX J_dQ =  1.0 / (2 * v["Ap"]) * Qnm;
        SX J_dC = 1.0/(2 * J_C) * SX::adj(v["C2d"]).T();//1.0/2 * SX::adj(v["C2d"]);
        res.push_back(Invariant{"J", J_v, J_C, J_dQ, J_dC});
        SX I3_v = SX::sq(J_v);
        SX I3_C = SX::det(v["C2d"]);//SX::sq(J_C);
        SX I3_dQ = J_v / v["Ap"] * Qnm;
        SX I3_dC = SX::adj(v["C2d"]).T();//2 * J_C * J_dC;
        res.push_back(Invariant{"I3", I3_v, I3_C, I3_dQ, I3_dC});

        req = true;
    }
    return res;
}

HH::vector<DefaultGenInputVar> HH::getDefaultAvailableInputVars() {
    vector<DefaultGenInputVar> res;
    DefaultGenInputVar Ap("Ap");
    struct ApRequier: public DataRequier{
        ConstProperty<F_ind, DReal> Ap;
        void registerObj(Object3D* obj) override {
            Ap = set_Ap(*obj);
        }
        void saveData(const DReal** target, F_ind f, std::vector<V_ind>& around) override {
            *target = &Ap[f];
        }
        [[nodiscard]] std::unique_ptr<DataRequier> copy() const override { return std::make_unique<ApRequier>(*this); }
    };
    Ap.setDataRequier(std::make_unique<ApRequier>());
    DefaultGenInputVar Aq("Aq");
    struct AqRequier: public DataRequier{
        UpdatableProperty<F_ind, Vector> S;
        DReal Aq = -1;
        void registerObj(Object3D* obj) override {
            S = set_S(*obj);
        }
        void saveData(const DReal** target, F_ind f, std::vector<V_ind>& around) override {
            Eigen::Matrix<DReal, 3, 1> s{S[f][0], S[f][1], S[f][2]};
            Aq = s.stableNorm();//sqrt(S[f].squared_length());
            *target = &Aq;
        }
        [[nodiscard]] std::unique_ptr<DataRequier> copy() const override { return std::make_unique<AqRequier>(*this); }
    };
    Aq.setDataRequier(std::make_unique<AqRequier>());

    DefaultGenInputVar n("n", 3, 1);
    struct nRequier: public DataRequier{
        Object3D* _obj = nullptr;
        void registerObj(Object3D* obj) override {
            _obj = obj;
        }
        void saveData(const DReal** target, F_ind f, std::vector<V_ind>& around) override {
            *target = &_obj->m_normal[f][0];
        }
        [[nodiscard]] std::unique_ptr<DataRequier> copy() const override { return std::make_unique<nRequier>(*this); }
    };
    n.setDataRequier(std::make_unique<nRequier>());

    DefaultGenInputVar np("np", 3, 1);
    struct npRequier: public DataRequier{
        Object3D* _obj = nullptr;
        std::array<DReal, 3> dat = {0};
        void registerObj(Object3D* obj) override {
            _obj = obj;
        }
        void saveData(const DReal** target, F_ind f, std::vector<V_ind>& around) override {
            auto& v = around;
            Vector np = CGAL::cross_product(_obj->m_x0[v[1]] - _obj->m_x0[v[0]], _obj->m_x0[v[2]] - _obj->m_x0[v[0]]);
            np /= sqrt(np.squared_length());
            dat = {np[0], np[1], np[2]};
            *target = dat.data();
        }
        [[nodiscard]] std::unique_ptr<DataRequier> copy() const override { return std::make_unique<npRequier>(*this); }
    };
    np.setDataRequier(std::make_unique<npRequier>());

    DefaultGenInputVar DD_in("DD_in", 6, 1);
    struct DDRequier: public DataRequier{
        ConstProperty<F_ind, std::array<DReal, 6>> DD;
        void registerObj(Object3D* obj) override {
            DD = set_DD_matrix(*obj);
        }
        void saveData(const DReal** target, F_ind f, std::vector<V_ind>& around) override {
            *target = &DD[f][0];
        }
        [[nodiscard]] std::unique_ptr<DataRequier> copy() const override { return std::make_unique<DDRequier>(*this); }
    };
    DD_in.setDataRequier(std::make_unique<DDRequier>());

    DefaultGenInputVar Pnt[6];
    for (int i = 0; i < 3; ++i){
        Pnt[i] = DefaultGenInputVar("Pnt" + std::to_string(i), 3, 1, std::make_unique<PntRequier>(i));
        Pnt[3+i] = DefaultGenInputVar("Pnt" + std::to_string(i) + "_init", 3, 1, std::make_unique<PntInitRequier>(i));
    }
    struct S2dRequier: public DataRequier{
        bool flat;
        ConstProperty<F_ind, std::array<DReal, 6>> S2dAll;
        std::array<DReal, 6> S2dFlat;
        void registerObj(Object3D* obj) override {
            flat = check_flat_initial_template(*obj);
            LocalCartesian S(*obj, flat);
            if (flat){
                Eigen::Matrix<DReal, 2, 3> q = S(*(obj->m_mesh.faces().begin()));
                for (int i = 0; i < 6; ++i)
                    S2dFlat[i] = q(i%2, i/2);
            } else {
                S2dAll = set_S2d(*obj, flat);
            }
        }
        void saveData(const DReal** target, F_ind f, std::vector<V_ind>& around) override {
            if (flat) *target = S2dFlat.data();
            else *target = &S2dAll[f][0];
        }
        [[nodiscard]] std::unique_ptr<DataRequier> copy() const override { return std::make_unique<S2dRequier>(*this); }
    };
    DefaultGenInputVar S2d("S2d", 2, 3, std::make_unique<S2dRequier>());
    DefaultGenInputVar L_in("L_in", 6, 1);
    struct LRequier: public DataRequier{
        ConstProperty<F_ind, std::array<DReal, 6>> L;
        void registerObj(Object3D* obj) override {
            L = set_L(*obj);
        }
        void saveData(const DReal** target, F_ind f, std::vector<V_ind>& around) override {
            *target = &L[f][0];
        }
        [[nodiscard]] std::unique_ptr<DataRequier> copy() const override { return std::make_unique<LRequier>(*this); }
    };
    L_in.setDataRequier(std::make_unique<LRequier>());
    struct DRequier: public DataRequier{
        ConstProperty<F_ind, std::array<Vector, 3>> D;
        std::array<DReal, 9> data;
        void registerObj(Object3D* obj) override {
            D = set_D_vecs(*obj);
        }
        void saveData(const DReal** target, F_ind f, std::vector<V_ind>& around) override {
            std::array<Vector, 3>& q = D[f];
            for (int i = 0; i < 9; ++i)
                data[i] = q[i/3][i%3];
            *target = data.data();
        }
        [[nodiscard]] std::unique_ptr<DataRequier> copy() const override { return std::make_unique<DRequier>(*this); }
    };
    DefaultGenInputVar D("D", 3, 3, std::make_unique<DRequier>());
    res.push_back(std::move(Ap));
    res.push_back(std::move(Aq));
    res.push_back(std::move(n));
    res.push_back(std::move(np));
    res.push_back(std::move(DD_in));
    for (int i = 0; i < 6; ++i)
        res.push_back(std::move(Pnt[i]));
    res.push_back(std::move(D));
    res.push_back(std::move(S2d));
    res.push_back(std::move(L_in));

    return std::move(res);
}

HH::vector<DefaultGenInputVar> HH::makeOrthotropicInvariantInput(
        string vecName, string vec2dName, string vecdtAName, string vectorTag) {
    DefaultGenInputVar A(vecName, 3, 1, std::make_unique<VectorRequier<3>>(vectorTag));
    struct SARequier: public DataRequier{
        ConstProperty<F_ind, std::array<DReal, 2>> SA;
        string tagName;
        void registerObj(Object3D* obj) override {
            auto A = obj->m_mesh.property_map<F_ind, std::array<DReal, 3>>(tagName);
            if (!A.second) throw std::runtime_error("Tag \"" +  tagName + "\" is not exists on the object \"" + obj->name + "\"");
            SA = obj->m_mesh.add_property_map<F_ind, std::array<DReal, 2>>(tagName+"::SA");
//                    if (!SA.second) return;
            bool flat = check_flat_initial_template(*obj);
            LocalCartesian S(*obj, flat);
            for (auto f: obj->m_mesh.faces()){
                Eigen::Matrix<DReal, 3, 1> a{A.first[f][0], A.first[f][1], A.first[f][2]};
                Eigen::Matrix<DReal, 2, 1> v = S(f) * a;
                SA[f][0] = v[0], SA[f][1] = v[1];
            }
            return;
        }
        void saveData(const DReal** target, F_ind f, std::vector<V_ind>& around) override {
            *target = &SA[f][0];
        }
        SARequier(string tagName): tagName{tagName} {}
        [[nodiscard]] std::unique_ptr<DataRequier> copy() const override { return std::make_unique<SARequier>(*this); }
    };
    DefaultGenInputVar SA(vec2dName, 2, 1, std::make_unique<SARequier>(vectorTag));
    struct dQRequier: public DataRequier{
        ConstProperty<F_ind, std::array<DReal, 3>> dQ;
        string tagName;
        void registerObj(Object3D* obj) override {
            auto A_it = obj->m_mesh.property_map<F_ind, std::array<DReal, 3>>(tagName);
            if (!A_it.second) throw std::runtime_error("Tag \"" +  tagName + "\" is not exists on the object \"" + obj->name + "\"");
            auto D_it = set_D_vecs(*obj);
            dQ = obj->m_mesh.add_property_map<F_ind, std::array<DReal, 3>>(tagName+"::DtA");
//                    if (!dQ.second) return;
            for (auto f: obj->m_mesh.faces()){
                auto& d = D_it[f];
                auto& a = A_it.first[f];
                Eigen::Matrix<DReal, 3, 3> D;
                D <<    d[0][0], d[1][0], d[2][0],
                        d[0][1], d[1][1], d[2][1],
                        d[0][2], d[1][2], d[2][2];
                Eigen::Matrix<DReal, 3, 1> A {a[0], a[1], a[2]};
                Eigen::Matrix<DReal, 3, 1> M = D.transpose() * A;
                for (int i = 0; i < 3; ++i)
                    dQ[f][i] = M[i];
            }
            if (D_it.second) obj->m_mesh.remove_property_map(D_it.first);
        }
        void saveData(const DReal** target, F_ind f, std::vector<V_ind>& around) override {
            *target = &dQ[f][0];
        }
        dQRequier(string tagName): tagName{tagName} {}
        [[nodiscard]] std::unique_ptr<DataRequier> copy() const override { return std::make_unique<dQRequier>(*this); }
    };
    DefaultGenInputVar dQ(vecdtAName, 3, 1, std::make_unique<dQRequier>(vectorTag));
    vector<DefaultGenInputVar> def;
    def.push_back(std::move(A));
    def.push_back(std::move(SA));
    def.push_back(std::move(dQ));
    return def;
}

HH::vector<DefaultGenInputVar> HH::makeOrthotropicInvariantInput(string vecName, string vec2dName, string vecdtAName, std::array<DReal, 3> Av) {
    DefaultGenInputVar A(vecName, 3, 1, std::make_unique<ConstVectorRequier<3>>(Av));
    struct SARequier: public DataRequier{
        ConstProperty<F_ind, std::array<DReal, 2>> SA;
        std::array<DReal, 3> A;
        std::vector<DReal> SAFlat;
        string prefix;
        void registerObj(Object3D* obj) override {
//                    if (!SA.second) return;
            bool flat = check_flat_initial_template(*obj);
            LocalCartesian S(*obj, flat);
            if (!flat) {
                SA = obj->m_mesh.add_property_map<F_ind, std::array<DReal, 2>>("f:" + prefix + "::SA");
                for (auto f: obj->m_mesh.faces()) {
                    Eigen::Matrix<DReal, 3, 1> a{A[0], A[1], A[2]};
                    Eigen::Matrix<DReal, 2, 1> v = S(f) * a;
                    SA[f][0] = v[0], SA[f][1] = v[1];
                }
            } else {
                SAFlat.resize(2);
                Eigen::Matrix<DReal, 3, 1> a{A[0], A[1], A[2]};
                Eigen::Matrix<DReal, 2, 1> v = S(*(obj->m_mesh.faces().begin())) * a;
                SAFlat[0] = v[0], SAFlat[1] = v[1];
            }
            return;
        }
        SARequier(std::array<DReal, 3> A, string prefix = ""): A{A}, prefix{prefix} {}
        void saveData(const DReal** target, F_ind f, std::vector<V_ind>& around) override {
            if (SAFlat.empty())
                *target = &SA[f][0];
            else
                *target = &SAFlat[0];
        }
        [[nodiscard]] std::unique_ptr<DataRequier> copy() const override { return std::make_unique<SARequier>(*this); }
    };
    DefaultGenInputVar SA(vec2dName, 2, 1, std::make_unique<SARequier>(Av, vec2dName));
    struct dQRequier: public DataRequier{
        ConstProperty<F_ind, std::array<DReal, 3>> dQ;
        std::array<DReal, 3> Av;
        string prefix;
        void registerObj(Object3D* obj) override {
            auto D_it = set_D_vecs(*obj);
            dQ = obj->m_mesh.add_property_map<F_ind, std::array<DReal, 3>>("f:" + prefix+"::DtA");
//                    if (!dQ.second) return;
            for (auto f: obj->m_mesh.faces()){
                auto& d = D_it[f];
                Eigen::Matrix<DReal, 3, 3> D;
                D <<    d[0][0], d[1][0], d[2][0],
                        d[0][1], d[1][1], d[2][1],
                        d[0][2], d[1][2], d[2][2];
                Eigen::Matrix<DReal, 3, 1> A {Av[0], Av[1], Av[2]};
                Eigen::Matrix<DReal, 3, 1> M = D.transpose() * A;
                for (int i = 0; i < 3; ++i)
                    dQ[f][i] = M[i];
            }
            if (D_it.second) obj->m_mesh.remove_property_map(D_it.first);
        }
        void saveData(const DReal** target, F_ind f, std::vector<V_ind>& around) override {
            *target = &dQ[f][0];
        }
        dQRequier(std::array<DReal, 3> A, string prefix = ""): Av{A}, prefix{prefix} {}
        [[nodiscard]] std::unique_ptr<DataRequier> copy() const override { return std::make_unique<dQRequier>(*this); }
    };
    DefaultGenInputVar dQ(vecdtAName, 3, 1, std::make_unique<dQRequier>(Av, vecdtAName));
    vector<DefaultGenInputVar> def;
    def.push_back(std::move(A));
    def.push_back(std::move(SA));
    def.push_back(std::move(dQ));
    return def;
}

ArgsDef HH::makeAnisotropicInvariant(string InvName, string vectorTag1, string vectorTag2){
    std::array<std::string, 2> 
        vecName = {InvName + "_A", InvName + "_B"},
        vec2dName = {InvName + "_A2d", InvName + "_B2d"},
        vecdQName = {InvName + "_DtA", InvName + "_DtB"},
        vectorTag = {vectorTag1, vectorTag2};
    auto def = makeOrthotropicInvariantInput(vecName[0], vec2dName[0], vecdQName[0], vectorTag[0]);
    auto def2 = makeOrthotropicInvariantInput(vecName[1], vec2dName[1], vecdQName[1], vectorTag[1]);
    for (std::size_t i = 0; i < def2.size(); ++i)
        def.push_back(std::move(def2[i]));
    def2.clear();

    auto li = makeAnisotropicInvariantVars(InvName, vecName, vec2dName, vecdQName);
    VarContainer<LocalVar>& loc = li.first;
    VarContainer<Invariant>& inv = li.second;

    return ArgsDef{std::move(def), std::move(loc), std::move(inv)};    
}

ArgsDef HH::makeOrthotropicInvariant(string InvName, string vectorTag) {
    string vecName = InvName + "_A";
    string vec2dName = InvName + "_A2d";
    string vecdQName = InvName + "_DtA";
    auto def = makeOrthotropicInvariantInput(vecName, vec2dName, vecdQName, vectorTag);
    auto li = makeOrthotropicInvariantVars(InvName, vecName, vec2dName, vecdQName);
    VarContainer<LocalVar>& loc = li.first;
    VarContainer<Invariant>& inv = li.second;

    return ArgsDef{std::move(def), std::move(loc), std::move(inv)};
}

ArgsDef HH::makeOrthotropicInvariant(string InvName, std::array<DReal, 3> A) {
    string vecName = InvName + "_A";
    string vec2dName = InvName + "_A2d";
    string vecdQName = InvName + "_DtA";
    auto def = makeOrthotropicInvariantInput(vecName, vec2dName, vecdQName, A);
    auto li = makeOrthotropicInvariantVars(InvName, vecName, vec2dName, vecdQName);
    VarContainer<LocalVar>& loc = li.first;
    VarContainer<Invariant>& inv = li.second;

    return ArgsDef{std::move(def), std::move(loc), std::move(inv)};
}

std::pair<VarContainer<LocalVar>, VarContainer<Invariant>> HH::makeOrthotropicInvariantVars(
        string InvName, string vecName, string vec2dName, string vecdtAName){
    LocalVar A_l(vecName, 3);
    LocalVar SA_l(vec2dName, 2);
    LocalVar dQ_l(vecdtAName, 3);
    VarContainer<LocalVar> loc;
    loc.push_back({A_l, SA_l, dQ_l});

    VarContainer<Invariant> inv;
    VarContainer<LocalVar> v = getDefaultLocalVars();
    SX I_v = SX::mtimes(SXVector{v["Q"], dQ_l.sym});
    I_v =  SX::sq(I_v(0, 0)) + SX::sq(I_v(1, 0)) + SX::sq(I_v(2, 0));
    SX I_C = SX::mtimes({SA_l.sym.T(), v["C2d"], SA_l.sym});
    SX I_dQ = 2 * SX::mtimes({v["Q"], dQ_l.sym, dQ_l.sym.T()});
    SX I_dC = SX::mtimes(SA_l.sym, SA_l.sym.T());
    inv.push_back(Invariant{InvName, I_v, I_C, I_dQ, I_dC});
    return std::pair<VarContainer<LocalVar>, VarContainer<Invariant>>{std::move(loc), std::move(inv)};
}

std::pair<VarContainer<LocalVar>, VarContainer<Invariant>> HH::makeAnisotropicInvariantVars(
        string InvName, std::array<string, 2> vecName, std::array<string, 2> vec2dName, std::array<string, 2> vecdtAName){
    LocalVar A_l(vecName[0], 3), B_l(vecName[1], 3);
    LocalVar SA_l(vec2dName[0], 2), SB_l(vec2dName[1], 2);
    LocalVar dQA_l(vecdtAName[0], 3), dQB_l(vecdtAName[1], 3);
    VarContainer<LocalVar> loc;
    loc.push_back({A_l, SA_l, dQA_l, B_l, SB_l, dQB_l});

    VarContainer<Invariant> inv;
    VarContainer<LocalVar> v = getDefaultLocalVars();
    SX Ff = SX::mtimes(SXVector{v["Q"], dQA_l.sym}), Fs = SX::mtimes(SXVector{v["Q"], dQB_l.sym});
    SX I_v = SX::mtimes(Ff.T(), Fs);
    SX I_C = SX::mtimes({SA_l.sym.T(), v["C2d"], SB_l.sym});
    SX I_dQ = SX::mtimes({v["Q"], dQA_l.sym, dQB_l.sym.T()}) + SX::mtimes({v["Q"], dQB_l.sym, dQA_l.sym.T()});
    SX I_dC = (SX::mtimes(SA_l.sym, SB_l.sym.T()) + SX::mtimes(SB_l.sym, SA_l.sym.T())) / 2;
    inv.push_back(Invariant{InvName, I_v, I_C, I_dQ, I_dC});
    return std::pair<VarContainer<LocalVar>, VarContainer<Invariant>>{std::move(loc), std::move(inv)};
}

void HH::correlateInputAndVariables(const vector<DefaultGenInputVar> &generatedInput,
                                                     VarContainer<LocalVar> &localVars) {
    for (const auto& i: generatedInput){
        auto lv = localVars.find(i.getName());
        if (lv != localVars.end())
            lv->val = i.var;
    }
}

SX HH::applyVariablesToExpr(SX expr, VarContainer<LocalVar> &loc_vars) {
    for (int i = loc_vars.size() - 1; i >= 0; --i){
        expr = SX::simplify(SX::substitute(expr, loc_vars.at(i).sym, loc_vars.at(i).val));
    }
    return expr;
}

void HH::correlateInputAndVariables(const vector<DefaultGenInputVar> &generatedInput, VarContainer<LocalVar> &localVars,
                                VarContainer<Invariant> &invariants) {
    for (const auto& i: generatedInput){
        auto lv = localVars.find(i.getName());
        if (lv != localVars.end())
            lv->val = i.var;
        auto it = std::find_if(invariants.seq.begin(), invariants.seq.end(), [name = i.getName()](const Invariant i) {return i.getName() == name;});
        if (it != invariants.seq.end()) {
            it->val = i.var;
        }
    }
}

GenFunction HH::makeJacobianFunction(string f_name, string save_dir, GenFunction func, SX var, bool regenerate, bool print_gen) {
    GenFunction genJacobian;
    genJacobian.generatedInput = func.generatedInput;

    SX jacob = func.func_expr;
    jacob = SX::reshape(jacob, jacob.size1()*jacob.size2(), 1);
    SX vars = SX::reshape(var, var.size1()*var.size2(), 1);
    jacob = SX::jacobian(jacob, vars);
    genJacobian.func_expr = jacob.T();

    genJacobian.update_input();
    if (regenerate){
        std::cout << "Generate jacobian for the function: " << func.func->name() << std::endl;
        if (print_gen) genJacobian.print_input();
        genJacobian.generate(f_name + "_Jacob", save_dir);
    }
    genJacobian.link_function(f_name + "_Jacob", save_dir);
    return genJacobian;
}

HH::ArgsDef HH::makeDefaultInits(const DefaultHyperElasticForce& f, const VarContainer<LocalVar>& local_vars){
    ArgsDef res;
    auto& gen_input = get<0>(res);
    auto& gen_loc = get<1>(res);
    auto& invariants = get<2>(res);
    gen_input = HH::getDefaultAvailableInputVars();
    vector<DefaultGenInputVar> input_vars = f.inputVars;
    for (int i = 0; i < input_vars.size(); ++i){
        auto it = std::find_if(gen_input.begin(), gen_input.end(),
                               [nm = input_vars[i].getName()] (DefaultGenInputVar& v){ return v.getName() == nm;});
        if (it == gen_input.end())
            gen_input.push_back(std::move(input_vars[i]));
        else
            *it = std::move(input_vars[i]);
    }

    gen_loc = local_vars;
    auto loc_vars = f.localVars;
    for (int i = 0; i < loc_vars.size(); ++i){
        auto it = gen_loc.find(loc_vars.at(i).getName());
        if (it == gen_loc.end())
            gen_loc.push_back(loc_vars.at(i));
        else
            *it = loc_vars.at(i);
    }
    invariants = f.invariants;
    return res;
}

DefaultHyperElasticForce::DefaultHyperElasticForce(const DefaultHyperElasticForce &f) {
    inputVars = f.inputVars;
    localVars = f.localVars;
    invariants = f.invariants;
    U = f.U;
    genForce = f.genForce;
    genJacobian = f.genJacobian;
    gen_dir = f.gen_dir;
    force_name = f.force_name;
    _obj = nullptr;
    force_verb = f.force_verb; jacob_verb = f.jacob_verb;
#ifdef USE_DEBUG_WATCHES
    m_watches = f.m_watches;
#endif
}

DefaultHyperElasticForce &DefaultHyperElasticForce::operator=(const DefaultHyperElasticForce &f) {
    if (this == &f) return *this;
    inputVars = f.inputVars;
    localVars = f.localVars;
    invariants = f.invariants;
    U = f.U;
    genForce = f.genForce;
    genJacobian = f.genJacobian;
    gen_dir = f.gen_dir;
    force_name = f.force_name;
    _obj = nullptr;
    force_verb = f.force_verb; jacob_verb = f.jacob_verb;
#ifdef USE_DEBUG_WATCHES
    m_watches = f.m_watches;
#endif

    return *this;
}

DefaultHyperElasticForce::DefaultHyperElasticForce(string f_name, string save_dir,
       vector<DefaultGenInputVar> input_vars, VarContainer<LocalVar> loc_vars, VarContainer<Invariant> invariants, SX U,
       bool regenerate) : localVars{loc_vars}, invariants{invariants}, gen_dir{save_dir}, force_name{f_name}, U{U}
{
    if (gen_dir != "" && gen_dir[gen_dir.size()-1] != '/')
        gen_dir += "/";
    inputVars = input_vars;
    genForce.generatedInput = std::move(input_vars);
    correlateInputAndVariables();
    genForce.func_expr = compute_force();
    genForce.update_input();
    if (regenerate)
        generate();
    link_function();
}

void DefaultHyperElasticForce::prepareJacobianFunction(bool regenerate) {
    auto get_input = [&gi = genForce.generatedInput](string name){
        auto it = std::find_if(gi.begin(), gi.end(),
                            [&name](DefaultGenInputVar& v) { return v.getName() == name;});
//        if(it == gi.end()) throw std::runtime_error("Can't find \"" + name + "\" in generated input [size=" + std::to_string(gi.size()) + "]");
        return it;
    };
    std::array<SX, 3> q = {get_input("Pnt0")->var, get_input("Pnt1")->var, get_input("Pnt2")->var};
    SX Q = SX::horzcat({q[0], q[1], q[2]});
    SX S2 = SX::cross(q[1] - q[0], q[2] - q[0]);
    SX S2n = SX::norm_2(S2);
    SX rforce = genForce.func_expr;
    {
        auto it = get_input("Aq");
        if (it != genForce.generatedInput.end()){
            SX Aq = it->var;
            genForce.func_expr = SX::substitute(genForce.func_expr, Aq, S2n/2);
        }
        it = get_input("n");
        if (it != genForce.generatedInput.end()){
            SX n = it->var;
            genForce.func_expr = SX::simplify(SX::substitute(genForce.func_expr, n, S2/S2n));
        }
    }
    genJacobian = HH::makeJacobianFunction(force_name, gen_dir, genForce, Q, regenerate, true);
    genForce.func_expr = rforce;
}

SX DefaultHyperElasticForce::compute_force() {
    SX force(3, 3);
    SX u = U;
    for (const auto& i: invariants){
        force -= SX::simplify((SX::jacobian(u, i.sym) * i.Qderiv));
        force = SX::simplify(SX::substitute(force, i.sym, i.val));
        u = SX::substitute(u, i.sym, i.val);
    }
    for (int i = localVars.size() - 1; i >= 0; --i){
        force = SX::simplify(SX::substitute(force, localVars.at(i).sym, localVars.at(i).val));
    }
    return force;
}

SX DefaultHyperElasticForce::compute_HS_tensor_expr(VarContainer<LocalVar>& lVars, VarContainer<Invariant>& lInvs, casadi::SX U){
    SX u = SX::simplify(U / lVars["Ap"]);
    SX S(2, 2);
    for (const auto& i: lInvs){
        S += SX::simplify((SX::jacobian(u, i.sym) * i.C2dderiv));
        S = SX::simplify(SX::substitute(S, i.sym, i.C2dval));
        u = SX::substitute(u, i.sym, i.C2dval);
    }
    S *= 2;
    for (int i = lVars.size() - 1; i >= 0; --i){
        S = SX::simplify(SX::substitute(S, lVars.at(i).sym, lVars.at(i).val));
    }
    return S;
}

SX DefaultHyperElasticForce::compute_S_tensor_expr(VarContainer<LocalVar>& lVars, VarContainer<Invariant>& lInvs, casadi::SX U){
    return SX::simplify(compute_HS_tensor_expr(lVars, lInvs, U) / lVars["H"]);
}

SX DefaultHyperElasticForce::compute_jacobian() {
    SX jacob = genForce.func_expr;
    jacob = SX::reshape(jacob, jacob.size1()*jacob.size2(), 1);
    SX Q = localVars["Q"];
    SX vars = SX::reshape(Q, Q.size1()*Q.size2(), 1);
    jacob = SX::jacobian(jacob, vars.T());
    return jacob;
}

void DefaultHyperElasticForce::generate() {
    std::cout << "Generate elastic force for potential: " << U << std::endl;
    genForce.print_input();
    genForce.generate(force_name, gen_dir);
}

void DefaultHyperElasticForce::generate_jacobian() {
    std::cout << "Generate jacobian for elastic force of potential: " << U << std::endl;
    genJacobian.print_input();
    genJacobian.generate(force_name + "_Jacob", gen_dir);
}

int DefaultHyperElasticForce::operator()(Object3D &obj) {
    if (&obj != _obj) registerObj(&obj);
    auto& m = obj.m_mesh;
    DReal force[9];
    DReal* pforce = &(force[0]);
    std::vector<V_ind>& v = _vdat; v.resize(3);
    int status = 0;
    for (auto f: m.faces()) {
        vert_around(m, f, v);
        if ((status = genForce(&pforce, f, v, force_verb-1)) < 0) {
            if (force_verb >= 1) std::cout << "Error(" << status << ") in generic hyperelastic force \"" << force_name << "\" at f = " << f << " on obj \"" << obj.name << "\"" << std::endl;
            return status;
        }
        for (auto i = 0; i < 3; ++i) {
            auto bc = obj.getBC(v[i]);
            obj.m_F[v[i]] += Vector((bc & 1) ?force[3 * i + 0] : 0 , (bc & 2) ?force[3 * i + 1] : 0, (bc & 4) ?force[3 * i + 2] : 0);
        }
    }
    return status;
}

int DefaultHyperElasticForce::element_matrix(Object3D &obj, ForceBase::LocMatrix &matr, F_ind element) {
    if (&obj != _jobj) registerjObj(&obj);
    auto& m = obj.m_mesh;
    matr.setZero();
    matr.resize(3*3, 3*3);
    DReal* pforce = matr.dat.data();
    std::vector<V_ind>& v = _vdat; v.resize(3);
    vert_around(m, element, v); //remove debugVar
    int genStat = genJacobian(&pforce, element, v, jacob_verb-1);
    if (genStat < 0){
        if (jacob_verb >= 1) std::cout << "Error(" << genStat << ") in jacobian of generic hyperelastic force \"" << force_name << "\" at f = " << element << " on obj \"" << obj.name << "\"" << std::endl;
        return genStat;
    }
    for (auto i = 0; i < 3; ++i) {
        auto bc = obj.getBC(v[i]);
        for (int l = 0; l < 3; ++l){
            if (!(bc & (1 << l)))
                matr.ApplyDirichlet(3 * i + l);
        }
    }
    return genStat;
}

int DefaultHyperElasticForce::add_watch(const std::function<SX(VarContainer<LocalVar>&, VarContainer<Invariant>&, casadi::SX U)>& watch_expr, bool quiet){
#ifdef USE_DEBUG_WATCHES
    GenFunction watch;
    watch.set_generatedInput(inputVars);
    watch.func_expr = watch_expr(localVars, invariants, U);
    watch.update_input();
    static int nwatch = 0;
    m_watches.emplace_back(std::move(watch));
    string watch_name = "HyperElastic_watch" + std::to_string(nwatch);
    if (!quiet){
        std::cout << "Generate watch \"" << watch_name << "\" on \"" << force_name << "\" with input: \n";
        m_watches.back().print_input();
    }
    m_watches.back().generate(watch_name, gen_dir);
    m_watches.back().link_function(watch_name, gen_dir);
    nwatch++;
#else
    throw std::runtime_error("This function is not avaliable in current build. Please add symbol -DUSE_DEBUG_WATCHES ");
#endif
    return m_watches.size() - 1;
}

#ifdef USE_DEBUG_WATCHES
Eigen::MatrixXd DefaultHyperElasticForce::eval_watch(Object3D& obj, F_ind f, GenFunction* watch, bool print){
    Eigen::MatrixXd mres;
    if (!watch) return mres;
    int w1s = watch->func->numel_out();
    long long ssr = watch->func->sparsity_out(0).rows(),
              ssc = watch->func->sparsity_out(0).columns();
    mres.resize(ssr, ssc);
    DReal* pforce = mres.data();
    _vdat.resize(3);
    auto& m = obj.m_mesh;
    auto vv = vert_around(m, f);
    std::copy(vv.begin(), vv.end(), _vdat.data());
    auto& v = _vdat;
    {
        watch->registerObj(&obj);
        watch->operator()(&pforce, f, v, print ? 10 : 0);
    }
    return mres;
}
Eigen::MatrixXd DefaultHyperElasticForce::eval_watch(Object3D& obj, F_ind f, int watch_num, bool print){
    return eval_watch(obj, f, m_watches.data() + watch_num, print);
}
std::vector<Eigen::MatrixXd> DefaultHyperElasticForce::eval_watches(Object3D& obj, F_ind f, bool print){
    std::vector<Eigen::MatrixXd> res(m_watches.size());
    for (int i = 0; i < m_watches.size(); ++i){
        if (print) std::cout << "Watch " << i << ":\n";
        res[i] = eval_watch(obj, f, &(m_watches[i]), print);
    }
    return res;
}
#endif

void HyperElasticForceBase::prepareJacobianFunction(bool regenerate) {
    f.prepareJacobianFunction(regenerate);
    _with_jacob = true;
}

std::unique_ptr<ForceBase> HyperElasticForceBase::clone() {
    auto res = copy();
    if (_with_jacob) res->prepareJacobianFunction(false);
    return res;
}