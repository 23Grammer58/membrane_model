//
// Created by alex on 02.02.2021.
//

#include "DataGatherer.h"
#include "../ForceAppliers/ForceAppliers.h"
using namespace World3d;
using namespace HH;

SX HEDataGatherer::compute_response_expression(const LocalVars &loc_vars, const Invariants &invariants) const{
    SX u = U;
    auto Ap = loc_vars["Ap"];
    u = SX::simplify(u / Ap);
    SX S(2, 2);

    for (const auto& i: invariants){
        S += SX::simplify((SX::jacobian(u, i.sym) * i.C2dderiv));
        S = SX::simplify(SX::substitute(S, i.sym, i.val));
        u = SX::substitute(u, i.sym, i.val);
    }
    S *= 2;
    auto H_it = loc_vars.find("H");
    if (H_it == loc_vars.end())
        throw std::runtime_error("To use this force you should provide parameter \"H\" that means height");
    SX H = H_it->sym;
    SX C2d = loc_vars["C2d"];
    SX F11 = SX::sqrt(C2d(0, 0));
    SX F12 = C2d(1, 0) / F11;
    SX F22 = sqrt(C2d(1, 1) - SX::sq(F12));
    SX F = SX::vertcat({SX::horzcat({F11, F12}),
                        SX::horzcat({0, F22})});

    SX Fr = loc_vars["F"];
    SX F2d = SX::mtimes({loc_vars["Q"], loc_vars["L"].T()});
    SX f1 = SX::horzsplit(F2d, 1)[0],  f2 = SX::horzsplit(F2d, 1)[1];
    SX e[3];
    e[0] = f1/SX::norm_2(f1);
    e[1] = f2 - SX::dot(f2, e[0])*e[0]; e[1] /= SX::norm_2(e[1]);
    e[2] = SX::cross(e[0], e[1]);
    SX S2d = loc_vars["S2d"];
    SX E[3];
    E[0] = SX::vertsplit(S2d, 1)[0].T();
    E[1] = SX::vertsplit(S2d, 1)[1].T();
    E[2] = SX::cross(E[0], E[1]);
    SX J = loc_vars["Aq"] / loc_vars["Ap"];
//    SX g[] = {loc_vars["Pnt1"] - loc_vars["Pnt0"], loc_vars["Pnt2"] - loc_vars["Pnt0"]};
//    g[0] /= SX::norm_2(g[0]); g[1] /= SX::norm_2(g[1]); g[1] = SX::cross(SX::cross(g[0], g[1]), g[1]);
//    g[1] /= SX::norm_2(g[1]);
//    SX R(2, 2);
//    for (int i = 0; i < 2; ++i)
//    for (int j = 0; j < 2; ++j)
//        R(i, j) = SX::dot(g[i], e[j]);

//    SX R = SX::horzcat({e[0], e[1]});
//    SX JT = SX::mtimes({F, S, F.T()});
//    SX Tr = JT / J; //SX::mtimes({R.T(), JT / J, R});
//    SX T3d = SX::mtimes({R, Tr, R.T()});
    SX em = SX::horzcat({e[0], e[1], e[2]}), Em = SX::horzcat({E[0], E[1], E[2]});
    SX Rq = SX::mtimes({em});//, Em.T()});
    SX Cq = SX::horzcat({SX::eye(2), SX(2, 1)});
    SX JT3d = SX::mtimes({F2d, S, F2d.T()});
    SX JT = SX::mtimes({Rq.T(), JT3d, Rq});//TODO: SX::mtimes({Rq, JT3d, Rq.T()});
    SX Tr = JT / J;
    SX T3d = JT3d / J;

    SX Tcomp = SX::vertcat({Tr(0,0), Tr(1, 1), Tr(1, 0)});
    SX T3dcomp = SX::vertcat({T3d(0, 0), T3d(1, 1), T3d(2, 2), T3d(0, 1), T3d(0, 2), T3d(1, 2)});

    SX response = SX::vertcat({JT(0, 0), JT(1, 1), F(0, 0) / F(1, 1) * JT(1, 0) });
    SX ksi = SX::vertcat({0.5 * SX::log(C2d(0, 0)),
                          0.5 * SX::log(C2d(1, 1) - SX::sq(C2d(1, 0))/C2d(0, 0)),
                          C2d(1, 0) / C2d(0, 0)});
    SX res = SX::vertcat({ksi, response, Tcomp, T3dcomp});

    for (int i = loc_vars.size() - 1; i >= 0; --i){
        res = SX::simplify(SX::substitute(res, loc_vars.at(i).sym, loc_vars.at(i).val));
    }
    return res;
}

void HEDataGatherer::compute_responses(HEDataGatherer::DefaultGenInputVars &generatedInput,
                                       HEDataGatherer::LocalVars &loc_vars,
                                       HEDataGatherer::Invariants &invariants) {
    response.set_generatedInput(generatedInput);
    response.func_expr = compute_response_expression(loc_vars, invariants);
    response.update_input();
}

HEDataGatherer::RS HEDataGatherer::gatherData(Object3D &obj) {
    RS res;
    gatherData(obj, res);
    return res;
}

HEDataGatherer::RS &HEDataGatherer::gatherData(Object3D &obj, HEDataGatherer::RS &rs) {
    if (&obj != _obj) registerObj(&obj);
    auto& m = obj.m_mesh;
    if (rs.stat.size() < m.num_faces()) rs.stat.resize(m.num_faces());
    DReal data[15] = {0}; DReal* pdat = &(data[0]);
    std::vector<V_ind>& v = _vdat; v.resize(3);
    RS::Elem el{};
    int f_i = 0;
    for (auto f: m.faces()) {
        vert_around(m, f, v);
        response(&pdat, f, v);
        auto& rsc = rs.stat[f_i++];
//        //TODO: instead responses collect tension
//        for (int i = 0; i < 3; ++i) el.ksi[i] = data[i], el.response[i] = data[6+i];
        for (int i = 0; i < 3; ++i) el.ksi[i] = data[i], el.response[i] = data[3+i];
        rsc.emplace_back(el);
    }
    return rs;
}

void HEDataGatherer::gatherData(Object3D &obj, std::function<void(const CompDat& dat)> my_filler) {
    if (&obj != _obj) registerObj(&obj);
    auto& m = obj.m_mesh;
    CompDat cd;
    DReal data[15]; DReal* pdat = &(data[0]);
    cd.ksi = data+0, cd.response = data+3, cd.Tee = data+6, cd.T = data+9;
    std::vector<V_ind>& v = _vdat; v.resize(3);
    RS::Elem el{};
    int f_i = 0;
    for (auto f: m.faces()) {
        cd.f = f;
        vert_around(m, f, v);
        response(&pdat, f, v);
        my_filler(cd);
    }
}

HEDataGatherer::HEDataGatherer(const DefaultHyperElasticForce &f, bool regenerate) :
        name{f.force_name + "_gather"}, gen_dir{f.gen_dir}, U{f.U} {
    auto t = HH::makeDefaultInits(f, getDefaultLocalVars());
    Init(std::move(get<0>(t)), get<1>(t), get<2>(t), regenerate);
}

void HEDataGatherer::Init(HEDataGatherer::DefaultGenInputVars input_vars,
                          HEDataGatherer::LocalVars &loc_vars,
                          HEDataGatherer::Invariants &invariants, bool regenerate) {
    if (gen_dir != "" && gen_dir[gen_dir.size()-1] != '/') gen_dir += "/";
    HH::correlateInputAndVariables(input_vars, loc_vars, invariants);
    compute_responses(input_vars, loc_vars, invariants);
    if (regenerate) generate();
    link_function();
}

HEDataGatherer::HEDataGatherer(std::string name,
                               std::string save_dir,
                               HEDataGatherer::DefaultGenInputVars input_vars,
                               HEDataGatherer::LocalVars &loc_vars,
                               HEDataGatherer::Invariants &invariants,
                               SX U, bool regenerate) : name{name}, gen_dir{save_dir}, U{U}
{
    Init(std::move(input_vars), loc_vars, invariants, regenerate);
}

std::array<double, 3> HEDataGatherer::responseAt(Object3D &obj, F_ind f) {
    if (&obj != _obj) registerObj(&obj);
    auto& m = obj.m_mesh;
    DReal data[6]; DReal* pdat = &(data[0]);
    std::vector<V_ind>& v = _vdat; v.resize(3);
    RS::Elem el{};
    vert_around(m, f, v);
    response(&pdat, f, v);
    for (int i = 0; i < 3; ++i) el.ksi[i] = data[i], el.response[i] = data[3+i];

    return el.response;
}

HEDataGatherer::LocalVars HEDataGatherer::getDefaultLocalVars() {
    static LocalVars res;
    static bool req = false;
    if (!req){
        req = true;
        res = DefaultHyperElasticForce::getDefaultLocalVars();
        HH::LocalVar H("H");
        res.push_back(H);
    }
    return res;
}

std::array<double, 3> _AnalyticChoose::choose(const std::array<double, 3>& ksi, ForceBase* dd, Object3D* obj, F_ind f){
    auto de = static_cast<DataDrivenAnalytic*>(dd);
    return de->dat.responseAt(*obj, f);
}

#include "Forces.h"

StaticDefiniteGather& StaticDefiniteGather::computeTensionedState(){
    if (!m_obj) return *this;
    m_w.clear();

    Object3D& obj = *m_obj;
    AniMesh am;
    am.vertices.resize(3*obj.m_mesh.num_vertices());
    am.faces.resize(3*obj.m_mesh.num_faces());
    am.vertex_label.resize(obj.m_mesh.num_vertices(), 0);
    for (auto v: obj.m_mesh.vertices()){
        for (int k = 0; k < 3; ++k)
            am.vertices[3*v.idx() + k] = obj.m_x[v][k];
    }
    for (auto f: obj.m_mesh.faces()){
        auto v = vert_around(obj.m_mesh, f);
        for (int k = 0; k < 3; ++k)
            am.faces[f.idx()*3 + k] = v[k].idx() + 1;
    }
    for (auto e: obj.m_mesh.edges()){
        auto f = face_around(obj.m_mesh, e);
        if (!f[1].second){
            auto v = vert_around(obj.m_mesh, e);
            am.vertex_label[v[0].idx()] = 1;
            am.vertex_label[v[1].idx()] = 1;
        }
    }
    Object3D obj1(convert_to_Mesh(am, "v:boundary_lbl"));

    auto is_movable = [] (Object3D& obj, const V_ind& v) -> bool { return !(obj.m_boundary[v] & 1); };
    obj1.name = "derived";
    obj1.set_is_movable_checker(is_movable);

    SimplePressureLoad Pr(m_t.m_P);
    auto id = m_w.addObject3D(move(obj1), 1);

    auto Ap_ = set_Ap(obj);
    auto S_ = set_S(obj);
    ConstProperty<F_ind, DReal> h_map = m_w.obj(id).m_mesh.add_property_map<F_ind, DReal>("f:Thickness");
    for (auto f: m_w.obj(id).m_mesh.faces()){ h_map[f] = m_t.m_H /*/ sqrt(S_[f].squared_length()) * Ap_[f]*/; }
    //TODO: можно написать версию силы, которая будет учитывать, что смещения очень маленькие
//    //TODO: wrong force, recomment this line!!!
//    auto aniso_s1 = m_w.obj(id).m_mesh.add_property_map<F_ind, std::array<DReal , 3>>("f:aniso_s").first;
//    for (auto f: m_w.obj(id).m_mesh.faces()){
//        auto v = vert_around(obj.m_mesh, f);
//        auto t = ((m_w.obj(id).m_x0[v[0]] - CGAL::ORIGIN) + (m_w.obj(id).m_x0[v[1]] - CGAL::ORIGIN) + (m_w.obj(id).m_x0[v[2]] - CGAL::ORIGIN)) / 3;
//        double phi = atan2(t[1], t[0]);
//        double teta = atan2(hypot(t[0], t[1]), t[2]);
//        aniso_s1[f] = std::array<DReal , 3>{cos(teta)*cos(phi), cos(teta)*sin(phi), -sin(teta)};
//    }
//    double alpha = 0.0245, beta = 0.0297;
//    Force elastic = MooneyRivlinModel(m_t.m_H, m_t.m_mu, alpha, beta, {NAN, NAN, NAN}, "../generated", true, "f:aniso_s");
//    elastic.target<HyperElasticForceBase>()->prepareJacobianFunction(true);
//    NeoGookModel1 elastic = NeoGookModel1(m_t.m_mu, "f:Thickness", "../generated", false);
//    elastic.prepareJacobianFunction(false);
    Force elastic = NeoGookModelOpt(m_t.m_mu, [h_map](Object3D& , F_ind f)->DReal{ return h_map.first[f]; });
    m_w.addForce(elastic, id);
    m_w.addForce(Pr, id);

    auto stopCondition = [](auto &freq, auto &err, auto &time, auto &maxits) {
        return [&freq, &err, &time, &maxits](StepSimInfo& info, World *w) -> bool {
            static double resid_init;
            int it = info.it;
            if (it == 1) resid_init = w->getWorldNorm("v:force");
            if (it > 0 && it % freq == 0) {
                double resid = w->getWorldNorm("v:force");
                double eps = resid / resid_init;
                std::cout << "it " << it << ": eps = " << eps << " abs = " << resid << " time = " << time.elapsed()
                          << "\n";
                if (eps < err) {
                    std::cout << "Algorithm is converged: \n";
                    return true;
                }
                if (it > maxits || std::isnan(eps)) {
                    std::cout << "Algorithm is diverged or reached maximum of time - iters: \n";
                    return true;
                }
            }
            return false;
        };
    };
    if (m_t.m_method == Traits::SIMPLE_RELAXATON) {
        m_w.setForceApplier(StaticForceApplier(m_t.m_delta), id);
    } else if (m_t.m_method == Traits::SIMPLE_NEWTON){
        m_w.setStepSimulationAlgo([&ssa = m_t.m_nssa](World* w) { return ssa(w); });
    }
    auto time = World3d::Timer();
    m_w.Simulation(stopCondition(m_t.m_freq, m_t.m_err, time, m_t.m_maxits));

    std::array<double, 3> shft = {0};
    for (auto v: m_w.obj(id).m_mesh.vertices()) {
        Vector d = m_w.obj(id).m_x[v] - m_w.obj(id).m_x0[v];
        for (int k = 0; k < 3; ++k)
            if (fabs(d[k]) > shft[k]) shft[k] = fabs(d[k]);
    }
    std::cout << "||X_inv - X||_inf = {" << shft[0] << ", " << shft[1] << ", " << shft[2] << "}" << std::endl;

    return *this;
}

StaticDefiniteGather& StaticDefiniteGather::gatherResponseData(RS& rs){
    if (m_w.objs().empty() || !m_obj) return *this;
    auto& obj1 = m_w.objs().begin()->second.first;
    auto& obj = *m_obj;
    auto D1_ = set_D_vecs(obj1);
//    auto DD1_ = set_DD_matrix(obj1);
    auto Ap1_ = set_Ap(obj1);
    auto S1_ = set_S(obj1);
    auto x1_ = obj1.m_x;
    auto p1_ = obj1.m_x0;
    auto h1 = obj1.m_mesh.property_map<F_ind, DReal>("f:Thickness").first;

    auto L_ = set_L(obj);
//    auto D_ = set_D_vecs(obj);
    auto x_ = obj.m_x;
    auto p_ = obj.m_x0;
    auto Ap_ = set_Ap(obj);
    auto S_ = set_S(obj);
    auto S2d_ = set_S2d(obj, check_flat_initial_template(obj));
    double H = m_t.m_H, mu = m_t.m_mu;

    rs.stat.resize(obj1.m_mesh.num_faces());

    for (auto f: obj1.m_mesh.faces()){
        auto v = vert_around(obj1.m_mesh, f);
        auto Ap1 = Ap1_[f];
        auto Aq1 = sqrt(S1_[f].squared_length());
        Eigen::Matrix<DReal, 3, 3> Q, Q1, P1, D1;
        Q <<    x_[v[0]][0], x_[v[1]][0], x_[v[2]][0],
                x_[v[0]][1], x_[v[1]][1], x_[v[2]][1],
                x_[v[0]][2], x_[v[1]][2], x_[v[2]][2];
        Q1 <<   x1_[v[0]][0], x1_[v[1]][0], x1_[v[2]][0],
                x1_[v[0]][1], x1_[v[1]][1], x1_[v[2]][1],
                x1_[v[0]][2], x1_[v[1]][2], x1_[v[2]][2];
        P1 <<   p1_[v[0]][0], p1_[v[1]][0], p1_[v[2]][0],
                p1_[v[0]][1], p1_[v[1]][1], p1_[v[2]][1],
                p1_[v[0]][2], p1_[v[1]][2], p1_[v[2]][2];
        D1 <<   D1_[f][0][0], D1_[f][1][0], D1_[f][2][0],
                D1_[f][0][1], D1_[f][1][1], D1_[f][2][1],
                D1_[f][0][2], D1_[f][1][2], D1_[f][2][2];
        auto dF1 = (Q1 - P1) * D1.transpose();
        DReal J1 = Aq1 / Ap1;

//                auto F1 = Q1 * D1.transpose();
//                auto T = h1[f]*mu / J1 * (F1 * F1.transpose() - Eigen::Matrix<DReal, 3, 3>::Identity() / (J1 * J1));
//        //TODO: здесь несогласованность!!! Не уничтожается лишняя компонента!
//        auto I = Eigen::Matrix<DReal, 3, 3>::Identity();
//        Eigen::Vector3d n1(obj1.m_normal[f][0], obj1.m_normal[f][1], obj1.m_normal[f][2]);
//        Eigen::Matrix<DReal, 3, 3> Qnm1, DD1;
//        Qnm1.col(0) = n1.cross(Q1.col(2) - Q1.col(1));
//        Qnm1.col(1) = n1.cross(Q1.col(0) - Q1.col(2));
//        Qnm1.col(2) = - (Qnm1.col(0) + Qnm1.col(1));
//        for (auto i = 0; i < 3; ++i)
//            for (auto j = i; j < 3; ++j)
//                DD1(j, i) = DD1(i, j) = DD1_[f][j + (i + 1) * (i > 0)];
//        auto T = h1[f] * mu / J1 * ((dF1 + I) * (2*Q1*DD1.transpose() - Qnm1 / (J1 * J1 * J1 * Ap1)) * (dF1 + I).transpose());
//        auto T = h1[f] * mu / J1 * (dF1 * dF1.transpose() + dF1 + dF1.transpose() + Eigen::Matrix<DReal, 3, 3>::Identity() * (J1 - 1) * (J1 + 1) / (J1 * J1));
//        Eigen::Matrix<DReal, 3, 2> G, g;
//        G.col(0) = P1.col(0) - P1.col(1); G.col(0) /= G.col(0).norm();
//        G.col(1) = P1.col(0) - P1.col(2) - (P1.col(0) - P1.col(2)).dot(G.col(0)) * G.col(0); G.col(1) /= G.col(1).norm();
//
//        g.col(0) = Q1.col(0) - Q1.col(1); g.col(0) /= g.col(0).norm();
//        g.col(1) = Q1.col(0) - Q1.col(2) - (Q1.col(0) - Q1.col(2)).dot(g.col(0)) * g.col(0); g.col(1) /= g.col(1).norm();
//        auto GG = G.transpose() * G, gg = g.transpose() * g;
//        auto T = h1[f] * mu / J1 * g * (GG - gg / (J1 * J1)) * g.transpose();
        auto L = Eigen::Map<Eigen::Matrix<DReal, 2, 3>>(L_[f].data());
        auto G = P1 * L.transpose(),  g = Q1 * L.transpose();
        auto GG = (G.transpose() * G).inverse(); auto gg = (g.transpose() * g).inverse();
        auto T = h1[f] * mu / J1 * g * (GG - gg / (J1 * J1)) * g.transpose();

        auto F = Q1 * L.transpose();
        Eigen::Matrix<DReal, 3, 2> R;
        R.col(0) = F.col(0) / F.col(0).norm();
        R.col(1) = F.col(1) - F.col(1).dot(R.col(0))*R.col(0); R.col(1) /= R.col(1).norm();
        auto Tee = R.transpose() * T * R;
        auto Fr = Q * L.transpose();
        auto C2d = Fr.transpose() * Fr;
//        auto F11 = sqrt(C2d(0, 0));
//        auto F12 = C2d(0, 1) / F11;
//        auto F22 = sqrt(C2d(1, 1) - C2d(0, 1) * C2d(0, 1) / C2d(0, 0));
        DReal J = sqrt(S_[f].squared_length()) / Ap_[f];
        auto sq = [](auto t) { return t * t; };
        RS::Elem e{0};

        e.response = {J * Tee(0, 0)            ,
                      J * Tee(1, 1)            ,
                      J * Tee(1, 0) * C2d(0, 0) / sqrt(C2d.determinant())};
//        e.response = {Tee(0, 0)            ,
//                      Tee(1, 1)            ,
//                      Tee(1, 0)};
        e.ksi = {   0.5 * log(C2d(0, 0)),
                    0.5 * log(C2d(1, 1) - sq(C2d(1, 0))/C2d(0, 0)),
                    C2d(1, 0) / C2d(0, 0)   };

        if (fabs(C2d(0, 0) - 1) < 1e-7 || fabs(C2d(1, 1) - 1) < 1e-7){
            auto S2d = Eigen::Map<Eigen::Matrix<DReal, 2, 3>>(S2d_[f].data());
            Eigen::Matrix<DReal, 3, 3> P;
            P <<    p_[v[0]][0], p_[v[1]][0], p_[v[2]][0],
                    p_[v[0]][1], p_[v[1]][1], p_[v[2]][1],
                    p_[v[0]][2], p_[v[1]][2], p_[v[2]][2];
            auto dFr = (Q - P) * L.transpose();
            auto SdF = S2d * dFr;
            auto dC2d = dFr.transpose() * dFr + SdF + SdF.transpose();
            if (fabs(C2d(0, 0) - 1) < 1e-7)
                e.ksi[0] = 0.5 * log1p(dC2d(0, 0));
            if (fabs(C2d(1, 1) - 1) < 1e-7)
                e.ksi[1] = 0.5 * (log1p(dC2d(1, 1)) + log1p(-sq(C2d(1, 0)) / (C2d(0, 0) * C2d(1, 1))));
        }
        rs.stat[f.idx()].push_back(e);
    }

    return *this;
}

StaticDefiniteGather &StaticDefiniteGather::gatherMyData(function<void(const ObjsData &)> my_filler) {
    if (m_w.objs().empty() || !m_obj) return *this;
    auto& obj1 = m_w.objs().begin()->second.first;
    auto& obj = *m_obj;
    auto D1_ = set_D_vecs(obj1);
//    auto DD1_ = set_DD_matrix(obj1);
    auto Ap1_ = set_Ap(obj1);
    auto S1_ = set_S(obj1);
    auto x1_ = obj1.m_x;
    auto p1_ = obj1.m_x0;
    auto h1 = obj1.m_mesh.property_map<F_ind, DReal>("f:Thickness").first;
    auto S2d1_ = set_S2d(obj1, false);

    auto L_ = set_L(obj);
    auto D_ = set_D_vecs(obj);
    auto x_ = obj.m_x;
    auto p_ = obj.m_x0;
    auto Ap_ = set_Ap(obj);
    auto S_ = set_S(obj);
    auto S2d_ = set_S2d(obj, check_flat_initial_template(obj));
    double mu = m_t.m_mu;

    ObjsData odat;
    odat.nfaces = obj1.m_mesh.num_faces();

    for (auto f: obj1.m_mesh.faces()){
        auto v = vert_around(obj1.m_mesh, f);
        odat.f = f;
        odat.Ap = Ap_[f];
        odat.Aq = sqrt(S_[f].squared_length());
        odat.Ap1 = Ap1_[f];
        odat.Aq1 = sqrt(S1_[f].squared_length());
        odat.Q << x_[v[0]][0], x_[v[1]][0], x_[v[2]][0],
                x_[v[0]][1], x_[v[1]][1], x_[v[2]][1],
                x_[v[0]][2], x_[v[1]][2], x_[v[2]][2];
        odat.P <<   p_[v[0]][0], p_[v[1]][0], p_[v[2]][0],
                p_[v[0]][1], p_[v[1]][1], p_[v[2]][1],
                p_[v[0]][2], p_[v[1]][2], p_[v[2]][2];
        odat.D <<   D_[f][0][0], D_[f][1][0], D_[f][2][0],
                D_[f][0][1], D_[f][1][1], D_[f][2][1],
                D_[f][0][2], D_[f][1][2], D_[f][2][2];
        odat.Q1 <<   x1_[v[0]][0], x1_[v[1]][0], x1_[v[2]][0],
                x1_[v[0]][1], x1_[v[1]][1], x1_[v[2]][1],
                x1_[v[0]][2], x1_[v[1]][2], x1_[v[2]][2];
        odat.P1 <<   p1_[v[0]][0], p1_[v[1]][0], p1_[v[2]][0],
                p1_[v[0]][1], p1_[v[1]][1], p1_[v[2]][1],
                p1_[v[0]][2], p1_[v[1]][2], p1_[v[2]][2];
        odat.D1 <<   D1_[f][0][0], D1_[f][1][0], D1_[f][2][0],
                D1_[f][0][1], D1_[f][1][1], D1_[f][2][1],
                D1_[f][0][2], D1_[f][1][2], D1_[f][2][2];
        odat.L = Eigen::Map<Eigen::Matrix<DReal, 2, 3>>(L_[f].data());
        odat.F2d = odat.Q1 * odat.L.transpose();
        odat.C2d = odat.F2d.transpose() * odat.F2d;

        Eigen::Matrix<DReal, 3, 3> U1 = (odat.Q1 - odat.P1);
        auto dF1 = (odat.Q1 - odat.P1) * odat.D1.transpose();
        DReal J1 = odat.Aq1 / odat.Ap1;
        DReal sqJ2 = S1_[f].squared_length() / (odat.Ap1 * odat.Ap1);

//        auto G = odat.P1 * odat.L.transpose(),  g = odat.Q1 * odat.L.transpose();
//        auto GG = (G.transpose() * G).inverse(); auto gg = (g.transpose() * g).inverse();
//        odat.T = h1[f] * mu / J1 * g * (GG - gg / (J1 * J1)) * g.transpose();
//        {
//            Eigen::Matrix<DReal, 3, 2> Fs = odat.P1 * odat.L.transpose(), Fr = odat.Q1 * odat.L.transpose();
//            Eigen::Matrix<DReal, 3, 3> FF = odat.P1 * odat.D1.transpose();
//            Eigen::Matrix<DReal, 3, 2> dF = U1 * odat.D1.transpose() * Fs,
//                                        F = FF * Fs;
//            Eigen::Matrix<DReal, 2, 2> C2d = F.transpose() * F,
//                dC2d = dF.transpose() * F + F.transpose() * dF + dF.transpose() * dF;
//            double J_ = C2d.determinant(), J1_ = (C2d + dC2d).determinant();
//            Eigen::Matrix<DReal, 2, 2> gss = (Fs.transpose() * Fs);
//            Eigen::Matrix<DReal, 2, 2> S = gss - C2d.adjoint() / (sqJ2 * sqJ2);
//            Eigen::Matrix<DReal, 2, 2> gs = (Fs.transpose() * Fs).inverse();
//            odat.T = h1[f] * mu / J1 * F * gs * S * gs * F.transpose();
//        }
        if (true){
            auto I = Eigen::Matrix<DReal, 2, 2>::Identity();
            auto S2d = Eigen::Map<Eigen::Matrix<DReal, 2, 3>>(S2d1_[f].data());
            Eigen::Matrix<DReal, 3, 3> FF = odat.P1 * odat.D1.transpose();
            Eigen::Matrix<DReal, 3, 2> dF = U1 * odat.D1.transpose() * S2d.transpose(),
                                        F = FF * S2d.transpose();
            Eigen::Matrix<DReal, 2, 2> C2d = F.transpose() * F,
                            dC2d = dF.transpose() * F + F.transpose() * dF + dF.transpose() * dF;
            double dsqJ = dC2d.trace() + dC2d.determinant();
            double J_ = C2d.determinant();
//            Eigen::Matrix<DReal, 2, 2> S = I - C2d.adjoint() / (sqJ2 * sqJ2);
            Eigen::Matrix<DReal, 2, 2> S = ((2*dsqJ + dsqJ*dsqJ)*I - dC2d.adjoint()) / (1 + dsqJ) / (1 + dsqJ);
            odat.T = h1[f] * mu / J1 * (F + dF) * S * (F + dF).transpose();
        }
//        if (f.idx() == 388){
//            std::cout << "h =\n" << h1[f] << std::endl;
//            std::cout << "L =\n" << odat.L << std::endl;
//            std::cout << "P1 =\n" << odat.P1 << std::endl;
//            std::cout << "Q1 =\n" << odat.Q1 << std::endl;
//            std::cout << "G =\n" << G <<std::endl;
//            std::cout << "g =\n" << g <<std::endl;
//            std::cout << "GG =\n" << GG <<std::endl;
//            std::cout << "gg =\n" << gg <<std::endl;
//            std::cout << "T =\n" << odat.T << std::endl;
//        }

//        auto I = Eigen::Matrix<DReal, 3, 3>::Identity();
//        Eigen::Vector3d n1(obj1.m_normal[f][0], obj1.m_normal[f][1], obj1.m_normal[f][2]);
//        Eigen::Matrix<DReal, 3, 3> Qnm1, DD1;
//        Qnm1.col(0) = n1.cross(odat.Q1.col(2) - odat.Q1.col(1));
//        Qnm1.col(1) = n1.cross(odat.Q1.col(0) - odat.Q1.col(2));
//        Qnm1.col(2) = - (Qnm1.col(0) + Qnm1.col(1));
//        for (auto i = 0; i < 3; ++i)
//            for (auto j = i; j < 3; ++j)
//                DD1(j, i) = DD1(i, j) = DD1_[f][j + (i + 1) * (i > 0)];
//        odat.T = h1[f] * mu / J1 * ((dF1 + I) * (2*odat.Q1*DD1.transpose() - Qnm1 / (J1 * J1 * J1 * odat.Ap1)) * (dF1 + I).transpose());
//        odat.T = h1[f] * mu / J1 * (dF1 * dF1.transpose() + dF1 + dF1.transpose() + Eigen::Matrix<DReal, 3, 3>::Identity() * (J1 - 1) * (J1 + 1) / (J1 * J1));

//        Eigen::Matrix<DReal, 3, 2> G, g;
//        auto& P1 = odat.P1, Q1 = odat.Q1;
//        G.col(0) = P1.col(0) - P1.col(1); G.col(0) /= G.col(0).norm();
//        auto G1 = P1.col(0) - P1.col(2) - (P1.col(0) - P1.col(2)).dot(G.col(0)) * G.col(0);
//        G.col(1) = G1 / G1.norm();
//
//        g.col(0) = Q1.col(0) - Q1.col(1); g.col(0) /= g.col(0).norm();
//        g.col(1) = Q1.col(0) - Q1.col(2) - (Q1.col(0) - Q1.col(2)).dot(g.col(0)) * g.col(0); g.col(1) /= g.col(1).norm();
//        Eigen::Matrix<DReal, 3, 2> dg;
//        auto dU01 = U1.col(0) - U1.col(1); auto dP01 = P1.col(0) - P1.col(1);
//        auto dU02 = U1.col(0) - U1.col(2); auto dP02 = P1.col(0) - P1.col(2);
//        dg.col(0) = (dU01 * dP01.norm() - dU01.dot(dP01) * G.col(0)) / dP01.squaredNorm();
//        auto ddg2 = -dU02.dot(G.col(0))*G.col(0)-dP02.dot(dg.col(0))*g.col(0) - dP02.dot(g.col(0))*dg.col(0) + dU02;
//        dg.col(1) = (ddg2 * G1.norm() - ddg2.dot(G1) * G.col(1)) / G1.squaredNorm();
//        dg = dg.eval();
//        auto dS = (dU01.cross(dP02) + dP01.cross(dU02) + dU01.cross(dU02)).eval();
//        auto Sp = (dP01.cross(dP02)).eval();
//        auto dJ = dS.dot(Sp) / Sp.squaredNorm();
//        auto dgg = dg.transpose() * G + G.transpose() * dg + dg.transpose() * dg;
//        auto GG = G.transpose() * G, gg = g.transpose() * g;
//        odat.T = h1[f] * mu / J1 * g * (GG - gg / (J1 * J1)) * g.transpose();
//        odat.T = h1[f] * mu / J1 * g * (2*dJ*GG - dgg) * g.transpose();
//        std::cout << odat.T << std::endl;


        auto& F = odat.F2d;

//        Eigen::Matrix<DReal, 3, 2> R;
//        R.col(0) = F.col(0) / F.col(0).norm();
//        R.col(1) = F.col(1) - F.col(1).dot(R.col(0))*R.col(0); R.col(1) /= R.col(1).norm();
//        odat.Tee = R.transpose() * odat.T * R;
        Eigen::Matrix<DReal, 3, 3> R;
        R.col(0) = F.col(0) / F.col(0).norm();
        R.col(1) = F.col(1) - F.col(1).dot(R.col(0))*R.col(0); R.col(1) /= R.col(1).norm();
        R.col(2) = R.col(0).cross(R.col(1));
        auto Tee3d = R.transpose() * odat.T * R;//TODO: R * odat.T * R.transpose();
        odat.Tee << Tee3d(0, 0), Tee3d(0, 1),
                    Tee3d(1, 0), Tee3d(1, 1);

        my_filler(odat);
    }

    return *this;
}
