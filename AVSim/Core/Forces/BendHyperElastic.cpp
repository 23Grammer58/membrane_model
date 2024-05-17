//
// Created by alex on 29.01.2021.
//

#include "BendHyperElastic.h"
using namespace World3d;

void BendingForce::_BoundaryData::set_boundary(Object3D &obj, const std::map<E_ind, Entry> &bdata) {
    _obj = &obj;
    std::map<E_ind, InternalEntry> data;
    bool is_flat = check_flat_initial_template(obj);
    LocalCartesian S(obj, is_flat);
    for (auto& it: bdata){
        auto v = vert_around(obj.m_mesh, it.first);
        Vector a0(obj.m_x0[v[0]], obj.m_x0[v[1]]);
        DReal l0 = sqrt(a0.squared_length());
        auto f = face_around(obj.m_mesh, it.first)[0].first;
        auto vv = vert_around(obj.m_mesh, f);
        V_ind o = vv[0];
        for (int i = 0; i < 3; ++i)
            if (vv[i] != v[0] && vv[i] != v[1]){
                o = vv[i];
                i = 3;
            }
        Vector a(obj.m_x0[v[0]], obj.m_x0[v[1]]), b(obj.m_x0[v[0]], obj.m_x0[o]);
        Vector normal = CGAL::cross_product(a, CGAL::cross_product(a, b));
        if (b * normal > 0)
            normal *= -1;

        Eigen::Matrix<DReal, 2, 1> nn = S(f)*Eigen::Matrix<DReal, 3, 1>(normal[0], normal[1], normal[2]);
        nn.normalize();
        nn.eval();
        data.insert({it.first, {it.second, {nn[0], nn[1]}, l0}});
    }
    boundary = [bdata = data](E_ind e) -> InternalEntry{
        static InternalEntry zero = {{0, {0, 0, 0}}, {0, 0}, 0};
        auto it = bdata.find(e);
        if (it != bdata.end())
            return it->second;
        else return zero;
    };
    prepared = true;
}

void BendingForce::_BoundaryData::generate_data(Object3D &obj) {
    std::cout << "Boundary data for " << obj.name << " is not presented! Will be used free boundary condition" << std::endl;
    auto conv = [&obj](V_ind v){
        return FREE;
    };
    generate_data(obj, conv);
}

int BendingForce::_BoundaryData::check_boundary(Object3D &obj) {
    for (auto e: obj.m_mesh.edges()) {
        auto v = vert_around(obj.m_mesh, e);
        auto ff = face_around(obj.m_mesh, e);
        if (ff[0].second && ff[1].second)
            continue;
        auto it = boundary(e);
        if (it.entry.btype == 0)
            return 1;
        if (it.entry.btype == CLAMPED && fabs(it.entry.normal[0]) + fabs(it.entry.normal[1]) + fabs(it.entry.normal[2]) < 1e-6)
            return 2;
    }
    return 0;
}

void BendingForce::_BoundaryData::repair_boundary(Object3D &obj) {
    std::map<E_ind, Entry> bdata;
    for (auto e: obj.m_mesh.edges()){
        auto v = vert_around(obj.m_mesh, e);
        auto ff = face_around(obj.m_mesh, e);
        if (ff[0].second && ff[1].second)
            continue;
        auto it = boundary(e);
        if (it.entry.btype != 0) {
            if (it.entry.btype == CLAMPED && fabs(it.entry.normal[0]) + fabs(it.entry.normal[1]) + fabs(it.entry.normal[2]) < 1e-6){
                auto h = CGAL::halfedge(e, obj.m_mesh);
                auto vv = CGAL::target(CGAL::next(h, obj.m_mesh), obj.m_mesh);
                Vector nn = CGAL::cross_product(obj.m_x[vv] - obj.m_x[v[0]], obj.m_x[vv] - obj.m_x[v[1]]);
                nn /= sqrt(nn.squared_length());
                nn = CGAL::cross_product(obj.m_x[v[1]] - obj.m_x[v[0]], nn);
                nn /= sqrt(nn.squared_length());
                it.entry = Entry(it.entry.btype, nn[0], nn[1], nn[2]);
            }
            bdata.insert({e, it.entry});
        }
        else
            bdata.insert({e, {static_cast<DReal>(FREE), {0, 0, 0}}});
    }
    set_boundary(obj, bdata);
}

void BendingForce::_BoundaryData::repair_if_necessary(Object3D &obj) {
    if (int check = check_boundary(obj); check) {
        if (check == 1)
            std::cout << "Bending boundary was set not completely, the remaining part of boundary will be set as free boundary" << std::endl;
        if (check == 2)
            std::cout << "Clamped part of bending boundary has broken normals, they will be set automatically" << std::endl;
        repair_boundary(obj);
    }
}

void BendingForce::_BoundaryData::generate_data(Object3D &obj, const std::function<int(V_ind)> &conv) {
    std::map<E_ind, Entry> bdata;
    for (auto e: obj.m_mesh.edges()){
        auto v = vert_around(obj.m_mesh, e);
        auto ff = face_around(obj.m_mesh, e);
        if (ff[0].second && ff[1].second)
            continue;
        auto f = ff[0].first;

        int types[2] = {conv(v[0]), conv(v[1])};
        int type = FREE;
        if (types[0] == CLAMPED && types[1] == CLAMPED)
            type = CLAMPED;
        if (type == CLAMPED){
            auto vv = vert_around(obj.m_mesh, f);
            V_ind o = vv[0];
            for (int i = 0; i < 3; ++i)
                if (vv[i] != v[0] && vv[i] != v[1]){
                    o = vv[i];
                    i = 3;
                }
            Vector a(obj.m_x[v[0]], obj.m_x[v[1]]), b(obj.m_x[v[0]], obj.m_x[o]);
            Vector normal = CGAL::cross_product(a, CGAL::cross_product(a, b));
            if (b * normal > 0) normal *= -1;
            normal /= sqrt(normal.squared_length());
            bdata.insert({e, {static_cast<DReal>(type), {normal[0], normal[1], normal[2]}}});
        } else bdata.insert({e, {static_cast<DReal>(type), {0, 0, 0}}});
    }
    set_boundary(obj, bdata);
}

const DReal *BendingForce::_BoundaryData::get_data(F_ind f, std::vector<V_ind> &v) {
    auto ee = edge_around(_obj->m_mesh, f);
    std::array<int, 3> remap = {1, 2, 0};
    for (int i = 0; i < 3; ++i) {
        data[remap[i]] = boundary(ee[i]);
    }
    return reinterpret_cast<const DReal*>(&data[0]);
}

void BendingForce::BoundaryRequier::registerObj(Object3D *obj) {
    if (!p->prepared || p->_obj != obj)
        p->generate_data(*obj);
    else p->repair_if_necessary(*obj);
}
#ifndef USE_BST
BendingForce::VarContainer<HH::LocalVar> BendingForce::getDefaultLocalVars() {
    static VarContainer<LocalVar> res;
    static bool req = false;
    if (!req){
        res = DefaultHyperElasticForce::getDefaultLocalVars();
        LocalVar H("H");
        LocalVar Pnt[6];
        for (int i = 3; i < 6; ++i){
            Pnt[i-3] = LocalVar("Pnt" + std::to_string(i), 3);
            Pnt[i] = LocalVar("Pnt" + std::to_string(i) + "_init", 3);
        }
        LocalVar N_in("N_in", 4*2*3);
        LocalVar JBnd_in("JBend", 2, 2*3);
        LocalVar k0("k0", 3);
        LocalVar Bnd_in("Bnd_in", 3*(1+3+2+1));
        res.push_back({H, Pnt[0], Pnt[1], Pnt[2], Pnt[3], Pnt[4], Pnt[5], N_in, JBnd_in, k0, Bnd_in});

        std::array<LocalVar, 3> N, b_lbl, b_phi0, b_n, l0;
        for (int k = 0; k < 3; ++k) {
            SX _N(6, 2);
            const bool use_N_in = false;
            if (use_N_in) {
                for (int j = 0; j < 2; ++j)
                    for (int i = 0; i < 6; ++i)
                        if (i < 3)
                            _N(i, j) = N_in.sym(i + 4 * j + 2 * 4 * k, 0);
                        else if (i == 3 + k)
                            _N(i, j) = N_in.sym(3 + 4 * j + 2 * 4 * k, 0);
                        else _N(i, j) = 0.0;
            } else {
                static const char _N_doubled[3][2][6] = {
                        {{-1, 1, -1, 1, 0,  0}, {-1, -1, 1, 1, 0, 0}},
                        {{-1, 1, 1,  0, -1, 0}, {-2, 0,  2, 0, 0, 0}},
                        {{-2, 2, 0,  0, 0,  0}, {-1, 1,  1, 0, 0, -1}}
                };
                SX NNd(2, 6);
                for (int i = 0; i < 2; ++i)
                    for (int j = 0; j < 6; ++j)
                        NNd(i, j) = _N_doubled[k][i][j];
                auto JJ = SX::horzsplit(JBnd_in.sym, 2);

                _N = SX::simplify(SX::mtimes(JJ[k],  NNd).T());
            }
            N[k] = LocalVar("N" + std::to_string(k), _N);

            int i = k;
            b_lbl[i] = LocalVar("b_lbl" + std::to_string(i), Bnd_in.sym(7*i+0, 0));
            b_phi0[i] = LocalVar("b_phi0" + std::to_string(i), SX::vertcat({Bnd_in.sym(7*i+1, 0), Bnd_in.sym(7*i+2, 0), Bnd_in.sym(7*i+3, 0)}));
            b_n[i] = LocalVar("b_n" + std::to_string(i), SX::vertcat({Bnd_in.sym(7*i+4, 0), Bnd_in.sym(7*i+5, 0)}));
            l0[i] = LocalVar("b_l0" + std::to_string(i), Bnd_in.sym(7*i+6, 0));
            res.push_back({N[k], b_lbl[k], b_phi0[k], b_n[k], l0[k]});
        }
        {
            auto Bb_func = [](bool is_bnd) {
                SX _Ap = res["Ap"], _Aq = res["Aq"], Q = res["Q"], L = res["L"];
                std::array<SX, 3> _N = {res["N0"], res["N1"], res["N2"]};
                std::array<SX, 6> _P;
                for (int i = 0; i < 6; ++i)
                    _P[i] = res["Pnt"+std::to_string(i)];
                std::array<SX, 3> b_lbl, b_phi0, b_n, l0;
                for (int i = 0; i < 3; ++i){
                    b_lbl[i] = res["b_lbl" + std::to_string(i)];
                    b_phi0[i] = res["b_phi0" + std::to_string(i)];
                    b_n[i] = res["b_n" + std::to_string(i)];
                    l0[i] = res["b_l0" + std::to_string(i)];
                }
                SX t = res["n"];
                SX lambda = _Ap / _Aq;
                SX phi_d_M = SX::mtimes(L, Q.T()).T();
                std::array<SX, 3> N_M_I;
                std::array<SX, 3> N_I_I;
                std::array<SX, 3> phi_d_I;
                std::array<SX, 3> loc_n_m;
                for (int i = 0; i < 3; ++i) {
                    auto tmp = SX::horzsplit(_N[i].T(), 1);
                    N_M_I[i] = SX::horzcat({tmp[0], tmp[1], tmp[2]});
                    N_I_I[i] = tmp[3 + i];

                    phi_d_I[i] = (SX::mtimes(N_M_I[i], Q.T()) + SX::mtimes(N_I_I[i], _P[i + 3].T())).T();
                    if (is_bnd) {
                        //здесь условие на не деление на ноль т.к. иначе при автогенерации производных это деление может почему-то выполняться приводя к nan!!!
                        SX phi_d_s = (_P[(i + 2) % 3] - _P[(i + 1) % 3]) / SX::if_else(b_lbl[i] == 2, l0[i], 1);
                        SX phi_d_n = b_phi0[i] / (lambda * SX::norm_2(phi_d_s));
                        loc_n_m[i] = SX::vertcat({SX::horzcat({b_n[i](0, 0), -b_n[i](1, 0)}),
                                                  SX::horzcat({b_n[i](1, 0), b_n[i](0, 0)})});
                        SX phi_d_B = SX::mtimes(loc_n_m[i], SX::vertcat({phi_d_n.T(), phi_d_s.T()})).T();
                        phi_d_I[i] = SX_switch(0.01,b_lbl[i], phi_d_I[i], SX_Case{1, SX::zeros(3, 2)/*phi_d_M*/}, SX_Case{2, phi_d_B});
                    }
                }
                auto L_I = [L](int i){
                    SX LI(3, 2);
                    LI(0, 0) = L(0, i);
                    LI(1, 1) = L(1, i);
                    LI(2, 0) = L(1, i);
                    LI(2, 1) = L(0, i);
                    return LI;
                };

                SX H = 2 * (SX::mtimes(L_I(0), phi_d_I[0].T()) + SX::mtimes(L_I(1), phi_d_I[1].T()) + SX::mtimes(L_I(2), phi_d_I[2].T()));
                SX k = SX::mtimes(H, t);
                k(2, 0) /= 2;

                SX phi_contr = lambda * SX::horzcat(
                        {SX::cross(SX::horzsplit(phi_d_M, 1)[1], t), -SX::cross(SX::horzsplit(phi_d_M, 1)[0], t)});
                SX Ro = SX::mtimes({H, phi_contr, L});

                SX Bb;
                if (true) {
                    SX dA = SX::sym("dA", 1, 6);
                    auto du = SX::horzsplit(dA, 1);
                    SX dQ = SX::horzcat({du[0], du[1], du[2]});
                    std::array<SX, 3> dH_I;
                    SX dk(3, 1);
                    if (!is_bnd){
                        for (int i = 0; i < 3; ++i) {
                            dH_I[i] = 2 * SX::mtimes(L_I(i),
                                                     SX::mtimes(N_M_I[i], dQ.T()) +
                                                     SX::mtimes(N_I_I[i], du[i + 3].T()));

                            dk += dH_I[i];
                        }
                        dk -= SX::mtimes({Ro, dQ.T()});
                    } else {
                        for (int i = 0; i < 3; ++i) {
                            dH_I[i] = 2 * SX::mtimes(L_I(i),
                                                     SX::mtimes(N_M_I[i], dQ.T()) +
                                                     SX::mtimes(N_I_I[i], du[i + 3].T()));
                            SX dH_bnd1 = SX::zeros(3, 1);//2 * SX::mtimes({L_I(i), L, dQ.T()});
                            SX dH_bnd2 = 2 / l0[i] * SX::mtimes({L_I(i), loc_n_m[i], SX::vertcat(
                                    {0, (du[(i + 2) % 3] - du[(i + 1) % 3]).T()})});
                            dH_I[i] = SX_switch(b_lbl[i], dH_I[i], SX_Case{1, dH_bnd1}, SX_Case{2, dH_bnd2});

                            dk += dH_I[i];
                        }
                        dk -= SX::mtimes({Ro, dQ.T()});
                    }
                    Bb = SX::jacobian(dk, dA);
                }
                return std::pair<SX, SX>{k, Bb};
            };
            auto _int = Bb_func(false);
            auto _bnd = Bb_func(true);
            LocalVar k_int("k_int", _int.first);
            LocalVar Bb_int("BbMatrix_int", _int.second);
            LocalVar k_bnd("k_bnd", _bnd.first);
            LocalVar Bb_bnd("BbMatrix_bnd", _bnd.second);
            res.push_back({Bb_int, k_int, Bb_bnd, k_bnd});
        }

        req = true;
    }
    return res;
}
#else
BendingForce::VarContainer<HH::LocalVar> BendingForce::getDefaultLocalVarsBST() {
    static VarContainer<LocalVar> res;
    static bool req = false;
    if (!req){
        res = DefaultHyperElasticForce::getDefaultLocalVars();
        LocalVar H("H");
        LocalVar Pnt[6];
        for (int i = 3; i < 6; ++i){
            Pnt[i-3] = LocalVar("Pnt" + std::to_string(i), 3);
            Pnt[i] = LocalVar("Pnt" + std::to_string(i) + "_init", 3);
        }
        LocalVar N_in("N_in", 4*2*3);
        LocalVar JBnd_in("JBend", 2, 2*3);
        LocalVar k0("k0", 3);
        LocalVar Bnd_in("Bnd_in", 3*(1+3+2+1));
        res.push_back({H, Pnt[0], Pnt[1], Pnt[2], Pnt[3], Pnt[4], Pnt[5], N_in, JBnd_in, k0, Bnd_in});

        std::array<LocalVar, 3> N, b_lbl, b_phi0, b_n, l0;
        for (int k = 0; k < 3; ++k) {
            SX _N(6, 2);
            const bool use_N_in = false;
            if (use_N_in) {
                for (int j = 0; j < 2; ++j)
                    for (int i = 0; i < 6; ++i)
                        if (i < 3)
                            _N(i, j) = N_in.sym(i + 4 * j + 2 * 4 * k, 0);
                        else if (i == 3 + k)
                            _N(i, j) = N_in.sym(3 + 4 * j + 2 * 4 * k, 0);
                        else _N(i, j) = 0.0;
            } else {
                static const char _N_doubled[3][2][6] = {
                        {{-1, 1, -1, 1, 0,  0}, {-1, -1, 1, 1, 0, 0}},
                        {{-1, 1, 1,  0, -1, 0}, {-2, 0,  2, 0, 0, 0}},
                        {{-2, 2, 0,  0, 0,  0}, {-1, 1,  1, 0, 0, -1}}
                };
                SX NNd(2, 6);
                for (int i = 0; i < 2; ++i)
                    for (int j = 0; j < 6; ++j)
                        NNd(i, j) = _N_doubled[k][i][j];
//                        std::cout << "JBnd_in = " << JBnd_in.sym <<  std::endl;
                auto JJ = SX::horzsplit(JBnd_in.sym, 2);
//                        std::cout << "JJ = \n" << JJ[0] << "\n" << JJ[1] << "\n" << JJ[2] <<  std::endl;

                _N = SX::simplify(SX::mtimes(JJ[k],  NNd).T());
            }
            N[k] = LocalVar("N" + std::to_string(k), _N);

            int i = k;
            b_lbl[i] = LocalVar("b_lbl" + std::to_string(i), Bnd_in.sym(7*i+0, 0));
            b_phi0[i] = LocalVar("b_phi0" + std::to_string(i), SX::vertcat({Bnd_in.sym(7*i+1, 0), Bnd_in.sym(7*i+2, 0), Bnd_in.sym(7*i+3, 0)}));
            b_n[i] = LocalVar("b_n" + std::to_string(i), SX::vertcat({Bnd_in.sym(7*i+4, 0), Bnd_in.sym(7*i+5, 0)}));
            l0[i] = LocalVar("b_l0" + std::to_string(i), Bnd_in.sym(7*i+6, 0));
            res.push_back({N[k], b_lbl[k], b_phi0[k], b_n[k], l0[k]});
        }
        {
            auto Bb_func = [](bool is_bnd) {
                SX _Ap = res["Ap"], _Aq = res["Aq"], Q = res["Q"];
                std::array<SX, 3> _N = {res["N0"], res["N1"], res["N2"]};
                std::array<SX, 6> _P, PP;
                for (int i = 0; i < 6; ++i) {
                    _P[i] = res["Pnt" + std::to_string(i)];
                    PP[i] = res["Pnt" + std::to_string(i) + "_init"];
                }
                SX LL_0;
                for (int i = 0; i < 3; ++i){
                    std::array<int, 3> q = {0, 1, 2};
                    SX Ap2 = SX::norm_2(SX::cross(PP[q[1]] - PP[q[0]], PP[q[2]] - PP[q[0]]));
                    SX n = SX::cross(PP[q[1]] - PP[q[0]], PP[q[2]] - PP[q[0]])/Ap2;
                    SX D[3];
                    for (int j = 0; j < 3; ++j)
                        D[j] = 1/Ap2 * SX::cross(n, PP[q[(j+2)%3]] - PP[q[(j+1)%3]]);
                    LL_0 = SX::horzcat({D[0], D[1], D[2]});
                    LL_0 = SX::vertcat({SX::vertsplit(LL_0, 1)[0], -SX::vertsplit(LL_0, 1)[1]});
                }
                SX L = res["L"];
                std::array<SX, 3> LL_I;
                for (int i = 0; i < 3; ++i){
                    std::array<int, 3> q = {(i+1)%3, i + 3, (i + 2)%3};
                    SX Ap2 = SX::norm_2(SX::cross(PP[q[1]] - PP[q[0]], PP[q[2]] - PP[q[0]]));
                    SX n = SX::cross(PP[q[1]] - PP[q[0]], PP[q[2]] - PP[q[0]])/Ap2;;
                    SX D[3];
                    for (int j = 0; j < 3; ++j)
                        D[j] = 1/Ap2 * SX::cross(n, PP[q[(j+2)%3]] - PP[q[(j+1)%3]]);
                    LL_I[i] = SX::horzcat({D[0], D[1], D[2]});
                    LL_I[i] = SX::vertcat({SX::vertsplit(LL_I[i], 1)[0], -SX::vertsplit(LL_I[i], 1)[1]});
                }
                std::array<SX, 3> b_lbl, b_phi0, b_n, l0;
                for (int i = 0; i < 3; ++i){
                    b_lbl[i] = res["b_lbl" + std::to_string(i)];
                    b_phi0[i] = res["b_phi0" + std::to_string(i)];
                    b_n[i] = res["b_n" + std::to_string(i)];
                    l0[i] = res["b_l0" + std::to_string(i)];
                }
                SX t = res["n"];
                SX lambda = _Ap / _Aq;
                SX phi_d_M = SX::mtimes(L, Q.T()).T();
                std::array<SX, 3> N_M_I;
                std::array<SX, 3> N_I_I;
                std::array<SX, 3> phi_d_I;
                std::array<SX, 3> loc_n_m;
                for (int i = 0; i < 3; ++i) {
                    std::array<int, 3> q = {(i+1)%3,  i + 3, (i + 2)%3};
                    SX lPP = SX::horzcat({_P[q[0]], _P[q[1]], _P[q[2]]});
                    phi_d_I[i] = (phi_d_M + SX::mtimes(LL_I[i], lPP.T()).T())/2;

                    if (is_bnd) {
                        SX phi_d_s = (_P[(i + 2) % 3] - _P[(i + 1) % 3]) / l0[i];
                        SX phi_d_n = b_phi0[i] / (lambda * SX::norm_2(phi_d_s));
                        loc_n_m[i] = SX::vertcat({SX::horzcat({b_n[i](0, 0), -b_n[i](1, 0)}),
                                                  SX::horzcat({b_n[i](1, 0), b_n[i](0, 0)})});
                        SX phi_d_B = SX::mtimes(loc_n_m[i], SX::vertcat({phi_d_n.T(), phi_d_s.T()})).T();
                        phi_d_I[i] = SX_switch(b_lbl[i], phi_d_I[i], SX_Case{1, SX::zeros(3, 2)}, SX_Case{2, phi_d_B});
                    }
                }
                auto L_I = [L](int i){
                    SX LI(3, 2);
                    LI(0, 0) = L(0, i);
                    LI(1, 1) = L(1, i);
                    LI(2, 0) = L(1, i);
                    LI(2, 1) = L(0, i);
                    return LI;
                };
                SX H = 2 * (SX::mtimes(L_I(0), phi_d_I[0].T()) + SX::mtimes(L_I(1), phi_d_I[1].T()) + SX::mtimes(L_I(2), phi_d_I[2].T()));
                SX k = SX::mtimes(H, t);
                k(2, 0) /= 2;

                SX phi_contr = lambda * SX::horzcat(
                        {SX::cross(SX::horzsplit(phi_d_M, 1)[1], t), -SX::cross(SX::horzsplit(phi_d_M, 1)[0], t)});
                SX Ro = SX::mtimes({H, phi_contr, L});

                SX Bb;
                if (true) {
                    SX dA = SX::sym("dA", 1, 6);
                    auto du = SX::horzsplit(dA, 1);
                    SX dQ = SX::horzcat({du[0], du[1], du[2]});
                    std::array<SX, 3> dH_I;
                    SX dk(3, 1);
                    if (!is_bnd){
                        for (int i = 0; i < 3; ++i) {
                            std::array<int, 3> q = {(i+1)%3,  i + 3, (i + 2)%3};
                            SX lPP = SX::horzcat({du[q[0]], du[q[1]], du[q[2]]});
                            dH_I[i] = 2 * SX::mtimes(L_I(i),
                                                     (SX::mtimes(L, dQ.T()) + SX::mtimes(LL_I[i], lPP.T())) / 2);
                            dk += dH_I[i];
                        }
                        dk -= SX::mtimes({Ro, dQ.T()});
                    } else {
                        for (int i = 0; i < 3; ++i) {
                            std::array<int, 3> q = {(i+1)%3,  i + 3, (i + 2)%3};
                            SX lPP = SX::horzcat({du[q[0]], du[q[1]], du[q[2]]});
                            dH_I[i] = 2 * SX::mtimes(L_I(i),
                                                     (SX::mtimes(L, dQ.T()) + SX::mtimes(LL_I[i], lPP.T())) / 2);
                            SX dH_bnd1 = 2 * SX::mtimes({L_I(i), L, dQ.T()});
                            SX dH_bnd2 = 2 / l0[i] * SX::mtimes({L_I(i), loc_n_m[i], SX::vertcat(
                                    {0, (du[(i + 2) % 3] - du[(i + 1) % 3]).T()})});
                            dH_I[i] = SX_switch(b_lbl[i], dH_I[i], SX_Case{1, dH_bnd1}, SX_Case{2, dH_bnd2});

                            dk += dH_I[i];
                        }
                        dk -= SX::mtimes({Ro, dQ.T()});
                    }
                    Bb = SX::jacobian(dk, dA);
                }
                return std::pair{k, Bb};
            };
            auto _int = Bb_func(false);
            auto _bnd = Bb_func(true);
            LocalVar k_int("k_int", _int.first);
            LocalVar Bb_int("BbMatrix_int", _int.second);
            LocalVar k_bnd("k_bnd", _bnd.first);
            LocalVar Bb_bnd("BbMatrix_bnd", _bnd.second);
            res.push_back({Bb_int, k_int, Bb_bnd, k_bnd});
        }

        req = true;
    }
    return res;
}
#endif

BendingForce::vector<HH::DefaultGenInputVar> BendingForce::getDefaultAvailableInputVarsSimple() {
    vector<DefaultGenInputVar> res = DefaultHyperElasticForce::getDefaultAvailableInputVars();
    DefaultGenInputVar N("N_in", 4*2*3, 1);
    struct NRequier: public DataRequier{
        ConstProperty<World3d::F_ind, std::array<DReal, 4 * 2 * 3>> N;
        void registerObj(Object3D* obj) override {
            N = set_N_vecs(*obj);
        }
        void saveData(const DReal** target, F_ind f, std::vector<V_ind>& around) override {
            *target = &N[f][0];
        }
        [[nodiscard]] std::unique_ptr<DataRequier> copy() const override { return std::make_unique<NRequier>(*this); }
    };
    N.setDataRequier(std::make_unique<NRequier>());
    DefaultGenInputVar JBend("JBend", 2, 2*3);
    struct JBendRequier: public DataRequier{
        ConstProperty<World3d::F_ind, std::array<DReal, 2 * 2 * 3>> JBend;
        void registerObj(Object3D* obj) override {
            JBend = set_JBend(*obj);
        }
        void saveData(const DReal** target, F_ind f, std::vector<V_ind>& around) override {
            *target = &JBend[f][0];
        }
        [[nodiscard]] std::unique_ptr<DataRequier> copy() const override { return std::make_unique<JBendRequier>(*this); }
    };
    JBend.setDataRequier(std::make_unique<JBendRequier>());
    DefaultGenInputVar Pnt[6];
    using PntRequier = DefaultHyperElasticForce::PntRequier;
    using PntInitRequier = DefaultHyperElasticForce::PntInitRequier;
    for (int i = 3; i < 6; ++i) {
        Pnt[i - 3] = DefaultGenInputVar("Pnt" + std::to_string(i), 3, 1, std::make_unique<PntRequier>(i));
        Pnt[i] = DefaultGenInputVar("Pnt" + std::to_string(i) + "_init", 3, 1, std::make_unique<PntInitRequier>(i));
    }
    DefaultGenInputVar Boundary("Bnd_in", 3*(1+3+2+1), 1);
    Boundary.setDataRequier(std::make_unique<BoundaryRequier>(nullptr));
    res.push_back(std::move(N));
    res.push_back(std::move(JBend));
    for (int i = 0; i < 6; ++i)
        res.push_back(std::move(Pnt[i]));
    res.push_back(std::move(Boundary));
#ifdef BEND_WATCH
    struct KExactRequier: public DataRequier{
                std::array<double, 3> k_exact = {0, 0, 0};
                Object3D* _obj;
                void registerObj(Object3D* obj) override {
                    _obj = obj;
                }
                void saveData(const DReal** target, F_ind f, std::vector<V_ind>& around) override {
                    double s = (_obj->m_x0[around[0]].x() + _obj->m_x0[around[1]].x() + _obj->m_x0[around[2]].x())/3;
                    g_Kproxy.set(s);
                    k_exact[0] = g_Kproxy.get();
                    *target = k_exact.data();
                }
                [[nodiscard]] std::unique_ptr<DataRequier> copy() const override { return std::make_unique<KExactRequier>(*this); }
            };
            DefaultGenInputVar K("k_exact", 3, 1, std::make_unique<KExactRequier>());
            res.push_back(K);
#endif
    return res;
}

void BendingForce::k0Requier::registerObj(Object3D *obj) {
    k0 = obj->m_mesh.add_property_map<F_ind, std::array<DReal, 3>>("f:BendingForce:k0");
    if (k0.second) {
        static bool regen = true;
        auto input = getDefaultAvailableInputVarsSimple();
        for (auto& i: input){
            if (i.name == "Bnd_in")
                static_cast<BoundaryRequier*>(i.require.get())->p = m_p;
        }

#ifdef USE_BST
        auto localVars = getDefaultLocalVarsBST();
#else
        auto localVars = getDefaultLocalVars();
#endif
        HH::correlateInputAndVariables(input, localVars);
        SX k0e = localVars["k_bnd"];
        for (int i = localVars.size() - 1; i >= 0; --i){
            k0e = SX::simplify(SX::substitute(k0e, localVars.at(i).sym, localVars.at(i).val));
        }
        for (int i = 0; i < 6; ++i)
            k0e = SX::substitute(k0e, localVars.at("Pnt"+std::to_string(i)).val, localVars.at("Pnt" + std::to_string(i) + "_init").val);
        k0e = SX::simplify(SX::substitute(k0e, localVars.at("Aq").val, localVars.at("Ap").val));
        k0e = SX::simplify(SX::substitute(k0e, localVars.at("n").val, localVars.at("np").val));
        for (int i = 0; i < 3; ++i){
            k0e = SX::simplify(SX::substitute(k0e, localVars.at("Bnd_in").val(7*i+0, 0),
                          SX::if_else(localVars.at("Bnd_in").val(7*i+0, 0) == 2, 1, localVars.at("Bnd_in").val(7*i+0, 0))));
        }

        GenFunction k_expr;
        k_expr.generatedInput = std::move(input);
        k_expr.func_expr = HH::applyVariablesToExpr(k0e, localVars);
        k_expr.update_input();
        if (regenerate && regen) {
            k_expr.generate("BendingForce_k0", gen_dir);
            regen = false;
        }
        k_expr.link_function("BendingForce_k0", gen_dir);
        k_expr.registerObj(obj);

        vector<V_ind> v(6);
        for (auto f: obj->m_mesh.faces()) {
            for (int i = 0; i < 3; ++i)
                k0.first[f][i] = 0.0;
            auto vv = vert_around(obj->m_mesh, f);
            auto vo = vert_opposite(obj->m_mesh, f);
            if (!vo[0].second && !vo[1].second && !vo[2].second)
                continue;

            for (int i = 0; i < 3; ++i) {
                v[i] = vv[i];
                v[i + 3] = vo[i].first;
            }

            auto pk = &k0.first[f][0];
            if (vo[0].second || vo[1].second || vo[2].second) {
                k_expr(&pk, f, v);
            }
        }
    }
}

BendingForce::vector<HH::DefaultGenInputVar> BendingForce::getDefaultAvailableInputVars() {
    auto res = getDefaultAvailableInputVarsSimple();
    DefaultGenInputVar k0("k0", 3, 1, std::make_unique<k0Requier>());
    res.push_back(std::move(k0));

    return res;
}

BendingForce::BendingForce(string f_name, string save_dir,
   vector<DefaultGenInputVar> input_vars, VarContainer<LocalVar> &loc_vars, VarContainer<Invariant> &invariants, SX U,
   bool regenerate) : force_name{f_name}, gen_dir{save_dir}, U{U}
{
    Init(std::move(input_vars), loc_vars, invariants, regenerate);
}

void BendingForce::compute_forces(
    vector<DefaultGenInputVar> &generatedInput, VarContainer<LocalVar> &loc_vars, VarContainer<Invariant> &invariants)
{
    force_int.set_generatedInput(generatedInput);
    force_int.func_expr = compute_force_expression(loc_vars, invariants, false);
    force_int.update_input();
    force_bnd.set_generatedInput(generatedInput);
    force_bnd.func_expr = compute_force_expression(loc_vars, invariants, true);
    force_bnd.update_input();
}

BendingForce& BendingForce::add_watch(const std::function<SX(VarContainer < LocalVar > &, VarContainer < Invariant > &, SX, bool)>& watch_expr){
#ifdef USE_DEBUG_WATCHES
    GenFunction watch_int, watch_bnd;
    watch_int.set_generatedInput(m_input_vars);
    watch_int.func_expr = watch_expr(m_loc_vars, m_invariants, U, false);
    watch_int.update_input();
    watch_bnd.set_generatedInput(m_input_vars);
    watch_bnd.func_expr = watch_expr(m_loc_vars, m_invariants, U, true);
    watch_bnd.update_input();
    std::array<GenFunction*, 2> fp = {&watch_int, &watch_bnd};
    for (int fi = 0; fi < fp.size(); ++fi){
        for (auto& i: fp[fi]->generatedInput){
            if (i.name == "Bnd_in")
                static_cast<BoundaryRequier*>(i.require.get())->p = &bdata;
            else if (i.name == "k0")
            static_cast<k0Requier*>(i.require.get())->m_p = &bdata;
        }
    }
    static int nwatch = 0;
    watches[0].emplace_back(std::move(watch_int));
    watches[1].emplace_back(std::move(watch_bnd));
    for (int i = 0; i < 2; ++i){
        string watch_name = "Bend_watch" + std::to_string(nwatch) + ((i == 0) ? "_int" : "_bnd");
        if (i == 1){
            std::cout << "Generate watch " << nwatch << " on \"" << force_name << "\" with input: \n";
            watches[i].back().print_input();
        }
        watches[i].back().generate(watch_name, gen_dir);
        watches[i].back().link_function(watch_name, gen_dir);
    }
    nwatch++;
#else
    throw std::runtime_error("This function is not avaliable in current build. Please add symbol -DUSE_DEBUG_WATCHES ");
#endif
    return *this;
}

HH::SX BendingForce::compute_current_curvature_expr(VarContainer<LocalVar> &loc_vars, VarContainer<Invariant> &invariants, SX u, bool is_bnd){
    using SX = casadi::SX;
    SX S(2, 2);
    for (const auto& i: invariants){
        S += SX::simplify((SX::jacobian(u, i.sym) * i.C2dderiv));
        S = SX::simplify(SX::substitute(S, i.sym, i.C2dval));
        u = SX::substitute(u, i.sym, i.C2dval);
    }
    S *= 2;
    auto H_it = loc_vars.find("H");
    if (H_it == loc_vars.end())
        throw std::runtime_error("To use this force you should provide parameter \"H\" that means height");
    SX H = H_it->sym;
    SX Ap = loc_vars["Ap"];
    SX Aq = loc_vars["Aq"];
    SX lambda = Ap / Aq;

    SX z = SX::sym("z");
    SX k = loc_vars["k_int"];
    SX k0 = loc_vars["k0"];
    SX nf = 0;
    SX B[3] = {SX(3, 6), SX(3, 6), SX(3, 6)};
    if (is_bnd) {
        k = loc_vars["k_bnd"];
        SX Bb = loc_vars["BbMatrix_bnd"];
        SX C2d = loc_vars["C2d"];
        SX SS = SX::simplify(S / (Ap * H));
        SX D = SX::vertcat({SS(0, 0), SS(1, 1), SS(0, 1), SS(1, 0)});
        D = SX::simplify(SX::jacobian(D, SX::horzcat({C2d(0, 0), C2d(1, 1), C2d(0, 1), C2d(1, 0)})));
        auto Ds = SX::horzsplit(D, 1);
        D = SX::horzcat({Ds[0], Ds[1], Ds[2] + Ds[3]});
        Ds = SX::vertsplit(D, 1);
        D = SX::vertcat({Ds[0], Ds[1], Ds[2] + Ds[3]});

        SX xi = SX::vertcat({(k - k0)(0, 0), (k - k0)(1, 0), (k - k0)(2, 0)});
        nf = SX::if_else_zero(loc_vars["b_lbl0"] == 1, 1) +
                SX::if_else_zero(loc_vars["b_lbl1"] == 1, 1) +
                SX::if_else_zero(loc_vars["b_lbl2"] == 1, 1);
        SX n = SX::if_else_zero(loc_vars["b_lbl0"] == 1, loc_vars["b_n0"]) +
                SX::if_else_zero(loc_vars["b_lbl1"] == 1, loc_vars["b_n1"]) +
                SX::if_else_zero(loc_vars["b_lbl2"] == 1, loc_vars["b_n2"]);
        SX t = SX::vertcat({-n(1, 0), n(0, 0)});
        SX NN = SX::vertcat({SX::sq(n(0, 0)), SX::sq(n(1, 0)), n(0, 0) * n(1, 0)});
        SX Tn = SX::vertcat(
                {n(0, 0) * t(0, 0), t(1, 0) * n(1, 0), (n(1, 0) * t(0, 0) + n(0, 0) * t(1, 0)) / 2});
        SX T = SX::vertcat({SX::sq(t(0, 0)), SX::sq(t(1, 0)), t(0, 0) * t(1, 0)});

        SX A = SX::vertcat({SX::mtimes(NN.T(), D), SX::mtimes(Tn.T(), D), T.T()});

        SX R = SX::vertcat({0, 0, T(0, 0) * xi(0, 0) + T(1, 0) * xi(1, 0) + 2 * T(2, 0) * xi(2, 0)});
        SX R_is_zero = (SX::abs(R(2, 0)) <= (10*(SX::norm_1(k) + SX::norm_1(k0)) * std::numeric_limits<double>::epsilon())); 
        R = SX::if_else_zero(nf < 1.5 && nf > 0.5, R);
        SX invA = SX::if_else_zero(nf < 1.5 && nf > 0.5, SX::inv(A));
        xi = SX::if_else_zero(!R_is_zero, SX::horzsplit(invA, 1)[2] * R(2, 0));
        k = SX::if_else(nf < 1.5 && nf > 0.5, xi + k0, k);
        // B[0] = Bb;
        // B[1] = SX::horzcat({SX::vertcat({C2d(0, 0), C2d(1, 1), C2d(0, 1)}), D, SX::vertcat({nf, t})});
        // B[2] = SX::horzcat({A, invA, SX::vertcat({SX::det(A), R(2, 0), R_is_zero})});
    }

    SX res = SX::horzcat({k});

    for (int i = loc_vars.size() - 1; i >= 0; --i) {
        res = SX::simplify(SX::substitute(res, loc_vars.at(i).sym, loc_vars.at(i).val));
    }
    return res;
} 

HH::SX BendingForce::compute_force_expression(VarContainer<LocalVar> &loc_vars, VarContainer<Invariant> &invariants, bool is_bnd) {
    SX u = U;
    SX S(2, 2);
    for (const auto& i: invariants){
        S += SX::simplify((SX::jacobian(u, i.sym) * i.C2dderiv));
        S = SX::simplify(SX::substitute(S, i.sym, i.C2dval));
        u = SX::substitute(u, i.sym, i.C2dval);
    }
    S *= 2;
    auto H_it = loc_vars.find("H");
    if (H_it == loc_vars.end())
        throw std::runtime_error("To use this force you should provide parameter \"H\" that means height");
    SX H = H_it->sym;
    SX Ap = loc_vars["Ap"];
    SX Aq = loc_vars["Aq"];
    SX lambda = Ap / Aq;

    SX z = SX::sym("z");
    SX k = loc_vars["k_int"];
    SX k0 = loc_vars["k0"];
    SX nf = 0;
    SX B[3] = {SX(3, 6), SX(3, 6), SX(3, 6)};
    if (is_bnd) {
        k = loc_vars["k_bnd"];
        SX Bb = loc_vars["BbMatrix_bnd"];
        SX C2d = loc_vars["C2d"];
        SX SS = SX::simplify(S / (Ap * H));
        SX D = SX::vertcat({SS(0, 0), SS(1, 1), SS(0, 1), SS(1, 0)});
        D = SX::simplify(SX::jacobian(D, SX::horzcat({C2d(0, 0), C2d(1, 1), C2d(0, 1), C2d(1, 0)})));
        auto Ds = SX::horzsplit(D, 1);
        D = SX::horzcat({Ds[0], Ds[1], Ds[2] + Ds[3]});
        Ds = SX::vertsplit(D, 1);
        D = SX::vertcat({Ds[0], Ds[1], Ds[2] + Ds[3]});

        SX xi = SX::vertcat({(k - k0)(0, 0), (k - k0)(1, 0), (k - k0)(2, 0)});
        nf = SX::if_else_zero(loc_vars["b_lbl0"] == 1, 1) +
             SX::if_else_zero(loc_vars["b_lbl1"] == 1, 1) +
             SX::if_else_zero(loc_vars["b_lbl2"] == 1, 1);
        SX n = SX::if_else_zero(loc_vars["b_lbl0"] == 1, loc_vars["b_n0"]) +
               SX::if_else_zero(loc_vars["b_lbl1"] == 1, loc_vars["b_n1"]) +
               SX::if_else_zero(loc_vars["b_lbl2"] == 1, loc_vars["b_n2"]);
        SX t = SX::vertcat({-n(1, 0), n(0, 0)});
        SX NN = SX::vertcat({SX::sq(n(0, 0)), SX::sq(n(1, 0)), n(0, 0) * n(1, 0)});
        SX Tn = SX::vertcat(
                {n(0, 0) * t(0, 0), t(1, 0) * n(1, 0), (n(1, 0) * t(0, 0) + n(0, 0) * t(1, 0)) / 2});
        SX T = SX::vertcat({SX::sq(t(0, 0)), SX::sq(t(1, 0)), t(0, 0) * t(1, 0)});

        SX A = SX::vertcat({SX::mtimes(NN.T(), D), SX::mtimes(Tn.T(), D), T.T()});

        SX R = SX::vertcat({0, 0, T(0, 0) * xi(0, 0) + T(1, 0) * xi(1, 0) + 2 * T(2, 0) * xi(2, 0)});
        R = SX::if_else_zero(nf < 1.5 && nf > 0.5, R);
        SX R_is_zero = (SX::abs(R(2, 0)) <= (10*(SX::norm_1(k) + SX::norm_1(k0)) * std::numeric_limits<double>::epsilon())); 
        SX invA = SX::if_else_zero(nf < 1.5 && nf > 0.5 && !R_is_zero, SX::inv(A));
        xi = SX::horzsplit(invA, 1)[2] * R(2, 0);
        k = SX::if_else(nf < 1.5 && nf > 0.5, xi + k0, k);

        SX A2 = SX::horzsplit(invA, 1)[2]; A2(2,0) *= 2;
        for (int i = 0; i < 3; ++i) {
            B[i] = SX::simplify(A2(i, 0)*SX::mtimes({loc_vars["n"], T.T(), Bb}));
            auto L = SX::vertsplit(loc_vars["L"], 1);
            SX Q = loc_vars["Q"];
            L[0] = L[0].T(), L[1] = L[1].T();
            SX dA = SX::jacobian(A2(i, 0), C2d(0, 0)) * SX::mtimes({2*Q, L[0], L[0].T()}) +
                    SX::jacobian(A2(i, 0), C2d(1, 1)) * SX::mtimes({2*Q, L[1], L[1].T()}) +
                    2 * SX::jacobian(A2(i, 0), C2d(0, 1)) * SX::mtimes({Q, SX::mtimes(L[1], L[0].T()) + SX::mtimes(L[0], L[1].T())});
            dA = SX::horzcat({dA, SX::zeros(3, 3)});
            B[i] += SX::simplify(dA * R(2, 0));
        }
    }

    SX dk = k - k0;
    SX xi = SX::vertcat({SX::horzcat({dk(0, 0), dk(2, 0)}), SX::horzcat({dk(2,0), dk(1, 0)})});

    SX C2d = loc_vars["C2d"] + 2 * z * xi;
    S = SX::substitute(S, loc_vars["C2d"], C2d );
    SX tmp(3 ,1);

    tmp(0, 0) = S(0, 0), tmp(1, 0) = S(1, 1), tmp(2, 0) = S(0, 1);

    SX M = SX::simplify(
            Num_integrator1D::integrate(tmp * z, z, -H * lambda/2.0, H * lambda/2.0, 9) / (H * lambda));
    SX N = SX::simplify(
            Num_integrator1D::integrate(tmp - SX::substitute(tmp, z, 0), z, -H * lambda/2.0, H * lambda/2.0, 9) / (H * lambda));
    SX Q = loc_vars["Q"];
    auto L = SX::horzsplit(loc_vars["L"].T(), 1);
    SX force_e_part = -SX::mtimes({Q, N(0,0)*SX::mtimes(L[0], L[0].T()) + N(1,0)*SX::mtimes(L[1], L[1].T()) +
                                      N(2,0)*(SX::mtimes(L[0], L[1].T()) + SX::mtimes(L[1], L[0].T()))});
    SX force_e = SX::horzcat({force_e_part, SX(3, 3)});

    SX xiM = (dk(0, 0) * M(0,0) + 2*dk(1, 0) * M(1,0) + dk(2, 0) * M(2,0))/lambda;
    SX dQ = SX::sym("dQ", 3, 3);
    auto dQQ = SX::horzsplit(dQ, 1);
    SX dl = (2 * Ap) / SX::norm_2(SX::cross(dQQ[1] - dQQ[0], dQQ[2] - dQQ[0]));
    SX G = xiM * SX::jacobian(dl, SX::horzcat({dQQ[0].T(), dQQ[1].T(), dQQ[2].T()}));
    auto GG = SX::horzsplit(G, 3);
    G = SX::horzcat({GG[0].T(), GG[1].T(), GG[2].T()});
    G = SX::substitute(G, dQ, Q);
    SX force_t = -SX::horzcat({G, SX(3, 3)});


    SX Bb = (is_bnd) ? loc_vars["BbMatrix_bnd"] : loc_vars["BbMatrix_int"];
    SX t = loc_vars["n"];
    SX force_b = -SX::mtimes({t, M.T(), Bb});
    if (is_bnd) {
        SX force_B = -(M(0, 0) * B[0] + M(1, 0) * B[1] + M(2, 0) * B[2]);
        force_b = SX::if_else(nf < 0.5, force_b, force_B);
    }
    SX force = force_b;// + force_e + force_t;

    for (int i = loc_vars.size() - 1; i >= 0; --i){
        force = SX::simplify(SX::substitute(force, loc_vars.at(i).sym, loc_vars.at(i).val));
    }

    return force;
}

void BendingForce::Init(vector<DefaultGenInputVar> input_vars, VarContainer<LocalVar> &loc_vars,
                        VarContainer<Invariant> &invariants, bool regenerate) {
    type = "BendingForce";
    if (gen_dir != "" && gen_dir[gen_dir.size()-1] != '/')
        gen_dir += "/";
    for (auto& i: input_vars){
        if (i.name == "k0") {
            static_cast<k0Requier *>(i.require.get())->regenerate = regenerate;
            static_cast<k0Requier *>(i.require.get())->gen_dir = gen_dir;
            break;
        }
    }
    HH::correlateInputAndVariables(input_vars, loc_vars, invariants);
    compute_forces(input_vars, loc_vars, invariants);
#ifdef USE_DEBUG_WATCHES
    m_input_vars = std::move(input_vars);
    m_loc_vars = std::move(loc_vars);
    m_invariants = std::move(invariants);
#endif
    actualize_internal_data();
    if (regenerate)
        generate();
    link_function();
}

int BendingForce::element_matrix(Object3D &obj, ForceBase::LocMatrix &matr, F_ind f) {
    if (&obj != _jobj) registerjObj(&obj);
    _vdat.resize(6);
    auto& m = obj.m_mesh;
    matr.setZero();
    matr.resize(2*3*3, 2*3*3);
    int status = 0;
    DReal* pforce = matr.dat.data();
    auto vv = vert_around(m, f);
    auto vo = vert_opposite(obj.m_mesh, f);
    if (!vo[0].second && !vo[1].second && !vo[2].second){
        for (int i = 0; i < 2*3*3; ++i)
            matr[i][i] = 1;
        return 0;
    }
    auto& v = _vdat;
    for (int i = 0; i < 3; ++i) {
        v[i] = vv[i];
        v[i+3] = vo[i].second ? vo[i].first : vv[0];
    }
    if (!vo[0].second || !vo[1].second || !vo[2].second){
        status = jacobian_bnd(&pforce, f, v, jacob_verb-1);
        if (status < 0){
            if (jacob_verb >= 1) std::cout << "Error(" << status << ") in jacobian of bending force \"" << force_name << "\" at f = " << f << " on boundary part of obj \"" << obj.name << "\"" << std::endl;
            return status;
        }
        for (auto i = 0; i < 3; ++i) {
            auto bc = obj.getBC(v[i]);
            for (int l = 0; l < 3; ++l)
                if (!(bc & (1 << l))) matr.ApplyDirichlet(3 * i + l);
        }
        for (auto i = 3; i < 6; ++i)
            if (!vo[i - 3].second){
                for (int l = 0; l < 3; ++l)
                    matr.ApplyNonExists(3*i+l);
            } else  {
                auto bc = obj.getBC(v[i]);
                for (int l = 0; l < 3; ++l)
                    if (!(bc & (1 << l))) matr.ApplyDirichlet(3 * i + l);
            }
    } else {
        status = jacobian_int(&pforce, f, v, jacob_verb-1);
        if (status < 0){
            if (jacob_verb >= 1) std::cout << "Error(" << status << ") in jacobian of bending force \"" << force_name << "\" at f = " << f << " on internal part of obj \"" << obj.name << "\"" << std::endl;
            return status;
        }
        for (auto i = 0; i < 6; ++i) {
            auto bc = obj.getBC(v[i]);
            for (int l = 0; l < 3; ++l)
                if (!(bc & (1 << l))) matr.ApplyDirichlet(3 * i + l);
        }
    }

    return 0;
}

#ifdef USE_DEBUG_WATCHES
Eigen::MatrixXd World3d::BendingForce::eval_watch(Object3D& obj, F_ind f, GenFunction* watch_int, GenFunction* watch_bnd, bool print) {
    Eigen::MatrixXd mres;
    if (!watch_int && !watch_bnd) return mres;
    int w1s = watch_int ? watch_int->func->numel_out() : 0;
    int w2s = watch_bnd ? watch_bnd->func->numel_out() : 0;
    long long ssr[2] = {watch_int ? watch_int->func->sparsity_out(0).rows() : 0, watch_bnd ? watch_bnd->func->sparsity_out(0).rows() : 0};
    long long ssc[2] = {watch_int ? watch_int->func->sparsity_out(0).columns() : 0, watch_bnd ? watch_bnd->func->sparsity_out(0).columns() : 0};
    mres.resize(std::max(ssr[0], ssr[1]), std::max(ssc[0], ssc[1]));
    DReal* pforce = mres.data();

    _vdat.resize(6);
    auto& m = obj.m_mesh;
    auto vv = vert_around(m, f);
    auto vo = vert_opposite(obj.m_mesh, f);
    if (!vo[0].second && !vo[1].second && !vo[2].second)
        throw std::runtime_error("Bending watches doesn't computable on isolated elements");
    auto& v = _vdat;
    for (int i = 0; i < 3; ++i) {
        v[i] = vv[i];
        v[i+3] = vo[i].second ? vo[i].first : vv[0];
    }
    if ((!vo[0].second || !vo[1].second || !vo[2].second) && watch_bnd){
        watch_bnd->registerObj(&obj);
        watch_bnd->operator()(&pforce, f, v, print ? 10 : 0);
    } else if (watch_int){
        watch_int->registerObj(&obj);
        watch_int->operator()(&pforce, f, v, print ? 10 : 0);
    }
    return mres;
}
std::vector<Eigen::MatrixXd> World3d::BendingForce::eval_watches(Object3D& obj, F_ind f, bool print){
    std::vector<Eigen::MatrixXd> res(watches[0].size());
    for (int i = 0; i < watches[0].size(); ++i){
        if (print) std::cout << "Watch " << i << ":\n";
        res[i] = eval_watch(obj, f, &watches[0][i], &watches[1][i], print);
    }
    return res;
}
#endif

int BendingForce::operator()(Object3D &obj) {
    _vdat.resize(6);
    if (&obj != _obj) registerObj(&obj);
    auto& m = obj.m_mesh;
    DReal force[9*2];
    DReal* pforce = &(force[0]);

    int status = 0;
    for (auto f: m.faces()) {
        auto vv = vert_around(m, f);
        auto vo = vert_opposite(obj.m_mesh, f);
        if (!vo[0].second && !vo[1].second && !vo[2].second)
            continue;

        auto& v = _vdat;
        for (int i = 0; i < 3; ++i) {
            v[i] = vv[i];
            v[i+3] = vo[i].second ? vo[i].first : vv[0];
        }

        if ((!vo[0].second || !vo[1].second || !vo[2].second)){
            status = force_bnd(&pforce, f, v, force_verb - 1);
            if (status < 0) {
                if (force_verb >= 1) std::cout << "Error(" << status << ") in bending force \"" << force_name << "\" at f = " << f << " on boundary part of obj \"" << obj.name << "\"" << std::endl;
                return status;
            }
            for (auto i = 0; i < 3; ++i)
                obj.m_F[v[i]] += obj.withBCmask(v[i], Vector(force[3*i + 0], force[3*i + 1], force[3*i + 2]));
            for (auto i = 3; i < 6; ++i)
                if (vo[i - 3].second)
                    obj.m_F[v[i]] += obj.withBCmask(v[i], Vector(force[3*i + 0], force[3*i + 1], force[3*i + 2]));
        } else {
            status = force_int(&pforce, f, v, force_verb - 1);
            if (status < 0) {
                if (force_verb >= 1) std::cout << "Error(" << status << ") in bending force \"" << force_name << "\" at f = " << f << " on internal part of obj \"" << obj.name << "\"" << std::endl;
                return status;
            }
            for (auto i = 0; i < 6; ++i)
            {
                obj.m_F[v[i]] += obj.withBCmask(v[i], Vector(force[3 * i + 0], force[3 * i + 1], force[3 * i + 2]));
            }
        }
    }
    return status;
}

void BendingForce::registerjObj(Object3D *obj) {
    _jobj = obj;
    jacobian_int.registerObj(obj);
    jacobian_bnd.registerObj(obj);
}

void BendingForce::registerObj(Object3D *obj) {
    _obj = obj;
    force_int.registerObj(obj);
    force_bnd.registerObj(obj);
}

void BendingForce::link_function() {
    force_int.link_function(force_name + "_int", gen_dir);
    force_bnd.link_function(force_name + "_bnd", gen_dir);
}

void BendingForce::generate() {
    std::cout << "Generate bending force for potential: " << U << std::endl;
    force_bnd.print_input();
    force_int.generate(force_name + "_int", gen_dir);
    force_bnd.generate(force_name + "_bnd", gen_dir);
}

BendingForce::BendingForce(const DefaultHyperElasticForce &f, bool regenerate) :
        force_name{f.force_name + "_bending"}, gen_dir{f.gen_dir}, U{f.U}
{
    vector<DefaultGenInputVar> gen_input = BendingForce::getDefaultAvailableInputVars();
    vector<DefaultGenInputVar> input_vars = f.inputVars;
    for (int i = 0; i < input_vars.size(); ++i){
        auto it = std::find_if(gen_input.begin(), gen_input.end(), [nm = input_vars[i].getName()](DefaultGenInputVar& v) { return v.getName() == nm;});
        if (it == gen_input.end())
            gen_input.push_back(std::move(input_vars[i]));
        else
            *it = std::move(input_vars[i]);
    }

#ifdef USE_BST
    auto gen_loc = BendingForce::getDefaultLocalVarsBST();
#else
    auto gen_loc = getDefaultLocalVars();
#endif
    auto loc_vars = f.localVars;
    for (int i = 0; i < loc_vars.size(); ++i){
        auto it = gen_loc.find(loc_vars.at(i).getName());
        if (it == gen_loc.end())
            gen_loc.push_back(loc_vars.at(i));
        else
            *it = loc_vars.at(i);
    }
    auto invariants = f.invariants;
    Init(std::move(gen_input), gen_loc, invariants, regenerate);
}

void BendingForce::prepareJacobianFunction(bool regenerate) {
    auto get_input = [](GenFunction& gi, string name){
        return std::find_if(gi.generatedInput.begin(), gi.generatedInput.end(),
                            [&name](DefaultGenInputVar& v) { return v.getName() == name;});
    };
    std::array<SX, 6> q;
    SX Aq, n, S2, S2n, Q, rforce;
    for (int i = 0; i < 6; ++i) q[i] = get_input(force_int, "Pnt" + std::to_string(i))->var;
    Aq = get_input(force_int, "Aq")->var, n = get_input(force_int, "n")->var;
    S2 = SX::cross(q[1] - q[0], q[2] - q[0]); S2n = SX::norm_2(S2);
    Q = SX::horzcat({q[0], q[1], q[2], q[3], q[4], q[5]});
    rforce = force_int.func_expr;
    force_int.func_expr = SX::substitute(force_int.func_expr, Aq, S2n/2);
    force_int.func_expr = SX::simplify(SX::substitute(force_int.func_expr, n, S2/S2n));
    jacobian_int = HH::makeJacobianFunction(force_name+"_int", gen_dir, force_int, Q, regenerate, false);
    force_int.func_expr = rforce;

    for (int i = 0; i < 6; ++i) q[i] = get_input(force_bnd, "Pnt" + std::to_string(i))->var;
    Aq = get_input(force_bnd, "Aq")->var, n = get_input(force_bnd, "n")->var;
    Q = SX::horzcat({q[0], q[1], q[2], q[3], q[4], q[5]});
    S2 = SX::cross(q[1] - q[0], q[2] - q[0]); S2n = SX::norm_2(S2);
    rforce = force_bnd.func_expr;
    force_bnd.func_expr = SX::substitute(force_bnd.func_expr, Aq, S2n/2);
    force_bnd.func_expr = SX::simplify(SX::substitute(force_bnd.func_expr, n, S2/S2n));
    jacobian_bnd = HH::makeJacobianFunction(force_name+"_bnd", gen_dir, force_bnd, Q, regenerate, true);
    force_bnd.func_expr = rforce;
}

BendingForce &BendingForce::operator=(BendingForce &&a) noexcept {
    if (this == &a) return *this;
    type = a.type;
    U = a.U;
    force_bnd = std::move(a.force_bnd);
    force_int = std::move(a.force_int);
    jacobian_bnd = std::move(a.jacobian_bnd);
    jacobian_int = std::move(a.jacobian_int);
#ifdef USE_DEBUG_WATCHES
    watches = std::move(a.watches);
    m_input_vars = std::move(a.m_input_vars);
    m_loc_vars = std::move(a.m_loc_vars);
    m_invariants = std::move(a.m_invariants);
#endif
    bdata = std::move(a.bdata);
    _obj = a._obj;
    _jobj = a._jobj;
    force_name = std::move(a.force_name);
    gen_dir = std::move(a.gen_dir);
    force_verb = a.force_verb; jacob_verb = a.jacob_verb;
    actualize_internal_data();
    return *this;
}

BendingForce &BendingForce::operator=(const BendingForce &a) {
    if (this == &a) return *this;
    type = a.type;
    U = a.U;
    force_bnd = a.force_bnd;
    force_int = a.force_int;
    jacobian_bnd = a.jacobian_bnd;
    jacobian_int = a.jacobian_int;
#ifdef USE_DEBUG_WATCHES
    watches = a.watches;
    m_input_vars = a.m_input_vars;
    m_loc_vars = a.m_loc_vars;
    m_invariants = a.m_invariants;
#endif
    bdata = a.bdata;
    _obj = a._obj;
    _jobj = a._jobj;
    force_name = a.force_name;
    gen_dir = a.gen_dir;
    force_verb = a.force_verb; jacob_verb = a.jacob_verb;
    actualize_internal_data();
    return *this;
}

void BendingForce::actualize_internal_data() {
    std::array<GenFunction*, 4> fp = {&force_int, &force_bnd, &jacobian_int, &jacobian_bnd};
    for (int fi = 0; fi < fp.size(); ++fi){
        for (auto& i: fp[fi]->generatedInput){
            if (i.name == "Bnd_in")
                static_cast<BoundaryRequier*>(i.require.get())->p = &bdata;
            else if (i.name == "k0")
                static_cast<k0Requier*>(i.require.get())->m_p = &bdata;
        }
    }
#ifdef USE_DEBUG_WATCHES
    for (int i = 0; i < 2; ++i)
    for (auto& f: watches[i]){
        for (auto& i: f.generatedInput){
            if (i.name == "Bnd_in")
                static_cast<BoundaryRequier*>(i.require.get())->p = &bdata;
            else if (i.name == "k0")
                static_cast<k0Requier*>(i.require.get())->m_p = &bdata;
        }
    }
#endif
}

