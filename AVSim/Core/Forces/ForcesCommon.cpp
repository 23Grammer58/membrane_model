//
// Created by alex on 29.01.2021.
//

#include "ForcesCommon.h"

World3d::LocalCartesian::LocalCartesian(Object3D &obj, bool flat) : obj(obj), flat(flat){
    if (flat){
        auto f0 = *obj.m_mesh.faces().begin();
        auto v0 = vert_around(obj.m_mesh, f0);
        Vector t1 = Vector(obj.m_x0[v0[0]], obj.m_x0[v0[1]]);
        t1 /= sqrt(t1.squared_length());
        auto t3 = CGAL::cross_product(Vector(obj.m_x0[v0[0]], obj.m_x0[v0[1]]),
                                      Vector(obj.m_x0[v0[0]], obj.m_x0[v0[2]]));
        t3 /= sqrt(t3.squared_length());
        Vector t2 = CGAL::cross_product(t3, t1);
        t2 /= sqrt(t2.squared_length());
        if (fabs(fabs(t3[2]) - 1) < 1e-10) {
            t1 = Vector(1, 0, 0);
            t2 = (fabs(t3[2] - 1) < 1e-10 ? 1 : -1) * Vector(0, 1, 0);
        }
        if (fabs(fabs(t3[1]) - 1) < 1e-10) {
            t1 = Vector(1, 0, 0);
            t2 = (fabs(t3[1] - 1) < 1e-10 ? 1 : -1) * Vector(0, 0, -1);
        }
        if (fabs(fabs(t3[0]) - 1) < 1e-10) {
            t1 = Vector(0, 1, 0);
            t2 = (fabs(t3[0] - 1) < 1e-10 ? 1 : -1) * Vector(0, 0, 1);
        }
        init << t1[0], t1[1], t1[2],
                t2[0], t2[1], t2[2];
    }
}

Eigen::Matrix<World3d::DReal, 2, 3> World3d::LocalCartesian::operator()(World3d::F_ind f) {
    if (flat)
        return init;
    else {
        auto v0 = vert_around(obj.m_mesh, f);
        Vector t1 = Vector(obj.m_x0[v0[0]], obj.m_x0[v0[1]]);
        t1 /= sqrt(t1.squared_length());
        auto t3 = CGAL::cross_product(Vector(obj.m_x0[v0[0]], obj.m_x0[v0[1]]),
                                      Vector(obj.m_x0[v0[0]], obj.m_x0[v0[2]]));
        t3 /= sqrt(t3.squared_length());
        Vector t2 = CGAL::cross_product(t3, t1);
        t2 /= sqrt(t2.squared_length());
        init << t1[0], t1[1], t1[2],
                t2[0], t2[1], t2[2];
        return init;
    }
}

World3d::ConstProperty<World3d::F_ind, std::array<World3d::DReal, 2 * 3>> World3d::set_S2d(Object3D& obj, bool flat, bool forced){
    auto it = obj.m_mesh.add_property_map<F_ind, std::array<World3d::DReal, 2 * 3>>("f:LocalBasis");
    if (!forced && !it.second) return it;
    LocalCartesian S(obj, flat);
    auto& S2d = it.first;
    for (auto f: obj.m_mesh.faces()){
        auto q = S(f);
        for (int i = 0; i < 6; ++i)
            S2d[f][i] = q(i%2, i/2);
    }

    return it;
}

World3d::ConstProperty <World3d::F_ind, std::array<World3d::DReal, 2 * 3>> World3d::set_L(Object3D &obj, bool forced) {
    auto it = obj.m_mesh.add_property_map<F_ind, std::array<World3d::DReal, 2 * 3>>("f:BendingForce::L");
    if (!forced && !it.second) return it;
    bool flat = check_flat_initial_template(obj);
    LocalCartesian S(obj, flat);
    auto& L = it.first;
    auto Dit = set_D_vecs(obj, forced);
    auto& D_V = Dit.first;
    for (auto f: obj.m_mesh.faces()){
        auto& D = D_V[f];
        Eigen::Matrix<World3d::DReal, 3, 3> Dm;
        Dm <<   D[0][0], D[1][0], D[2][0],
                D[0][1], D[1][1], D[2][1],
                D[0][2], D[1][2], D[2][2];
        auto Lm = (S(f) * Dm).eval();
        for (int i = 0; i < 6; ++i) {
            L[f][i] = Lm(i % 2, i / 2);
        }
    }
    if (Dit.second) obj.m_mesh.remove_property_map(Dit.first);

    return it;
}

World3d::ConstProperty<World3d::F_ind, std::array<World3d::DReal, 2 * 2 * 3>> World3d::set_JBend(Object3D& obj, bool forced) {
    auto it = obj.m_mesh.add_property_map<F_ind, std::array<World3d::DReal, 2 * 2 * 3>>("f:JBend");
    if (!forced && !it.second) return it;
    auto &JBend = it.first;
    static const char _N_doubled[3][2][6] = {
            {{-1, 1, -1, 1, 0,  0}, {-1, -1, 1, 1, 0, 0}},
            {{-1, 1, 1,  0, -1, 0}, {-2, 0,  2, 0, 0, 0}},
            {{-2, 2, 0,  0, 0,  0}, {-1, 1,  1, 0, 0, -1}}
    };
    bool flat = check_flat_initial_template(obj);
    LocalCartesian S(obj, flat);

    for (auto f: obj.m_mesh.faces()) {
        std::fill(JBend[f].begin(), JBend[f].end(), 0.0);
        auto vv = vert_around(obj.m_mesh, f);
        auto vo = vert_opposite(obj.m_mesh, f);
        if (!vo[0].second && !vo[1].second && !vo[2].second)
            continue;
        std::array<V_ind, 6> v;
        for (int i = 0; i < 3; ++i) {
            v[i] = vv[i];
            v[i + 3] = vo[i].first;
        }

        Eigen::Matrix<World3d::DReal, 2, 6> __N[3];
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 2; ++j)
                for (int k = 0; k < 6; ++k)
                    __N[i](j, k) = 1.0/2 * _N_doubled[i][j][k];
        Eigen::Matrix<World3d::DReal, 3, 6> P__;
        for (auto k = 0; k < 6; ++k)
            for (int l = 0; l < 3; ++l)
                P__(l, k) = obj.m_x0[v[k]][l] - obj.m_x0[v[0]][l];
        Eigen::Matrix<World3d::DReal, 2, 2> phi_d_eta_[3];
        for (int i = 0; i < 3; ++i)
            phi_d_eta_[i] = __N[i] * (S(f) * P__).transpose();
        std::vector<World3d::DReal> Nvec(6*2*3);
        for (int k = 0; k < 3; ++k){
            if (!vo[k].second) continue;
            Eigen::Matrix<DReal, 2, 2> p = phi_d_eta_[k].inverse();
            for (int i = 0; i < 4; ++i)
                JBend[f][i + 4 * k] = p(i%2, i/2)/2;
        }
    }

    return it;
}

World3d::ConstProperty<World3d::F_ind, std::array<World3d::DReal, 4 * 2 * 3>> World3d::set_N_vecs(Object3D& obj, bool forced){
    auto it = obj.m_mesh.add_property_map<F_ind, std::array<World3d::DReal, 4 * 2 * 3>>("f:N_vecs");
    if (!forced && !it.second) return it;
    auto& N_vec = it.first;
    static const char _N_doubled[3][2][6] = {
            {{-1, 1, -1, 1,  0, 0}, {-1, -1, 1, 1, 0,  0} },
            {{-1, 1,  1, 0, -1, 0}, {-2,  0, 2, 0, 0,  0} },
            {{-2, 2,  0, 0,  0, 0}, {-1,  1, 1, 0, 0, -1} }
    };
    bool flat = check_flat_initial_template(obj);
    LocalCartesian S(obj, flat);

    for (auto f: obj.m_mesh.faces()){
        auto vv = vert_around(obj.m_mesh, f);
        auto vo = vert_opposite(obj.m_mesh, f);
        if (!vo[0].second && !vo[1].second && !vo[2].second)
            continue;
        std::array<V_ind, 6> v;
        for (int i = 0; i < 3; ++i) {
            v[i] = vv[i];
            v[i+3] = vo[i].first;
        }

//        Eigen::Matrix<double, 2, 6> __N[3];
//        for (int i = 0; i < 3; ++i)
//            for (int j = 0; j < 2; ++j)
//                for (int k = 0; k < 6; ++k)
//                    __N[i](j, k) = 1.0/2 * _N_doubled[i][j][k];
//        Eigen::Matrix<double, 3, 6> P__;
//        for (auto k = 0; k < 6; ++k)
//            for (int l = 0; l < 3; ++l)
//                P__(l, k) = obj.m_x0[v[k]][l] - obj.m_x0[v[0]][l];
//        Eigen::Matrix<double, 2, 2> phi_d_eta_[3];
//        for (int i = 0; i < 3; ++i)
//            phi_d_eta_[i] = __N[i] * (S(f) * P__).transpose();
//        std::vector<double> Nvec(6*2*3);
//        for (int k = 0; k < 3; ++k){
//            auto p = phi_d_eta_[k].inverse();
//            auto R = p * __N[k];
//            for (int i = 0; i < 6; ++i)
//                for (int j = 0; j < 2; ++j)
//                    Nvec[i + 6*j + 2*6*k] = R(j, i);
//        }

        Vector_2 p[6];
        for (auto i = 0; i < 6; ++i) {
            Vector q = obj.m_x0[v[i]] - obj.m_x0[v[0]];
            auto pp = (S(f) * Eigen::Matrix<DReal, 3, 1>(q[0], q[1], q[2])).eval();
            p[i] = {pp[0], pp[1]};
        }
        Vector_2 phi_d_eta[3][2];
        for (auto i = 0; i < 3; ++i)
            for (auto j = 0; j < 2; ++j) {
                phi_d_eta[i][j] = CGAL::NULL_VECTOR;
                for (auto k = 0; k < 6; ++k)
                    if (_N_doubled[i][j][k] != 0)
                        phi_d_eta[i][j] += _N_doubled[i][j][k] / 2.0 * p[k];
            }
        Eigen::Matrix<DReal, 2, 2> m[3];
        for (auto i = 0; i < 3; ++i)
            for (auto j = 0; j < 2; ++j) {
                m[i](0, j) = phi_d_eta[i][j][0];
                m[i](1, j) = phi_d_eta[i][j][1];
            }
        for (int k = 0; k < 3; ++k){
            auto p = m[k].inverse().transpose();
            for (int i = 0; i < 4; ++i){
                Eigen::Matrix<DReal, 2, 1> _N;
                if (i < 3)
                    _N = (p * Eigen::Matrix<DReal, 2, 1>(_N_doubled[k][0][i], _N_doubled[k][1][i]) / 2).eval();
                else
                    _N = (p * Eigen::Matrix<DReal, 2, 1>(_N_doubled[k][0][3+k], _N_doubled[k][1][3+k]) / 2).eval();
                for (int j = 0; j < 2; ++j)
                    N_vec[f][i + 4 * j + 2 * 4 * k] = _N[j];
            }
        }
    }

    return it;
}
