//
// Created by alex on 29.01.2021.
//

#ifndef AORTIC_VALVE_NEOGOOK_H
#define AORTIC_VALVE_NEOGOOK_H

#include "HyperElastic.h"
namespace World3d{
    class NeoGookModel: public HyperElasticForceBase{
    public:
        DReal _mu;
        DReal _h;
        NeoGookModel(double mu, double h, std::string gen_dir, bool regenerate = true):
                _mu{mu}, _h{h}
        {
            type = "NeoGookModel";
            SX Mu = SX::sym("mu");
            SX H  = SX::sym("H" );
            auto lv = DefaultHyperElasticForce::getDefaultLocalVars();
            SX Ap = lv["Ap"];
            auto I = DefaultHyperElasticForce::getDefaultInvariants();
            SX I1 = I["I1"];
            SX I3 = I["I3"];
            SX U = Ap * H * Mu / 2 * (I1 + 1.0 / I3 - 3);
            auto iv = DefaultHyperElasticForce::getDefaultAvailableInputVars();
            iv.emplace_back("mu", Mu, std::make_unique<ConstDataRequier>(_mu));
            iv.emplace_back("H", H, std::make_unique<ConstDataRequier>(_h));
            f = DefaultHyperElasticForce("NeoGook", gen_dir, std::move(iv), lv, I, U, regenerate);
        }
        std::unique_ptr<HyperElasticForceBase> copy() override{
            return std::make_unique<NeoGookModel>(_mu, _h, f.gen_dir, false);
        }
    };

    class NeoGookModel1: public HyperElasticForceBase{
    public:
        DReal _mu;
        std::string h_tag;
        NeoGookModel1(double mu, std::string h_tag, std::string gen_dir, bool regenerate = true):
                _mu{mu}, h_tag{h_tag}
        {
            type = "NeoGookModel";
            SX Mu = SX::sym("mu");
            SX H  = SX::sym("H" );
            auto lv = DefaultHyperElasticForce::getDefaultLocalVars();
            SX Ap = lv["Ap"];
            auto I = DefaultHyperElasticForce::getDefaultInvariants();
            SX I1 = I["I1"];
            SX I3 = I["I3"];
            SX U = Ap * H * Mu / 2 * (I1 + 1.0 / I3 - 3);
            auto iv = DefaultHyperElasticForce::getDefaultAvailableInputVars();
            iv.emplace_back("mu", Mu, std::make_unique<ConstDataRequier>(_mu));
            iv.emplace_back("H", H, std::make_unique<ConstDataTagRequier<double>>(h_tag));
            f = DefaultHyperElasticForce("NeoGook", gen_dir, std::move(iv), lv, I, U, regenerate);
        }
        std::unique_ptr<HyperElasticForceBase> copy() override{
            return std::make_unique<NeoGookModel1>(_mu, h_tag, f.gen_dir, false);
        }
    };

    class NeoGookModelOpt: public ForceBase{
    public:
        DReal mu;
        std::function<DReal(Object3D& obj, F_ind f)> h;
        NeoGookModelOpt(double mu, double h):
                mu{mu}, h{[h](Object3D& obj, F_ind f) -> DReal {return h;}} { type = "NeoGookModelOpt"; }
        NeoGookModelOpt(double mu, std::function<DReal(Object3D& obj, F_ind f)> h):
                mu{mu}, h{h} { type = "NeoGookModelOpt"; }
        std::unique_ptr<ForceBase> clone() override{
            return std::make_unique<NeoGookModelOpt>(mu, h);
        }
        int operator()(Object3D& obj) override{
            auto DD_ = set_DD_matrix(obj);
            auto Ap_ = set_Ap(obj);
            auto S_ = set_S(obj);
            auto x_ = obj.m_x;

            auto n_ = obj.m_normal;
            for (auto f: obj.m_mesh.faces()){
                auto v = vert_around(obj.m_mesh, f);

                Eigen::Vector3d n(n_[f][0], n_[f][1], n_[f][2]);
                auto Ap = Ap_[f];
                auto Aq = sqrt(S_[f].squared_length());
                Eigen::Matrix<DReal, 3, 3> Q, Qnm, DD;
                Q <<    x_[v[0]][0], x_[v[1]][0], x_[v[2]][0],
                        x_[v[0]][1], x_[v[1]][1], x_[v[2]][1],
                        x_[v[0]][2], x_[v[1]][2], x_[v[2]][2];
                Qnm.col(0) = n.cross(Q.col(2) - Q.col(1));
                Qnm.col(1) = n.cross(Q.col(0) - Q.col(2));
                Qnm.col(2) = - (Qnm.col(0) + Qnm.col(1));
                for (auto i = 0; i < 3; ++i)
                    for (auto j = i; j < 3; ++j)
                        DD(j, i) = DD(i, j) = DD_[f][j + (i + 1) * (i > 0)];

                DReal J = Aq / Ap;
                Eigen::Matrix<DReal, 3, 3> F = - mu * h(obj, f) * (Ap * Q * DD.transpose() - Qnm / (2 * J * J * J));
                for (auto i = 0; i < 3; ++i) {
                    auto bc = obj.getBC(v[i]);
                    assert(!std::isnan((bc & 1) ? F(0, i) : 0 + (bc & 2) ? F(1, i) : 0 + (bc & 4) ? F(2, i) : 0));
                    obj.m_F[v[i]] += Vector((bc & 1) ? F(0, i) : 0 , (bc & 2) ? F(1, i) : 0, (bc & 4) ? F(2, i) : 0);
                }
            }

            return 0;
        };

        int fill_matrix(Object3D& obj, SparseMatrix::IndexFrom& m) override{
            auto DD_ = set_DD_matrix(obj);
            auto Ap_ = set_Ap(obj);
            auto S_ = set_S(obj);
            auto x_ = obj.m_x;

            LocMatrix matr;
            matr.resize(3*3, 3*3);
            Eigen::Matrix<DReal, 3, 3> Q, DD, d, Qd, QdQT, QdTQd, QdQTQd;
            d <<     0,  1, -1,
                    -1,  0,  1,
                     1, -1,  0;

            for (auto f: obj.m_mesh.faces()){
                matr.setZero();
                auto v = vert_around(obj.m_mesh, f);
                Q <<    x_[v[0]][0], x_[v[1]][0], x_[v[2]][0],
                        x_[v[0]][1], x_[v[1]][1], x_[v[2]][1],
                        x_[v[0]][2], x_[v[1]][2], x_[v[2]][2];
                for (auto i = 0; i < 3; ++i)
                    for (auto j = i; j < 3; ++j)
                        DD(j, i) = DD(i, j) = DD_[f][j + (i + 1) * (i > 0)];
                auto Ap = Ap_[f];
                auto Aq = sqrt(S_[f].squared_length());
                DReal J = Aq / Ap; DReal J2 = J * J;
                DReal Aq2 = Aq * Aq; DReal Aq4 = Aq2 * Aq2;
                Qd = Q * d;
                QdQT = Qd * Q.transpose();
                QdTQd = Qd.transpose() * Qd;
                QdQTQd = QdQT * Qd;
                DReal muHAp = mu * h(obj, f) * Ap;
                DReal Aq2J2m4 = 4 * Aq2 * J2, Aq4J2m4 = 4 * Aq4 * J2;
                for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j)
                    for (int y = 0; y < 3; ++y)
                    for (int k = 0; k < 3; ++k)
                        matr[i + j*3][y + k*3] = - muHAp * (((i == y) ? DD(j, k) : 0) -
                                    (((i == y) ? QdTQd(k, j) : 0) - Qd(i, k)*Qd(y, j) - QdQT(i, y)*d(k, j)) / Aq2J2m4 +
                                    QdQTQd(i, j) * QdQTQd(y, k) / Aq4J2m4);

                for (auto i = 0; i < 3; ++i) {
                    auto bc = obj.getBC(v[i]);
                    for (int l = 0; l < 3; ++l){
                        if (!(bc & (1 << l)))
                            matr.ApplyDirichlet(3 * i + l);
                    }
                }

                for (int i = 0; i < 9; ++i)
                    for (int j = 0; j < 9; ++j)
                        if (matr[i][j] != 0.0)
                            m(v[i/3].idx()*3+i%3, v[j/3].idx()*3+j%3) += matr[i][j];
            }

            return 0;
        }
    };
}


#endif //AORTIC_VALVE_NEOGOOK_H
