//
// Created by alex on 23.01.2022.
//

#include "LinearElastic.h"

void World3d::LinearElasticModel::registerObj(Object3D *obj) {
    _obj = obj;
    _Ap = set_Ap(*obj);
    _Dv = set_D_vecs(*obj);
}

int World3d::LinearElasticModel::operator()(Object3D &obj) {
    if (_obj != &obj)
        registerObj(&obj);
    auto& x = obj.m_x;
    auto& x0 = obj.m_x0;
    for (auto f: obj.m_mesh.faces()){
        auto v = vert_around(obj.m_mesh, f);
        Eigen::Matrix<DReal, 3, 3> Q, P, D;
        Q <<     x[v[0]][0],  x[v[1]][0],  x[v[2]][0],
                 x[v[0]][1],  x[v[1]][1],  x[v[2]][1],
                 x[v[0]][2],  x[v[1]][2],  x[v[2]][2];
        P <<    x0[v[0]][0], x0[v[1]][0], x0[v[2]][0],
                x0[v[0]][1], x0[v[1]][1], x0[v[2]][1],
                x0[v[0]][2], x0[v[1]][2], x0[v[2]][2];
        D <<   _Dv[f][0][0],_Dv[f][1][0],_Dv[f][2][0],
               _Dv[f][0][1],_Dv[f][1][1],_Dv[f][2][1],
               _Dv[f][0][2],_Dv[f][1][2],_Dv[f][2][2];
        auto I = Eigen::Matrix<DReal, 3, 3>::Identity();
        Eigen::Matrix<DReal, 3, 3> er = (Q-P)*D.transpose();
        auto e = (er + er.transpose())/2;
        Eigen::Matrix<DReal, 3, 3> F = -_Ap[f]*_H*(_lambda*(e.trace())*I + 2*_mu*e)*D;
        for (auto i = 0; i < 3; ++i) {
            auto bc = obj.getBC(v[i]);
            assert(!std::isnan((bc & 1) ? F(0, i) : 0 + (bc & 2) ? F(1, i) : 0 + (bc & 4) ? F(2, i) : 0));
            obj.m_F[v[i]] += Vector((bc & 1) ? F(0, i) : 0 , (bc & 2) ? F(1, i) : 0, (bc & 4) ? F(2, i) : 0);
        }
    }
    return 0;
}

int World3d::LinearElasticModel::fill_matrix(Object3D &obj, SparseMatrix::IndexFrom &m) {
    throw std::runtime_error("Jacobian is not prepared! Method is not correct");
    LocMatrix matr;
    matr.resize(3 * 3, 3 * 3);
    Eigen::Matrix<DReal, 3, 3> D, DD;
    for (auto f: obj.m_mesh.faces()) {
        matr.setZero();
        D <<    _Dv[f][0][0],_Dv[f][1][0],_Dv[f][2][0],
                _Dv[f][0][1],_Dv[f][1][1],_Dv[f][2][1],
                _Dv[f][0][2],_Dv[f][1][2],_Dv[f][2][2];
        auto v = vert_around(obj.m_mesh, f);
        DD = D.transpose()*D;
        auto ApH = _Ap[f]*_H;
        for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k)
            for (int l = 0; l < 3; ++l)//TODO: should be repaired
                matr[i + j * 3][k + l * 3] = -ApH*(_lambda*D(i,j)*D(k,l) + _mu*(D(k,j)*D(i,l) + ((i==l) ? DD(j,l): 0)));

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


World3d::LinearElasticModel1::LinearElasticModel1(double h, double lambda, double mu, std::string gen_dir,
                                                  bool regenerate) : _H{h}, _lambda{lambda}, _mu{mu} {
    type = "LinearElastic1";

    using HH = DefaultHyperElasticForce;
    SX L = SX::sym("lambda");
    SX Mu = SX::sym("mu");
    SX H  = SX::sym("H" );
    auto lv = HH::getDefaultLocalVars();
    SX Ap = lv["Ap"];
    HH::VarContainer<HH::Invariant> I;
    SX e_v = (lv["F"] + lv["F"].T()) / 2;// - SX::eye(lv["F2d"].size1());
    SX tre_v = SX::trace(e_v) - 2;
    SX tre2_v = SX::trace(SX::mtimes(e_v, e_v.T())) - 2*SX::trace(e_v) + 2;
    SX ltre_v = tre_v; lv.at("F").substituteIn(ltre_v);
    SX ltre2_v = tre2_v; lv.at("F").substituteIn(ltre2_v);
    SX tre_dQ = SX::gradient(ltre_v, lv["Q"]); //lv["D"];//SX::gradient(ltre_v, lv["Q"]);
    SX tre2_dQ = SX::gradient(ltre2_v, lv["Q"]);//2*SX::mtimes(e_v, lv["D"]);//SX::gradient(ltre2_v, lv["Q"]);
//    std::cout << "tre_dQ = " << tre_dQ << std::endl;
//    std::cout << "tre2_dQ = " << tre2_dQ << std::endl;
    SX dF3d = SX::mtimes(lv["Q"] - lv["P"], lv["D"].T());
    SX e3d = (dF3d + dF3d.T()) / 2;
    SX e2d = SX::mtimes({lv["S2d"], e3d, lv["S2d"].T()});
    lv.at("C2d").val = 2*e2d + SX::eye(2);
    SX e_c2d = (lv["C2d"] - SX::eye(2)) / 2;
    SX tre_C = SX::trace(e_c2d);//not implemented
    SX tre2_C = SX::trace(SX::mtimes(e_c2d, e_c2d));//not implemented
    SX tre_dC = SX::eye(2)/2; //not implemented
    SX tre2_dC = e_c2d; //not implemented

    I.push_back(HH::Invariant{"tr_e", tre_v, tre_C, tre_dQ, tre_dC});
    I.push_back(HH::Invariant{"tr_e2", tre2_v, tre2_C, tre2_dQ, tre2_dC});
    auto iv = HH::getDefaultAvailableInputVars();
    iv.emplace_back("lambda", L, std::make_unique<ConstDataRequier>(_lambda));
    iv.emplace_back("mu", Mu, std::make_unique<ConstDataRequier>(_mu));
    iv.emplace_back("H", H, std::make_unique<ConstDataRequier>(_H));

    SX U = Ap * H * (L/2 * SX::sq(I["tr_e"]) + Mu * I["tr_e2"]);
    f = DefaultHyperElasticForce("LinearElastic1", gen_dir, std::move(iv), lv, I, U, regenerate);
}
