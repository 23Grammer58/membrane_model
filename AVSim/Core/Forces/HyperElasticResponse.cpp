#include "HyperElasticResponse.h"
#include "ForcesCommon.h"

using namespace World3d;

int HEForce::operator()(Object3D &obj){
    if (_obj != &obj) registerObj(&obj);
    auto _L = set_L(obj);
    auto _Ap = set_Ap(obj);
    for (auto f: obj.m_mesh.faces()){
        auto vv = vert_around(obj.m_mesh, f);
        auto Ldat = _L[f];
        auto Ap = _Ap[f];
        DReal Ht = m_ht->operator()(f);
        auto L = Eigen::Map<Eigen::Matrix<World3d::DReal, 2, 3>>(Ldat.data());
        Eigen::Matrix<World3d::DReal, 3, 3> Q, P;
        Q <<    obj.m_x[vv[0]][0], obj.m_x[vv[1]][0], obj.m_x[vv[2]][0],
                obj.m_x[vv[0]][1], obj.m_x[vv[1]][1], obj.m_x[vv[2]][1],
                obj.m_x[vv[0]][2], obj.m_x[vv[1]][2], obj.m_x[vv[2]][2];
        // Eigen::Matrix<World3d::DReal, 2, 3> R = L * Q.transpose();
        // // Eigen::Matrix<World3d::DReal, 2, 2> C2d = R * R.transpose();
        // // std::array<double, 3> C{C2d(0,0), C2d(1,1), C2d(0,1)}; 
        // std::array<double, 3> E{
        //     (R.row(0).dot(R.row(0).transpose()) - 1)/2,
        //     (R.row(1).dot(R.row(1).transpose()) - 1)/2,
        //     (R.row(0).dot(R.row(1).transpose())) / 2
        // };       
        P <<    obj.m_x0[vv[0]][0], obj.m_x0[vv[1]][0], obj.m_x0[vv[2]][0],
                obj.m_x0[vv[0]][1], obj.m_x0[vv[1]][1], obj.m_x0[vv[2]][1],
                obj.m_x0[vv[0]][2], obj.m_x0[vv[1]][2], obj.m_x0[vv[2]][2];  
        Eigen::Matrix<World3d::DReal, 3, 3> U = Q - P; 
        Eigen::Matrix<World3d::DReal, 3, 2> grU = U * L.transpose(), I2d = P * L.transpose();
        Eigen::Matrix<World3d::DReal, 2, 2> IgrU = I2d.transpose() * grU;
        Eigen::Matrix<World3d::DReal, 2, 2> E = (IgrU + IgrU.transpose() +  grU.transpose() * grU) / 2; 
        Eigen::Matrix<World3d::DReal, 2, 3> R = (I2d + grU).transpose();     
        auto Sd = m_mat->PK2_tensor(std::array<double, 3>{E(0, 0), E(1, 1), E(0, 1)}, f);
        Eigen::Matrix<World3d::DReal, 3, 1> S;
        S << Sd[0], Sd[1], Sd[2];
        // Eigen::Matrix<World3d::DReal, 2, 2> S;
        // S <<    Sd[0], Sd[2],
        //         Sd[2], Sd[1];
        for (int k = 0; k < 3; ++k) {
            auto v = vv[k];
            Eigen::Matrix<World3d::DReal, 3, 3> dEdQ;
            dEdQ.col(0) = R.row(0).transpose() * L(0, k);
            dEdQ.col(1) = R.row(1).transpose() * L(1, k);
            dEdQ.col(2) = R.row(0).transpose() * L(1, k) + R.row(1).transpose() * L(0, k); // /2

            Eigen::Matrix<World3d::DReal, 3, 1> F = -Ap * Ht * dEdQ * S;
            if (obj.is_movable(v))
                obj.m_F[v] += obj.withBCmask(v, Vector(F(0, 0), F(1, 0), F(2, 0)));
        }
    } 

    return 0;
}

int HEForce::fill_matrix(Object3D& obj, SparseMatrix::IndexFrom& m){
    if (_obj != &obj) registerObj(&obj);
    auto _L = set_L(obj);
    auto _Ap = set_Ap(obj);
    LocMatrix matr;
    matr.resize(3*3, 3*3);

    for (auto f: obj.m_mesh.faces()){
        matr.setZero();
        auto vv = vert_around(obj.m_mesh, f);
        auto Ldat = _L[f];
        auto Ap = _Ap[f];
        DReal Ht = m_ht->operator()(f);
        auto L = Eigen::Map<Eigen::Matrix<World3d::DReal, 2, 3>>(Ldat.data());
        Eigen::Matrix<World3d::DReal, 3, 3> Q;
        Q <<    obj.m_x[vv[0]][0], obj.m_x[vv[1]][0], obj.m_x[vv[2]][0],
                obj.m_x[vv[0]][1], obj.m_x[vv[1]][1], obj.m_x[vv[2]][1],
                obj.m_x[vv[0]][2], obj.m_x[vv[1]][2], obj.m_x[vv[2]][2];
        Eigen::Matrix<World3d::DReal, 2, 3> R = L * Q.transpose();
        std::array<double, 3> E{
            (R.row(0).dot(R.row(0).transpose())-1)/2,
            (R.row(1).dot(R.row(1).transpose())-1)/2,
            (R.row(0).dot(R.row(1).transpose()))/2
        };    
        auto SdS = m_mat->PK2_dPK2_tensor(E, f); 
        auto dS = SdS.second; 
        Eigen::Matrix<World3d::DReal, 3, 1> S;
        S << SdS.first[0], SdS.first[1], SdS.first[2]; 
        Eigen::Matrix<World3d::DReal, 3, 3> dSm;
        dSm <<  dS[0], dS[3], dS[4],
                dS[3], dS[1], dS[5],
                dS[4], dS[5], dS[2];
        for (int vi = 0; vi < 3; ++vi)
        for (int ki = 0; ki < 3; ++ki)
            for (int vj = 0; vj < 3; ++vj){
                int kj = ki;
                Eigen::Matrix<World3d::DReal, 3, 1> A;
                A << L(0,vi)*L(0,vj), L(1,vi)*L(1,vj), (L(0,vi)*L(1,vj) + L(1,vi)*L(0,vj));
                matr[ki+3*vi][kj+3*vj] += S.dot(A);
                // matr[ki+3*vi][kj+3*vj] += 2*S[0]*L(0,vi)*L(0,vj) + 2*S[1]*L(1,vi)*L(1,vj) + 2*S[2]*(L(0,vi)*L(1,vj) + L(1,vi)*L(0,vj)); 
            }
        for (int vi = 0; vi < 3; ++vi)
        for (int ki = 0; ki < 3; ++ki)
            for (int vj = 0; vj < 3; ++vj)
            for (int kj = 0; kj < 3; ++kj) {
                double A0i = L(0,vi)*R(0,ki), A1i = L(1,vi)*R(1,ki), A2i = L(0,vi)*R(1,ki) + L(1,vi)*R(0,ki);
                double A0j = L(0,vj)*R(0,kj), A1j = L(1,vj)*R(1,kj), A2j = L(0,vj)*R(1,kj) + L(1,vj)*R(0,kj);
                Eigen::Matrix<World3d::DReal, 3, 1> Ai, Aj;
                Ai << A0i, A1i, A2i;   Aj << A0j, A1j, A2j; 

                matr[ki+3*vi][kj+3*vj] += Ai.transpose() * dSm * Aj; 
                    // 4*dS[0]*A0i*A0j + 
                    // 4*dS[1]*A1i*A1j + 
                    // 4*dS[2]*A2i*A2j + 
                    // 4*dS[3]*(A0i*A1j + A1i*A0i) + 
                    // 4*dS[4]*(A0i*A2j + A2i*A0j) + 
                    // 4*dS[5]*(A1i*A2j + A2i*A1j);
                // matr[ki+3*vi][kj+3*vj] += 
                //     4*dS[0]*L(0,vi)*R(0,ki)*L(0,vj)*R(0,kj)    //S00'00
                // +   4*dS[1]*L(1,vi)*R(1,ki)*L(1,vj)*R(1,kj)    //S11'11
                // +   4*dS[2]*(L(0,vi)*R(1,ki) + L(1,vi)*R(0,ki))*(L(0,vj)*R(1,kj) + L(1,vj)*R(0,kj))   //S01'01
                // +   4*dS[3]*(L(0,vi)*R(0,ki)*L(1,vj)*R(1,kj) + L(1,vi)*R(1,ki)*L(0,vj)*R(0,kj)) //S00'11
                // +   4*dS[4]*(L(0,vi)*R(0,ki)*(L(0,vj)*R(1,kj) + L(1,vj)*R(0,kj)) + L(0,vj)*R(0,kj)*(L(0,vi)*R(1,ki) + L(1,vi)*R(0,ki))) //S00'01
                // +   4*dS[5]*(L(1,vi)*R(1,ki)*(L(1,vj)*R(0,kj) + L(0,vj)*R(1,kj)) + L(1,vj)*R(1,kj)*(L(1,vi)*R(0,ki) + L(0,vi)*R(1,ki))); //S11'01
            }
        double scl = -Ap * Ht;
        for (int vi = 0; vi < 3; ++vi)
        for (int ki = 0; ki < 3; ++ki)
            for (int vj = 0; vj < 3; ++vj)
            for (int kj = 0; kj < 3; ++kj)
                matr[ki+3*vi][kj+3*vj] *= scl;

        for (auto i = 0; i < 3; ++i) {
            auto bc = obj.getBC(vv[i]);
            for (int l = 0; l < 3; ++l){
                if (!(bc & (1 << l)))
                    matr.ApplyDirichlet(3 * i + l);
            }
        }

        for (int i = 0; i < 9; ++i)
            for (int j = 0; j < 9; ++j) if (matr[i][j] != 0.0)
                m(vv[i/3].idx()*3+i%3, vv[j/3].idx()*3+j%3) += matr[i][j];        
    }

    return 0;
}