//
// Created by Liogky Alexey on 19.05.2022.
//

#include "ForceInterface.h"
#include "Object3D.h"

using namespace World3d;
void World3d::ForceBase::LocMatrix::resize(int N, int M){
    if (LocMatrix::N == N && LocMatrix::M == M) return;
    if (LocMatrix::N < 2 || LocMatrix::M == M){
        dat.resize(N * M, 0.0);
        LocMatrix::N = N;
        LocMatrix::M = M;
    } else {
        std::vector<double> tmp = dat;
        dat.resize(N * M);
        std::fill(dat.begin(), dat.end(), 0.0);
        for (int i = 0; i < LocMatrix::N; ++i)
            for (int j = 0; j < LocMatrix::M; ++j)
                dat[i * M + j] = tmp[i * LocMatrix::M + j];
    }
}

//#define APPLYDIRICHLETINFILL
int World3d::ForceBase::fill_matrix(Object3D& obj, SparseMatrix::IndexFrom& m){
    //deleted timers
//        Timer tt, tt1;
//        auto print_T = [&tt](std::string s = ""){ std::cout << "\t\ttime = " << tt.elapsed() << " " << s << std::endl; };
//        print_T("Start fill matrix by the force");
//        double t_comp = 0, t_fill = 0;

//        tt1.reset();
    LocMatrix matr;
    int stat = 0;
    for (auto f: obj.m_mesh.faces()){
        matr.setZero();
//            t_fill += tt1.elapsed(); tt1.reset();
        if ((stat = element_matrix(obj, matr, f)) < 0) return stat;
//            t_comp += tt1.elapsed(); tt1.reset();
        if (matr.M == matr.N && matr.M == 3*3){
            auto v = vert_around(obj.m_mesh, f);
#ifdef APPLYDIRICHLETINFILL
            for (int k = 0; k < 3; ++k){
                auto bc = obj.getBC(v[k]);
                for (int l = 0; l < 3; ++l){
                    if (!(bc & (1 << l)))
                        matr.ApplyDirichlet(3 * k + l);
                }
            }
#endif
            for (int i = 0; i < 9; ++i)
                for (int j = 0; j < 9; ++j)
                    if (matr[i][j] != 0.0)
                        m(v[i/3].idx()*3+i%3, v[j/3].idx()*3+j%3) += matr[i][j];
        } else if (matr.M == matr.N && matr.M == 2*3*3){
            auto vv = vert_around(obj.m_mesh, f);
            auto vo = vert_opposite(obj.m_mesh, f);
            std::array<int, 6> v;
            for (int i = 0; i < 3; ++i) {
                v[i] = vv[i].idx();
                v[i+3] = vo[i].second ? vo[i].first.idx() : -1;
            }
#ifdef APPLYDIRICHLETINFILL
            for (int k = 0; k < 6; ++k){
                if (v[k] == -1) continue;
                auto bc = obj.getBC(V_ind(v[k]));
                for (int l = 0; l < 3; ++l){
                    if (!(bc & (1 << l)))
                        matr.ApplyDirichlet(3 * k + l);
                }
            }
#endif
            for (int i = 0; i < 2*9; ++i)
                for (int j = 0; j < 2*9; ++j)
                    if (v[i/3] != -1 && v[j/3] != -1 && matr[i][j] != 0.0)
                        m(v[i/3]*3+i%3, v[j/3]*3+j%3) += matr[i][j];
        } else {
            throw std::runtime_error("Default fill_matrix for this case is not implemented");
        }
    }
//        t_fill += tt1.elapsed();
//        std::cout << "\t\ttime element evaluation = " << t_comp << " time filling result = " << t_fill << std::endl;
//        print_T("End fill matrix by the force");

    return stat;
}