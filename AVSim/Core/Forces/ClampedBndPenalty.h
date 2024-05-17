//
// Created by alex on 30.04.2021.
//

#ifndef AORTIC_VALVE_CLAMPEDBNDPENALTY_H
#define AORTIC_VALVE_CLAMPEDBNDPENALTY_H
#include "../Object3D.h"

namespace World3d {
    class ClampedBndPenalty : public ForceBase {
        UpdatableProperty<F_ind, Vector> _S;
        Object3D *_obj = nullptr;
        std::function<double(double, bool)> m_penalty = [](double cosPhi, bool posSinPhi) -> double{
            if (posSinPhi) return 1 - cosPhi;
            if (cosPhi <= 0) return 3 + cosPhi;
            else return cosPhi - 1;
        };
        std::function<double(double, bool)> m_penalty_deriv = [](double cosPhi, bool posSinPhi){ return (posSinPhi) ? -1.0 : 1.0; };
        bool m_updated_data = false;
    public:
        struct ClampData{
            E_ind e = E_ind(-1);
            Vector n_ext{0, 0, 0};
            F_ind f = F_ind(-1);
            std::array<V_ind, 3> fv = {V_ind(-1), V_ind(-1), V_ind(-1)};
        };
        std::vector<ClampData> m_clmp_dirs;
        double m_K = 0;

        ClampedBndPenalty(const ClampedBndPenalty& ) = default;
        ClampedBndPenalty(ClampedBndPenalty&& ) = default;
        ClampedBndPenalty(double K = 0, const std::map<E_ind, Vector>& data = std::map<E_ind, Vector>()): m_K{K} {
            if (!data.empty()) setClampedBnd(data);
        }
        //clamped_boundary is map between edges and wish external normal to element on the edge
        ClampedBndPenalty& setClampedBnd(const std::map<E_ind, Vector>& clamped_boundary) {
            m_clmp_dirs.reserve(clamped_boundary.size());
            for (auto& i: clamped_boundary){
                ClampData cd;
                cd.e = i.first, cd.n_ext = i.second;
                m_clmp_dirs.push_back(cd);
            }
            m_updated_data = true;
            return *this;
        }
        ClampedBndPenalty& setAnglePenaltyFunction(std::function<double(double cosPhi, bool sinPhiPositive)> penalty,
                                                   std::function<double(double cosPhi, bool sinPhiPositive)> penalty_deriv){
            return m_penalty = std::move(penalty), m_penalty_deriv = std::move(penalty_deriv), *this;
        }
        ClampedBndPenalty& setStiffnes(double K) { return m_K = K, *this; }
        void registerObj(Object3D *obj){
            _obj = obj;
            _S = set_S(*obj);
            for (auto& cd: m_clmp_dirs){
                auto f = face_around(obj->m_mesh, cd.e);
                cd.f = f[0].first;
                auto vv = vert_around(obj->m_mesh, cd.f);
                auto vo = vert_opposite(obj->m_mesh, cd.e).first[0];
                int j = 0;
                for (; j < 3; ++j) if (vv[j]==vo) break;
                if (j == 1) cd.fv = std::array<V_ind, 3>{vv[1], vv[2], vv[0]};
                else if (j == 2) cd.fv = std::array<V_ind, 3>{vv[2], vv[0], vv[1]};
                else cd.fv = vv;
                if (abs(cd.n_ext.squared_length() - 1) > 1e-5) cd.n_ext /= sqrt(cd.n_ext.squared_length());
            }
            m_updated_data = false;
        }
        int operator()(Object3D &obj) override{
            if (&obj != _obj || m_updated_data) registerObj(&obj);
            for (const auto& cd: m_clmp_dirs){
                Vector r = obj.m_x[cd.fv[2]] - obj.m_x[cd.fv[1]];
                double sq_r = r.squared_length();
                double e_len = sqrt(sq_r);
                Vector Q = obj.m_x[cd.fv[0]] - obj.m_x[cd.fv[1]];
                Vector q = Q - (Q*r)*r / sq_r; q /= sqrt(q.squared_length());
                auto& n = obj.m_normal.first[cd.f];
                auto& n_ext = cd.n_ext;
                obj.m_F[cd.fv[0]] += obj.withBCmask(cd.fv[0], m_K * m_penalty(-n_ext*q, m_K*(n*(-n_ext)) > 0) * e_len * n);
            }
            return 0;
        }
        std::unique_ptr<ForceBase> clone() override{ return std::make_unique<ClampedBndPenalty>(*this); }
        int fill_matrix(Object3D& obj, SparseMatrix::IndexFrom& m) override{
            for (const auto& cd: m_clmp_dirs){
                double Aq = sqrt(_S.first[cd.f].squared_length());
                Vector r = obj.m_x[cd.fv[2]] - obj.m_x[cd.fv[1]];
                double sq_r = r.squared_length();
                double e_len = sqrt(sq_r);
                Vector rnm = r / e_len;
                Vector Q = obj.m_x[cd.fv[0]] - obj.m_x[cd.fv[1]];
                Vector q = Q - (Q*r)*r / sq_r;
                double h = sqrt(q.squared_length());
                q /= h;
                auto& n = obj.m_normal.first[cd.f];
                auto& n_ext = cd.n_ext;
                Vector l = -n_ext;
                double lq = l*q; bool b = m_K*n*l > 0;
                double f_x = m_penalty(lq, b), df_x = m_penalty_deriv(lq, b);

                double rl = rnm * l, qr = q*rnm;
                Vector dxdQ = (l - (rl - lq*qr)*rnm - lq*q) / h;
                Vector nxr = CGAL::cross_product(n, r);
                std::array<double, 3> K = {r[2], -r[1], r[0]};
                auto K_mat = [K](int i, int j) { return ((i != j) ? ((i > j) ? 1 : -1) * K[i+j-1] : 0); };

                auto bc = obj.getBC(cd.fv[0]);
                for (int i = 0; i < 3; ++i) if (bc & (1 << i))
                for (int j = 0; j < 3; ++j) if (bc & (1 << j)) {
                    m(cd.fv[0].idx()*3+i, cd.fv[0].idx()*3+j) +=
                            m_K * e_len * ( df_x * n[i] * dxdQ[j] + f_x / (2*Aq) * (n[i] * nxr[j] + K_mat(i, j)) );
                }
            }

            return 0;
        }

    };
}


#endif //AORTIC_VALVE_CLAMPEDBNDPENALTY_H
