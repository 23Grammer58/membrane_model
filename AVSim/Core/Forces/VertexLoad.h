//
// Created by Liogky Alexey on 25.05.2022.
//

#ifndef AVSIM_VERTEXLOAD_H
#define AVSIM_VERTEXLOAD_H
#include <functional>
#include "../Object3D.h"


namespace World3d {
    template<typename VFPARAM = void>
    class VertexLoad : public ForceBase{
    public:
        struct VertParam{
            double weight;
            VFPARAM vfparam;
        };
        ///@param[in] x [3] - point in space
        ///@param[out] F [3] - force F(x) in the point x
        ///@param[in] prm - some parameters specific for the vertex
        ///@return true if force field is defined in x
        std::function<bool(double *x/*[3]*/, double *F/*[3]*/, VFPARAM& prm)> m_F;
        ///@param[in] x [3] - point in space
        ///@param[out] J [3][3] - J[3*i + j] = J[i][j] = ∇_j F_i(x) gradient of force field
        ///@param[in] prm - some parameters specific for the vertex
        ///@return true if gradient of force field is defined in x
        std::function<bool(double *x/*[3]*/, double *J/*[3][3]*/, VFPARAM& prm)> m_F_deriv;
        ///vector of acting vertices and weights
        std::map<V_ind, VertParam> m_vertices;

        VertexLoad<VFPARAM>& insertVertices(const std::vector<std::pair<V_ind, VertParam>>& vertices){
            for (auto& i: vertices) m_vertices.insert(i);
            return *this;
        }
        VertexLoad<VFPARAM>& setForceField(std::function<bool(double *x/*[3]*/, double *F/*[3]*/, VFPARAM& prm)> force, std::function<bool(double *x/*[3]*/, double *J/*[3][3]*/, VFPARAM& prm)> dforce){
            m_F = force;
            m_F_deriv = dforce;
            return *this;
        }
        VertexLoad<VFPARAM>& setForceField(std::function<bool(double *x/*[3]*/, double *F/*[3]*/)> force, std::function<bool(double *x/*[3]*/, double *J/*[3][3]*/)> dforce){
            m_F = [force](double *x/*[3]*/, double *F/*[3]*/, VFPARAM& prm){ return force(x, F); };
            m_F_deriv = [dforce](double *x/*[3]*/, double *J/*[3][3]*/, VFPARAM& prm){ return dforce(x, J); };
            return *this;
        }

        VertexLoad(){ type = std::string("VertexLoad<") + typeid(VFPARAM).name() + ">"; }
        VertexLoad(const VertexLoad& a) = default;
        int operator()(Object3D &obj) override{
            std::array<double, 3> x, F;
            for (auto i: m_vertices) {
                auto v = i.first;
                auto w = i.second;
                for (int k = 0; k < 3; ++k) x[k] = obj.m_x[v][k];
                m_F(x.data(), F.data(), w.vfparam);
                obj.m_F[v] += obj.withBCmask(v, w.weight * Vector(F[0], F[1], F[2]));
            }
            return 0;
        }
        std::unique_ptr<ForceBase> clone() override{ return std::make_unique<VertexLoad>(*this); }
        int fill_matrix(Object3D &obj, SparseMatrix::IndexFrom &m) override{
            std::array<double, 3> x;
            std::array<double, 9> J;
            for (auto i: m_vertices){
                auto v = i.first;
                auto w = i.second;
                m_F_deriv(x.data(), J.data(), w.vfparam);
                std::array<bool, 3> bc = { obj.is_movable(v, 0), obj.is_movable(v, 1), obj.is_movable(v, 2) };
                for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j)
                    if (bc[i] && bc[j]) m(v*3 + i, v*3 + j) += w.weight * J[3*i + j];
            }
            return 0;
        }
    };

    template<>
    class VertexLoad<void> : public ForceBase{
    public:
        ///@param[in] x [3] - point in space
        ///@param[out] F [3] - force F(x) in the point x
        ///@return true if force field is defined in x
        std::function<bool(double *x/*[3]*/, double *F/*[3]*/)> m_F;
        ///@param[in] x [3] - point in space
        ///@param[out] J [3][3] - J[3*i + j] = J[i][j] = ∇_j F_i(x) gradient of force field
        ///@return true if gradient of force field is defined in x
        std::function<bool(double *x/*[3]*/, double *J/*[3][3]*/)> m_F_deriv;
        ///vector of acting vertices and weights
        std::map<V_ind, double> m_vertices;
        using VertParam = double;

        VertexLoad& insertVertices(const std::vector<std::pair<V_ind, VertParam>>& vertices);
        VertexLoad& insertVertices(const std::vector<V_ind>& vertices);
        VertexLoad& setForceField(std::function<bool(double *x/*[3]*/, double *F/*[3]*/)> force, std::function<bool(double *x/*[3]*/, double *J/*[3][3]*/)> dforce);

        VertexLoad(){ type = "VertexLoad<void>"; }
        VertexLoad(const VertexLoad& a) = default;
        int operator()(Object3D &obj) override;
        std::unique_ptr<ForceBase> clone() override{ return std::make_unique<VertexLoad>(*this); }
        int fill_matrix(Object3D &obj, SparseMatrix::IndexFrom &m) override;
    };
}


#endif //AVSIM_VERTEXLOAD_H
