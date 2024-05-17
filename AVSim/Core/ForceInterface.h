//
// Created by Liogky Alexey on 19.05.2022.
//

#ifndef AVSYM_FORCEINTERFACE_H
#define AVSYM_FORCEINTERFACE_H
#include "AVMesh.h"
#include "AVSim/Solvers/LinearSolver.h"

namespace World3d {
    class Object3D;
    class ObjectID;

    struct ForceBase {
        struct LocMatrix {
            std::vector<DReal> dat;
            int N, M;

            LocMatrix(int N = 0, int M = 0) : N{N}, M{M}, dat(N * M) {}
            void resize(int N, int M);
            void setZero() { std::fill(dat.begin(), dat.end(), 0.0); }
            void setEye() {
                for (int i = 0; i < N; ++i)
                    for (int j = 0; j < M; ++j)
                        dat[j + M * i] = (i = j);
            }
            double *operator[](int i) { return dat.data() + M * i; }
            void ApplyDirichlet(int i) {
                for (int j = 0; j < N; ++j)
                    dat[M * j + i] = 0;
                for (int j = 0; j < M; ++j)
                    dat[M * i + j] = 0;
                dat[M * i + i] = 1;
            }

            void ApplyNonExists(int i) {
                for (int j = 0; j < N; ++j)
                    dat[M * j + i] = 0;
                for (int j = 0; j < M; ++j)
                    dat[M * i + j] = 0;
            }

            struct DoubleIndex {
                int N1, N2;
                DoubleIndex(int N1, int N2) : N1{N1}, N2{N2} {}
                int operator()(int i, int j) { return i * N2 + j; }
            };
        };

        std::string type;

        virtual int operator()(Object3D &obj) = 0;
        virtual int element_matrix(Object3D &obj, LocMatrix &matr, F_ind element) {
            throw std::runtime_error("This method is not implemented");
            return -1;
        }
        virtual int fill_matrix(Object3D &obj, SparseMatrix::IndexFrom &m);
        virtual std::unique_ptr<ForceBase> clone() = 0;
    };

    struct Force {
        std::unique_ptr<ForceBase> m_invoker;
        Force() : m_invoker{} {}

        template<typename ForceT>
        Force(const ForceT& f, typename std::enable_if<std::is_base_of<ForceBase, ForceT>::value>::type* = 0): m_invoker{new ForceT(f)} {}
        template<typename ForceT>
        Force(ForceT&& f, typename std::enable_if<std::is_base_of<ForceBase, ForceT>::value>::type* = 0): m_invoker{new ForceT(std::move(f))} {}
        Force(const Force &f) : m_invoker{f.m_invoker->clone()} {}

        Force &operator=(const Force &f) {
            m_invoker = f.m_invoker->clone();
            return *this;
        }

        Force(Force &&f) = default;

        Force &operator=(Force &&f) = default;

        int operator()(Object3D &obj) { return m_invoker->operator()(obj); }

        ForceBase *operator->() { return m_invoker.get(); }

        template<typename ForceT>
        ForceT *target() { return static_cast<ForceT *>(m_invoker.get()); }
    };

    class World;
    struct WorldForceBase {
        std::string type;

        virtual bool addObj(Object3D& obj, const ObjectID& id, bool is_dynamic) { return false; }
        virtual void removeObj(const ObjectID& id) {}
        virtual bool registerWorld(World& w) = 0;
        virtual int operator()() = 0;
        virtual int fill_matrix(SparseMatrix *sm) = 0;
        virtual std::unique_ptr<WorldForceBase> clone() = 0;
    };
    struct WorldForce {
        std::unique_ptr<WorldForceBase> m_invoker;
        WorldForce() : m_invoker{} {}

        template<typename ForceT>
        WorldForce(const ForceT& f, typename std::enable_if<std::is_base_of<WorldForceBase, ForceT>::value>::type* = 0): m_invoker{new ForceT(f)} {}
        template<typename ForceT>
        WorldForce(ForceT&& f, typename std::enable_if<std::is_base_of<WorldForceBase, ForceT>::value>::type* = 0): m_invoker{new ForceT(std::move(f))} {}

        WorldForce(const WorldForce &f) : m_invoker{f.m_invoker->clone()} {}

        WorldForce &operator=(const WorldForce &f) {
            m_invoker = f.m_invoker->clone();
            return *this;
        }

        WorldForce(WorldForce &&f) = default;
        WorldForce &operator=(WorldForce &&f) = default;
        int operator()() { return m_invoker->operator()(); }
        WorldForceBase *operator->() { return m_invoker.get(); }
        template<typename ForceT>
        ForceT *target() { return static_cast<ForceT *>(m_invoker.get()); }
    };
};


#endif //AVSYM_FORCEINTERFACE_H
