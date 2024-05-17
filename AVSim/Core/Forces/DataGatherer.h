//
// Created by alex on 02.02.2021.
//

#ifndef AORTIC_VALVE_DATAGATHERER_H
#define AORTIC_VALVE_DATAGATHERER_H

#include "../TriangularMeshHelpers.h"
#include "../World.h"
#include "HyperElastic.h"
#include "DataDriven/ResponseStatistics.h"
#include "DataDriven.h"
#include "NeoGook.h"
#include "Pressure.h"

namespace World3d {
    class HEDataGatherer {
    public:
        using DefaultGenInputVars = std::vector<HH::DefaultGenInputVar>;
        using LocalVars = HH::VarContainer<HH::LocalVar>;
        using Invariants = HH::VarContainer<HH::Invariant>;
        using RS = ResponseStatistics;

        std::string gen_dir, name;
        HH::SX U;
        Object3D *_obj = nullptr;
        HH::GenFunction response;

        HEDataGatherer(const HEDataGatherer& b) = default;
        HEDataGatherer(HEDataGatherer&& b) noexcept = default;
        HEDataGatherer& operator=(const HEDataGatherer& a) = default;
        HEDataGatherer& operator=(HEDataGatherer&& a) noexcept = default;
        HEDataGatherer(std::string name, std::string save_dir,
                       DefaultGenInputVars input_vars, LocalVars& loc_vars, Invariants& invariants,
                       HH::SX U, bool regenerate = true);
        explicit HEDataGatherer(const DefaultHyperElasticForce& f, bool regenerate = true);
        void Init(DefaultGenInputVars input_vars, LocalVars& loc_vars, Invariants& invariants, bool regenerate);

        struct CompDat{
            F_ind f;
            DReal* ksi;
            DReal* response;
            DReal* Tee;
            DReal* T;
        };

        RS& gatherData(Object3D& obj, RS& rs);
        RS gatherData(Object3D& obj);
        void gatherData(Object3D& obj, std::function<void(const CompDat& dat)> my_filler);
        std::array<double, 3> responseAt(Object3D& obj, F_ind f);

    private:
        std::vector<V_ind> _vdat;
        void generate() { response.generate(name, gen_dir); }
        void link_function() { response.link_function(name, gen_dir); }
        void registerObj(Object3D *obj) { response.registerObj(_obj = obj); }
        void compute_responses(DefaultGenInputVars &generatedInput, LocalVars &loc_vars, Invariants &invariants);
        HH::SX compute_response_expression(const LocalVars &loc_vars, const Invariants &invariants) const;
        static LocalVars getDefaultLocalVars();
    };

    struct _AnalyticChoose{
        static std::array<double, 3> choose(const std::array<double, 3>& ksi, ForceBase* dd, Object3D* obj, F_ind f);
    };

    class DataDrivenAnalytic: public DataDrivenElasticBase<_AnalyticChoose>{
    public:
        HEDataGatherer dat;

        DataDrivenAnalytic(const DataDrivenAnalytic&) = default;
        DataDrivenAnalytic(DataDrivenAnalytic&&) = default;
        DataDrivenAnalytic& operator=(const DataDrivenAnalytic&) = default;
        DataDrivenAnalytic& operator=(DataDrivenAnalytic&&) = default;
        DataDrivenAnalytic(const DefaultHyperElasticForce& f, bool regenerate = true): dat(f, regenerate)
            { type = "DataDrivenAnalytic"; }

        std::unique_ptr<ForceBase> clone() override{
            return std::make_unique<DataDrivenAnalytic>(*this);
        }
    };

    class StaticDefiniteGather{
    public:
        using RS = ResponseStatistics;
        struct Traits{
            //forces parameters
            double m_P = 5.0e3;             //pressure
            double m_H = 5e-2, m_mu = 1e11; //thickness and neogook stiffness term
            //solver parameters
            double m_err = 1e-5;
            int m_freq = 1000, m_maxits = 100000;
            enum SolveMethod{
                SIMPLE_RELAXATON,
                SIMPLE_NEWTON,
            };
            SolveMethod m_method = SIMPLE_RELAXATON;
            //simple relaxation solver parameter
            double m_delta = 1.5e-11;
            //newton solver parameters
            NewtonSolverStepAlgo m_nssa;

            Traits() = default;
            Traits(double P, double H, double mu = 1e11, double delta = 1.5e-11, double err = 1e-5, int freq = 1000, int maxits = 100000):
                m_P{P}, m_H{H}, m_mu{mu}, m_delta{delta}, m_err{err}, m_freq{freq}, m_maxits{maxits} {}
            Traits& setP(double P) { return m_P = P, *this; }
            Traits& setH(double H) { return m_H = H, *this; }
            Traits& setMu(double mu) { return m_mu = mu, *this; }
            Traits& setDelta(double delta) { return m_delta = delta, *this; }
            Traits& setErr(double err) { return m_err = err, *this; }
            Traits& setFreq(int freq) { return m_freq = freq, *this; }
            Traits& setMaxIters(int maxits) { return m_maxits = maxits, *this; }
            Traits& setSolveMethod(SolveMethod method) { return m_method = method, *this; }
            Traits& setLinearSolver(const LinearSolver& ls){ return m_nssa._solver.solver = ls, *this; }
            Traits& setNewtonAlgo(const std::function<int(NewtonSolverStepAlgo& ns, World* w)>& algo) { return m_nssa.newton_algo = algo, *this; }
        };
        Traits m_t;
        Object3D* m_obj = nullptr;
        World m_w;

        StaticDefiniteGather() = default;
        StaticDefiniteGather(Object3D& obj, Traits t): m_t{t}, m_obj{&obj} {}
        StaticDefiniteGather& setTraits(Traits t) { return m_t = t, *this; }
        StaticDefiniteGather& setObject(Object3D& obj) { return m_obj = &obj, *this; }
        StaticDefiniteGather& computeTensionedState();
        StaticDefiniteGather& gatherResponseData(RS& rs);
        template<class VECTOR>
        StaticDefiniteGather& gatherTension(VECTOR& vec){
            if (m_w.objs().empty() || !m_obj) return *this;
            auto& obj1 = m_w.objs().begin()->second.first;
            auto& obj = *m_obj;
            auto D1_ = set_D_vecs(obj1);
            auto Ap1_ = set_Ap(obj1);
            auto S1_ = set_S(obj1);
            auto x1_ = obj1.m_x;
            auto p1_ = obj1.m_x0;
            auto h1 = obj1.m_mesh.property_map<F_ind, DReal>("f:Thickness").first;

            auto L_ = set_L(obj);
//            auto Ap_ = set_Ap(obj);
//            auto S_ = set_S(obj);
            double mu = m_t.m_mu;

            for (auto f: obj1.m_mesh.faces()){
                auto v = vert_around(obj1.m_mesh, f);
                auto Ap1 = Ap1_[f];
                auto Aq1 = sqrt(S1_[f].squared_length());
                Eigen::Matrix<DReal, 3, 3>  Q1, P1, D1;
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
//                auto Ap = Ap_[f], Aq = sqrt(S_[f].squared_length());
//                auto J = Aq / Ap;
//                auto H = m_t.m_H;

                auto T = h1[f] * mu / J1 * (dF1 * dF1.transpose() + dF1 + dF1.transpose() + Eigen::Matrix<DReal, 3, 3>::Identity() * (J1 - 1) * (J1 + 1) / (J1 * J1));
                auto L = Eigen::Map<Eigen::Matrix<DReal, 2, 3>>(L_[f].data());
                auto F = Q1 * L.transpose();

                auto C2d = F.transpose() * F;
                auto S = F * C2d.inverse();
                Eigen::Matrix<DReal, 3, 2> R;
                R.col(0) = F.col(0) / F.col(0).norm();
                R.col(1) = F.col(1) - F.col(1).dot(R.col(0))*R.col(0); R.col(1) /= R.col(1).norm();
                auto Tee = /*S.transpose() * T * S;*/R.transpose() * T * R;
                vec[f.idx() + 0*obj1.m_mesh.num_faces()] = Tee(0, 0);
                vec[f.idx() + 1*obj1.m_mesh.num_faces()] = Tee(1, 1);
                vec[f.idx() + 2*obj1.m_mesh.num_faces()] = Tee(1, 0);
            }

            return *this;
        }
        struct ObjsData{
            int nfaces;
            F_ind f;
            Eigen::Matrix<DReal, 3, 3>  Q, P, D, Q1, P1, D1;
            double Ap, Aq, Ap1, Aq1;
            Eigen::Matrix<DReal, 3, 3> T;
            Eigen::Matrix<DReal, 2, 3> L;
            Eigen::Matrix<DReal, 3, 2> F2d;
            Eigen::Matrix<DReal, 2, 2> C2d;
            Eigen::Matrix<DReal, 2, 2> Tee;
        };

        StaticDefiniteGather& gatherMyData(std::function<void(const ObjsData& dat)> my_filler);
    };
};


#endif //AORTIC_VALVE_DATAGATHERER_H
