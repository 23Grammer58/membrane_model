//
// Created by alex on 01.02.2021.
//

#ifndef AORTIC_VALVE_DATADRIVEN_H
#define AORTIC_VALVE_DATADRIVEN_H

#include "DataDriven/InterpolatorBase.h"
#include "DataDriven/KNearestSearcher.h"
#include "DataDriven/ResponseStatistics.h"
#include "DataDriven/ConjugatePairSearcher.h"
#include "DataDriven/LinearZeroSeacher.h"
#include "../Object3D.h"
#include "ForcesCommon.h"

namespace World3d {
    //DerivsFunctor - std::function(std::array<double, 3>(const std::array<double, 3>& ksi, DataDrivenElasticBase* dd, Object3D* obj, F_IND f))
    template<class DerivsFunctor>
    class DataDrivenElasticBase: public ForceBase{
    public:
        Object3D* _obj = nullptr;
        ConstProperty<F_ind, std::array<DReal, 2 * 3>> _L;
        ConstProperty<F_ind, DReal> _Ap;

        DataDrivenElasticBase() { type = "DataDrivenElasticBase"; };
        DataDrivenElasticBase(const DataDrivenElasticBase&) = default;
        DataDrivenElasticBase(DataDrivenElasticBase&&) = default;
        DataDrivenElasticBase& operator=(const DataDrivenElasticBase&) = default;
        DataDrivenElasticBase& operator=(DataDrivenElasticBase&&) = default;

        void registerObj(Object3D *obj) {
            _obj = obj;
            _L = set_L(*obj);
            _Ap = set_Ap(*obj);
        }
        int operator()(Object3D &obj) override{
            if (_obj != &obj) registerObj(&obj);
            for (auto f: obj.m_mesh.faces()) {
                auto vv = vert_around(obj.m_mesh, f);
                auto Ldat = _L[f];
                auto Ap = _Ap[f];

                auto L = Eigen::Map<Eigen::Matrix<World3d::DReal, 2, 3>>(Ldat.data());
                Eigen::Matrix<World3d::DReal, 3, 3> Q;
                Q <<    obj.m_x[vv[0]][0], obj.m_x[vv[1]][0], obj.m_x[vv[2]][0],
                        obj.m_x[vv[0]][1], obj.m_x[vv[1]][1], obj.m_x[vv[2]][1],
                        obj.m_x[vv[0]][2], obj.m_x[vv[1]][2], obj.m_x[vv[2]][2];
                Eigen::Matrix<World3d::DReal, 2, 3> R = L * Q.transpose();
                World3d::DReal C[] = {
                        R.row(0).dot(R.row(0).transpose()),
                        R.row(0).dot(R.row(1).transpose()),
                        R.row(1).dot(R.row(1).transpose())
                };
                std::array<double, 3> ksi;
                ksi[0] = 0.5 * std::log(C[0]);                      //ln (F11)
                ksi[1] = 0.5 * std::log(C[2] - C[1] * C[1] / C[0]); //ln (F22)
                ksi[2] = C[1] / C[0];
                std::array<double, 3> derivs = DerivsFunctor::choose(ksi, this, &obj, f);
                auto response = Eigen::Map<Eigen::Matrix<World3d::DReal, 3, 1>>(derivs.data());
                double k1 = C[0] / (2 * (C[0] * C[2] - C[1] * C[1])), k2 = C[1] / C[0], k3 = 1/C[0];
                double k4 = k3 / 2, k5 = (-2*k2), k6 = k2*k2;
                Eigen::Matrix<World3d::DReal, 3, 3> Ksi;
                std::array<Eigen::Matrix<World3d::DReal, 3, 1>, 3> dCdQ;

                for (int k = 0; k < 3; ++k) {
                    auto v = vv[k];
                    dCdQ[0] = 2 * R.row(0).transpose() * L(0, k);
                    dCdQ[1] = R.row(0).transpose() * L(1, k) + R.row(1).transpose() * L(0, k);
                    dCdQ[2] = 2 * R.row(1).transpose() * L(1, k);
                    Ksi.col(0) = k4 * dCdQ[0];
                    Ksi.col(1) = k1 * (dCdQ[2] + k5 * dCdQ[1] + k6*dCdQ[0]);
                    Ksi.col(2) = k3 * (dCdQ[1] - k2 * dCdQ[0]);
                    Eigen::Matrix<World3d::DReal, 3, 1> F = -Ap * Ksi * response;
                    //TODO: добавлено ограничение движения
                    if (obj.is_movable(v))
                        obj.m_F[v] += obj.withBCmask(v, Vector(F(0, 0), F(1, 0), F(2, 0)));
                }
            }

            return 0;
        }
        std::unique_ptr<ForceBase> clone() override{
            return std::make_unique<DataDrivenElasticBase>(*this);
        }
    };

    struct _SimpleChoose{
        static std::array<double, 3> choose(const std::array<double, 3>& ksi, ForceBase* dd, Object3D* obj, F_ind f);
    };

    class DataDrivenElastic: public DataDrivenElasticBase<_SimpleChoose>{
    public:
        using ELEM = ResponseStatistics::Elem;
        using DATA = ResponseStatistics::Rdat;
        using Interp = InterpolantBase3D;

        std::shared_ptr<DATA> data;
        std::shared_ptr<Interp> interp;

        DataDrivenElastic() { type = "DataDrivenElasticSimple"; }
        DataDrivenElastic(const DataDrivenElastic&) = default;
        DataDrivenElastic(DataDrivenElastic&&) = default;
        DataDrivenElastic& operator=(const DataDrivenElastic&) = default;
        DataDrivenElastic& operator=(DataDrivenElastic&&) = default;

        DataDrivenElastic(DATA data): DataDrivenElastic(std::move(data), KNearestSearcher3D(1)) {}

        DataDrivenElastic(DATA data, std::shared_ptr<Interp> interp):
            data{std::make_shared<DATA>(std::move(data))}, interp{std::move(interp)} {
            type = "DataDrivenElasticSimple";
        }
        template<typename INTERP>
        DataDrivenElastic(DATA data, INTERP interp): data{std::make_shared<DATA>(std::move(data))}, interp{ interp.clone() } {
            type = "DataDrivenElasticSimple";
        }
        DataDrivenElastic& setData(DATA dat) { return data = std::make_shared<DATA>(std::move(dat)), *this; }
        DataDrivenElastic& setInterpolant(std::shared_ptr<Interp> interpolant) { return interp = std::move(interp), *this; }
        template<typename INTERP>
        DataDrivenElastic& setInterpolant(INTERP interpolant) { return interp.reset(interpolant.clone()), *this; }
        DataDrivenElastic& init() { return interp.get()->init(data.get()), *this; }

        std::unique_ptr<ForceBase> clone() override{
            return std::make_unique<DataDrivenElastic>(*this);
        }
    };

    struct _RegionChoose{
        static std::array<double, 3> choose(const std::array<double, 3>& ksi, ForceBase* dd, Object3D* obj, F_ind f);
    };

    class DataDrivenElasticRegion: public DataDrivenElasticBase<_RegionChoose>{
    public:
        using ELEM = ResponseStatistics::Elem;
        using DATA = ResponseStatistics;
        using Interp = InterpolantBase3D;
        using Initializer = std::function<void(std::shared_ptr<DATA>& data, std::vector<std::shared_ptr<Interp>>& interps)>;
        using InterpolantChooser = std::function<int (DataDrivenElasticRegion* dd, Object3D* obj, F_ind f)>;

        std::shared_ptr<DATA> data;
        std::vector<std::shared_ptr<Interp>> interps;
        Initializer initializer =
                [](std::shared_ptr<DATA>& data, std::vector<std::shared_ptr<Interp>>& interps) -> void {
                    interps[0].get()->init(&(data.get()->stat[0]));
        };
        InterpolantChooser interpolantChooser =
                [](DataDrivenElasticRegion* dd, Object3D* obj, F_ind f){ return 0; };

        DataDrivenElasticRegion() { type = "DataDrivenElasticRegion"; }
        DataDrivenElasticRegion(const DataDrivenElasticRegion&) = default;
        DataDrivenElasticRegion(DataDrivenElasticRegion&&) = default;
        DataDrivenElasticRegion& operator=(const DataDrivenElasticRegion&) = default;
        DataDrivenElasticRegion& operator=(DataDrivenElasticRegion&&) = default;

        DataDrivenElasticRegion(DATA dat, std::vector<std::shared_ptr<Interp>>interp):
                data{std::make_shared<DATA>(std::move(dat))}, interps{std::move(interp)} {
            type = "DataDrivenElasticRegion";
        }
        DataDrivenElasticRegion& setData(DATA dat) { return data = std::make_shared<DATA>(std::move(dat)), *this; }
        DataDrivenElasticRegion& setInterpolant(std::vector<std::shared_ptr<Interp>> interp) { return interps = std::move(interp), *this; }
        DataDrivenElasticRegion& setInitializer(Initializer init) { return initializer = std::move(init), *this;}
        DataDrivenElasticRegion& setInterpolantChooser(InterpolantChooser chooser) {
            return interpolantChooser = std::move(chooser), *this;
        }
        DataDrivenElasticRegion& init() { return initializer(data, interps), *this; }

        std::unique_ptr<ForceBase> clone() override{
            return std::make_unique<DataDrivenElasticRegion>(*this);
        }
    };


};

static void testExample(){
    using DAT = ResponseStatistics::Elem;
    using DATTraits = ResponseStatistics::Elem::Traits;
    using PNT = DATTraits::PNT;
    using VAL = DATTraits::VAL;

    ResponseStatistics::Rdat dat;
    auto f = makeInterpolantFactory3D(KNearestSearcher3D(1));
    World3d::DataDrivenElastic force(dat, f);

}

#endif //AORTIC_VALVE_DATADRIVEN_H
