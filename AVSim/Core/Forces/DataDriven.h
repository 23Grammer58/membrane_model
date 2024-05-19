//
// Created by alex on 01.02.2021.
//

#ifndef AORTIC_VALVE_DATADRIVEN_H
#define AORTIC_VALVE_DATADRIVEN_H

#include "HyperElasticResponse.h"
#include "DataDriven/InterpolatorBase.h"
#include "DataDriven/KNearestInterpolant.h"
#include "DataDriven/MaterialResponse.h"
#include "DataDriven/ConjugatePairSearcher.h"
#include "DataDriven/LinearElasticInterpolant.h"

namespace World3d{

struct DataDrivenMaterial: public HyperElasticMaterial{
    using Interp = ResponseInterpolant;
    using InterpolantChooser = std::function<int (DataDrivenMaterial* dd, Object3D* obj, F_ind f)>;

    std::vector<std::unique_ptr<Interp>> interps;
    InterpolantChooser interpolantChooser = [](DataDrivenMaterial* dd, Object3D* obj, F_ind f){ return 0; };
    Object3D* _obj = nullptr;

    DataDrivenMaterial() = default;
    DataDrivenMaterial(const DataDrivenMaterial& a): interpolantChooser{a.interpolantChooser}{
        interps.resize(a.interps.size());
        for (std::size_t i = 0; i < a.interps.size(); ++i)
            interps[i] = a.interps[i]->clone();
    }
    
    DataDrivenMaterial& setInterpolant(std::unique_ptr<Interp> interp) { interps.resize(1); return interps[0] = std::move(interp), *this; }
    DataDrivenMaterial& setInterpolant(std::vector<std::unique_ptr<Interp>> interp) { return interps = std::move(interp), *this; }
    template <typename INTERP>
    typename std::enable_if<std::is_base_of<Interp, INTERP>::value, DataDrivenMaterial&>::type
        setInterpolant(INTERP interp) { interps.resize(1); return interps[0] = std::make_unique<INTERP>(std::move(interp)), *this; }
    DataDrivenMaterial& setInterpolantChooser(InterpolantChooser chooser) { return interpolantChooser = std::move(chooser), *this; }
    void registerObj(Object3D* obj) override { _obj = obj; }
    std::array<double, 3> PK2_tensor(std::array<double, 3> E2d, F_ind f) override;
    std::pair<std::array<double, 3>, std::array<double, 6>> PK2_dPK2_tensor(std::array<double, 3> E2d, F_ind f) override;
    std::unique_ptr<HyperElasticMaterial> clone() const override { return std::make_unique<DataDrivenMaterial>(*this); }
};

}

#endif //AORTIC_VALVE_DATADRIVEN_H