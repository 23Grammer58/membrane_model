//
// Created by alex on 01.02.2021.
//

#include "DataDriven.h"

std::array<double, 3> World3d::_SimpleChoose::choose(const std::array<double, 3>& ksi, ForceBase* dd, Object3D* obj, F_ind f){
    auto de = static_cast<DataDrivenElastic*>(dd);
    return de->interp.get()->operator()(ksi);
}

std::array<double, 3> World3d::_RegionChoose::choose(const std::array<double, 3>& ksi, ForceBase* dd, Object3D* obj, F_ind f){
    auto de = static_cast<DataDrivenElasticRegion*>(dd);
    int id = de->interpolantChooser(de, obj, f);
    return de->interps[id].get()->operator()(ksi);
}
