//
// Created by alex on 01.12.2023.
//

#ifndef AORTIC_VALVE_MATERIAL_THICKNESS_HELPER_H
#define AORTIC_VALVE_MATERIAL_THICKNESS_HELPER_H

#include "HyperElasticResponse.h"
#include "HyperElastic.h"

namespace World3d{

struct ThicknessFromDataRequier: public MaterialThickness{
    using DataRequier = HyperElasticHelpers::DataRequier;

    std::unique_ptr<DataRequier> m_ht;
    Object3D* _obj = nullptr;

    ThicknessFromDataRequier() = default;
    ThicknessFromDataRequier(std::unique_ptr<DataRequier> ht): m_ht{std::move(ht)} {}
    ThicknessFromDataRequier& setSource(std::unique_ptr<DataRequier> ht) { return m_ht = std::move(ht), *this; }

    DReal operator()(F_ind f) override { 
        const DReal* res = 0;
        auto v = vert_around(_obj->m_mesh, f);
        std::vector<V_ind> vv(v.begin(), v.end());
        m_ht->saveData(&res, f, vv);
        return *res;
    }
    void registerObj(Object3D* obj) { _obj = obj, m_ht->registerObj(obj); }
    std::unique_ptr<MaterialThickness> clone() const { auto res = std::make_unique<ThicknessFromDataRequier>(); res->setSource(m_ht->copy()); return res; }
};

}

#endif //AORTIC_VALVE_MATERIAL_THICKNESS_HELPER_H