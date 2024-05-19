//
// Created by alex on 01.12.2023.
//

#ifndef AORTIC_VALVE_HYPER_ELASTIC_RESPONSE_H
#define AORTIC_VALVE_HYPER_ELASTIC_RESPONSE_H

#include <array>
#include <functional>
#include <memory>
#include "AVSim/Core/Object3D.h"

namespace World3d{

struct HyperElasticMaterial{
    virtual void registerObj(Object3D* obj) {}
    /// Return Second Piola-Kirchgoff tensor S_2d by Cauchy-Green deformation tensor C_2d
    /// @param E2d_tensor is (E00, E11, E01)
    /// @return S2d_tensor which is (S00, S11, S01)
    virtual std::array<double, 3> PK2_tensor(std::array<double, 3> E2d_tensor, F_ind f) = 0;
    /// Return S_2d and d S_2d / d E_2d tensor, required for newton-like solver methods
    /// @param E2d_tensor is (E00, E11, E01)
    /// @return (S00, S11, S01) and (S00'00, S11'11, S01'01, S00'11, S00'01, S11'01), Sij'kl = d S_ij / d E_kl
    virtual std::pair<std::array<double, 3>, std::array<double, 6>> PK2_dPK2_tensor(std::array<double, 3> E2d_tensor, F_ind f) = 0;
    virtual std::unique_ptr<HyperElasticMaterial> clone() const = 0;
};

struct MaterialThickness{
    virtual DReal operator()(F_ind f) = 0;
    virtual void registerObj(Object3D* obj) {}
    virtual std::unique_ptr<MaterialThickness> clone() const = 0;
};

struct ConstThickness: public MaterialThickness{
    double H = 0.5;
    explicit ConstThickness(double H = 0.5): H{H} {}
    DReal operator()(F_ind f) override { return H; }
    std::unique_ptr<MaterialThickness> clone() const override { return std::make_unique<ConstThickness>(*this); }
};
struct TagThickness: public MaterialThickness{
    ConstProperty<F_ind, DReal> H;
    std::string tag_name = "f:thickness"; 

    explicit TagThickness(std::string tag_name = "f:thickness"): tag_name{tag_name} {}
    void registerObj(Object3D* obj){ 
        H = obj->m_mesh.property_map<F_ind, DReal>(tag_name);
        if (!H.second) throw std::runtime_error("Tag \"" +  tag_name + "\" is not exists on the object \"" + obj->name + "\"");
    }
    DReal operator()(F_ind f) override { return H[f]; }
    std::unique_ptr<MaterialThickness> clone() const override { return std::make_unique<TagThickness>(*this); }
};
struct TagOrValThickness: public MaterialThickness{
    ConstProperty<F_ind, DReal> H;
    DReal Hdef = 0.5;
    std::string tag_name = "f:thickness"; 

    TagOrValThickness(double Hdef = 0.5, std::string tag_name = "f:thickness"): Hdef{Hdef}, tag_name{tag_name} {}
    void registerObj(Object3D* obj) override { H = obj->m_mesh.property_map<F_ind, DReal>(tag_name); }
    DReal operator()(F_ind f) override { return H.second ? H[f] : Hdef; }
    std::unique_ptr<MaterialThickness> clone() const override { return std::make_unique<TagOrValThickness>(*this); }
};

struct HEForce: public ForceBase{
    std::unique_ptr<HyperElasticMaterial> m_mat;
    std::unique_ptr<MaterialThickness> m_ht;
    Object3D* _obj = nullptr;

    HEForce(){ type = "HEForce"; }
    template<typename HEMatT>
    typename std::enable_if<std::is_base_of<HyperElasticMaterial, HEMatT>::value, HEForce&>::type 
        setMaterial(HEMatT f) { return m_mat = std::make_unique<HEMatT>(std::move(f)), *this; }
    HEForce& setMaterial(const HyperElasticMaterial& f)  { return m_mat = f.clone(), *this; }
    template<typename MatHT>
    typename std::enable_if<std::is_base_of<MaterialThickness, MatHT>::value, HEForce&>::type 
        setThickness(MatHT f) { return m_ht = std::make_unique<MatHT>(std::move(f)), *this; }
    HEForce& setThickness(const MaterialThickness& f)  { return m_ht = f.clone(), *this; }    

    HEForce(const HEForce &f) : m_mat{f.m_mat->clone()}, m_ht{f.m_ht->clone()} { type = f.type; }
    HEForce& operator=(const HEForce &f) = default;
    HEForce(HEForce &&f) = default;
    HEForce &operator=(HEForce &&f) = default;

    void registerObj(Object3D* obj) { _obj = obj, m_mat->registerObj(obj), m_ht->registerObj(obj); }
    std::unique_ptr<ForceBase> clone() override { return std::make_unique<HEForce>(*this); }
    int operator()(Object3D &obj) override;
    int fill_matrix(Object3D& obj, SparseMatrix::IndexFrom& m) override;
};

}

#endif //AORTIC_VALVE_HYPER_ELASTIC_RESPONSE_H