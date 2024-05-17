//
// Created by alex on 11.10.2023.
//

#ifndef AORTIC_VALVE_HYPER_ELASTIC_MODEL_H
#define AORTIC_VALVE_HYPER_ELASTIC_MODEL_H

#include "HyperElastic.h"
namespace World3d{
    class HyperElasticModel: public HyperElasticForceBase{
    public:
        struct ParamSource{
            enum Type{
                VAL_T,
                TAG_T,
                TAG_VAL_T,  //prefer tag but if not found get value
            };
            Type type = VAL_T;
            std::string tag_name;
            DReal val = 0;
            ParamSource() = default;
            ParamSource(Type type, std::string tag_name, double val): type{type}, tag_name{tag_name}, val{val} {}
            std::unique_ptr<HH::DataRequier> makeDataRequier() const {
                switch (type){
                    case VAL_T: return std::make_unique<ConstDataRequier>(val);
                    case TAG_T: return std::make_unique<ConstDataTagRequier<DReal>>(tag_name);
                    default: return std::make_unique<HH::ConstDataTagOrValRequier<DReal>>(tag_name, val);
                }
            } 
        };
        struct InvariantDescr{
            enum Type{
                I1_T,
                I2_T,
                J_T,
                If_T,
                Ifs_T,
            };
            Type type = I1_T;
            std::string f1, f2;
            InvariantDescr() = default;
            InvariantDescr(Type type, const std::string& f1 = "", const std::string& f2 = ""): type{type}, f1{f1}, f2{f2} {}
            HH::ArgsDef createInvariant(){
                HH::ArgsDef res;
                switch (type){
                    case I1_T: get<2>(res).push_back(HH::getDefaultInvariants().at("I1")); return res;
                    case I2_T: get<2>(res).push_back(HH::getDefaultInvariants().at("I3")); return res;
                    case J_T: get<2>(res).push_back(HH::getDefaultInvariants().at("J")); return res;
                    case If_T: return DefaultHyperElasticForce::makeOrthotropicInvariant("I4[" + f1 + "]", f1);
                    default: return HH::makeAnisotropicInvariant("I8[" + f1 + "," + f2 + "]", f1, f2);
                }
            } 
        };

        ParamSource thickness_src;
        std::vector<std::pair<SX, ParamSource>> potential_params;
        std::vector<std::pair<SX, InvariantDescr>> potential_invariants;
        SX potential;
        std::string name;

        HyperElasticModel(
            ParamSource h, std::vector<std::pair<SX, ParamSource>> prms, std::vector<std::pair<SX, InvariantDescr>> invs, SX potential,
            std::string gen_dir, bool regenerate = true, std::string model_name = "HyperElastic"):
            thickness_src{h}, potential_params{prms}, potential_invariants{invs}, potential{potential}, name{model_name}{
            type = "HyperElasticModel";
            auto lv = DefaultHyperElasticForce::getDefaultLocalVars();
            lv.push_back(HH::LocalVar("H"));    //TODO: this may be wrong
            SX H  = lv["H"]; //SX::sym("H");
            SX U = lv["Ap"] * H * potential;
            auto iv = DefaultHyperElasticForce::getDefaultAvailableInputVars();
            iv.emplace_back("H", H, thickness_src.makeDataRequier());
            for (auto& p: potential_params)
                iv.emplace_back("&" + p.first.name(), p.first, p.second.makeDataRequier());
            
            HH::ArgsDef args(std::move(iv), std::move(lv), HH::VarContainer<HH::Invariant>());    
            for (auto& p: potential_invariants){
                DefaultHyperElasticForce::ArgsDef I = p.second.createInvariant();
                U = SX::substitute(U, p.first, get<2>(I)[0]);
                for (auto& i: get<0>(I)) get<0>(args).push_back(std::move(i));
                for (auto& i: get<1>(I)) get<1>(args).push_back(std::move(i));
                for (auto& i: get<2>(I)) get<2>(args).push_back(std::move(i));
            }
            
            f = DefaultHyperElasticForce(name, gen_dir, std::move(args), U, regenerate);
        }
        std::unique_ptr<HyperElasticForceBase> copy() override{
            return std::make_unique<HyperElasticModel>(thickness_src, potential_params, potential_invariants, potential, f.gen_dir, false, name);
        }
    };
}

#endif  //HYPER_ELASTIC_MODEL