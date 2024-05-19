#include "model.h"
#include <filesystem>
#include "AVSim/Core/NSWorldWrapper.h"
#include "AVSim/Solvers/NonLinearSolverKinsol.h"
#include "AVSim/Core/NonLinearSolverCustom.h"
#include <bitset>

#ifdef WITH_DATADRIVEN
#include "AVSim/Core/Forces/DataDriven.h"
#endif

#if __cplusplus >= 201703L
using namespace std::filesystem;
#else
using namespace std::experimental::filesystem;
#endif

static int start_debug_gui(int argc, char** argv){
#ifdef USE_MAGNUM_GUI
    GuiApplication app({argc, argv});
    return app.exec();
#else
    return 0;
#endif
}

std::string Model::getObjectTraitTypeName(ObjectTraitType t){
    static const char* res[] = {
        "INITIAL_CONFIGURATION",
        "START_CONFIGURATION",
        "BOUNDARY",
        "PRESSURE",
        "LAGRANGE_PRESSURE",
        "BODY_FORCE",
        "EDGE_FORCE",
        "CLAMP_DIRECTION",
        "THICKNESS",
    };
    return res[static_cast<unsigned>(t)];
}
const std::map<std::string, Model::ObjectTraitType>& Model::getObjectTraitTypeToNameMap(){
    static std::map<std::string, ObjectTraitType> rmap{
        {"INITIAL_CONFIGURATION", INIT_X_VECTOR_DTAG},
        {"START_CONFIGURATION", START_X_VECTOR_DTAG},
        {"BOUNDARY", BOUNDARY_ITAG},
        {"PRESSURE", PRESSURE_DTAG},
        {"LAGRANGE_PRESSURE", LAGRANGE_PRESSURE_DTAG},
        {"BODY_FORCE", BODY_FORCE_VECTOR_DTAG},
        {"EDGE_FORCE", EDGE_FORCE_VECTOR_DTAG},
        {"CLAMP_DIRECTION", CLAMPED_EDGE_VECTOR_DTAG},
        {"THICKNESS", THICKNESS_DTAG},
    };
    return rmap;
}
Model::ObjectTraitType Model::getObjectTraitTypeByName(const std::string& s){
    const auto& rmap = getObjectTraitTypeToNameMap();
    auto it = rmap.find(s);
    if (it != rmap.end()) return it->second;
    else throw std::runtime_error("Can't find object trait with name = \"" + s + "\"");
}
std::string Model::getPotentialTypeName(PotentialType t){
    static const char* res[] = {
        "ANALYTICAL",
        //"EXTERNAL",
        "DATADRIVEN"
    };
    return res[static_cast<unsigned>(t)];
}
std::map<std::string, Model::PotentialType>& Model::getPotentialTypeToNameMap(){
    static std::map<std::string, Model::PotentialType> rmap{
        {"ANALYTICAL", ANALYTICAL},
        //{"EXTERNAL", EXTERNAL},
        {"DATADRIVEN", DATADRIVEN}
    };
    return rmap;
}
Model::PotentialType Model::getPotentialTypeByName(const std::string& s){
    const auto& rmap = getPotentialTypeToNameMap();
    auto it = rmap.find(s);
    if (it != rmap.end()) return it->second;
    else throw std::runtime_error("Can't find potential type with name = \"" + s + "\"");
}

std::ostream& Model::TraitValue::print(std::ostream& out) const {
    out << (m_real ? "DVAL" : "IVAL") << "[" << m_dim << "] ";
    out << "on "<< MetaTriMesh::get_sparsity_name(m_etype) << " ";
    if (m_tag_name.empty() && std::isnan(m_def_value[0])){
        out << "NOT SPECIFIED";
        return out;
    } 
    if (!m_tag_name.empty())
        out << "TAG = \"" << m_tag_name << "\" ";
    if (m_dim > 0 && !std::isnan(m_def_value[0])){
        out << "VAL = ";
        if (m_dim > 1) out << "{ ";
        out << m_def_value[0];
        for (std::size_t i = 1; i < m_dim; ++i)
            out << ", " << m_def_value[i];
        if (m_dim > 1) out << " }";
    }    
    return out;
}
std::ostream& Model::SaveTrait::print(std::ostream& out) const {
    out << "DVAL[" << m_dim << "] on " << MetaTriMesh::get_sparsity_name(m_etype) << " ";
    out << (m_tag_name.empty() ? std::string("NOT SPECIFIED") : ("TAG = \"" + m_tag_name + "\" "));
    return out; 
}
std::map<Model::ObjectTraitType, Model::TraitValue> Model::getDefaultTraits(){
    using ET = MetaTriMesh::ElementSparsity;
    return std::map<ObjectTraitType, TraitValue>{
        {INIT_X_VECTOR_DTAG, TraitValue{ET::NODE, 3, "v:point|self:point"}},//INIT_X_VECTOR_DTAG
        {START_X_VECTOR_DTAG, {ET::NODE, 3, "v:x|self:point"}},             //START_X_VECTOR_DTAG
        {BOUNDARY_ITAG, {ET::NODE, 1, "v:bnd", {NAN, NAN, NAN}, false}},    //BOUNDARY_ITAG
        {PRESSURE_DTAG, {ET::FACE, 1, "f:pressure"}},                       //PRESSURE_DTAG
        {LAGRANGE_PRESSURE_DTAG, {ET::FACE, 1, "f:lagrange_pressure"}},     //LAGRANGE_PRESSURE_DTAG
        {BODY_FORCE_VECTOR_DTAG, {ET::NODE, 3, "v:body_force"}},            //BODY_FORCE_VECTOR_DTAG
        {EDGE_FORCE_VECTOR_DTAG, {ET::EDGE, 3, "e:edge_force"}},            //EDGE_FORCE_VECTOR_DTAG
        {CLAMPED_EDGE_VECTOR_DTAG, {ET::EDGE, 3, "be:clamped_direction"}},  //CLAMPED_EDGE_VECTOR_DTAG
        {THICKNESS_DTAG, {ET::FACE, 1, "f:thickness"}},                     //THICKNESS_DTAG
    };
}
std::map<unsigned int, Model::SaveTrait> Model::getDefaultSaveTraits(){
    using ET = MetaTriMesh::ElementSparsity;
    return std::map<unsigned int, SaveTrait>{
        std::pair<unsigned int, SaveTrait>{SAVE_FINAL_STATE, SaveTrait{ET::NODE, 3, "v:x"}},
        {SAVE_RESIDUAL, {ET::NODE, 3, ""}},
        {SAVE_LOCAL_BASIS, {ET::FACE, 3, ""}},
        {SAVE_LOCAL_BASIS_EXT, {ET::FACE, 3, ""}},
        {SAVE_PK2_STRESS, {ET::FACE, 3, ""}},
        {SAVE_DEF_GRAD, {ET::FACE, 6, ""}},
        {SAVE_CAUCHY_STRAIN, {ET::FACE, 3, ""}},
        {SAVE_CURVATURE, {ET::FACE, 3, ""}},
        {SAVE_CAUCHY_STRESS3D, {ET::FACE, 6, ""}},
        {SAVE_CAUCHY_TENSION3D, {ET::FACE, 6, ""}},
    };
}
const std::map<std::string, unsigned int>& Model::getSaveTraitToNameMap(){
    static std::map<std::string, unsigned int> rmap{
        {"FINAL_CONFIGURATION", SAVE_FINAL_STATE},
        {"RESIDUAL", SAVE_RESIDUAL},
        {"LOCAL_BASIS", SAVE_LOCAL_BASIS|SAVE_LOCAL_BASIS_EXT},
        {"PK2_STRESS", SAVE_PK2_STRESS},
        {"DEFORMATION_GRADIENT", SAVE_DEF_GRAD},
        {"CAUCHY_STRAIN", SAVE_CAUCHY_STRAIN},
        {"FINAL_CURVATURE", SAVE_CURVATURE},
        {"CAUCHY_STRESS_3D", SAVE_CAUCHY_STRESS3D},
        {"CAUCHY_TENSION_3D", SAVE_CAUCHY_TENSION3D},
    };
    return rmap;
}
std::vector<std::pair<unsigned int, Model::SaveTrait*>> Model::getSaveTraitByName(std::map<unsigned int, SaveTrait>& save_traits, const std::string& name){
    std::vector<std::pair<unsigned int, SaveTrait*>> res;
    auto& rmap = getSaveTraitToNameMap();
    auto it = rmap.find(name);
    if (it == rmap.end())
        throw std::runtime_error("Faced unknown postcomputation value = \"" + name + "\"");
    auto stat = it->second;
    unsigned long comp_stat = 1;
    while(comp_stat <= MAX_SAVE_FLAG){
        if (stat & comp_stat){
            auto jt = save_traits.find(comp_stat);
            if (jt != save_traits.end()){
                res.push_back({comp_stat, &(jt->second)});
            }
        }
        comp_stat = (comp_stat << 1);
    }
    return res;    
}
std::string Model::getSaveTraitName(unsigned int trait){
    auto& rmap = getSaveTraitToNameMap();
    for (auto& v: rmap) if (v.second & trait) 
        return v.first;
    return std::string();
}
using MDDIR = Model::DataDrivenPotentialInfo::InterpRegion;
std::string MDDIR::typeName(Type tp){
    static const char* types[] = {"KNEAREST", "LINEAR"};
    return types[static_cast<int>(tp)];
}
Model::DataDrivenPotentialInfo::Type MDDIR::typeByName(const std::string& val){
    static std::map<std::string, Type> rmap{
        {"KNEAREST", KNEAREST},
        {"LINEAR", LINEAR},
    };
    return rmap.at(val);
}
std::ostream& MDDIR::print(std::ostream& out) const{
    switch (m_tp){
        case KNEAREST:{
            auto a = knearest();
            return out << typeName(KNEAREST) << " " << a.R_trust << " " << a.k << " " << a.metric_pow;
        }
        case LINEAR:{
            auto a = linear();
            return out << typeName(LINEAR) << " " << a.R_trust << " " << a.R_fit << " " << std::bitset<8>(a.octant_mask);
        }
    }
    return out;
}
std::string Model::CustomSolverIteration::solverName(int sol_type){
    switch(sol_type) {
        case 0: return "RELAX";
        case 1: return "NEWTON";
        case 2: return "KINSOL";
        default: return "KINSOL";
    }
} 
int Model::CustomSolverIteration::solverTypeByName(const std::string& val){
    static std::map<std::string, int> smap{
        {"RELAX", 0},
        {"NEWTON", 1},
        {"KINSOL", 2}
    };
    auto it = smap.find(val);
    if (it == smap.end())
        throw std::runtime_error("Faced unknown nonlinear solver type = \"" + val + "\"");
    return it->second;    
}
std::ostream& Model::CustomSolverIteration::print(std::ostream& out) const {
    switch(ns_type) {
        case 0: out << "RELAX " << delta_full_step << " " << relax_its << " " << relax_pfreq; break;
        case 1: out << "NEWTON " << max_ns_its << " " << delta_full_step; break;
        default: out << "KINSOL " << max_ns_its; break;
    }
    return out;
}
Model& Model::parse(int argc, const char* argv[], int ignore_argv_cnt){
    std::vector<std::string> args(argc - ignore_argv_cnt);
    for (std::size_t i = 0; i < argc - ignore_argv_cnt; ++i)
        args[i] = argv[ignore_argv_cnt + i];
    return parse(args);    
}

Model& Model::parse(std::vector<std::string> args){
    MetaTriMesh mtm;    bool specific_mesh = false;
    Energy_Parser energy_expr; bool specific_energy = false;
    std::map<ObjectTraitType, TraitValue> traits = getDefaultTraits();
    std::map<ObjectTraitType, bool> trait_associated;
    std::map<unsigned int, SaveTrait> save_traits = getDefaultSaveTraits();
    std::map<MetaTriMesh::ElementSparsity, std::map<std::string, SaveInTag>> save_intag;
    bool print_final_state = true;
    bool end_after_parsing = false;

    auto set_mesh = [&](){
        path p(fmesh); 
        mtm.clear();
        bool status = true;
        if (p.extension() == ".tmh")
            status = TMH_File().read(fmesh, mtm);
        else if (p.extension() == ".vtk")
            VTK_File().read(fmesh, mtm);
        else 
            throw std::runtime_error("Mesh file have unexpected extension \"" + p.extension().string() + "\" try .vtk or .tmh files");    
        if (!status) 
            throw std::runtime_error("Can't read data from file " + fmesh);
    };
    auto associate_mesh_data_with_energy_params = [&](){
        if (mtm.empty()) return;
        auto& prm = energy_expr.getParameters();
        auto& fib = energy_expr.getFibers();
        if (prm.empty() && fib.empty()) return;
        for (auto& p: prm){
            if (!p.def_tag.empty()){
                auto it = mtm.cell_data.find(p.def_tag);
                if (it != mtm.cell_data.end() && it->second.m_values.dim() == 1) p.def_val = NAN;
                else p.def_tag.resize(0);
            }
        }
        for (auto& p: fib){
            if (!p.def_tag.empty()){
                auto it = mtm.cell_data.find(p.def_tag);
                if (it == mtm.cell_data.end() || it->second.m_values.dim() != 3)
                    p.def_tag.resize(0);
            }
        }
    };
    auto associate_mesh_with_traits = [&](){
        if (mtm.empty()){
            for(auto& i: traits) trait_associated[i.first] = false;
            return;
        }
        {
            auto& v = traits[INIT_X_VECTOR_DTAG];
            if (v.m_tag_name == "v:point|self:point" || !mtm.value_names<V_ind>().count(v.m_tag_name) || mtm.point_data[v.m_tag_name].m_values.dim() != 3){
                if (mtm.value_names<V_ind>().count("v:point") && mtm.point_data["v:point"].m_values.dim() == 3){
                    v.m_tag_name = "v:point";
                }else
                    v.m_tag_name = "self:point";
            }
            v.m_def_value = std::array<double, 3>{NAN, NAN, NAN};
            trait_associated[INIT_X_VECTOR_DTAG] = true;
        }
        {
            auto& v = traits[START_X_VECTOR_DTAG];
            if (v.m_tag_name == "v:x|self:point"|| !mtm.value_names<V_ind>().count(v.m_tag_name) || mtm.point_data[v.m_tag_name].m_values.dim() != 3){
                if (mtm.value_names<V_ind>().count("v:x") && mtm.point_data["v:x"].m_values.dim() == 3){
                    v.m_tag_name = "v:x";
                }else
                    v.m_tag_name = "self:point";
            }
            v.m_def_value = std::array<double, 3>{NAN, NAN, NAN};
            trait_associated[START_X_VECTOR_DTAG] = true;
        }
        for (auto& v: traits){
            if (v.first == INIT_X_VECTOR_DTAG || v.first == START_X_VECTOR_DTAG) continue;
            auto& r = v.second;
            if (!r.m_tag_name.empty()){
                auto& d = mtm.data_array(r.m_etype);
                auto it = d.find(r.m_tag_name);
                if (it == d.end() || it->second.m_values.m_dim != r.m_dim)
                    r.m_tag_name.resize(0);
                else   
                    r.m_def_value = std::array<double, 3>{NAN, NAN, NAN};    
            }
            trait_associated[v.first] = (!r.m_tag_name.empty() || !std::isnan(r.m_def_value[0]));
        }
    };
    auto make_same_strings = [](std::vector<std::string>& words){
        std::size_t msz = 0;
        for (auto& s: words) msz = std::max(msz, s.size());
        for (auto& s: words) s.resize(msz, ' ');
    };
    auto show_traits = [&](std::ostream& out = std::cout, const std::string& prefix = ""){
        std::vector<std::string> type_names;
        std::vector<std::string> info_str;
        std::vector<bool> assoc;
        for (auto& v: traits){
            type_names.push_back(getObjectTraitTypeName(v.first));
            std::stringstream oss;
            v.second.print(oss);
            info_str.push_back(oss.str());
            assoc.push_back(trait_associated[v.first]);
        }
        make_same_strings(type_names); 
        make_same_strings(info_str);
        for (std::size_t i = 0; i < type_names.size(); ++i){
            out << prefix << type_names[i] << ": " << info_str[i] << " -> " << (assoc[i] ? "ASSOCIATED" : "NOT ASSOCIATED") << "\n";
        }
    };
    auto show_save_traits = [&](std::ostream& out = std::cout, const std::string& prefix = ""){
        std::vector<std::string> names, info;
        {
            names.push_back(getSaveTraitName(SAVE_LOCAL_BASIS | SAVE_LOCAL_BASIS_EXT));
            std::stringstream oss;
            auto &v0 = save_traits[SAVE_LOCAL_BASIS], &v1 = save_traits[SAVE_LOCAL_BASIS_EXT];
            oss << "DVAL[" << v0.m_dim << "] on " << MetaTriMesh::get_sparsity_name(v0.m_etype) << " ";
            if (v0.m_tag_name.empty() || v1.m_tag_name.empty()){
                oss << "NOT SPECIFIED";
            } else {
                oss << "TAGS = \"" << v0.m_tag_name + "\", \"" << v1.m_tag_name + "\"";
            }
            info.push_back(oss.str());
        }
        for (auto& v: save_traits){
            if (v.first & (SAVE_LOCAL_BASIS | SAVE_LOCAL_BASIS_EXT)) continue;
            names.push_back(getSaveTraitName(v.first));
            std::stringstream oss; 
            v.second.print(oss);
            info.push_back(oss.str());
        }
        make_same_strings(names); 
        make_same_strings(info);
        for (std::size_t i = 0; i < names.size(); ++i)
            out << prefix << names[i] << ": " << info[i] << "\n";
    };

    auto is_double = [](const std::string& s) -> std::pair<bool, double>{
        std::istringstream iss(s);
        double f = NAN;
        iss >> f; 
        return {iss.eof() && !iss.fail(), f};    
    }; 
    auto is_int = [](const std::string& s) -> std::pair<bool, int>{
        std::istringstream iss(s);
        int f = 0;
        iss >> f; 
        return {iss.eof() && !iss.fail(), f};    
    };     
    auto to_bool = [](const std::string& s) -> bool{
        std::string s0 = s;
        std::transform(s0.begin(), s0.end(), s0.begin(), [](char c){ return std::tolower(c); });
        if (s0 == "true") return true;
        else if (s0 == "false") return false;
        int r = stoi(s0);
        if (r == 0) return false;
        if (r == 1) return true;
        throw std::runtime_error("Can't convert \"" + s + "\" to boolean");
        return false;
    };
    
    auto compare = [&args](std::size_t j, const char* long_nm, const char* short_nm){
        return (long_nm && !strcmp(long_nm, args[j].c_str())) || (short_nm && !strcmp(short_nm, args[j].c_str()));
    }; 
    for (std::size_t i = 0; i < args.size(); ++i){
        auto expect_at_least_n_args = [&i, &args](std::size_t n) { 
            if (i + n >= args.size()) 
                throw std::runtime_error("Expected at least " + std::to_string(n) + " arguments after " + args[i]);
        };
        auto update_save_flag = [&](unsigned char state){
            expect_at_least_n_args(1);
            bool save = to_bool(args[++i]);
            if (save) save_flag |= state;
            else save_flag &= ~state;
        };
        if (compare(i, "-cf", "--config")){
            expect_at_least_n_args(1);
            std::string config_file = args[++i];
            parse_config_file(config_file, args, i-1, 2);
            i -= 2;
            continue;
        } else if (compare(i, "-m", "--mesh")){
            expect_at_least_n_args(1);
            fmesh = args[++i];
            set_mesh();
            specific_mesh = true;
            continue;    
        } else if (compare(i, "-e", "--energy")){
            expect_at_least_n_args(1);
            m_potential.setType(ANALYTICAL);
            auto& a = m_potential.analytical();
            a.fpotential = args[++i];
            energy_expr = Energy_Parser().setCodeFromFile(a.fpotential).LexicalAnalysis().SyntaxAnalysis().Simplify();
            specific_energy = true;
            // energy_expr.print_parsed_state(std::cout) << "\n";
            continue;
        } 
#ifdef WITH_DATADRIVEN
        else if (compare(i, "-dd", "--data_tabs")){
            expect_at_least_n_args(1);
            m_potential.setType(DATADRIVEN);
            auto& d = m_potential.datadriven();
            do{
                d.data_files.push_back(args[++i]);
            } while (i+1 < args.size() && !args[i+1].empty() && args[i+1][0] != '-');
            continue;
        }
#endif 
        else if (compare(i, "-t", "--target")){
            expect_at_least_n_args(1);
            save_dir = args[++i];
            continue; 
        }  else if (compare(i, "-gd", "--gen_dir")){
            expect_at_least_n_args(1);
            to_gen = args[++i];
            continue;
        } else if (compare(i, "-nm", "--name")){
            expect_at_least_n_args(1);
            save_prefix = args[++i];
            continue;
        } else if (compare(i, "-sm", "--show_mesh")){
            mtm.print_content_info() << "\n";
            continue;
        } else if (compare(i, "-se", "--show_elast")){
            if (m_potential.m_type == ANALYTICAL){
                associate_mesh_data_with_energy_params();
                energy_expr.print_parsed_state(std::cout) << "\n";
            } else {
                std::cout << "Warning: elastic parameters available only for " << getPotentialTypeName(ANALYTICAL) << " potential, "
                          << "but current potential type is " << getPotentialTypeName(m_potential.m_type) << "\n" << std::endl;
            }
            continue;
        } else if (compare(i, "-sl", "--show_lin")){
            auto slvs = LinearSolver::getAvailableSolvers();
            std::cout << "Available solvers:\n";
            for (auto& v: slvs)
                std::cout << "    " << v << "\n";
            std::cout << "\n";     
            continue;
        } else if (compare(i, "-mt", "--model_type")){
            expect_at_least_n_args(1);
            auto val = args[++i];
            if (strcmp(val.c_str(), "MEMBRANE") == 0)
                is_membrane_approx = true;
            else if (strcmp(val.c_str(), "SHELL") == 0)  
                is_membrane_approx = false;
            else 
                throw std::runtime_error("Expected name of model theory (MEMBRANE or SHELL) but faced \"" + std::string(val) + "\"");
            continue;          
        } else if (compare(i, "-ep", "--elast_prm")){
            if (m_potential.m_type != ANALYTICAL)
                throw std::runtime_error("Specification of elastic parameters available only for " + getPotentialTypeName(ANALYTICAL) + 
                    " potential, but current potential type is " + getPotentialTypeName(m_potential.m_type));
            expect_at_least_n_args(2);
            auto name = args[++i];
            auto sval = args[++i];
            auto& fib = energy_expr.getFibers();
            auto fit = std::find_if(fib.begin(), fib.end(), [&name](auto& v){ return v.name == name; });
            if (fit != fib.end()) {
                fit->def_tag = sval;
                continue;
            }
            auto& prm = energy_expr.getParameters();
            auto pit = std::find_if(prm.begin(), prm.end(), [&name](auto& v){ return v.name == name; });
            if (pit != prm.end()){
                auto val = is_double(sval);
                if (val.first){
                    pit->def_val = val.second;
                    pit->def_tag = "";
                } else {
                    pit->def_val = NAN;
                    pit->def_tag = sval;
                }
            }
            continue;
        } else if (compare(i, "-st", "--show_trait")){
            associate_mesh_with_traits();
            show_traits();
            continue;
        } else if (compare(i, "-tr", "--trait_from")){
            expect_at_least_n_args(2);
            auto name = args[++i];
            auto sval = args[++i];
            auto tp = getObjectTraitTypeByName(name);
            auto& v = traits[tp];
            auto r = is_double(sval);
            if (!r.first) {
                v.m_tag_name = sval;
                continue;
            }
            v.m_def_value[0] = r.second;
            if (v.m_dim > 1){
                expect_at_least_n_args(v.m_dim - 1);
                for (std::size_t k = 1; k < v.m_dim; ++k){
                    auto sval1 = args[++i];
                    auto r1 = is_double(sval1);
                    if (!r1.first)
                        throw std::runtime_error("Expected real value for trait \"" + std::string(name) + "\" but found \"" + sval1 + "\"");
                    v.m_def_value[k] = r1.second;
                }
            }
            v.m_tag_name.resize(0);
            continue;
        } else if (compare(i, "-vw", "--view")){
            expect_at_least_n_args(1);
            m_view = to_bool(args[++i]);
            continue;
        } else if (compare(i, "-rg", "--regen_expr")){
            expect_at_least_n_args(1);
            regen_potential = to_bool(args[++i]);
            continue;
        } 
#ifdef WITH_DATADRIVEN        
        else if (compare(i, "-di", "--dat_interp")){
            if (m_potential.m_type != DATADRIVEN)
                throw std::runtime_error("Specification of data interpolation available only for " + getPotentialTypeName(DATADRIVEN) + 
                    " potential type, but current potential type is " + getPotentialTypeName(m_potential.m_type));
            expect_at_least_n_args(3);
            do {
                auto interp_name = args[++i];
                using DDI = DataDrivenPotentialInfo;
                DDI::InterpRegion ir;
                ir.setType(DDI::InterpRegion::typeByName(interp_name));
                int expect_args = 0;
                switch (ir.m_tp){
                    case DDI::KNEAREST: expect_args = 3; break;
                    case DDI::LINEAR: expect_args = 3; break;
                }
                expect_at_least_n_args(expect_args);
                auto r0 = is_double(args[++i]);
                if (!r0.first)
                    throw std::runtime_error("Expected real parameter TRUSTED_REGION, but faced = \"" + args[i] + "\"");
                std::string v1 = args[++i];
                std::string v2 = args[++i];
                switch (ir.m_tp){
                    case DDI::KNEAREST:{
                        auto i1 = is_int(v1);
                        if (!i1.first)
                            throw std::runtime_error("Expected integer parameter k_neighbour, but faced = \"" + v1 + "\"");
                        auto r2 = is_double(v2);
                        if (!r2.first)
                            throw std::runtime_error("Expected real parameter inverse_distance_metric_power, but faced = \"" + v2 + "\"");    
                        auto& a = ir.knearest();
                        a.R_trust = r0.second; a.k = i1.second; a.metric_pow = r2.second;
                        break;
                    }
                    case DDI::LINEAR: {
                        auto r1 = is_double(v1);
                        if (!r1.first)
                            throw std::runtime_error("Expected real parameter fit_radius, but faced = \"" + v1 + "\""); 
                        auto i2 = is_double(v2);
                        if (!i2.first)
                            throw std::runtime_error("Expected uchar parameter octant_mask, but faced = \"" + v2 + "\"");
                        auto& a = ir.linear();
                        a.R_trust = r0.second; a.R_fit = r1.second; a.octant_mask = i2.second;
                        break;  
                    }
                }  
                m_potential.datadriven().interp.push_back(ir); 
            } while(i+1 < args.size() && !args[i+1].empty() && args[i+1][0] != '-');
            continue;
        } 
#endif        
        else if (compare(i, "-ns", "--nlin_sol")){
            expect_at_least_n_args(2);
            do{
                auto nslv_name = args[++i];
                CustomSolverIteration csi;
                csi.ns_type = CustomSolverIteration::solverTypeByName(nslv_name);
                int expect_args = 0;
                switch(csi.ns_type){
                    case 0: expect_args = 3; break;
                    case 1: expect_args = 2; break;
                    default: expect_args = 1; break;
                }
                expect_at_least_n_args(expect_args);
                std::string v0 = args[++i], v1, v2;
                if (csi.ns_type != 2) v1 = args[++i];
                if (csi.ns_type == 0) v2 = args[++i];
                switch(csi.ns_type){
                    case 0: {
                        auto r0 = is_double(v0);
                        if (!r0.first)
                            throw std::runtime_error("Expected real parameter for RELAX solver, but faced = \"" + v0 + "\"");
                        auto r1 = is_int(v1), r2 = is_int(v2);
                        if (!r1.first || !r2.first)
                            throw std::runtime_error("Expected integer parameter for RELAX solver, but faced = \"" + (!r1.first ? v1 : v2) + "\"");
                        csi.delta_full_step = r0.second;
                        csi.relax_its = r1.second, csi.relax_pfreq = r2.second;    
                        break;
                    }
                    case 1: {
                        auto r1 = is_int(v0); auto r2 = is_double(v1);
                        if (!r1.first || !r2.first)
                            throw std::runtime_error("Expected parameter for NEWTON solver, but faced = \"" + (!r1.first ? v0 : v1) + "\"");
                        csi.max_ns_its = r1.second, csi.delta_full_step = r2.second;    
                    }
                    default:{
                        auto r1 = is_int(v0);
                        if (!r1.first)
                            throw std::runtime_error("Expected maximal number of iterations for KINSOL solver, but faced = \"" + v0 + "\"");
                        csi.max_ns_its = r1.second;    
                    }
                }
                m_solverItData.push_back(csi);
            } while(i+1 < args.size() && !args[i+1].empty() && args[i+1][0] != '-');
            continue;
        } else if (compare(i, "-ls", "--lin_sol")){
            expect_at_least_n_args(1);
            m_lin_sol_name = args[++i];
            auto slvs = LinearSolver::getAvailableSolvers();
            auto it = std::find(slvs.begin(), slvs.end(), m_lin_sol_name);
            if (it == slvs.end())
                throw std::runtime_error("Solver \"" + m_lin_sol_name + "\" is not available");
            continue;
        }  else if (compare(i, "-lp", "--lin_prm")){
            expect_at_least_n_args(2);
            auto prm_name = args[++i];
            auto prm_val = args[++i];
            m_lin_sol_prms.push_back(std::pair<std::string, std::string>{prm_name, prm_val});
            continue;
        }  else if (compare(i, "-nt", "--nlin_tol")){
            expect_at_least_n_args(2);
            auto abs_tol_str = args[++i];
            auto rel_tol_str = args[++i];
            auto r1 = is_double(abs_tol_str);
            if (!r1.first)  throw std::runtime_error("Expected absolute tolerance value but faced \"" + std::string(abs_tol_str) + "\"" );
            auto r2 = is_double(rel_tol_str);
            if (!r2.first)  throw std::runtime_error("Expected relative tolerance value but faced \"" + std::string(rel_tol_str) + "\"" );
            m_abs_tol = r1.second; m_rel_tol = r2.second;
            continue;
        } else if (compare(i, "-q", "--quit")){
            exit(0);
            continue;
        } else if (compare(i, "-h", "--help")){
            printArgsHelpMessage();
            if (i == args.size() - 1) exit(0);
            continue;
        } else if (compare(i, nullptr, "--bin_vtk")){
            update_save_flag(SAVE_BIN_VTK);
            continue;
        } else if (compare(i, nullptr, "--txt_vtk")){
            update_save_flag(SAVE_TXT_VTK);
            continue;
        } else if (compare(i, nullptr, "--bin_tmh")){
            update_save_flag(SAVE_BIN_TMH);
            continue;
        } else if (compare(i, nullptr, "--txt_tmh")){
            update_save_flag(SAVE_TXT_TMH);
            continue;
        } else if (compare(i, "-mv", "--rename_tag")){
            expect_at_least_n_args(3);
            auto sp = MetaTriMesh::get_sparsity_by_name(args[++i]);
            std::string old_name = args[++i];
            std::string new_name = args[++i];
            auto& dt = mtm.data_array(sp);
            auto it = dt.find(old_name);
            if (it == dt.end()) continue;
            auto jt = dt.find(new_name);
            if (jt != dt.end())
                throw std::runtime_error("Can't renme tag \"" + old_name + "\" into \"" + new_name + "\" because tag with name \"" + new_name + "\" exists");
            dt.emplace(std::move(new_name), std::move(it->second));
            dt.erase(it);
            continue;
        } else if (compare(i, "-rm", "--remove_tag")){
            expect_at_least_n_args(2);
            auto sp = MetaTriMesh::get_sparsity_by_name(args[++i]);
            std::string tag_name = args[++i];
            auto& dt = mtm.data_array(sp);
            dt.erase(tag_name);
            continue;
        } else if (compare(i, "-sc", "--show_comp")){
            show_save_traits();
            continue;
        } else if (compare(i, "-ac", "--save_comp")){
            expect_at_least_n_args(2);
            auto trait_name = args[++i];
            auto tag_name = args[++i];
            auto v = getSaveTraitByName(save_traits, trait_name);
            if (v.size() == 0)
                throw std::runtime_error("Expected name of postcomputation trait but faced unknown trait =\"" + std::string(trait_name) + "\"");
            if (v.size() == 1){
                v[0].second->m_tag_name = tag_name;
            }else{
                for (std::size_t i = 0; i < v.size(); ++i)
                    v[i].second->m_tag_name = std::string(tag_name) + "_" + std::to_string(i);
            }    
            continue;
        } else if (compare(i, nullptr, "--show_save_intag")){
            using ET = MetaTriMesh::ElementSparsity;
            std::vector<std::string> names, info, status;
            for (std::size_t i = ET::NODE; i <= ET::HALFEDGE; ++i){
                const auto& vd = mtm.data_array(static_cast<ET>(i));
                auto& sd = save_intag[static_cast<ET>(i)];
                for (auto& v: vd){
                    names.push_back("\"" + v.first + "\"");
                    std::stringstream oss;
                    oss << MetaTriMesh::get_type_name(v.second.m_values.dataType()) 
                    << "[" << v.second.m_values.dim() << "] on " 
                    << MetaTriMesh::get_sparsity_name(static_cast<ET>(i));
                    info.push_back(oss.str());
                    auto jt = sd.find(v.first);
                    if (jt == sd.end()){
                        status.push_back("NOT SAVE");
                    }  else {
                        if (!jt->second.m_new_tag_name.empty()) {
                            status.push_back("SAVE in \"" + jt->second.m_new_tag_name + "\"");
                        }
                        else
                            status.push_back("SAVE");
                    } 
                }
            }
            make_same_strings(names);
            make_same_strings(info);
            for (std::size_t k = 0; k < names.size(); ++k)
                std::cout << names[k] << ": " << info[k] << " " << status[k] << "\n";
            continue;
        } else if (compare(i, nullptr, "--save_intag")){
            expect_at_least_n_args(2);
            auto etype = MetaTriMesh::get_sparsity_by_name(args[++i]);
            std::string tag_name = args[++i];
            auto& vd = mtm.data_array(etype);
            auto it = vd.find(tag_name);
            if (it != vd.end())
                save_intag[etype][tag_name] = SaveInTag(etype, tag_name);
            continue;
        } else if (compare(i, nullptr, "--save_move_intag")){
            expect_at_least_n_args(3);
            auto etype = MetaTriMesh::get_sparsity_by_name(args[++i]);
            std::string tag_name = args[++i];
            std::string new_tag_name = args[++i];
            auto& vd = mtm.data_array(etype);
            auto it = vd.find(tag_name);
            if (it != vd.end())
                save_intag[etype][tag_name] = SaveInTag(etype, tag_name, new_tag_name);
            continue;
        } else if (compare(i, nullptr, "--save_all_intag")){
            using ET = MetaTriMesh::ElementSparsity;
            for (std::size_t i = ET::NODE; i <= ET::HALFEDGE; ++i){
                auto et = static_cast<ET>(i);
                const auto& vd = mtm.data_array(et);
                auto& sd = save_intag[et];
                for (auto& v: vd){
                    auto jt = sd.find(v.first);
                    if (jt == sd.end()){
                        sd.insert({v.first, SaveInTag(et, v.first)});
                    }
                }
            }
            continue;
        } else if (compare(i, nullptr, "--save_used_intag")){
            using ET = MetaTriMesh::ElementSparsity;
            associate_mesh_with_traits();
            associate_mesh_data_with_energy_params();
            for (auto& v: traits){
                if (!trait_associated[v.first]) continue;
                if (v.second.m_tag_name.empty()) continue;
                if ((v.first == INIT_X_VECTOR_DTAG || v.first == START_X_VECTOR_DTAG) && v.second.m_tag_name == "self:point") continue;
                auto& sd = save_intag[v.second.m_etype];
                auto jt = sd.find(v.second.m_tag_name);
                if (jt == sd.end()) 
                    sd.insert({v.second.m_tag_name, SaveInTag(v.second.m_etype, v.second.m_tag_name)});
            }
            auto& fd = save_intag[ET::FACE];
            auto& prm = energy_expr.getParameters();
            auto& fib = energy_expr.getFibers();
            for (auto& p: prm) if (!p.def_tag.empty()){
                auto jt = fd.find(p.def_tag);
                if (jt == fd.end())
                    fd.insert({p.def_tag, SaveInTag(ET::FACE, p.def_tag)});
            }
            for (auto& p: fib) if (!p.def_tag.empty()){
                auto jt = fd.find(p.def_tag);
                if (jt == fd.end())
                    fd.insert({p.def_tag, SaveInTag(ET::FACE, p.def_tag)});
            }
            continue;
        } else if (compare(i, nullptr, "--show_end_parsed")){
            expect_at_least_n_args(1);
            print_final_state = to_bool(args[++i]);
            continue;
        } else if (compare(i, nullptr, "--end_after_parse")){
            end_after_parsing = true;
            continue;
        } else if (compare(i, nullptr, "--perform_solve")){
            expect_at_least_n_args(1);
            m_perform_solve = to_bool(args[++i]);
            continue;
        } else if (compare(i, nullptr, "--perform_pcomp")){
            expect_at_least_n_args(1);
            m_perform_postcomps = to_bool(args[++i]);
            continue;
        }
        throw std::runtime_error("Faced unparsed command line argument \"" + std::string(args[i]) + "\"");
    }
    if (!specific_mesh) set_mesh();
    if (!specific_energy && m_potential.m_type == ANALYTICAL) {
        auto& a = m_potential.analytical();
        energy_expr = Energy_Parser().setCodeFromFile(a.fpotential).LexicalAnalysis().SyntaxAnalysis().Simplify();
    } else if (m_potential.m_type == DATADRIVEN){
        auto& a = m_potential.datadriven();
        if (a.interp.empty())
            a.interp.push_back(DataDrivenPotentialInfo::InterpRegion());
    }
    if (m_solverItData.empty()){
        CustomSolverIteration csi;
        csi.ns_type = 2;
        csi.max_ns_its = 200;
        m_solverItData.emplace_back(std::move(csi));
    }
    associate_mesh_with_traits();
    for (auto t: trait_associated) if (!t.second){
        traits[t.first].m_def_value[0] = NAN;
        traits[t.first].m_tag_name.resize(0);
    }
    if (m_potential.m_type == ANALYTICAL)
        associate_mesh_data_with_energy_params();
    if (!save_dir.empty() && save_dir.back() != '/') save_dir += "/";
    #ifndef USE_MAGNUM_GUI
    m_view = false;
    #endif

    if (print_final_state){
        std::cout << "\n";
        std::cout << "Mesh from file \"" + fmesh + "\"\n";
        std::cout << "#P = " << mtm.points.size() << " #E = " << mtm.edges.size() << " #C = " << mtm.cells.size() << "\n";
        using ET = MetaTriMesh::ElementSparsity;
        std::vector<std::string> names, storage_type, data_type, size_info, dim_val, status, sparsity;
        for (std::size_t i = ET::NODE; i <= ET::HALFEDGE; ++i){
            auto et = static_cast<ET>(i);
            const auto& vd = mtm.data_array(et);
            auto& sd = save_intag[et];
            for (auto& v: vd){
                sparsity.push_back(MetaTriMesh::get_sparsity_name(et));
                names.push_back("\"" + v.first + "\"");
                storage_type.push_back(MetaTriMesh::MeshValue::getTypeName(v.second.m_type));
                data_type.push_back(MetaTriMesh::get_type_name(v.second.m_values.dataType()));
                switch(v.second.m_type){
                    case MetaTriMesh::MeshValue::CONSTANT: size_info.push_back(""); break;
                    case MetaTriMesh::MeshValue::CONTIGUOUS: size_info.push_back("size = " + std::to_string(v.second.m_values.size())); break;
                    case MetaTriMesh::MeshValue::SPARSE: size_info.push_back("sparse_size = " + std::to_string(v.second.m_sparse_indexes.size())); break;
                }
                dim_val.push_back(std::to_string(v.second.m_values.dim()));
                auto jt = sd.find(v.first);
                if (jt == sd.end()){
                    status.push_back("WILL NOT SAVE");
                }  else {
                    if (!jt->second.m_new_tag_name.empty()) {
                        status.push_back("WILL SAVE in \"" + jt->second.m_new_tag_name + "\"");
                    }
                    else
                        status.push_back("WILL SAVE");
                }
            }
        }
        make_same_strings(names);
        make_same_strings(sparsity);
        make_same_strings(storage_type);
        make_same_strings(data_type);
        make_same_strings(size_info);
        make_same_strings(dim_val);
        make_same_strings(status);
        for (std::size_t i = 0; i < names.size(); ++i){
            std::cout << "  " << names[i] << ": " << storage_type[i] << " " << size_info[i] << " " << data_type[i] << "[" << dim_val[i] << "] on " << sparsity[i] << " " << status[i] << "\n";
        }
        std::cout << "\n";
        if (m_potential.m_type == ANALYTICAL){
            auto& a = m_potential.analytical();
            std::cout << getPotentialTypeName(m_potential.m_type) << " potential from file \"" << a.fpotential << "\"" << (regen_potential ? " will be regenerated" : " will NOT be regenerated") << " in directory \"" << to_gen << "\"\n";
            energy_expr.print_parsed_state(std::cout, "  ");
        } else if (m_potential.m_type == DATADRIVEN){
            auto& a = m_potential.datadriven();
            std::cout << getPotentialTypeName(m_potential.m_type) << " potential with data from files: \n";
            for (auto f: a.data_files)
                std::cout << "  \"" << f << "\"\n";
            std::cout << "Data interpolation pipeline: ";
            for (auto& v: a.interp){
                v.print(std::cout) << " "; 
            }
            std::cout << "\n"; 
        }
        std::cout << "Mechanical model: " << (is_membrane_approx ? "MEMBRANE" : "SHELL") << "\n\n";
        std::cout << "Input traits:\n";
        show_traits(std::cout, "  ");
        std::cout << "\n";
        std::cout << "Output postcomputation values:\n";
        show_save_traits(std::cout, "  ");
        std::cout << "\n";
        #ifdef USE_MAGNUM_GUI
        std::cout << "With GUI window = " << (m_view ? "true" : "false") << "\n";
        #else
        std::cout << "With GUI window = " << "false" << "\n";
        #endif
        std::cout << "Save directory = \"" << save_dir << "\"\n";
        std::cout << " Save files: ";
        for (int i = 0; i < 4; ++i){
            if (save_flag & (1 << i))
                std::cout << (i%2 == 0 ? "binary" : "ascii") << " \"" << save_prefix << (i/2 == 0 ? ".vtk" : ".tmh") << "\" ";
        }
        std::cout << "\n";
        std::cout << "Linear solver \"" << m_lin_sol_name << "\"\n  with params {";
        for (auto& v: m_lin_sol_prms)
            std::cout << "{ \"" << v.first << "\", \"" << v.second << "\"}, ";
        std::cout << "}\n";
        std::cout << "Nonlinear solver:\n"; 
        std::cout << "  Stop conditions: relative_tolerance = " << m_rel_tol << " absolute_tolerance = " << m_abs_tol << "\n";
        std::cout << "  Pipeline: ";
        for (auto& v: m_solverItData){
            v.print(std::cout) << " ";
        }
        std::cout << "\n\n";
    }

    static std::vector<ObjectTraitType> required_traits{
        INIT_X_VECTOR_DTAG,
        START_X_VECTOR_DTAG,
        BOUNDARY_ITAG,
        THICKNESS_DTAG,
    };
    for (std::size_t i = 0; i < required_traits.size(); ++i){
        auto& t = traits[required_traits[i]];
        if (std::isnan(t.m_def_value[0]) && t.m_tag_name.empty())
            throw std::runtime_error("Can't start computations without " + getObjectTraitTypeName(required_traits[i]) + " trait");
    }
    if (m_potential.m_type == ANALYTICAL){
        auto& a = m_potential.analytical();
        for (auto& p: energy_expr.getParameters()) 
            if (p.def_tag.empty() && std::isnan(p.def_val))
                throw std::runtime_error("Can't start computations without association any value with potential parameter \"" + p.name + "\"");
        for (auto& p: energy_expr.getFibers()) 
            if (p.def_tag.empty())
                throw std::runtime_error("Can't start computations without association vector tag on mesh with potential fiber \"" + p.name + "\"");
        a.m_energy_expr = std::move(energy_expr);
    } else if (m_potential.m_type == DATADRIVEN)
        if (m_potential.datadriven().data_files.empty())
            throw std::runtime_error("Can't start computation with DATADRIVEN potential if data tables is not specified");
    
    m_mtm = std::move(mtm);
    m_input_traits = std::move(traits);
    m_out_traits = std::move(save_traits);
    for (auto& sd: save_intag)
        for (auto& v: sd.second)
            m_save_intag[sd.first][v.first] = v.second.m_new_tag_name;
    if (end_after_parsing)
        exit(0);        
    return *this;        
}
void Model::parse_config_file(const std::string& config_file, std::vector<std::string>& args, std::size_t pos, std::size_t delete_next){
    path p(config_file);
    path cur_dir = p;
    cur_dir.remove_filename();
    if (!exists(p))
        throw std::runtime_error("Configuration file \"" + std::string(absolute(p)) + "\" doesn't exists");
    std::ifstream f(p);
    if (!f.is_open())
        throw std::runtime_error("Can't open configuration file \"" + config_file + "\"");
    struct TokenData{
        unsigned id = std::numeric_limits<unsigned>::max();
        unsigned num_args = 0;
        bool is_path_arg = false;
        TokenData() = default;
        TokenData(unsigned id, unsigned num_args = 0, bool is_path_arg = false): id{id}, num_args{num_args}, is_path_arg{is_path_arg} {}
    }; 
#define SETTOKEN(X, Y, NUM, Z) std::pair<std::string, TokenData>{#X, {NUM, Z}}, {#Y, {NUM, Z}}
#define SETTOKENF(X, Y, NUM, F, Z) std::pair<std::string, TokenData>{#X, {NUM, Z, F}}, {#Y, {NUM, Z, F}}
    static std::map<std::string, TokenData> tokens{
        SETTOKENF(-cf, --config, 0, 1, true),
        SETTOKENF(-m , --mesh, 1, 1, true),
        SETTOKENF(-e , --energy, 2, 1, true),
        SETTOKENF(-t , --target, 3, 1, true),
        SETTOKEN(-rm, --remove_tag, 4, 2),
        SETTOKEN(-mv, --rename_tag, 5, 2),
        SETTOKEN(-nm, --name, 6, 1),
        SETTOKEN(-sm, --show_mesh, 7, 0),
        SETTOKEN(-se, --show_elast, 8, 0),
        SETTOKEN(-st, --show_trait, 9, 0),
        SETTOKEN(-sl, --show_lin, 10, 0),
        SETTOKEN(-ep, --elast_prm, 11, 2),
        SETTOKEN(-tr, --trait_from, 12, 2),
        SETTOKEN(-ls, --lin_sol, 13, 1),
        SETTOKEN(-lp, --lin_prm, 14, 2),
        SETTOKEN(-nt, --nlin_tol, 15, 2),
        SETTOKEN(-mt, --model_type, 16, 1),
        SETTOKEN(-vw, --view, 17, 1),
        SETTOKEN(-rg, --regen_expr, 18, 1),
        {"--bin_vtk", {19, 1}}, {"--txt_vtk", {20, 1}}, {"--bin_tmh", {21, 1}}, {"--txt_tmh", {22, 1}},
        SETTOKEN(-sc, --show_comp, 23, 0),
        SETTOKEN(-ac, --save_comp, 24, 2),
        {"--show_save_intag",{25, 0}}, {"--save_intag", {26, 2}}, {"--save_move_intag", {27, 3}},
        {"--save_all_intag", {28, 0}}, {"--save_used_intag", {29, 0}}, {"--show_end_parsed", {30, 1}},
        // SETTOKEN(-q , --quit, 31, 0),
        SETTOKEN(-h , --help, 32, 0),
        SETTOKEN(-ns, --nlin_sol, 33, 0), {"KINSOL", {34, 1}}, {"NEWTON", {35, 2}}, {"RELAX", {36, 3}}, 
        SETTOKENF(-gd, --gen_dir, 37, 1, true),
        // {"--end_after_parse", {38, 0}},
        {"--perform_solve", {39, 1}},
        {"--perform_pcomp", {40, 1}},
    };
#undef SETTOKENF    
#undef SETTOKEN
    std::vector<std::string> new_strs;
    auto read_token = [&]()->std::string{
        if (f.peek() == std::ifstream::traits_type::eof())
            throw std::runtime_error("Expected token but not found in config file = \"" + config_file + "\"");
        f >> std::ws;
        int st_sym = -1;
        while((st_sym = f.peek()) == '#'){
            f.ignore(numeric_limits<streamsize>::max(),'\n');
            f >> std::ws;
        }
        if (st_sym == std::ifstream::traits_type::eof())
            throw std::runtime_error("Expected token but not found in config file = \"" + config_file + "\"");
        std::string val;
        if (st_sym == '\"'){
            char c = '\0';
            f >> c;
            c = '\0';
            while (f >> c && !f.eof() && (c != '\"' || (!val.empty() && val.back() == '\\'))){
                if (c != '\"')
                    val.push_back(c);
                else
                    val.back() = c;   
            }
        } else {
            f >> val;
        }
        return val;
    };
    while ((f >> std::ws) && f.peek() != std::ifstream::traits_type::eof()){
        auto token = read_token();
        new_strs.push_back(token);
        auto it = tokens.find(token);
        if (it == tokens.end())
            throw std::runtime_error("Can't parse token \"" + token + "\" in cofig file = \"" + config_file + "\"");
        auto info = it->second;
        for (unsigned i = 0; i < info.num_args; ++i)
            new_strs.push_back(read_token());
        if (info.num_args > 0 && info.is_path_arg)
            new_strs[new_strs.size() - info.num_args] = (cur_dir / path(new_strs[new_strs.size() - info.num_args])).c_str();
        if (info.id == 0)
            parse_config_file(new_strs.back(), new_strs, new_strs.size()-2, 2);
    }
    if (delete_next > 0)
        args.erase(args.begin() + pos, args.begin() + pos + delete_next);
    args.insert((pos == 0 ? args.begin() :  (args.begin() + pos)), new_strs.begin(), new_strs.end());
}
void Model::printArgsHelpMessage(const std::string& prefix, std::ostream& out){
    out << prefix << "This program is designed to simulate the behavior of a thin shell in a membrane or shell approximation.\n"
        << prefix << "Copyright (C) 2023, the Marchuk Institute of Numerical Mathematics of the Russian Academy of Sciences.\n\n"
        << prefix << " Command line options: " << "\n"
        << prefix << "  -cf, --config     FILE    <Configuration file>\n"
        << prefix << "  -m , --mesh       FILE    <Input mesh file, default=\"mesh.vtk\">" << "\n"
        << prefix << "  -rm, --remove_tag STR STR <Remove tag in input mesh: first STR should be NODE|EDGE|FACE|HALFEDGE, second STR is the name of tag to be deleted>\n"
        << prefix << "  -mv, --rename_tag STR STR STR <Rename tag in input mesh: first STR should be NODE|EDGE|FACE|HALFEDGE, second STR is current name of tag, \n"
        << prefix << "                                   third STR is new name of tag>\n"
        << prefix << "  -e , --energy     FILE    <Specify input energy file, default=\"potential.energy\">" << "\n"
#ifdef WITH_DATADRIVEN
        << prefix << "  -dd, --data_tabs  FILE[+] <Specify energy xi-response data table files, default=\"\">" << "\n"
#endif        
        << prefix << "  -t , --target     PATH    <Directory to save results, default=\"\">" << "\n"
        << prefix << "  -nm, --name       STR     <Prefix for saved results, default=\"res\">" << "\n"
        << prefix << "  -sm, --show_mesh          <Print list of mesh tag names and common mesh information>" << "\n"
        << prefix << "  -se, --show_elast         <Print current elastic potential state>" << "\n"
        << prefix << "  -st, --show_trait         <Print list of supported mesh traits and it's current expected values>" << "\n"
        << prefix << "  -sl, --show_lin           <Print list of available linear solvers>" << "\n"
        << prefix << "  -mt, --model_type STR     <Set type of theory model: MEMBRANE or SHELL, default=MEMBRANE>\n"
        << prefix << "  -ep, --elast_prm  STR VAL <"
#ifdef WITH_DATADRIVEN          
        << "For ANALYTICAL potential only\n" 
        << prefix << "                             " 
#endif          
                                                << "Associate elastic parameter or fiber with name STR to tag name VAL or raw DVAL[] value VAL depending on VAL type>" << "\n"
#ifdef WITH_DATADRIVEN        
        << prefix << "  -di, --dat_interp STR[+]  <For DATADRIVEN potential only\n"
        << prefix << "                             Add specific interpolation rules followed after previously added rules parameters\n"
        << prefix << "                             Syntax is follow: TYPE TRUSTED_REGION TYPE_PARAMETERS\n"
        << prefix << "                             TYPE: KNEAREST | LINEAR; TRUSTED_REGION (DVAL) is radius of domain in xi-domain\n"
        << prefix << "                             TYPE_PARAMETERS: for KNEAREST <k_neighbour inverse_distance_metric_power>,\n"
        << prefix << "                                              for LINEAR <fit_radius octant_mask>, default=KNEAREST 1e20 1 1>" << "\n" 
#endif        
        << prefix << "  -tr, --trait_from STR VAL <Associate mesh trait with name STR to tag name VAL or raw DVAL[*] value VAL depending on VAL type>" << "\n"
        << prefix << "  -vw, --view       BVAL    <Create ImGui debug view window if available>" << "\n"
        << prefix << "  -ns, --nlin_sol   STR[+]  <Add specific solver parameters followed after previous added specific parameters\n"
        << prefix << "                             Syntax is follow: NONLINEAR_SLV_TYPE NONLINEAR_SLV_TYPE_PARAMETERS\n"
        << prefix << "                             NONLINEAR_SLV_TYPE: KINSOL | NEWTON | RELAX\n"
        << prefix << "                             NONLINEAR_SLV_TYPE_PARAMETERS: for KINSOL <maxnits>,\n"
        << prefix << "                                                            for NEWTON <maxnits full_step>\n"
        << prefix << "                                                            for RELAX: <delta relax_its relax_pfreq>, default=KINSOL 200>" << "\n"
        << prefix << "  -ls, --lin_sol    STR     <Set type of linear solver, default=\"inner_mptiluc\">" << "\n"
        << prefix << "  -lp, --lin_prm    STR STR <Set some parameter of linear solver: name_of_parameter value_of_parameter>" << "\n"
        << prefix << "  -nt, --nlin_tol   DVAL[2] <Set nonlinear absolute tolerance and relative tolerance, default=1e-7 1e-7>\n"
        << prefix << "       --bin_vtk    BVAL    <Should we save results in binary vtk file, default=false>\n"
        << prefix << "       --txt_vtk    BVAL    <Should we save results in ascii vtk file, default=false>\n"
        << prefix << "       --bin_tmh    BVAL    <Should we save results in binary tmh file, default=false>\n"
        << prefix << "       --txt_tmh    BVAL    <Should we save results in ascii tmh file, default=false>\n"
        << prefix << "  -sc, --show_comp          <Print list of postcomputation values to save>\n"
        << prefix << "  -ac, --save_comp  STR STR <Add specific postcomputation trait to save>\n"
        << prefix << "  --show_save_intag         <Print list of input tags and their save-status>\n"
        << prefix << "      --save_intag  STR STR <Add input tag to save-list of tags, first STR should be NODE|EDGE|FACE|HALFEDGE, second STR is the name of tag to be saved>\n"
        << prefix << "  --save_move_intag STR STR STR <Add input tag to save-list of tags, first STR should be NODE|EDGE|FACE|HALFEDGE, second STR is the name of tag to be saved, \n"
        << prefix << "                                 third name of the tag to be used on final mesh>\n"
        << prefix << "  --save_all_intag          <Add all input tags to save-list of tags>\n"
        << prefix << "  --save_used_intag         <Add all currently associated input tags to save-list of tags>\n"
        << prefix << "  --show_end_parsed BVAL    <Should we print final parsed state after all argument processing, default=true>\n"
        << prefix << "  -rg, --regen_expr BVAL    <Should program regenerate shared libraries with potential derivatives, default=true>\n"
        << prefix << "  -gd, --gen_dir    PATH    <Directory to save generated objects, default=\"generated\">\n"
        << prefix << "  --perform_solve   BVAL    <Should program run solution of nonlinear system, default=true>\n"
        << prefix << "  --perform_pcomp   BVAL    <Should program perform postcomputations and save results, default=true>\n"
        << prefix << "  -q , --quit               <Exit after parsing this argument>" << "\n"
        << prefix << "  --end_after_parse         <Exit after constructing final parsed state>\n"
        << prefix << "  -h , --help               <Print this message and exit>" << "\n\n"
        << prefix << "Please send messages about detected bugs to the email: al.liogky@yandex.ru\n";
}
int Model::execute(int argc, char** argv){
    using namespace World3d;

    // read mesh geometry
    if (m_mtm.empty_mesh()) throw std::runtime_error("Empty input mesh");
    Mesh _m = m_mtm.convertToMesh();

    // read mesh input required data
    TraitValue trait = m_input_traits[INIT_X_VECTOR_DTAG];
    if (!trait.m_tag_name.empty() && trait.m_tag_name != "self:point"){
         if (!m_mtm.createTagOnMesh<V_ind, Point>(&_m, trait.m_tag_name, "v:point", true))
            throw std::runtime_error("Can't create \"v:point\" tag from meta mesh tag \"" + trait.m_tag_name + "\"");
    } else if (trait.m_tag_name.empty())
        throw std::runtime_error("Trait " + getObjectTraitTypeName(INIT_X_VECTOR_DTAG) + " isn't specified");        
    
    trait = m_input_traits[START_X_VECTOR_DTAG];
    if (trait.m_tag_name == "self:point"){
        if (!m_mtm.createTagOnMeshFromSelfPoints<Point>(&_m, "v:x", false))
            throw std::runtime_error("Can't overwrite \"v:x\" tag from meta mesh self points");
    } else if (!trait.m_tag_name.empty()){
        if (!m_mtm.createTagOnMesh<V_ind, Point>(&_m, trait.m_tag_name, "v:x", true))
            throw std::runtime_error("Can't create \"v:x\" tag from meta mesh tag \"" + trait.m_tag_name + "\"");
    } else 
        throw std::runtime_error("Trait " + getObjectTraitTypeName(START_X_VECTOR_DTAG) + " isn't specified");
    {
        auto next_x = _m.add_property_map<V_ind, Point>("v:next_x").first;
        auto x = _m.property_map<V_ind, Point>("v:x").first;
        for (auto v: _m.vertices()) next_x[v] = x[v];
    }
    trait = m_input_traits[BOUNDARY_ITAG];
    if (!trait.m_tag_name.empty()){
        if (!m_mtm.createTagOnMesh<V_ind, int>(&_m, trait.m_tag_name, "v:boundary_lbl", true))
            throw std::runtime_error("Can't create \"v:boundary_lbl\" tag from meta mesh tag \"" + trait.m_tag_name + "\"");
    } else if (!std::isnan(trait.m_def_value[0])){
        if (!_m.add_property_map<V_ind, int>("v:boundary_lbl", static_cast<int>(trait.m_def_value[0])).second)
            throw std::runtime_error("Can't create new \"v:boundary_lbl\" tag");
    } else 
        throw std::runtime_error("Trait " + getObjectTraitTypeName(BOUNDARY_ITAG) + " isn't specified");

    // create object
    Object3D _obj{std::move(_m)};
    _obj.setName(path(fmesh).stem().string());
    World w;
    ObjectID oid = w.addObject3D(std::move(_obj), 1);
    Object3D& obj = w.obj(oid);
    Mesh& m = obj.m_mesh;
    obj.setDirichletCond([](Object3D& obj, const V_ind& v) -> unsigned char{ return (obj.m_boundary[v] < 0) ? 7 : (7 & (~obj.m_boundary[v])); });

    // set hyperelastic force
    using HEM = HyperElasticModel;
    HEM::ParamSource h;
    trait = m_input_traits[THICKNESS_DTAG];
    if (!trait.m_tag_name.empty()){
        if (!m_mtm.createTagOnMesh<F_ind, DReal>(&m, trait.m_tag_name, "f:thickness", true))
            throw std::runtime_error("Can't create \"f:thickness\" tag from meta mesh tag \"" + trait.m_tag_name + "\"");
        h = HEM::ParamSource(HEM::ParamSource::TAG_T, "f:thickness", NAN);
    } else if (!std::isnan(trait.m_def_value[0]))   
        h = HEM::ParamSource(HEM::ParamSource::VAL_T, "", trait.m_def_value[0]); 
    else
        throw std::runtime_error("Trait " + getObjectTraitTypeName(THICKNESS_DTAG) + " isn't specified");    
    
    Force_ID ef_id = -1, bf_id = -1;
    if (m_potential.m_type == ANALYTICAL){
        auto a = m_potential.analytical();
        std::vector<std::pair<casadi::SX, HEM::ParamSource>> eprm;
        for (auto& x: a.m_energy_expr.getParameters()){
            std::pair<casadi::SX, HEM::ParamSource> p;
            p.first = casadi::SX::sym(x.name);
            if (!x.def_tag.empty()){
                if (!m_mtm.createTagOnMesh<F_ind, DReal>(&m, x.def_tag, "f:prm:" + x.name, true))
                    throw std::runtime_error("Can't create \"" + x.name + "\" tag from meta mesh tag \"" + x.def_tag + "\"");
                p.second = HEM::ParamSource(HEM::ParamSource::TAG_T, "f:prm:" + x.name, NAN);
            } else if (!std::isnan(x.def_val)){
                p.second = HEM::ParamSource(HEM::ParamSource::VAL_T, "", x.def_val);
            } else  
                throw std::runtime_error("Elastic parameter \"" + x.name + "\" isn't specified");
            eprm.emplace_back(std::move(p));    
        }
        const auto& efibs = a.m_energy_expr.getFibers();
        if (!efibs.empty()){
            // specify fiber data
            auto _n0 = set_n0(obj);
            for (auto f: efibs){
                if (!f.def_tag.empty()){
                    if (!m_mtm.createTagOnMesh<F_ind, std::array<DReal, 3>>(&m, f.def_tag, "f:fib:" + f.name, true))
                        throw std::runtime_error("Can't create \"" + ("f:fib:" + f.name) + "\" tag from meta mesh tag \"" + f.def_tag + "\""); 
                    auto p = m.property_map<F_ind, std::array<DReal, 3>>("f:fib:" + f.name).first;
                    bool have_out_surface_vec = false;
                    for (auto f: m.faces()){
                        auto d = p[f];
                        if (abs(d[0]) < 1e-14 && abs(d[1]) < 1e-14 && abs(d[2]) < 1e-14){
                            p[f] = std::array<DReal, 3>{0, 0, 0};
                            continue;
                        }
                        auto d_len = std::hypot(d[0], d[1], d[2]);
                        for (int k = 0; k < 3; ++k) d[k] /= d_len;
                        auto scl = Vector(d[0], d[1], d[2]) * _n0[f];
                        // project data if required
                        if (abs(scl) > 100 * std::numeric_limits<DReal>::epsilon())
                            have_out_surface_vec = true;   
                        if (1 - abs(scl) < 100 * std::numeric_limits<DReal>::epsilon()){
                            p[f] = std::array<DReal, 3>{0, 0, 0};
                        } else {
                            Vector r = Vector(d[0], d[1], d[2]) - scl * _n0[f];
                            r /= sqrt(r.squared_length());
                            p[f] = std::array<DReal, 3>{r[0], r[1], r[2]};
                        }
                    }
                    if (have_out_surface_vec)
                        std::cout << "Warning: Fiber from meta tag \"" << f.def_tag << "\" had not in-initial-surface directions and was projected\n";       
                } else 
                    throw std::runtime_error("Elastic fiber \"" + f.name + "\" isn't specified"); 
            }
            if (_n0.second) m.remove_property_map(_n0.first);
        }
        std::vector<std::pair<casadi::SX, HEM::InvariantDescr>> einvs;
        for (auto x: a.m_energy_expr.getInvariants()){
            switch(x.type){
                case Energy_Parser::Invariants::I1_T:{
                    einvs.push_back({casadi::SX::sym("I[1]"), {HEM::InvariantDescr::I1_T}}); break;
                }
                case Energy_Parser::Invariants::I2_T:{
                    einvs.push_back({casadi::SX::sym("I[2]"), {HEM::InvariantDescr::I2_T}}); break;
                }
                case Energy_Parser::Invariants::J_T:{
                    einvs.push_back({casadi::SX::sym("J[]"), {HEM::InvariantDescr::J_T}}); break;
                }
                case Energy_Parser::Invariants::If_T:{
                    einvs.push_back({casadi::SX::sym("I[" + efibs[x.v1].name + "]"), {HEM::InvariantDescr::If_T, "f:fib:" + efibs[x.v1].name}}); break;
                }
                default:{ //Ifs_T,
                    einvs.push_back({casadi::SX::sym("I[" + efibs[x.v1].name + "," + efibs[x.v2].name + "]"), {HEM::InvariantDescr::Ifs_T, "f:fib:" + efibs[x.v1].name, "f:fib:" + efibs[x.v2].name}}); break;
                }
            }
        }
        casadi::SX epotential;
        {
            std::vector<casadi::SX> prm(eprm.size()), inv(einvs.size());
            for (std::size_t i = 0; i < eprm.size(); ++i) prm[i] = eprm[i].first;
            for (std::size_t i = 0; i < einvs.size(); ++i) inv[i] = einvs[i].first;
            epotential = casadi::SX::simplify(a.m_energy_expr.Evaluate(prm, inv));
        }
        
        Force elast_f = HyperElasticModel(h, eprm, einvs, epotential, to_gen, regen_potential, path(a.fpotential).stem().c_str());
        elast_f.target<HyperElasticForceBase>()->prepareJacobianFunction(regen_potential); elast_f.target<HyperElasticForceBase>()->f.setVerbosity(2);
        ef_id = w.addForce(std::move(elast_f), oid);
        if (!is_membrane_approx){
            // prepare bending force
            Force bend_f = BendingForce(obj.m_forces[ef_id].target<HyperElasticForceBase>()->f, regen_potential);
            bend_f.target<BendingForce>()->prepareJacobianFunction(regen_potential); bend_f.target<BendingForce>()->setVerbosity(4);
            trait = m_input_traits[CLAMPED_EDGE_VECTOR_DTAG];
            //prepare shell boundary marks (FREE or CLAMPED)
            std::map<E_ind, BendingForce::BoundaryData::Entry> boundary;
            if (!trait.m_tag_name.empty()){
                auto rit = m_mtm.edge_data.find(trait.m_tag_name);
                if (rit == m_mtm.edge_data.end()) throw std::runtime_error("Can't find specified trait " + getObjectTraitTypeName(CLAMPED_EDGE_VECTOR_DTAG));
                auto& r = rit->second;
                for (auto e: m.edges()){
                    auto f = face_around(m, e);
                    if (f[1].second) continue;
                    auto v = vert_around(m, e);
                    if (!((obj.m_boundary[v[0]] & 0xF) && (obj.m_boundary[v[1]] & 0xF))){
                        boundary[e] = BendingForce::BoundaryData::Entry(BendingForce::BoundaryData::FREE);
                    } else {
                        std::array<double, 3> res;
                        r.get_element<std::array<double, 3>>(e.idx(), res);
                        auto len = sqrt(res[0]*res[0] + res[1]*res[1] + res[2]*res[2]);
                        if (len > 1e-5){
                            for (int k = 0; k < 3; ++k) res[k] /= len;
                            boundary[e] = BendingForce::BoundaryData::Entry(BendingForce::BoundaryData::CLAMPED, res);
                        } else
                            boundary[e] = BendingForce::BoundaryData::Entry(BendingForce::BoundaryData::FREE);
                    } 
                }
            } else if (!std::isnan(trait.m_def_value[0]) && !std::isnan(trait.m_def_value[1]) && !std::isnan(trait.m_def_value[2])){
                std::array<double, 3> res;
                auto len = sqrt(res[0]*res[0] + res[1]*res[1] + res[2]*res[2]);
                bool is_zero = (len < 1e-5);
                if (!is_zero) for (int k = 0; k < 3; ++k) res[k] /= len;
                for (auto e: m.edges()){
                    auto f = face_around(m, e);
                    if (f[1].second) continue;
                    auto v = vert_around(m, e);
                    if (!((obj.m_boundary[v[0]] & 0xF) && (obj.m_boundary[v[1]] & 0xF))){
                        boundary[e] = BendingForce::BoundaryData::Entry(BendingForce::BoundaryData::FREE);
                    } else if (is_zero){
                            boundary[e] = BendingForce::BoundaryData::Entry(BendingForce::BoundaryData::CLAMPED, res);
                    } else
                        boundary[e] = BendingForce::BoundaryData::Entry(BendingForce::BoundaryData::FREE);
                }
            } 
            // else {
            //     for (auto e: m.edges()){
            //         auto f = face_around(m, e);
            //         if (f[1].second) continue;
            //         boundary[e] = BendingForce::BoundaryData::Entry(BendingForce::BoundaryData::FREE);
            //     }
            // }
            if (!boundary.empty())
                bend_f.target<BendingForce>()->set_boundary_data(obj, boundary);
            bf_id = w.addForce(std::move(bend_f), oid);
        }
    } else if (m_potential.m_type == DATADRIVEN){
#ifdef WITH_DATADRIVEN        
        auto a = m_potential.datadriven();
        //read data from tables
        RegionalResponseTable com_tab;
        for (auto f: a.data_files){
            RegionalResponseTable tmp;
            tmp.read(f);
            tmp.squashAll();
            com_tab.append(tmp);
        }
        com_tab.squashAll();
        ResponseInterpFactory interps;
        std::cout << "Dataset contains " << com_tab.cregion(0).size() << " unique points" << std::endl;
        using DDI = DataDrivenPotentialInfo;
        // create interpolation pipeline
        for (std::size_t i = 0; i < a.interp.size(); ++i){
            auto& p = a.interp[i];
            switch(p.m_tp){
                case DDI::KNEAREST:{
                    auto kd = p.knearest();
                    ResponseKNearest kns(kd.k);
                    kns.setLocalInterpolator(std::make_unique<ResponseLocKInvDist>(kd.metric_pow));
                    auto shell = makeRespConstraintShell(kns, makeRespSphericalConstraint(kd.R_trust));
                    interps.append(std::make_unique<decltype(shell)>(std::move(shell)));
                    interps.m_factory.back()->fit(com_tab.cregion(0));
                    break;
                }
                case DDI::LINEAR:{
                    auto kd = p.linear();
                    LinearElasticResponse le;
                    auto constr = [oct = kd.octant_mask, r = kd.R_trust](const std::array<double, 3>& x)->bool{
                        char mask = (x[0] > 0 ? 0 : 1) + (x[1] > 0 ? 0 : 2) + (x[2] > 0 ? 0 : 4);
                        if (oct & (1 << mask)) return true; 
                        return (x[0]*x[0] + x[1]*x[1] + x[2]*x[2] < r*r);
                    };
                    auto shell = makeRespConstraintShell(std::move(le), constr);
                    interps.append(std::make_unique<decltype(shell)>(std::move(shell)));
                    interps.m_factory.back()->fit(RegionalResponseTable::filter(com_tab, makeRespSphericalConstraint(kd.R_fit)).cregion(0));
                    break;
                }
                default:
                    throw std::runtime_error("Faced unknown DATADRIVEN interpolant type");
            }
        }
        // create elastic force from interpolant
        DataDrivenMaterial ddm;
        ddm.setInterpolant(std::make_unique<decltype(interps)>(std::move(interps)));
        HEForce elast_f;
        elast_f.setMaterial(std::move(ddm));
        auto H_ = h.makeDataRequier();
        H_->registerObj(&obj);
        elast_f.setThickness(ThicknessFromDataRequier(std::move(H_)));
        ef_id = w.addForce(std::move(elast_f), oid);
        if (!m_perform_solve)
            obj.m_forces[ef_id].target<HEForce>()->m_mat->registerObj(&obj);
        if (!is_membrane_approx)
            throw std::runtime_error("DATADRIVEN potential currently does not support shell mechanical model");    
#else
        throw std::runtime_error("Atteption to use DATADRIVEN potential in version of model without DATADRIVEN potential support");
#endif        
    }

    // set pressure force
    trait = m_input_traits[PRESSURE_DTAG];
    Force_ID pr_id = -1;
    if (!trait.m_tag_name.empty()){
        if (!m_mtm.createTagOnMesh<F_ind, DReal>(&m, trait.m_tag_name, "f:pressure", true))
            throw std::runtime_error("Can't create \"f:pressure\" tag from meta mesh tag \"" + trait.m_tag_name + "\"");
        auto ptag = m.property_map<F_ind, DReal>("f:pressure").first;
        pr_id = w.addForce(SimplePressureLoad([ptag](Object3D& obj, F_ind f)->DReal{ return ptag[f]; }), oid);
    } else if (!std::isnan(trait.m_def_value[0]))
        pr_id = w.addForce(SimplePressureLoad(trait.m_def_value[0]), oid);
    
    // set pressure force acting on initial configuration
    trait = m_input_traits[LAGRANGE_PRESSURE_DTAG];
    Force_ID lpr_id = -1;
    if (!trait.m_tag_name.empty()){
        if (!m_mtm.createTagOnMesh<F_ind, DReal>(&m, trait.m_tag_name, "f:lagr_pressure", true))
            throw std::runtime_error("Can't create \"f:lagr_pressure\" tag from meta mesh tag \"" + trait.m_tag_name + "\"");
        auto ptag = m.property_map<F_ind, DReal>("f:lagr_pressure").first;
        lpr_id = w.addForce(SimpleLagrangePressureLoad([ptag](Object3D& obj, F_ind f)->DReal{ return ptag[f]; }), oid);
    } else if (!std::isnan(trait.m_def_value[0]))
        lpr_id = w.addForce(SimpleLagrangePressureLoad(trait.m_def_value[0]), oid);

    // set body force
    trait = m_input_traits[BODY_FORCE_VECTOR_DTAG];
    Force_ID body_id = -1;
    if (!trait.m_tag_name.empty()){
        if (!m_mtm.createTagOnMesh<V_ind, Vector>(&m, trait.m_tag_name, "v:body_force_density", true))
            throw std::runtime_error("Can't create \"v:body_force_density\" tag from meta mesh tag \"" + trait.m_tag_name + "\"");
        auto ptag = m.property_map<V_ind, Vector>("v:body_force_density").first;
        std::function<bool(double *x/*[3]*/, double *F/*[3]*/, V_ind& id)> body_vector = [ptag](double *x, double *F, V_ind& id){ 
            F[0] = ptag[id][0], F[1] = ptag[id][1], F[2] = ptag[id][2];
            return true; 
        };
        std::function<bool(double *x/*[3]*/, double *J/*[3][3]*/, V_ind& id)> body_deriv  = [](double *x, double *J, V_ind& id){
            std::fill(J, J+9, 0.0);
            return true;
        };
        std::vector<std::pair<V_ind, VertexLoad<V_ind>::VertParam>> verts;
        auto _Ap = set_Ap(obj);
        for (auto v: m.vertices()) if (ptag[v].squared_length() != 0.0){
            auto f = face_around(m, v);
            double A = 0;
            for (auto ff: f)
                A += _Ap[ff];
            A /= 3;    
            verts.push_back({v, VertexLoad<V_ind>::VertParam{A, v}});  
        } 
        if (_Ap.second) obj.m_mesh.remove_property_map(_Ap.first); 
        body_id = w.addForce(VertexLoad<V_ind>().setForceField(body_vector, body_deriv).insertVertices(verts), oid);
    } else if ( !std::isnan(trait.m_def_value[0]) && !std::isnan(trait.m_def_value[1]) && !std::isnan(trait.m_def_value[2]) 
            && (!(trait.m_def_value[0] == 0.0 && trait.m_def_value[1] == 0.0 && trait.m_def_value[2] == 0.0)) ){
        std::function<bool(double *x/*[3]*/, double *F/*[3]*/)> body_vector = [v = trait.m_def_value](double *x, double *F){ return std::copy(v.begin(), v.end(), F), true; };
        std::function<bool(double *x/*[3]*/, double *J/*[3][3]*/)> body_deriv  = [](double *x, double *J){ return std::fill(J, J+9, 0.0), true; };
        std::vector<std::pair<V_ind, VertexLoad<>::VertParam>> verts;
        auto _Ap = set_Ap(obj);
        for (auto v: m.vertices()){
            auto f = face_around(m, v);
            double A = 0;
            for (auto ff: f)
                A += _Ap[ff];
            verts.push_back({v, VertexLoad<>::VertParam{A/3}});  
        } 
        if (_Ap.second) obj.m_mesh.remove_property_map(_Ap.first); 
        body_id = w.addForce(VertexLoad<>().setForceField(body_vector, body_deriv).insertVertices(verts), oid);
    }

    // set edge force
    trait = m_input_traits[EDGE_FORCE_VECTOR_DTAG];
    Force_ID edge_body_id = -1;
    if (!trait.m_tag_name.empty()){
        if (!m_mtm.createTagOnMesh<E_ind, Vector>(&m, trait.m_tag_name, "e:edge_body_density", true))
            throw std::runtime_error("Can't create \"e:edge_body_density\" tag from meta mesh tag \"" + trait.m_tag_name + "\"");
        auto ptag = m.property_map<E_ind, Vector>("e:edge_body_density").first;
        edge_body_id = w.addForce(BodyEdgeLoad([ptag](Object3D& obj, E_ind e)->Vector{ return ptag[e]; }), oid);
    } else if (!std::isnan(trait.m_def_value[0]) && !std::isnan(trait.m_def_value[1]) && !std::isnan(trait.m_def_value[2])
             && (!(trait.m_def_value[0] == 0.0 && trait.m_def_value[1] == 0.0 && trait.m_def_value[2] == 0.0)) ){
        edge_body_id = w.addForce(BodyEdgeLoad(Vector(trait.m_def_value[0], trait.m_def_value[1], trait.m_def_value[1])), oid);
    }

    //create gui vindow if available
    std::thread* view_window = nullptr;
    if (m_view && m_perform_solve){
        view_window = new std::thread(start_debug_gui, argc, argv);
        w.setRenderer(std::make_unique<World3d::DefaultRenderer>());
    }

    int status = 0;
    obj.apply_updaters();
    if (m_perform_solve){
        //init solvers
        w.UpdateForces();
        double init_abs_residual = w.getWorldNorm("v:force");
        double abs_tol = std::max(m_abs_tol, init_abs_residual*m_rel_tol);

        NSWorldWrapper nsww(w);
        NLProblem nlp = NSWorldWrapper::makeProblem(nsww);
        LinearSolver ls(m_lin_sol_name);
        ls.setVerbosityLevel(1);
        for (auto& v: m_lin_sol_prms)
            ls.SetParameter(v.first, v.second);
        ls.SetParameterReal("relative_tolerance", 1e-20);
        auto monitor_fcn = [&nsww](const double *x) {
            nsww.setX(x);
            nsww.RenderIteration();
            return 0;
        };

        //initialize KINSOL solver
        NonLinearSolverKinsol nlsp(nlp);
        nlsp.SetLinearSolver(&ls);
        static_cast<SUNLinearSolverContent_Com>(static_cast<SUNLinSolCustom>(nlsp.GetLinearSolver()->content)->content)->verbosity = 1;
        nlsp.SetVerbosityLevel(1);
        nlsp.SetInitialGuess(nsww.getCurrentX());
        nlsp.SetFuncNormTol(abs_tol);
        nlsp.SetMaxNewtonStep(1e20);
        nlsp.SetParameterInt("MaxSetupCalls", 1);
        nlsp.SetParameterInt("MaxSubSetupCalls", 1);
        nlsp.SetScaledStepTol(1e-12);
        nlsp.SetMaxBetaFails(40);
        nlsp.SetSolveStrategy(NonLinearSolverKinsol::LINESEARCH);
        World3d::Timer sym_time, com_time;
        com_time.reset();
        nlsp.SetInfoHandlerFn([&time = sym_time, &com_time](const char *module, const char *function, char *msg) {
            std::cout << "[" << module << "] " << function << "\n   " << msg << "\n" << "time = " << time.elapsed()
                        << " com_time = " << com_time.elapsed() << std::endl;
        });
        nlsp.setInterIterationMonitorFunc(monitor_fcn);
        auto result_str_nlsp = [](NonLinearSolverKinsol &nlsp) {
            std::stringstream sstr;
            sstr << "\t#NonLinIts = " << nlsp.GetNumNolinSolvIters() << " Residual = " << nlsp.GetResidualNorm()
                    << "\n\t#linIts = " << nlsp.GetNumLinIters() << " #funcEvals = " << nlsp.GetNumFuncEvals()
                    << " #jacEvals = " << nlsp.GetNumJacEvals()
                    << "\n\t#convFails = " << nlsp.GetNumLinConvFails() << " #betaCondFails = "
                    << nlsp.GetNumBetaCondFails() << " #backtrackOps = " << nlsp.GetNumBacktrackOps() << "\n";
            sstr << "Reason: " << nlsp.GetReason();
            return sstr.str();
        };

        //initialize custom full_step NEWTON solver
        NonLinearSolverCustom nlsc;
        nlsc.setProblem(nlp);
        nlsc.setInterIterationMonitorFunc(monitor_fcn);
        nlsc.SetLinearSolver(ls);
        nlsc.SetInitialGuess(nsww.getCurrentX());
        nlsc.SetNSWorldWrapper(nsww);
        double nlsc_step = 1.0;
        int nlsc_maxnits = 50;
        nlsc.tau_algo = [&nlsc_step](int it, double rel_resid, double abs_resid, NonLinearSolverCustom& ns)->double{
            std::cout << "\tit = " << it << ": norm = " << abs_resid << " rel = " << rel_resid << std::endl;
            return nlsc_step;
        };
        nlsc.stop_cond = [&nlsc_maxnits, abs_tol](int it, double rel, double abs, int& reason, NonLinearSolverCustom& ns){
            if (abs < abs_tol) {
                return reason = 0, true;
            } else if (abs > 1e20 || std::isnan(abs)) {
                return reason = -1, true;
            }
            if (it >= nlsc_maxnits) {
                return reason = 2, true;
            }
            return reason = 1, false;
        };
        auto result_str_nlsc = [](NonLinearSolverCustom &nlsc) {
            std::stringstream sstr;
            sstr << "\t#NonLinIts = " << nlsc.GetNumNolinSolvIters() << " Residual = " << nlsc.GetResidualNorm()
                    << "\n\t#linIts = " << nlsc.GetNumLinIters() << " #funcEvals = " << nlsc.GetNumFuncEvals()
                    << " #jacEvals = " << nlsc.GetNumJacEvals() <<  "\n";
            sstr << "Reason: " << nlsc.GetReason();
            return sstr.str();
        };

        //initialize relaxation solver
        double Ht = 1e40;
        for (auto e: m.edges()){
            auto v = vert_around(m, e);
            auto len = (obj.m_x0[v[0]] - obj.m_x0[v[1]]).squared_length();
            if (Ht > len) Ht = len;
        }
        Ht = sqrt(Ht);
        auto relaxate_solver = [Ht, &w, &sym_time, oid, abs_tol](double &delta, int nits, int pfreq)->std::pair<int, int> {
            if (nits <= 0) return {false, 0};
            // std::cout << "delta = " << delta << "\n";
            static double diverge_lim = 1e15;
            double err = abs_tol;
            auto stopCond = World3d::StaticStopCondition(nullptr, &diverge_lim, &nits, &pfreq, &sym_time, nullptr, true);
            stopCond.setAbsErr(err);
            StaticForceApplierP sfa(delta);
            sfa.addDamper(Damper(0.05 * Ht, 0.9, 15, 1000000).setInitDxScale(3.0)
            .setReInitDxParams(-1, -1).setMinDeltaFactor(1e-5));
            w.setForceApplier(sfa);
            sym_time.reset();
            int status = w.Simulation([&stopCond](StepSimInfo& s, World* w)->bool{ return stopCond(s, w); });
            // double factor = w.obj(oid)._next_x_updater.target<StaticForceApplierP>()->_dampers[0].target<Damper>()->getUsedDeltaFactor();
            // delta *= factor;
            // if (1 - factor > 1e-2) std::cout << "New relaxation delta = " << delta << std::endl;
            return {status, stopCond.last_ssi.it};
        };
        int func_evals = 0, jac_evals = 0, accum_func_evals = 0, accum_jac_evals = 0;
        auto print_iterate_info = [&fe = func_evals, &je = jac_evals, &afe = accum_func_evals, &aje = accum_jac_evals, &com_time](int it){
            std::cout << "ITERATE INFO: it: " << it << " com_time = " << com_time.elapsed() <<
                        " func_evals = " << fe << " jac_evals = " << je <<
                        " accum_func_evals = " << afe << " accum_jac_evals = " << aje << std::endl;
        };
        //perform solution of nonlinear system
        bool slvFlag = (init_abs_residual <= abs_tol);
        for (std::size_t it = 0; it < m_solverItData.size() && !slvFlag; ++it){
            auto csi = m_solverItData[it];
            switch(csi.ns_type) {
                case 0: { //RELAX
                    auto s = relaxate_solver(csi.delta_full_step, csi.relax_its, csi.relax_pfreq);
                    func_evals += s.second;
                    if (s.first == 0) slvFlag = true;
                    if (s.first <= -10) it = m_solverItData.size();
                    break;
                }
                case 1: { //NEWTON
                    bool update_ls_abs_tol = false;
                    for (auto& v: m_lin_sol_prms)
                        if (v.first == "absolute_tolerance"){
                            ls.SetParameter(v.first, v.second);
                            update_ls_abs_tol = true;
                        }
                    if (!update_ls_abs_tol)  
                        ls.SetParameterReal("absolute_tolerance", 1e-9);
                    nlsc_step = csi.delta_full_step;
                    nlsc_maxnits = csi.max_ns_its;
                    nlsc.SetInitialGuess(nsww.getCurrentX());
                    slvFlag = nlsc.Solve();
                    slvFlag  = (nlsc.getReasonFlag() == 0);
                    std::cout << result_str_nlsc(nlsc) << std::endl;
                    func_evals += nlsc.GetNumFuncEvals(), jac_evals += nlsc.GetNumJacEvals();
                    break;
                }
                default: { //KINSOL
                    nlsp.SetInitialGuess(nsww.getCurrentX());
                    nlsp.SetNumMaxIters(csi.max_ns_its);
                    slvFlag = nlsp.Solve();
                    std::cout << result_str_nlsp(nlsp) << std::endl;
                    func_evals += nlsp.GetNumFuncEvals(), jac_evals += nlsp.GetNumJacEvals();
                    break;
                }
            }
            accum_func_evals += func_evals, accum_jac_evals += jac_evals;
            print_iterate_info(it);
            func_evals = jac_evals = 0;
        }
        if (!slvFlag) status = -1;

        std::cout << (slvFlag ? "Convergence is reached\n" : "Convergence is NOT reached\n") << std::endl;
    }
    
    if ( m_perform_postcomps && (save_flag & (SAVE_BIN_VTK|SAVE_TXT_VTK|SAVE_BIN_TMH|SAVE_TXT_TMH)) ){
        //save intags
        using ET = MetaTriMesh::ElementSparsity;
        for (std::size_t i = ET::NODE; i <= ET::HALFEDGE; ++i){
            auto& vd = m_mtm.data_array(static_cast<ET>(i));
            if (vd.empty()) continue;
            auto jt = m_save_intag.find(static_cast<ET>(i));
            if (jt == m_save_intag.end()){
                vd.clear();
                continue;
            }
            const auto& sd = jt->second;
            for (auto it = vd.begin(); it != vd.end(); ){
                auto kt = sd.find(it->first);
                if (kt == sd.end())
                    it = vd.erase(it);
                else 
                    ++it;    
            }
            std::map<std::string, MetaTriMesh::MeshValue> vals;
            for (auto it = sd.begin(); it != sd.end(); ++it){
                auto kt = vd.find(it->first);
                if (kt == vd.end()) continue;
                vals.insert({it->second.empty() ? it->first : it->second, std::move(kt->second)});
                vd.erase(kt);
            }
            vd = std::move(vals);
        }

        auto strait = m_out_traits[SAVE_FINAL_STATE];
        if (!strait.m_tag_name.empty())
            m_mtm.readTagDataFromMesh<V_ind>(m, "v:x", strait.m_tag_name);
        strait = m_out_traits[SAVE_RESIDUAL];
        if (!strait.m_tag_name.empty()){
            if (!m_perform_solve) w.UpdateForces();
            m_mtm.readTagDataFromMesh<V_ind>(m, "v:force", strait.m_tag_name);
        }
        strait = m_out_traits[SAVE_LOCAL_BASIS]; 
        auto strait1 = m_out_traits[SAVE_LOCAL_BASIS_EXT];  
        if (!strait.m_tag_name.empty() || !strait1.m_tag_name.empty()){
            bool flat = check_flat_initial_template(obj);
            LocalCartesian S(obj, flat);
            auto r1 = m.add_property_map<F_ind, Vector>("f:basis_1").first;
            auto r2 = m.add_property_map<F_ind, Vector>("f:basis_2").first;
            for (auto f: m.faces()){
                auto q = S(f);
                r1[f] = Vector(q(0,0), q(0,1), q(0,2));
                r2[f] = Vector(q(1,0), q(1,1), q(1,2));
            }
            if (!strait.m_tag_name.empty())
                m_mtm.readTagDataFromMesh<F_ind>(m, "f:basis_1", strait.m_tag_name);
            if (!strait1.m_tag_name.empty())
                m_mtm.readTagDataFromMesh<F_ind>(m, "f:basis_2", strait1.m_tag_name);    
        }
        strait = m_out_traits[SAVE_DEF_GRAD];
        strait1 = m_out_traits[SAVE_CAUCHY_STRAIN];
        auto strait2 = m_out_traits[SAVE_CAUCHY_STRESS3D];
        auto strait3 = m_out_traits[SAVE_CAUCHY_TENSION3D];
        bool comp_stress_3d = !strait2.m_tag_name.empty() || !strait3.m_tag_name.empty();
        bool prepare_strains = !strait.m_tag_name.empty() || !strait1.m_tag_name.empty() || comp_stress_3d || 
                    (m_potential.m_type == DATADRIVEN && !m_out_traits[SAVE_PK2_STRESS].m_tag_name.empty());
        if (prepare_strains){
            auto F_ = m.add_property_map<F_ind, std::array<double, 6>>("f:tensor:F").first;
            auto C2d_ = m.add_property_map<F_ind, std::array<double, 3>>("f:tensor:C2d").first;

            auto L_ = set_L(obj);
            auto x_ = obj.m_x;
            auto p_ = obj.m_x0;
            for (auto f: m.faces()){
                auto v = vert_around(m, f);
                Eigen::Matrix<DReal, 3, 3> P, Q;
                Q <<    x_[v[0]][0], x_[v[1]][0], x_[v[2]][0],
                        x_[v[0]][1], x_[v[1]][1], x_[v[2]][1],
                        x_[v[0]][2], x_[v[1]][2], x_[v[2]][2];
                P <<   p_[v[0]][0], p_[v[1]][0], p_[v[2]][0],
                        p_[v[0]][1], p_[v[1]][1], p_[v[2]][1],
                        p_[v[0]][2], p_[v[1]][2], p_[v[2]][2];
                auto L = Eigen::Map<Eigen::Matrix<DReal, 2, 3>>(L_[f].data());
                Eigen::Matrix<DReal, 3, 2> F2d = Q * L.transpose();
                F_[f] = std::array<double, 6>{F2d(0,0), F2d(1,0), F2d(2,0), F2d(0,1), F2d(1,1), F2d(2,1)};
                Eigen::Matrix<DReal, 2, 2> C2d = F2d.transpose() * F2d;     
                C2d_[f] = std::array<double, 3>{C2d(0,0), C2d(1,1), C2d(0,1)};
            }
            if (L_.second) m.remove_property_map(L_.first);
            if (!strait.m_tag_name.empty())
                m_mtm.readTagDataFromMesh<F_ind>(m, "f:tensor:F", strait.m_tag_name);
            if (!strait1.m_tag_name.empty())
                m_mtm.readTagDataFromMesh<F_ind>(m, "f:tensor:C2d", strait1.m_tag_name);    
        }
        strait = m_out_traits[SAVE_PK2_STRESS];
        bool prepare_pk2_stress = !strait.m_tag_name.empty() || comp_stress_3d;
        if (prepare_pk2_stress){
            auto S_ = m.add_property_map<F_ind, std::array<double, 3>>("f:tensor:S").first;
            if (m_potential.m_type == ANALYTICAL){
                auto watch_num = obj.m_forces[ef_id].target<HyperElasticForceBase>()->f.add_watch(DefaultHyperElasticForce::compute_S_tensor_expr, false, true);
                for (auto f: m.faces()){
                    auto Stensor = obj.m_forces[ef_id].target<HyperElasticForceBase>()->f.eval_watch(obj, f, watch_num);
                    S_[f] = std::array<double, 3>{Stensor(0,0), Stensor(1,1), Stensor(0,1)};
                }
                auto gen_dir = obj.m_forces[ef_id].target<HyperElasticForceBase>()->f.gen_dir;
                auto watch_name = obj.m_forces[ef_id].target<HyperElasticForceBase>()->f.m_watches.back().func->name();
                remove(path(gen_dir)/path(watch_name).concat(".c"));
                remove(path(gen_dir)/path(watch_name).concat(".so"));
            } 
#ifdef WITH_DATADRIVEN
            else if (m_potential.m_type == DATADRIVEN){
                auto& he = *(obj.m_forces[ef_id].target<HEForce>()->m_mat);
                auto C2d_ = m.property_map<F_ind, std::array<double, 3>>("f:tensor:C2d").first;
                for (auto f: obj.m_mesh.faces())
                    S_[f] = he.PK2_tensor({(C2d_[f][0]-1)/2, (C2d_[f][1]-1)/2, (C2d_[f][2])/2}, f);
            }
#endif
            if (!strait.m_tag_name.empty())
                m_mtm.readTagDataFromMesh<F_ind>(m, "f:tensor:S", strait.m_tag_name);
        }
        if (comp_stress_3d){
            auto s_ = m.add_property_map<F_ind, std::array<double, 6>>("f:tensor:sigma").first;
            auto S_ = m.property_map<F_ind, std::array<double, 3>>("f:tensor:S").first;
            auto F_ = m.property_map<F_ind, std::array<double, 6>>("f:tensor:F").first;
            for (auto f: m.faces()){
                auto F = F_[f];
                auto S = S_[f];
                Eigen::Matrix<double, 3, 2> F2d;
                F2d(0,0)=F[0], F2d(1,0)=F[1], F2d(2,0)=F[2],
                F2d(0,1)=F[3], F2d(1,1)=F[4], F2d(2,1)=F[5];
                Eigen::Matrix<double, 2, 2> S2d;
                S2d(0,0)=S[0], S2d(1,0)=S2d(0,1)=S[2], S2d(1,1)=S[1];
                Eigen::Matrix<double, 3, 3> s = F2d * S2d * F2d.transpose();  
                s_[f] = std::array<double, 6>{s(0,0), s(1,1), s(2,2), s(0,1), s(0,2), s(1,2)};
            }
            if (!strait2.m_tag_name.empty())
                m_mtm.readTagDataFromMesh<F_ind>(m, "f:tensor:sigma", strait2.m_tag_name);
        }
        if (!strait3.m_tag_name.empty()){
            auto T_ = m.add_property_map<F_ind, std::array<double, 6>>("f:tensor:T").first;
            auto s_ = m.property_map<F_ind, std::array<double, 6>>("f:tensor:sigma").first;    
            auto Sq_ = set_S(obj);
            auto Ap_ = set_Ap(obj);
            auto H_ = h.makeDataRequier();
            H_->registerObj(&obj);
            for (auto f: m.faces()){
                double Aq = sqrt(Sq_[f].squared_length()), Ap = Ap_[f];
                double Jinv = Ap / Aq;   
                std::array<double, 6> s = s_[f];
                const double* Hp = nullptr; 
                {   //read initial thickness
                    std::vector<V_ind> tmp;
                    H_->saveData(&Hp, f, tmp);
                }
                double scl = Jinv * (*Hp);
                for (auto& v: s) v *= scl;
                T_[f] = s;
            }
            if (!strait3.m_tag_name.empty())
                m_mtm.readTagDataFromMesh<F_ind>(m, "f:tensor:T", strait3.m_tag_name);
        }
        strait = m_out_traits[SAVE_CURVATURE];
        if (!strait.m_tag_name.empty()){
            if (bf_id < 0){
                std::cout << "Warning: Curvature storing supported only for shell approximation and ANALYTICAL potential\n";
            } else {
                auto k_ = m.add_property_map<F_ind, std::array<double, 3>>("f:tensor:curvature").first;
                obj.m_forces[bf_id].target<BendingForce>()->add_watch(BendingForce::compute_current_curvature_expr);
                for (auto f: obj.m_mesh.faces()) {
                    auto dat = obj.m_forces[bf_id].target<BendingForce>()->eval_watches(obj, f)[0];
                    k_[f] = std::array<double, 3>{dat(0, 0), dat(1, 0), dat(2, 0)};
                }
                m_mtm.readTagDataFromMesh<F_ind>(m, "f:tensor:curvature", strait.m_tag_name);
            }
        }
        
        //save results in files
        if (save_flag & SAVE_BIN_VTK){
            std::string res_file = (path(save_dir)/save_prefix).string() + "_bin.vtk"; 
            VTK_File().write(res_file, m_mtm, true);
            std::cout << "Save results in \"" << res_file << "\"\n";
        } 
        if (save_flag & SAVE_TXT_VTK){
            std::string res_file = (path(save_dir)/save_prefix).string() + "_txt.vtk";
            VTK_File().write(res_file, m_mtm, false);
            std::cout << "Save results in \"" << res_file << "\"\n";
        } 
        if (save_flag & SAVE_BIN_TMH){
            std::string res_file = (path(save_dir)/save_prefix).string() + "_bin.tmh";
            TMH_File().write_bin(res_file, m_mtm);
            std::cout << "Save results in \"" << res_file << "\"\n";
        }
        if (save_flag & SAVE_TXT_TMH){
            std::string res_file = (path(save_dir)/save_prefix).string() + "_txt.tmh";
            TMH_File().write_txt(res_file, m_mtm);  
            std::cout << "Save results in \"" << res_file << "\"\n"; 
        }    
    }

    if (m_view && m_perform_solve){
        view_window->join();
        delete view_window;
    }

    return status;
}