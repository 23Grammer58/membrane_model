//
// Created by alex on 01.10.2023.
//

#ifndef GENERAL_MODEL_H
#define GENERAL_MODEL_H

#include <iostream>
#include <thread>
#include "AVSim/Core/Renderers/GuiApplication.h"
#include "AVSim/Core/Object3D.h"
#include "AVSim/Core/TriangularMeshHelpers.h"
#include "AVSim/DomainCollection/AniMeshCollection.h"
#include "AVSim/Core/Renderers/GuiBackInterconnection.h"
#include "AVSim/Core/World.h"
#include "AVSim/Core/ForceAppliers/ForceAppliers.h"
#include "AVSim/Core/Forces/Forces.h"
#include "AVSim/Core/Renderers/Renderer.h"
#include "AVSim/Core/IO/VTKFormat.h"
#include "AVSim/Core/IO/MetaTriFormat.h"

#include "AVSim/Core/IO/HyperElasticFormat.h"

struct Model{
    static constexpr unsigned int SAVE_BIN_VTK = 0x1;
    static constexpr unsigned int SAVE_TXT_VTK = 0x2;
    static constexpr unsigned int SAVE_BIN_TMH = 0x4;
    static constexpr unsigned int SAVE_TXT_TMH = 0x8;

    static constexpr unsigned int SAVE_FINAL_STATE = 0x10;
    static constexpr unsigned int SAVE_RESIDUAL = 0x20;
    static constexpr unsigned int SAVE_LOCAL_BASIS = 0x40;
    static constexpr unsigned int SAVE_LOCAL_BASIS_EXT = 0x400;
    static constexpr unsigned int SAVE_PK2_STRESS = 0x80;
    static constexpr unsigned int SAVE_DEF_GRAD = 0x100;
    static constexpr unsigned int SAVE_CAUCHY_STRAIN = 0x200;
    static constexpr unsigned int SAVE_CURVATURE = 0x800;
    static constexpr unsigned int SAVE_CAUCHY_STRESS3D = 0x1000;
    static constexpr unsigned int SAVE_CAUCHY_TENSION3D = 0x2000;
    static constexpr unsigned int MAX_SAVE_FLAG = SAVE_CAUCHY_TENSION3D;

    enum ObjectTraitType{
        INIT_X_VECTOR_DTAG = 0,
        START_X_VECTOR_DTAG = 1,
        BOUNDARY_ITAG = 2,
        PRESSURE_DTAG = 3,
        LAGRANGE_PRESSURE_DTAG = 4,
        BODY_FORCE_VECTOR_DTAG = 5,
        EDGE_FORCE_VECTOR_DTAG = 6,
        CLAMPED_EDGE_VECTOR_DTAG = 7,
        THICKNESS_DTAG = 8,
        // END_X_VECTOR_DTAG,
    };
    static std::string getObjectTraitTypeName(ObjectTraitType t);
    static const std::map<std::string, ObjectTraitType>& getObjectTraitTypeToNameMap();
    static ObjectTraitType getObjectTraitTypeByName(const std::string& s);

    struct TraitValue{
        MetaTriMesh::ElementSparsity m_etype = MetaTriMesh::ElementSparsity::NODE;
        std::string m_tag_name = "";
        std::array<double, 3> m_def_value = {NAN, NAN, NAN};
        std::size_t m_dim = 0;
        bool m_real = true;
        TraitValue() = default;
        TraitValue(MetaTriMesh::ElementSparsity etype, std::size_t dim, std::string tag_name = "", std::array<double, 3> def_value = {NAN, NAN, NAN}, bool is_real = true): 
            m_etype{etype}, m_dim{dim}, m_tag_name{tag_name}, m_def_value{def_value}, m_real{is_real} {}
        
        std::ostream& print(std::ostream& out = std::cout) const;   
    };
    static std::map<ObjectTraitType, TraitValue> getDefaultTraits();

    struct SaveTrait{
        std::string m_tag_name;
        MetaTriMesh::ElementSparsity m_etype = MetaTriMesh::ElementSparsity::NODE;
        std::size_t m_dim = 0;
        SaveTrait() = default;
        SaveTrait(MetaTriMesh::ElementSparsity etype, std::size_t dim, std::string tag_name = ""): m_etype{etype}, m_dim{dim}, m_tag_name{tag_name} {} 

        std::ostream& print(std::ostream& out = std::cout) const;
    };
    static std::map<unsigned int, SaveTrait> getDefaultSaveTraits();
    static const std::map<std::string, unsigned int>& getSaveTraitToNameMap();
    static std::vector<std::pair<unsigned int, SaveTrait*>> getSaveTraitByName(std::map<unsigned int, SaveTrait>& save_traits, const std::string& name);
    static std::string getSaveTraitName(unsigned int trait);

    struct SaveInTag{
        MetaTriMesh::ElementSparsity m_etype = MetaTriMesh::ElementSparsity::NODE;
        std::string m_tag_name;
        std::string m_new_tag_name;
        SaveInTag() = default;
        SaveInTag(MetaTriMesh::ElementSparsity etype, const std::string& tag_name, const std::string& new_tag_name = ""): m_etype{etype}, m_tag_name{tag_name}, m_new_tag_name{new_tag_name} {}
    };

    struct CustomSolverIteration{
        int relax_its = 1000, relax_pfreq = 50;
        int ns_type = 0;
        double delta_full_step = 0.1;
        int max_ns_its = 100;
        
        static std::string solverName(int sol_type); 
        static int solverTypeByName(const std::string& val);
        std::ostream& print(std::ostream& out  = std::cout) const;
    };

    Model& parse(int argc, const char* argv[], int ignore_argv_cnt = 1);
    Model& parse(std::vector<std::string> args);
    static void parse_config_file(const std::string& config_file, std::vector<std::string>& args, std::size_t pos, std::size_t delete_next);
    static void printArgsHelpMessage(const std::string& prefix = "", std::ostream& out = std::cout);

    int execute(int argc = 0, char** argv = nullptr);

    std::string fmesh = "mesh.vtk", 
                fpotential = "potental.energy", 
                save_dir = "", 
                save_prefix = "res",
                to_gen = "generated";
    bool regen_potential = true; 
    bool is_membrane_approx = true; 

    bool m_perform_solve = true; 
    std::vector<CustomSolverIteration> m_solverItData;         
    std::string m_lin_sol_name = "inner_mptiluc"; 
    std::vector<std::pair<std::string, std::string>> m_lin_sol_prms;            
    double m_abs_tol = 1e-7, m_rel_tol = 1e-7;
    
    bool m_view = false;
    
    unsigned int save_flag = 0x0;
    bool m_perform_postcomps = true; 
    
    MetaTriMesh m_mtm;
    Energy_Parser m_energy_expr;
    std::map<ObjectTraitType, TraitValue> m_input_traits;
    std::map<unsigned int, SaveTrait> m_out_traits;
    std::map<MetaTriMesh::ElementSparsity, std::map<std::string, std::string>> m_save_intag;
};


#endif //GENERAL_MODEL_H

