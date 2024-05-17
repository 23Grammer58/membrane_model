//
// Created by Liogky Alexey on 01.08.2022.
//

#ifndef AVSYM_VANLOONBENCH_H
#define AVSYM_VANLOONBENCH_H

#include "../BenchCommon.h"

#include "AVSim/AorticSimulator/Messurer.h"
#include "AVSim/Core/IO/Obj3dVTKSaver.h"
#include "AVSim/Core/NSWorldWrapper.h"
#include "AVSim/Solvers/NonLinearSolverKinsol.h"
#include "AVSim/Core/NonLinearSolverCustom.h"
#include <Eigen/Dense>
#include "AVSim/Core/MeshGen/Helper.h"
#include "AVSim/Core/MeshGen/StretchFigure.h"
#include "AVSim/LeafTemplates/TemplatesCollection.h"

static bool _choose_edge_default(Object3D& leaf, E_ind e, V_ind* v, V_ind vop){
    const int FIXED = 2, FREE = 1|4|8;
    return (leaf.m_boundary[v[0]] & FIXED) && (leaf.m_boundary[v[1]] & FIXED) && !(leaf.m_boundary[vop] & FREE);
}
struct ComputeInitMesh{
    double R = 12;
    double phi = 2*M_PI/3 - 0.4, Hc = 1.5, Hbc = 1, Omega = -40.0/180.0 * M_PI;//-30.0; +8.0;
    double mesh_h = 0.45; //0.3;
    double P = 90.0_mmHg / 1e3, mu = 900, Ht = 0.4; //Ht = 0.8
    double c_delta = 0.05, err = 1e-8, abs_err = 1e-4, diverge_lim = 1e10;
    int maxits = 1e2, pfreq = 5e2;
    bool use_clamped = true;
    bool use_true_bc = false;
    bool use_bending = true;
    string case_name = "1";
    string leaf_name = "ModelLeaf";
    string dir = "../../../result/VanLoonBench/temp1/";
    double bpln_Pc = 2, bpln_Htc = 0.1, bpln_bilc = -1.3;
    double fpln_Pc = 3, fpln_Htc = 0.01, fpln_posc = /*0.35**/0.2;
    double spln_Pc = 4, spln_Htc = 0.01, spln_posc = 0 /*0.03*/, spln_phi = M_PI/3;

    bool evalute_quality_of_res = true;
    bool with_view = true;
    int argc = 0; 
    char** argv = nullptr;
    
    ComputeInitMesh& set_args(int _argc, char** _argv){ argc = _argc, argv = _argv; return *this; }
};

bool ObSaveTag(const Object3D& obj, string filename, string tag);
bool ObReadTag(Object3D& obj, string filename, string tag);
StretchFigure prepareInitialLeaf(double R, double phi, double Hc, double Hb, double Omega);
Object3D computeInitialLeaf(double R, double phi, double Hc, double Hb, double Omega, double mesh_h);
std::map<E_ind, BendingForce::BoundaryData::Entry>
    computeBndDataOnCilAligned(Object3D& leaf, 
        std::function<bool(Object3D&, E_ind, V_ind*/*[2]*/, V_ind)> choose_edge = _choose_edge_default, 
        double angle_align = 0);
Mesh::Property_map<V_ind, Vector> addBdataTag(Object3D& leaf, const map<E_ind, BendingForce::BoundaryData::Entry>& bdata);
int computeInitialMesh(ComputeInitMesh prm);

// TODO: need to fix
// int save_ModelLeaf_to_mat(string postfix = "_4_3");
// int SaveDisplaceMap();
// int ModelLeafBench1(int argc, char* argv[]);
// int ModelLeafBench(int argc, char* argv[]);
// int VanLoonBench(int argc, char* argv[]);



#endif //AVSYM_VANLOONBENCH_H