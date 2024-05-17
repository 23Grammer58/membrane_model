//
// Created by alex on 17.11.2020.
//

#ifndef AORTIC_VALVE_CONFIGURATION_H
#define AORTIC_VALVE_CONFIGURATION_H

#include <string>
#include <array>
#include <vector>
#include <map>
#include <nlohmann/json.hpp>
#include "AVSim/Core/Forces/Forces.h"

class Configuration {
    using string = std::string;
    template<typename T>
    using vector = std::vector<T>;
    template<typename T, int n>
    using array = std::array<T, n>;
    template<typename T1, typename T2>
    using map = std::map<T1, T2>;
public:
    struct LeafletInput{
        bool generative = true;
        //if true then leaflet will be created according to model
        //else leaflet will be loaded from leaf_file
        string model = "ozaki";
        nlohmann::json model_params = {
                {"mesh_h", 1.2},
                {"template", 17}
        };

        string leaflet_file = "leaflet.txt";
        array<double, 3> main_axis_direction = {0, 1, 0};
        //this information is necessary for sewing with right orientation
        array<double, 3> second_line_symmetry= {0, 0, 1};
        /*
         * The plane formed by the main_axis_direction and second_line_symmetry
         * should divide the pattern so that the extreme points of the
         * pattern are extreme relative to this plane
         */
    };

    struct SewEnergyParams{
        /*The weights with which these restrictions are taken in the energy functional*/
        double plane_w = 1.0;
        double sqr_length_w = 1.0;
        double digedral_angle_w = 0.7;
        //Additional factor increasing energy from concavities
        //The more factor, the less concavities
        double convexity_w = 1;
        double force_w = 0.0;
        double shift = 0.5;
        struct EnergyMinimizerParams {
            //solver params - sp
            int freq = 50;
            double step_sz = 1.e-4;
            double tol = 1.e-4;
            double epsabs = 1.e-2;
            int maxits = 301;
            double time/*secs*/ = 3;
        };
        EnergyMinimizerParams sp;
    };

    struct SolverParams{
        string method = "simple_relaxation";
        nlohmann::json method_vals = {
                {"delta", 5.5e-7},  //coefficient converting force into shift
                {"eps", 0.001},     //stop condition
                {"max_possible_shift_scale", 200},
                {"max_possible_recommend_shift_scale", 200},
                {"maxits", 72000}
        };
    };

    struct ObjectTraits{
        string elastic_model_type = "NeoGookModel";
        nlohmann::json elastic_model_params = { //coefficients of model
                {"H", 0.5},
                {"mu", 1.0e6/3.0}
        };
        bool use_bending = false;
        double pressure = World3d::operator""_mmHg(80.0);                //pressure acted on object
        double collision_margin = 0.05;                      //The distance at which collision detection is turned on
    };
    struct SaveResults{
        string result_dir = "../result/";
        string postfix = "";
        int its_save_freq = 5000;
        bool save_aorta = true;
        bool save_cutoff_aorta = true;
        bool save_leafs = true;
        bool save_sewed_leafs = true;
        bool try_reload_last_state = true;
    };

    string bs_dir = "../data/";
    vector<LeafletInput> l_ins;
    string leafletInputFile = "";

    SewEnergyParams sep;
    string sewEnergyParamsInputFile = "";

    std::vector<ObjectTraits> ots;
    string objectTraitsInputFile = "";

    SolverParams sp;
    string solverParamsInputFile = "";

    int scenario = 0;
    SaveResults svr;
    Configuration() {}
    Configuration(int argc, const char* argv[]) { ConfigurationInit(argc, argv); }
    void ConfigurationInit(int argc, const char* argv[]);

    void processConfigFile(const char* filename);
    void processConfigFile_leaflets(string filename);
    void processConfigFile_sew_energy(string filename);
    void processConfigFile_solver(string filename);
};


#endif //AORTIC_VALVE_CONFIGURATION_H
