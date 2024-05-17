//
// Created by alex on 17.11.2020.
//

#include "Configuration.h"

#include <iostream>
#include <cstring>
#include <fstream>
#include <string>

using namespace std;

void Configuration::ConfigurationInit(int argc, const char **argv) {
    int i = 0;
    if (argc == 1) goto helpMessage;
    for (i = 1; i < argc; i++) {
        //Print help message and exit
        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            helpMessage:
            std::cout << "Help message: " << std::endl;
            std::cout << "Command line options: " << std::endl;
            std::cout << "-c, --config <Main configuration parameters file name>" << std::endl;
            std::cout << "-lc, --leafscnfg <Leaflets configuration parameters file name>" << std::endl;
            std::cout << "-sc, --solvercnfg <Solver configuration parameters file name>" << std::endl;
            exit(0);
        }
        //Parameters file name found with -d or --database options
        if (strcmp(argv[i], "-c") == 0 || strcmp(argv[i], "--config") == 0) {
            std::cout << "Main configuration file found: " << argv[i + 1] << std::endl;
            processConfigFile(argv[i + 1]);
            i++;
            continue;
        }
        if (strcmp(argv[i], "-lc") == 0 || strcmp(argv[i], "--leafscnfg") == 0) {
            std::cout << "Leaflets configuration file externally set to : " << argv[i + 1] << std::endl;
            processConfigFile_leaflets(argv[i + 1]);
            i++;
            continue;
        }
        if (strcmp(argv[i], "-sc") == 0 || strcmp(argv[i], "--solvercnfg") == 0) {
            std::cout << "Solver configuration file externally set to : " << argv[i + 1] << std::endl;
            processConfigFile_solver(argv[i + 1]);
            i++;
            continue;
        }
    }
}

using json = nlohmann::json;

json get_json(string filename){
    std::ifstream js;
    js.open(filename);
    if (! js.good()) {
        throw std::runtime_error("Unsuccessful attempt to open file \"" + filename + "\"");
    }
    json conf;
    js >> conf;
    js.close();

    return std::move(conf);
}

void no_field(string field){
    std::cout << "Warning: Can't find the field \"" + field + "\". For this will be used defaults.\n";
}

array<double, 3> from_vec(vector<double> v){
    array<double, 3> res;
    for (int i = 0; i < 3 && i < v.size(); ++i)
        res[i] = v[i];
    return std::move(res);
}

#define SET(X, Y, Z) if (Y.count(Z)) X = Y[Z]

#define SETW(X, Y, Z, W) if (Y.count(Z)) X = Y[Z]; else no_field(W)

string get_directory(string filepath){
    return filepath.substr(0, filepath.find_last_of("/")) + "/";
}

void Configuration::processConfigFile_leaflets(string filename) {
    if (leafletInputFile.size() && filename != leafletInputFile)
        cout << "Warning: The field \"leaflets\" in file " + leafletInputFile + " will be ignored.\n";
    leafletInputFile = filename;
    objectTraitsInputFile = filename;
    json conf = get_json(filename);

    if (!conf.count("leaflets")){
        no_field("leaflets");
        return;
    }
    conf.size();
    int n = conf["leaflets"].size();
    l_ins.resize(n);
    ots.resize(n);
    string direct = get_directory(filename);
    for (int i = 0; i < n; ++i){
        json leaf = conf["leaflets"][i];
        if (!leaf.count("mesh")){
            no_field("leaflets[" + to_string(i) + "]:mesh");
        } else {
            json mesh = leaf["mesh"];
            if (!mesh.count("generative")) {
                no_field("leaflets[" + to_string(i) + "]:mesh:generative");
            }
            else
                l_ins[i].generative = mesh["generative"];
            SET(l_ins[i].model, mesh, "model");
            if (mesh.count("model_params")) {
                l_ins[i].model_params = mesh["model_params"];
            }
            if (mesh.count("mesh_file"))
                l_ins[i].leaflet_file = direct + mesh["mesh_file"].get<string>();
            else
            if (!l_ins[i].generative)
                no_field("leaflets[" + to_string(i) + "]:mesh:mesh_file");
            if (mesh.count("main_axis_direction"))
                l_ins[i].main_axis_direction = from_vec(mesh["main_axis_direction"]);
            if (mesh.count("second_line_symmetry"))
                l_ins[i].second_line_symmetry = from_vec(mesh["second_line_symmetry"]);
        }
        if (!leaf.count("traits")){
            no_field("leaflets[" + to_string(i) + "]:traits");
        } else {
            json traits = leaf["traits"];
            bool has_model = true;
            if (!traits.count("elastic_model")) {
                no_field("leaflets[" + to_string(i) + "]:traits:elastic_model");
                has_model = false;
            } else ots[i].elastic_model_type = traits["elastic_model"];
            if (traits.count("elastic_model_params")){
                if (!has_model)
                    cout << "Warning: in leaflets[\" + to_string(i) + \"]:traits: \"elastic_model\" not setted but \"elastic_params\" setted.\n";
                ots[i].elastic_model_params = traits["elastic_model_params"];
            }
            SET(ots[i].pressure, traits, "pressure");
            SET(ots[i].use_bending, traits, "use_bending");
            SET(ots[i].collision_margin, traits, "collision_margin");
        }
    }
}

void Configuration::processConfigFile_sew_energy(string filename) {
    if (sewEnergyParamsInputFile.size() && filename != sewEnergyParamsInputFile)
        cout << "Warning: The field \"sew_energy_params\" in file " + sewEnergyParamsInputFile + " will be ignored.\n";
    sewEnergyParamsInputFile = filename;
    json conf = get_json(filename);

    if (!conf.count("sew_energy_params")){
        no_field("sew_energy_params");
        return;
    }
    conf = conf["sew_energy_params"];
    SET(sep.plane_w, conf, "plane_constraint_weight");
    SET(sep.sqr_length_w, conf, "mesh_edge_constraint_weight");
    SET(sep.digedral_angle_w, conf, "mesh_digedral_element_weight");
    SET(sep.convexity_w, conf, "mesh_convexity_factor");
    SET(sep.force_w, conf, "external_force_weight");
    SET(sep.shift, conf, "axis_shift_factor");
    if (conf.count("energy_minimizer_params")){
        json emp = conf["energy_minimizer_params"];
        SET(sep.sp.freq, emp, "freq");
        SET(sep.sp.step_sz, emp, "step_size");
        SET(sep.sp.tol, emp, "tol");
        SET(sep.sp.epsabs, emp, "epsabs");
        SET(sep.sp.maxits, emp, "maxits");
        SET(sep.sp.time, emp, "max_time_per_leaflet");
    }
}

void Configuration::processConfigFile_solver(string filename) {
    if (solverParamsInputFile.size() && filename != solverParamsInputFile)
        cout << "Warning: The field \"solver_params\" in file " + solverParamsInputFile + " will be ignored.\n";
    solverParamsInputFile = filename;
    json conf = get_json(filename);

    if (!conf.count("solver_params")){
        no_field("solver_params");
        return;
    }
    conf = conf["solver_params"];
    SET(sp.method, conf, "method");
    SET(sp.method_vals, conf, "method_params");
}

void Configuration::processConfigFile(const char* filename){
    json conf = get_json(filename);

    std::vector<std::pair<string, string>> conf_types = {
            {"aorta_file", "aorta"},
            {"leaflets_file", "leaflets"},
            {"sew_energy_file", "sew_energy_params"},
            {"solver_file", "solver_params"}
    };
    string thisfile = filename;
    int __pos = thisfile.find_last_of('/');
    __pos = (__pos == thisfile.size() - 1) ? 0 : __pos;
    string dir = thisfile.substr(0, __pos + ((thisfile.size() > 1) ? 1 : 0));

    bs_dir = dir;
    if (conf.count("leaflets_file")) {
        processConfigFile_leaflets(dir + string(conf["leaflets_file"]));
        if (conf.count("leaflets"))
            cout << "Warning: The field \"leaflets\" in file " + thisfile + " will be ignored.\n";
    }
    else processConfigFile_leaflets(thisfile);
    if (conf.count("sew_energy_file")) {
        processConfigFile_sew_energy(dir + string(conf["sew_energy_file"]));
        if (conf.count("sew_energy_params"))
            cout << "Warning: The field \"sew_energy_params\" in file " + thisfile + " will be ignored.\n";
    }
    else processConfigFile_sew_energy(thisfile);
    if (conf.count("solver_file")) {
        processConfigFile_solver(dir + string(conf["solver_file"]));
        if (conf.count("solver_params"))
            cout << "Warning: The field \"solver_params\" in file " + thisfile + " will be ignored.\n";
    }
    else processConfigFile_solver(thisfile);
    SET(scenario, conf, "scenario");
}