#include "VirtualSuturer.h"
#include <boost/program_options.hpp>

struct MeshD{
    double mesh_size; //mm
    double D = 0, L_f = 0, H = 0, H_f = 0; //mm
    ///The error code shows which variable needs to be changed to get rid of the error
    enum ErrCode{
        OK = 0,
        ERR_H_f = 1,
        ERR_D = 2,
        ERR_L_f = 3,
        ERR_H = 4,
        ERR_UNRECOVERED = 5
    };

    OzakiTempl getOzaki() const;
    std::pair<bool, double> check_Hf_condition() const;
    std::pair<bool, double> check_D_condition() const;
    ///return err code and approximate change to be applyeid to recover problem to unacceptable variable
    std::pair<ErrCode, double> check_state() const;
    double get_sutured_boundary_length() const;
    double get_free_boundary_length() const;
    friend std::ostream& operator<<(std::ostream& out, const MeshD& t);
};
std::ostream& operator<<(std::ostream& out, const MeshD& t);
struct PhysProperties{
    double Ht;
    double E;
    double P;
};
std::ostream& operator<<(std::ostream& out, const PhysProperties& t);

struct LeafMessures{
    std::array<Messurer::CoaptationDistrib, 2> half_distrib;
    Messurer::CoaptationDistrib union_distrib;
    double central_coaptation = -1;
    double length_coaptation = -1;
    double billowing = -1;
    Messurer::ColArea coaptation_area;
    bool form_closed_valve = false; 
    double closure_degree = -1;
    double effective_height = -1;
    double free_edge_init_length = -1, free_edge_res_length = -1;
    void print(std::ostream& out = std::cout);
};

struct TestVirtualSutureParams{
    std::array<MeshD, 3> mhd;
    PhysProperties props;
    int nleafs = 3;
    std::array<DReal, 3> phi_relations;
    bool with_gui_window = false;
    std::string root_directory = "";
    std::string result_directory = "";
    double R = 12;
    double dist_suture = 4;
    double Hc = 1.5, Hb = 12;
    std::array<int, 3> bnd_types{1, 1, 1};
    bool save_steps = false;
    bool with_bezier_contact = false;
    bool with_aorta = false;
    bool with_bending_forces = false;
    bool with_clamped_bc = false;
    bool use_quasi_optimal_leaf_width = false;
    std::vector<int> ileafs = {0};
    bool save_col_shapes = false;
    int verbose_level = 2;
    bool regenerate_expressions = false;
    int argc = 0;
    char** argv = nullptr;
    double collision_Ht_coef = 1.0, collision_Ht_shift = 0;

    TestVirtualSutureParams(); 
    std::array<double, 6> getLeafCommisureAngles() const;
    std::array<Point, 3> getCommissurePoints() const;
    void setQuasiOptWidth();
    void parseInputArgs(int argc, char* argv[]);
};

std::map<int, LeafMessures> run_suture_simulation(TestVirtualSutureParams p);