#include <array>
#include <limits>
#include "TestVirtualSuturer.h"

struct SimplifiedTemplateGeometry{
    double H = NAN;       ///< value in [5, 25], start H = 11 
    double H_f = NAN;     ///< value in [-H/2, 6], start H_f = 2
    double D = NAN;       ///< value in [9, 34], start D = 19
    
    double mesh_step = 0.7;///< mesh resolution
};
struct TestConstantParameters{
    double R = NAN;              ///< radius of aortic root
    double Hc = NAN, Hb = NAN;   ///< geometric heights of aortic suturing lines
    std::array<double, 3> phi_relations = {1, 1, 1}; ///<The ratio between the angular sectors of the aortic ring falling on the corresponding cusp pattern
};
struct EvaluatedValues{
    using _nl = std::numeric_limits<double>;
    bool form_closed_valve = false;    ///< true if valve closed well
    double closure_degree = _nl::max();///< if the valve is not closed, then the lower this positive parameter, the closer the valve is to the closed state
    double central_coaptation = -1;    ///< central coaptation, optimal >1.5
    double effective_height = -1;      ///< effective height, optimal 9-11
    double length_coaptation = -1;     ///< length of coaptation, optimal > 6   
    double billowing = -1;             ///< optimal <1.4
    double coaptation_area = -1;       ///< if the other parameters are in a good area, then the smaller the better
};
static std::ostream& operator<<(std::ostream& out, const EvaluatedValues& p){
    return out 
        << "form_closed_valve  = " << (p.form_closed_valve ? "true" : "false") << " ( closure_degree = " << p.closure_degree << " )\n"
        << "central_coaptation = " << p.central_coaptation << "\n"
        << "effective_height   = " << p.effective_height << "\n"
        << "length_coaptation  = " << p.length_coaptation << "\n"
        << "billowing          = " << p.billowing << "\n"
        << "coaptation_area    = " << p.coaptation_area << "\n";
}
struct SimResult: public EvaluatedValues{
    MeshD::ErrCode merr = MeshD::OK; ///< geometric error
    double add_value = 0;            ///< if merr is equal to ERR_H_f or ERR_D or ERR_H then add_value is expected change to be applied to the corresponding parameter to overcome error
};

static SimResult performOneLeafSimulation(SimplifiedTemplateGeometry g, TestConstantParameters q = TestConstantParameters(), int argc = 0, char* argv[] = nullptr){
    SimResult res;
    TestVirtualSutureParams p;
    p.parseInputArgs(argc, argv);
    if (!std::isnan(q.R) && !std::isnan(q.Hc) && !std::isnan(q.Hb)){
        p.R = q.R; p.Hc = q.Hc, p.Hb = q.Hb, p.phi_relations = q.phi_relations;
    }
    p.ileafs.resize(1); p.ileafs[0] = 0;
    if (!std::isnan(g.D) && !std::isnan(g.H) && !std::isnan(g.H_f)){
        p.mhd[0].D = g.D, p.mhd[0].H = g.H, p.mhd[0].H_f = g.H_f;
        p.mhd[0].mesh_size = g.mesh_step;
    }
    for (int i = 1; i < 3; ++i) p.mhd[i].D = p.mhd[i].H = p.mhd[i].H_f = -1;
    p.setQuasiOptWidth();
    auto state = p.mhd[0].check_state();
    if (state.first != MeshD::OK) {
        res.merr = state.first;
        res.add_value = state.second;
        return res;
    }
    auto lmm = run_suture_simulation(p);
    auto lm = lmm[0];
    res.form_closed_valve = lm.form_closed_valve;
    res.closure_degree = lm.closure_degree;
    res.central_coaptation = lm.central_coaptation;
    res.effective_height = lm.effective_height;
    res.length_coaptation = lm.length_coaptation;
    res.billowing = lm.billowing;
    res.coaptation_area = lm.coaptation_area.colArea;
    return res;
}