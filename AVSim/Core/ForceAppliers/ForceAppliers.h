//
// Created by alex on 30.06.2020.
//

#ifndef AORTIC_VALVE_FORCEAPPLIERS_H
#define AORTIC_VALVE_FORCEAPPLIERS_H

#include "../Object3D.h"

namespace World3d {
    class StaticForceApplier{
    public:
        double _delta;
        StaticForceApplier(double delta): _delta{delta} {}
        int operator() (Object3D& obj);
    };

    class StaticForceApplierP{
    public:
        using Damper = std::function<double(Object3D& obj, double cur_delta)>;
        double& _delta;
        std::vector<Damper> _dampers;

        void addDamper(Damper dmp){ _dampers.push_back(std::move(dmp)); }
        StaticForceApplierP(double& delta): _delta{delta} {}
        int operator() (Object3D& obj);
    };

    struct Damper{
    private:
        double factor = 1, init_dx = -1;
        int ok_its = 0, bad_its = 0, increase_its = 0, m_it = 0;
        double reinit_val = 0;
    public:
        double min_factor = 0.01;   //set minimal availiable factor, after achieving that factor will not decrease
        double init_incr = 3;       //set bound |dx| < init_incr * init_dx
        //if during re_init_nit iterations reinit_val = max(|dx|) < reinit_min_factor * |allow_dx| than init_dx = 1.5 * reinit_val
        double reinit_min_factor = 0.2; int re_init_nit = 5000;
        //if there are damp_it optional consecutive exceeding the permissible bounds (|dx| > |allow_dx|)
        //  than delta decrease in scl times (factor multiplies by scl < 1)
        //if there are no_damp_it consecutive iterations without exceeding the permissible bounds
        //  than delta increase in scl times (factor divides by scl < 1)
        int damp_it = 100, no_damp_it = 200000;
        double scl = 0.9;
        //manually setted bound: if initial delta was chosen lucky, than this value may be ignored
        double max_dx;

        Damper(double max_dx, double scl = 0.9, int damp_it = 100, int no_damp_it = 200000):
                max_dx{max_dx}, scl{scl}, damp_it{damp_it}, no_damp_it{no_damp_it} {}
        Damper& setMaxDx(double maxdx){ max_dx = maxdx; return *this; }
        Damper& setScaleFactor(double scale) { scl = scale; return *this; }
        Damper& setInitDxScale(double init_scl) { init_incr = init_scl; return *this; }
        Damper& setMinDeltaFactor(double min_fact) { min_factor = min_fact; return *this; }
        Damper& setReInitDxParams(int nit, double reinit_fact_bnd) { re_init_nit = nit, reinit_min_factor = reinit_fact_bnd; return *this; }
        Damper& setAvailBndErrNumIts(int errNumIts){ damp_it = errNumIts; return *this; }
        Damper& setIncrDeltaNumIt(int incrNumIt) { no_damp_it = incrNumIt; return *this; }
        double getUsedDeltaFactor() { return factor; }
        double operator()(Object3D& obj, double cur_delta);
    };

    class DynamicForceApplier{
    public:
        typedef Mesh::Property_map<V_ind, Vector> Velocity;
        typedef Mesh::Property_map<V_ind, double> Mass;
        typedef std::tuple<Velocity, Mass> MeshData;

        double _dt;
        std::map<Object3D*, MeshData> data;
        DynamicForceApplier(double dt): _dt{dt} {}
        auto registerObj(Object3D* obj);
        int operator() (Object3D& obj);
    };
}


#endif //AORTIC_VALVE_FORCEAPPLIERS_H
