//
// Created by alex on 17.06.2021.
//

#ifndef AORTIC_VALVE_NSWORLDWRAPPER_H
#define AORTIC_VALVE_NSWORLDWRAPPER_H

#include "World.h"
#include "AVSim/Solvers/NonLinearSolverInterface.h"

namespace World3d {
    struct NSWorldWrapper {
        Timer it_time;
        double render_freq = 1.0 / 10;
        World *pw;

        NSWorldWrapper() = default;
        NSWorldWrapper(World &w) : pw{&w} { it_time.reset(); }
        
        void resetWorld(World &w){ pw = &w; it_time.reset(); }
        int getNodfs();
        int setX(const double *X);
        int getCurrentX(double *X);
        std::vector<double> getCurrentX();
        //R resized
        int Residual(const double *X, double *R);
        //sm cleared and resized! (add jacobian to sm)
        int Jacobian(const double *X, SparseMatrix *sm);
        int System(const double *X, SparseMatrix *sm, double *R);
        int RenderIteration();
        int RenderVector(const double *X);
        int compResidual(double *R);
        int compJacobian(SparseMatrix *sm);
        double compCorrentForceNorm(unsigned int Nobj = 2, unsigned int Nmesh = 2, unsigned int Nnode = 2);
        static double VectorNorm(double* residual, long size, unsigned int N1, unsigned int N2);

        static NLProblem makeProblem(NSWorldWrapper& nsww);
    };
}


#endif //AORTIC_VALVE_NSWORLDWRAPPER_H
