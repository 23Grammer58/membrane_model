//
// Created by alex on 16.10.2023.
//

#ifndef GENERAL_MODEL_EXAMPLES_H
#define GENERAL_MODEL_EXAMPLES_H

#include <iostream>
#include "AVSim/Core/TriangularMeshHelpers.h"
#include "AVSim/DomainCollection/AniMeshCollection.h"
#include "AVSim/Core/IO/VTKFormat.h"
#include "AVSim/Core/IO/MetaTriFormat.h"

struct stretch_rectangle{
    double  a = 10, b = 10, ht = 0.4, mesh_h = 0.5, dx = 1, dy = 1, theta_f = 0, theta_s = M_PI_2; 
    bool is_2d = true, is_fixed = false;
    std::ostream& print_param_list(std::ostream& out = std::cout) const;
    std::ostream& print_state(std::ostream& out = std::cout) const;
    stretch_rectangle& init(int argc, char* argv[]);
    int generate() const;
};

struct stretch_cross_piece{
    double a = 20, ac = 5, r = 7.5, ht = 0.4, mesh_h = 0.5;
    double dx0 = 0.75, dy0 = 0.75, dx1 = 0, dy1 = 0;
    double theta_f = 0, theta_s = M_PI_2;
    int is_fixed = 1;
    double r1 = -1, r2 = -1, r3 = -1;

    std::ostream& print_param_list(std::ostream& out = std::cout) const;
    std::ostream& print_state(std::ostream& out = std::cout) const;
    stretch_cross_piece& init(int argc, char* argv[]);
    int generate() const;
};

struct inflation_circle{
    double R = 10, ht = 0.4, mesh_h = 0.5, theta_f = 0, theta_s = M_PI_2, R_c = -1, K = 0;
    std::ostream& print_param_list(std::ostream& out = std::cout) const;
    std::ostream& print_state(std::ostream& out = std::cout) const;
    inflation_circle& init(int argc, char* argv[]);
    int generate() const;
};


#endif //GENERAL_MODEL_EXAMPLES_H