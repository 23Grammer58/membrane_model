//
// Created by alex on 30.06.2020.
//

#include <stdexcept>
#include "Forces.h"

#ifdef BEND_WATCH
World3d::Mesh::Property_map<World3d::F_ind, std::array<double, 3*5>> K_comp;
int glob_k_it = 0;
KProxy g_Kproxy;
#endif




