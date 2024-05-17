//
// Created by Liogky Alexey on 01.08.2022.
//
#include "VanLoonBench.h"

int main(int argc, char* argv[]){
    auto p = ComputeInitMesh().set_args(argc, argv);
    p.Omega = -50.0/180.0 * M_PI;
    p.fpln_posc = 0.00;
    p.case_name = "2";

    return computeInitialMesh(p);
}