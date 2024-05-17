//
// Created by Liogky Alexey on 19.05.2022.
//

#ifndef AVSYM_TEMPLATESCOLLECTION_H
#define AVSYM_TEMPLATESCOLLECTION_H

#include "AVSim/Core/TriangularMeshHelpers.h"
#include "AVSim/Core/Object3D.h"

struct ThubricarValvePrms{
    double R1;      //radius of top base
    double R2;      //radius of bottom base
    double H;       //height of cutted conus
    double dh = 0;      //leaf incline
    double alpha = 0;   //leaf slope
    int nleafs = 3;     //count of leafs in valve
    double ht = 0;      //thickness
};

struct OzakiTempl{
    double D = 25;                              //size
    double h = 11.0 + erf(2*(25.0-16.0));        //height of linear part
    double w = 3.0;                             //ozaki circle parameter
    double alpha = 1.0;                         //template lugs parameter
    double beta = 2.5 + 0.5*erf(2*(25.0-24.0)); //upper shift
    double s = 3;                               //suture_width
    double m = 2.0;                             //margins_width
    OzakiTempl(){}
    OzakiTempl(double size): D{size}, h{11.0 + erf(2*(size-16.0))}, beta{2.5 + 0.5*erf(2*(size-24.0))} {}
    OzakiTempl(double size, double h): D{size}, h{h}, beta{2.5 + 0.5*erf(2*(size-24.0))} {}
    OzakiTempl(double size, double h, double beta): D{size}, h{h}, beta{beta} {}
};

AniMesh generate_old_ozaki(int ozaki_size, double size);
AniMesh generate_ozaki(OzakiTempl p, double size);
AniMesh generateThubricarFlatLeaf(ThubricarValvePrms tvp, double size);
std::vector<Object3D> generateThubricarOpenValve(ThubricarValvePrms tvp, double size);
std::vector<Object3D> generateThubricarCloseValve(ThubricarValvePrms tvp, double size);

#endif //AVSYM_TEMPLATESCOLLECTION_H
