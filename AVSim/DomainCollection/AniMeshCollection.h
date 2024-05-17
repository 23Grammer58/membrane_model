//
// Created by alex on 01.04.2021.
//

#ifndef AORTIC_VALVE_ANIMESHGENERATOR_H
#define AORTIC_VALVE_ANIMESHGENERATOR_H

#include "AVSim/Core/TriangularMeshHelpers.h"

AniMesh generate_eigth_sphere_part(double R0, double size);
AniMesh generate_half_sphere(double R, double size);
AniMesh generate_EllipsoidPart(double a, double b, double d, double size);
AniMesh generate_cilinder_part(double R0, double l, double phi, double size);
AniMesh generate_circle(double R, double size);
AniMesh generate_circle_with_coaxial(double R, double theta, double size, double K = 0);
AniMesh generate_half_circle(double R, double size);
AniMesh generate_annulus(double R1, double R2, double size);
AniMesh generate_cutted_annulus(double R1, double R2, double size);
AniMesh generate_holey_rectangle(double a, double b, double size, const std::vector<std::pair<std::array<double, 2>, double>>& holes);
AniMesh generate_half_splited_rectangle(double a, double b, double size);
AniMesh generate_rectangle(double a, double b, double size, int fragments = 1);
AniMesh generate_dual_aligned_rectangle(double a, double b, double size);
AniMesh generate_HalfEllipseWithRectangle(double H1, double H2, double l, double size);
AniMesh generate_ideal_mesh(double R, int depth);
AniMesh generate_cross_piece(double a, double ac, double r, double size);
AniMesh generate_cross_piece(std::array<double, 2> a, std::array<double, 2> ac, std::array<double, 4> r, double size);

#endif //AORTIC_VALVE_ANIMESHGENERATOR_H
