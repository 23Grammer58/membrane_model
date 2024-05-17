//
// Created by alex on 01.04.2021.
//

#ifndef AORTIC_VALVE_HELPER_H
#define AORTIC_VALVE_HELPER_H
#include "../TriangularMeshHelpers.h"
#include "ForAni3dFrtPrm.h"

int mark_vertices(int nV, int *vertexmaterial, int nE, int *edge, int *edgematerial, int defcolor = 0);
int mark_vertices(AniMesh& am, int defcolor = 0);

#endif //AORTIC_VALVE_HELPER_H
