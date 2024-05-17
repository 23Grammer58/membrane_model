//
// Created by alex on 01.04.2021.
//

#include "Helper.h"

int mark_vertices(AniMesh &am, int defcolor) {
    am.vertex_label.resize(am.vertices.size() / 3);
    return mark_vertices(am.vertices.size() / 3, am.vertex_label.data(), am.edges.size()/2, am.edges.data(), am.edge_label.data(), defcolor);
}

int mark_vertices(int nV, int *vertexmaterial, int nE, int *edge, int *edgematerial, int defcolor) {
    int i;
    for (i=0; i<nV; i++)
        vertexmaterial[i] = defcolor;
    for (i=0; i<nE; i++) {
        vertexmaterial[edge[2*i+0]-1] |= edgematerial[i];
        vertexmaterial[edge[2*i+1]-1] |= edgematerial[i];
    }
    return 0;
}

