//
// Created by alex on 01.04.2021.
//

#include "AVSim/Core/TriangularMeshHelpers.h"
#include <valarray>

static int mark_vertices(int nV, int *vertexmaterial, int nE, int *edge, int *edgematerial, int defcolor) {
    int i;
    for (i=0; i<nV; i++)
        vertexmaterial[i] = defcolor;
    for (i=0; i<nE; i++) {
        vertexmaterial[edge[2*i+0]-1] |= edgematerial[i];
        vertexmaterial[edge[2*i+1]-1] |= edgematerial[i];
    }
    return 0;
}

AniMesh generate_ideal_mesh(double R, int depth){
    using namespace std;
    vector<array<int,    3>> faceE[2];
    vector<array<int,    2>> edge[2];
    std::vector<valarray<double>> vertex;
    vector<int> facematerial, edgemateraial[2], vertexmaterial;

    int cur = 0;

    {
        int coef = pow(2, depth);
        int coef_sq = coef * coef;
        edge[depth % 2].reserve(3 * coef + 9 * coef_sq), edge[(depth + 1) % 2].reserve(3 * coef / 2 + 9 * coef_sq / 4);
        faceE[depth % 2].reserve(coef_sq * 6), faceE[(depth + 1) % 2].reserve(coef_sq * 6 / 2);
        vertex.reserve(3 * (coef + coef_sq) + 1);
        vertex.push_back(valarray<double>{0, 0, 0});
        for (int i = 0; i < 6; ++i)
            vertex.push_back({cos(M_PI / 6 + M_PI / 3 * i), sin(M_PI / 6 + M_PI / 3 * i), 0});
        for (int i = 1; i <= 6; ++i)
            edge[0].push_back({0, i}), edgemateraial[0].push_back(0);
        for (int i = 1; i <= 5; ++i) {
            edge[0].push_back({i, i + 1}), edgemateraial[0].push_back(1);
            faceE[0].push_back({i - 1, i, i + 5});
        }
        edge[0].push_back({6, 1}), edgemateraial[0].push_back(1);
        faceE[0].push_back({5, 0, 11});
        for (cur = 0; cur < depth; ++cur) {
            int off = vertex.size();
            cout << "\t       nV = " << vertex.size() << ", nF = " << faceE[cur % 2].size() << ", nE = "
                 << edge[cur % 2].size() << endl;
            edge[(cur + 1) % 2].clear(), edgemateraial[(cur + 1) % 2].clear();
            faceE[(cur + 1) % 2].clear();
            for (unsigned i = 0; i < edge[cur % 2].size(); ++i) {
                int id = vertex.size();
                vertex.push_back((vertex[edge[cur % 2][i][0]] + vertex[edge[cur % 2][i][1]]) / 2);
                edge[(cur + 1) % 2].push_back({edge[cur % 2][i][0], id}), edgemateraial[(cur + 1) % 2].push_back(
                        edgemateraial[cur % 2][i]);
                edge[(cur + 1) % 2].push_back({edge[cur % 2][i][1], id}), edgemateraial[(cur + 1) % 2].push_back(
                        edgemateraial[cur % 2][i]);
            }
            for (auto &i: faceE[cur % 2]) {
                int qq = edge[(cur + 1) % 2].size();
                for (int k = 0; k < 3; ++k) {
                    int com = edge[cur % 2][i[k]][0];
                    if (com != edge[cur % 2][i[(k + 1) % 3]][0] && com != edge[cur % 2][i[(k + 1) % 3]][1])
                        com = edge[cur % 2][i[k]][1];
                    array<int, 2> sv{2 * i[k], 2 * i[(k + 1) % 3]};
                    for (int l = 0; l < 2; ++l)
                        if (edge[(cur + 1) % 2][sv[l]][0] != com && edge[(cur + 1) % 2][sv[l]][1] != com)
                            ++sv[l];
                    faceE[(cur + 1) % 2].push_back({sv[0], sv[1], static_cast<int> (edge[(cur + 1) % 2].size())});
                    edge[(cur + 1) % 2].push_back({off + i[k], off + i[(k + 1) % 3]}), edgemateraial[(cur + 1) %
                                                                                                     2].push_back(0);
                }
                faceE[(cur + 1) % 2].push_back({qq, qq + 1, qq + 2});
            }
        }
    }

    vector<int> face;
    face.reserve(faceE[cur%2].size()*3);
    for (const auto& i: faceE[cur%2]){
        array<int, 3> vis {edge[cur%2][i[0]][0], edge[cur%2][i[0]][1], -1};
        int com = edge[cur%2][i[1]][0];
        if (com == vis[0] || com == vis[1]) com = edge[cur%2][i[1]][1];
        vis[2] = com;
        valarray<double> p1 = vertex[vis[1]] - vertex[vis[0]],
                p2 = vertex[vis[2]] - vertex[vis[0]];
        if (p1[0]*p2[1] - p1[1]*p2[0] < 0)
            swap(vis[1], vis[0]);
        for (int k = 0; k < 3; ++k)
            face.push_back(vis[k]+1);
    }
//    for (int i = 0; i < faceE[cur%2].size(); ++i){
//        valarray p1 = vertex[face[3*i+1]-1] - vertex[face[3*i]-1],
//                 p2 = vertex[face[3*i+2]-1] - vertex[face[3*i]-1];
//        if (p1[0]*p2[1] - p1[1]*p2[0] < 0)
//            swap(face[3*i+1], face[3*i+2]);
//    }
    vector<double> _vertx(3 * vertex.size());
    for (unsigned i = 0; i < vertex.size(); ++i)
        for (int k = 0; k < 3; ++k)
            _vertx[i*3 + k] = R*vertex[i][k];
    facematerial.resize(face.size());
    vertexmaterial.resize(vertex.size());
    vector<int> _edge;
    _edge.reserve(2 * edge[cur%2].size());
    for (auto& i: edge[cur%2]) {
        _edge.push_back(i[0]+1);
        _edge.push_back(i[1]+1);
    }
    int nV = vertex.size(), nF = faceE[cur%2].size(), nE = edge[cur%2].size();
    (void) nF;
    mark_vertices(nV, vertexmaterial.data(), nE, _edge.data(), edgemateraial[cur%2].data(), 0);

    AniMesh am;
    am.vertices = _vertx;
    am.faces = face;
    am.edges = _edge;
    am.vertex_label = vertexmaterial;

    return am;
}
