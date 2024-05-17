//
// Created by Liogky Alexey on 19.05.2022.
//

#ifndef AVSYM_AVMESH_H
#define AVSYM_AVMESH_H

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh/IO.h>

namespace World3d {
    using DReal = double;

    typedef CGAL::Simple_cartesian<DReal>::Point_3 Point;
    typedef CGAL::Simple_cartesian<DReal>::Point_2 Point_2;
    typedef CGAL::Simple_cartesian<DReal>::Vector_3 Vector;
    typedef CGAL::Simple_cartesian<DReal>::Vector_2 Vector_2;
    typedef CGAL::Surface_mesh <Point> Mesh;
    typedef Mesh::Vertex_index V_ind;
    typedef Mesh::Face_index F_ind;
    typedef Mesh::Edge_index E_ind;
    typedef Mesh::Halfedge_index HE_ind;
};

#endif //AVSYM_AVMESH_H
