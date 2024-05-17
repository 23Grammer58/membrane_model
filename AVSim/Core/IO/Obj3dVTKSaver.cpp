//
// Created by alex on 29.04.2021.
//

#include "Obj3dVTKSaver.h"

using namespace std;

int World3d::VTKSaver::save(string file){
    const static int prec = std::numeric_limits<double>::digits10 + 1;
    std::ofstream f(file);
    if (!f.is_open()) {
        std::cout << "Can't open file \"" << file << "\"" << std::endl;
        return -2;
    }
    if (!m_obj) {
        std::cout << "Object doesn't specified. Please correct this using method setObj(...)" << std::endl;
        return -1;
    }
    if (!m_x) setX(m_obj->m_x);

    f << "# vtk DataFile Version 3.0\n"
         " \n"
         "ASCII\n"
         "DATASET UNSTRUCTURED_GRID\n"
         "POINTS " << m_obj->m_mesh.num_vertices() << " double\n";

    auto vit = std::max_element(m_obj->m_mesh.vertices().begin(), m_obj->m_mesh.vertices().end());
    std::vector<int> loc_ind(*vit+1, -1);
    {
        int i = 0;
        for (auto v: m_obj->m_mesh.vertices()) loc_ind[v] = i++;
    }
    for (auto v: m_obj->m_mesh.vertices()){
        auto p = m_x(v);
        f << std::setprecision(prec) << std::fixed << p[0] << " " << p[1] << " " << p[2] << "\n";
    }
    f << "CELLS " <<  m_obj->m_mesh.num_faces() << " " << 4 * m_obj->m_mesh.num_faces() << "\n";
    for (auto fc: m_obj->m_mesh.faces()){
        auto v = vert_around(m_obj->m_mesh, fc);
        f << "3 " << loc_ind[v[0]] << " " << loc_ind[v[1]] << " " << loc_ind[v[2]] << "\n";
    }
    f << "CELL_TYPES " << m_obj->m_mesh.num_faces() << "\n";
    const int VTK_TRIANGLE  = 5;
    for (int i = 0; i < m_obj->m_mesh.num_faces(); ++i) f << VTK_TRIANGLE << "\n";
    bool has_cell_data = false;
    bool has_point_data = false;
    for (auto it = m_tags.begin(); it != m_tags.end() && (!has_cell_data || !has_point_data); ++it){
        switch (it->sp){
            case FACE: has_cell_data = true; break;
            case NODE: has_point_data = true; break;
        }
    }
    auto printData = [&tags = m_tags, obj = m_obj](Sparsity sparsity, std::ostream& out){
        auto getType = [](DType dtype)->string{
            switch (dtype) {
                case REAL_TP: return "double";
                case INT_TP: return "int";
                default: return "void";
            }
        };
        auto value_saver = [](Sparsity sparsity, Object3D* obj) -> std::function<void(TagWrapper& tg, std::ostream& out)>{
            switch (sparsity){
                case NODE:
                    return [obj](TagWrapper& tg, std::ostream& out){
                        if (tg.dim <= 0) return;
                        if (tg.sp != NODE) return;
                        if (tg.tp == REAL_TP)
                            for (auto v: obj->m_mesh.vertices()){
                                auto res = tg.evaluate<double>(v);
                                out << std::setprecision(prec)  << std::fixed << res[0];
                                for (int i = 1; i < tg.dim; ++i)
                                    out << " " << std::setprecision(prec)  << std::fixed << res[i];
                                out << "\n";
                            }
                        else if (tg.tp == INT_TP)
                            for (auto v: obj->m_mesh.vertices()){
                                auto res = tg.evaluate<int>(v);
                                out << res[0];
                                for (int i = 1; i < tg.dim; ++i)
                                    out << " " << res[i];
                                out << "\n";
                            }
                    };
                case EDGE:
                    return [obj](TagWrapper& tg, std::ostream& out){
                        if (tg.dim <= 0) return;
                        if (tg.sp != EDGE) return;
                        if (tg.tp == REAL_TP)
                            for (auto v: obj->m_mesh.edges()){
                                auto res = tg.evaluate<double>(v);
                                out << std::setprecision(prec)  << std::fixed << res[0];
                                for (int i = 1; i < tg.dim; ++i)
                                    out << " " << std::setprecision(prec)  << std::fixed << res[i];
                                out << "\n";
                            }
                        else if (tg.tp == INT_TP)
                            for (auto v: obj->m_mesh.edges()){
                                auto res = tg.evaluate<int>(v);
                                out << res[0];
                                for (int i = 1; i < tg.dim; ++i)
                                    out << " " << res[i];
                                out << "\n";
                            }
                    };
                case FACE:
                    return [obj](TagWrapper& tg, std::ostream& out){
                        if (tg.dim <= 0) return;
                        if (tg.sp != FACE) return;
                        if (tg.tp == REAL_TP)
                            for (auto v: obj->m_mesh.faces()){
                                auto res = tg.evaluate<double>(v);
                                out << std::setprecision(prec)  << std::fixed << res[0];
                                for (int i = 1; i < tg.dim; ++i)
                                    out << " " << std::setprecision(prec)  << std::fixed << res[i];
                                out << "\n";
                            }
                        else if (tg.tp == INT_TP)
                            for (auto v: obj->m_mesh.faces()){
                                auto res = tg.evaluate<int>(v);
                                out << res[0];
                                for (int i = 1; i < tg.dim; ++i)
                                    out << " " << res[i];
                                out << "\n";
                            }
                    };
                case HALFEDGE:
                    return [obj](TagWrapper& tg, std::ostream& out){
                        if (tg.dim <= 0) return;
                        if (tg.sp != HALFEDGE) return;
                        if (tg.tp == REAL_TP)
                            for (auto v: obj->m_mesh.halfedges()){
                                auto res = tg.evaluate<double>(v);
                                out << std::setprecision(prec)  << std::fixed << res[0];
                                for (int i = 1; i < tg.dim; ++i)
                                    out << " " << std::setprecision(prec)  << std::fixed << res[i];
                                out << "\n";
                            }
                        else if (tg.tp == INT_TP)
                            for (auto v: obj->m_mesh.halfedges()){
                                auto res = tg.evaluate<int>(v);
                                out << res[0];
                                for (int i = 1; i < tg.dim; ++i)
                                    out << " " << res[i];
                                out << "\n";
                            }
                    };
                default: throw std::runtime_error("Reached unreacheable code");
            }
        };
        auto print_internal = value_saver(sparsity, obj);
        for (auto& t: tags){
            if (t.sp != sparsity) continue;
            out << "SCALARS " << t.name << " " << getType(t.tp) << " " << t.dim << "\n"
                << "LOOKUP_TABLE default\n";
            print_internal(t, out);
        }
    };
    if (has_cell_data){
        f << "CELL_DATA " << m_obj->m_mesh.num_faces() << "\n";
        printData(FACE, f);
    }
    if (has_point_data){
        f << "POINT_DATA " << m_obj->m_mesh.num_vertices() << "\n";
        printData(NODE, f);
    }

    f.close();

    return 0;
}