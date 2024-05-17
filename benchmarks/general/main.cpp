#include "model.h"
#include "examples.h"

#include <filesystem>

#if __cplusplus >= 201703L
using namespace std::filesystem;
#else
using namespace std::experimental::filesystem;
#endif

int main_test_io(int argc, char* argv[]){
    using namespace World3d;
    double R = 1, mesh_h = 0.2;
    int test_id = 1;
    if (test_id == 0){
        AniMesh am = generate_half_sphere(R, mesh_h);
        Mesh m = convert_to_Mesh(am, "v:boundary_lbl");
        std::cout << "#P = " << m.number_of_vertices() << " #E = " << m.number_of_edges() << " #F = " << m.number_of_faces() << std::endl;
        Object3D obj{m};
        obj.name = "half_sphere";

        MetaTriMesh mtm(obj.m_mesh);
        mtm.print_content_info();
        mtm.readAllTagsFromMesh(obj.m_mesh, {"v:point"});
        VTK_File().write("sphere.vtk", mtm);
        obj.save("sphere_0.vtk");
    } else if (test_id == 1) {
        MetaTriMesh mtm;
        VTK_File().read("sphere.vtk", mtm);
        mtm.print_content_info();
        auto m = mtm.convertToMesh();
        std::cout << "#P = " << m.number_of_vertices() << " #E = " << m.number_of_edges() << " #F = " << m.number_of_faces() << std::endl;
        Object3D obj{m};

        bool status = true;
        auto check_status = [&status]() { if (!status) std::cout << "Error" << std::endl; };
        status = mtm.createTagOnMesh<F_ind, Vector>(&obj.m_mesh, "f:normal"); check_status();
        status = mtm.createTagOnMesh<V_ind, int>(&obj.m_mesh, "v:boundary_lbl"); check_status();
        status = mtm.createTagOnMesh<V_ind, Vector>(&obj.m_mesh, "v:force"); check_status();
        status = mtm.createTagOnMesh<V_ind, Point>(&obj.m_mesh, "v:next_x"); check_status();
        status = mtm.createTagOnMesh<V_ind, Point>(&obj.m_mesh, "v:x"); check_status();

        auto c0_it = obj.m_mesh.add_property_map<F_ind, double>("f:c0", 5.95e6);
        auto c1_it = obj.m_mesh.add_property_map<F_ind, double>("f:c1", 1.48e-3);
        auto c2_it = obj.m_mesh.add_property_map<F_ind, double>("f:c2", 1.);
        auto f_it = obj.m_mesh.add_property_map<F_ind, Vector>("f:fiber_f").first;
        auto s_it = obj.m_mesh.add_property_map<F_ind, Vector>("f:fiber_s").first;
        for (auto f: obj.m_mesh.faces()){
            auto v = vert_around(obj.m_mesh, f);
            auto c = ((obj.m_x[v[0]] - CGAL::ORIGIN) + (obj.m_x[v[1]] - CGAL::ORIGIN) + (obj.m_x[v[2]] - CGAL::ORIGIN)) / 3;
            auto n = obj.m_normal[f];
            auto f_vec = Vector(c[2], 0, -c[0]);
            f_vec = f_vec - (f_vec * n)*n;
            f_vec /= sqrt(f_vec.squared_length());
            auto s_vec = CGAL::cross_product(f_vec, n);
            f_it[f] = f_vec; s_it[f] = s_vec;
        }
        // status = mtm.createTagOnMesh<V_ind, Point>(&obj.m_mesh, "v:unknown_data"); check_status();
        
        VTK_File().write("half_sphere_txt.vtk", MetaTriMesh(obj.m_mesh).readAllTagsFromMesh(obj.m_mesh, {"v:point"}));
        VTK_File().write("half_sphere_bin.vtk", MetaTriMesh(obj.m_mesh).readAllTagsFromMesh(obj.m_mesh, {"v:point"}), true);
        TMH_File().write_txt("half_sphere_txt.tmh", MetaTriMesh(obj.m_mesh).readAllTagsFromMesh(obj.m_mesh, {"v:point"}));
        TMH_File().write_bin("half_sphere_bin.tmh", MetaTriMesh(obj.m_mesh).readAllTagsFromMesh(obj.m_mesh, {"v:point"}));
    }
    return 0;
}

int main_test_energy_io(int argc, char* argv[]){
    using namespace World3d;

    if (argc <= 1) {
        std::cerr << "Usage: " << argv[0] << " path_to_energy_file" << std::endl; 
        return -1;
    }
    path energy_file(argv[1]);
    if (energy_file.extension().string() != ".energy"){
        std::cerr << "Energy file should have .energy extension" << std::endl;
        return -1;
    }
    Energy_Parser p;
    p.setCodeFromFile(energy_file.string());
    p.LexicalAnalysis();
    //p.debud_print_lexical_buffer(); std::cout << "\n";
    p.SyntaxAnalysis();
    p.print_parsed_state(std::cout) << "\n";
    p.Simplify();
    p.print_parsed_state(std::cout) << "\n";
    std::vector<double> par = {5.95e6, 1.48e-3, 6.75e-4};
    std::vector<double> inv = {2.0, 1.0, 1.0};
    using namespace casadi;
    SXVector par_s = {SX::sym("c0"), SX::sym("c1"), SX::sym("c2")};
    SXVector inv_s = {SX::sym("I[1]"), SX::sym("I[2]"), SX::sym("I[f]")};
    std::cout << "W(par, inv) = " << p.Evaluate(par, inv) << std::endl; 
    std::cout << "symbol W(par, inv) = " << SX::simplify(p.Evaluate(par_s, inv_s)) << std::endl;

    return 0;
}

int main(int argc, char* argv[]){
    if (argc == 1 || argv[1][0] == '-')
        return Model().parse(argc, const_cast<const char**>(argv)).execute(argc, argv);
    else if (!strcmp(argv[1], "stretch_rectangle")){
        auto s = stretch_rectangle().init(argc-2, argv+2);
        std::cout << "Generate " << argv[1] << " with params { "; s.print_state(std::cout) << " }" << std::endl; 
        return s.generate();
    } else if (!strcmp(argv[1], "inflation_circle")) { 
        auto s = inflation_circle().init(argc-2, argv+2);
        std::cout << "Generate " << argv[1] << " with params { "; s.print_state(std::cout) << " }" << std::endl; 
        return s.generate();
    } else if (!strcmp(argv[1], "stretch_cross_piece")){
        auto s = stretch_cross_piece().init(argc-2, argv+2);
        std::cout << "Generate " << argv[1] << " with params { "; s.print_state(std::cout) << " }" << std::endl; 
        return s.generate();
    } else{ 
        throw std::runtime_error("Faced unknow token = \"" + std::string(argv[1]) + "\". Try restart with --help option");
        return -1;
    }    
}