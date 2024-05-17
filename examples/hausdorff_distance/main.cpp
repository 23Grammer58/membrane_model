//
// Created by alex on 31.01.2022.
//

#include "AVSim/Core/Object3D.h"
#include "AVSim/Core/IO/Obj3dVTKSaver.h"

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

typedef CGAL::AABB_face_graph_triangle_primitive<World3d::Mesh> Primitive;
typedef CGAL::AABB_traits<CGAL::Simple_cartesian<World3d::DReal>, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;

static bool ObReadTag(World3d::Object3D& obj, std::string filename, std::string tag){
    std::ifstream ob(filename, std::ios::binary);
    if (!ob) return false;
    auto m_x = obj.m_mesh.add_property_map<World3d::V_ind, World3d::Point>(tag);
    for (auto& v: obj.m_mesh.vertices()){
        std::array<double, 3> P;
        ob.read(reinterpret_cast<char *> ((P.data())), 3 * sizeof(double));
        m_x.first[v] = World3d::Point(P[0], P[1], P[2]);
    }
    ob.close();
    return true;
}

int main(int argc, const char* argv[]){
    std::string mesh_file1, mesh_file2;
    std::string out_mesh = "out.vtk", tag_name = "v:hausdorff";
    std::string x1_tag;
    std::vector<std::pair<std::string, std::string>> extra_point_tags;
    bool quiet = false;
//    std::string x2_tag;

    struct CmdInfo{
        std::string shrt_cmd;
        std::string lng_cmd;
        std::string discr;
        std::function<int(const char**, int maxcnt)> handler; //return cnt of read element
        bool optional = true;
        bool was_set = false;
        int operator()(const char** s, int maxcnt){
            auto n = handler(s, maxcnt);
            if (n >= 0) was_set = true;
            return n;
        }
        CmdInfo(std::string shrt_cmd, std::string lng_cmd, std::string discr, std::function<int(const char**, int maxcnt)> handler, bool optional = true):
                shrt_cmd{shrt_cmd}, lng_cmd{lng_cmd}, discr{discr}, optional{optional}, handler{handler}
        {}
    };
    std::vector<CmdInfo> vcmd;
    std::map<std::string, std::size_t> cmd_handler;
    auto add_cmd = [&cmd_handler, &vcmd](const CmdInfo& x) {
        vcmd.push_back(x);
        auto it = &vcmd.back();
        if (!x.lng_cmd.empty())
            cmd_handler.insert({x.lng_cmd, it - vcmd.data()});
        if (!x.shrt_cmd.empty())
            cmd_handler.insert({x.shrt_cmd, it - vcmd.data()});
    };
    auto print_help_msg = [&vcmd](std::string prefix = "", bool over = false, int over_val = 0){
        std::cout << prefix << " Command line options: \n";
        int max_shrt_sz = 0, max_lng_sz = 0;
        for (auto& i: vcmd)
            max_shrt_sz = std::max(max_shrt_sz, static_cast<int>(i.shrt_cmd.size())), max_lng_sz = std::max(max_lng_sz, static_cast<int>(i.lng_cmd.size()));
        max_shrt_sz += (max_shrt_sz > 0);
        max_lng_sz += (max_lng_sz > 0);

        for (auto& i: vcmd) {
            std::cout << prefix << "    " << i.shrt_cmd << ((!i.shrt_cmd.empty() && !i.lng_cmd.empty()) ? "," : "");
            for (int j = i.shrt_cmd.size(); j < max_shrt_sz; ++j) std::cout << " ";
            std::cout << i.lng_cmd;
            for (int j = i.lng_cmd.size(); j < max_lng_sz; ++j) std::cout << " ";
            std::cout << "<" << (i.optional ? "[Optional] " : "") <<  i.discr << ">\n";
        }
        if (over) exit(over_val);
    };
    auto read_str_func = [](std::string& s){ return [&arg = s](const char** s, int maxcnt){
        if (maxcnt < 1) throw std::runtime_error("Wrong count of args");
        arg = s[0];
        return 1;
    }; };
    add_cmd(CmdInfo{"-s", "--AVSim", "path to mesh on that distance will be set", read_str_func(mesh_file1), false});
    add_cmd(CmdInfo{"-t", "--target", "path to mesh to that distance will be computed", read_str_func(mesh_file2), false});
    add_cmd(CmdInfo{"-o", "--out", "name of out mesh [default=\"out.vtk\"]", read_str_func(out_mesh)});
    add_cmd(CmdInfo{"-g", "--tag", "name of tag to be saved [default=\"v:hausdorff\"]", read_str_func(tag_name)});
    add_cmd(CmdInfo{"-x1", "--x1_tag", "name of tag to be considered as coordinate tag on AVSim mesh [default: considered from mesh format]", read_str_func(x1_tag)});
    add_cmd(CmdInfo{"-h", "--help", "print help message and exit program", [&print_help_msg](const char** s, int maxcnt) { print_help_msg("", true); return 0; }} );
    add_cmd(CmdInfo{"-e", "--extra_tag", "read and save extra point tag path and name", [&t = extra_point_tags](const char** s, int maxcnt){
        if (maxcnt < 1) throw std::runtime_error("Wrong count of args");
        if (maxcnt >= 2 && s[1][0] != '-')
            return t.emplace_back(s[0], s[1]), 2;
        else return t.emplace_back(s[0], "tag_" + std::to_string(t.size())), 1;
    } });
    add_cmd(CmdInfo{"-q", "--quiet", "do not print information", [&quiet](const char** s, int maxcnt){ quiet = true; return 0; }} );
//    add_cmd(CmdInfo{"-x2", "--x2_tag", "name of tag to be considered as coordinate tag on target mesh [default: considered from mesh format]", read_str_func(x2_tag)});

    std::string run_cmd = argv[0];
    std::string env_shift = run_cmd.substr(0, run_cmd.find_last_of("/\\")+1);

    for (int i = 1; i < argc; i++) {
        auto _it = cmd_handler.find(argv[i]);
        if (_it != cmd_handler.end()){
            auto& it = vcmd[_it->second];
            auto n = it(argv + i + 1, argc - (i + 1));
            if (n < 0) {
                std::cerr << "Can't process \"" << std::string(argv[i]) << "\"\n";
                print_help_msg();
                throw std::runtime_error(
                        "Processor of command \"" + std::string(argv[i]) + "\" return error (" + std::to_string(n) +
                        ")");
            }
            else
                i += n;
        } else {
            std::cerr << "Can't find command \"" << std::string(argv[i]) << "\"\n";
            print_help_msg("", true);
        }
    }
    for (auto i: vcmd){
        if (!i.optional && !i.was_set) throw std::runtime_error("Value \"" + i.lng_cmd + "\" must be set");
    }

//    if (mesh_file1[0] != '/')
//        mesh_file1 = env_shift + mesh_file1;
//    if (mesh_file2[0] != '/')
//        mesh_file2 = env_shift + mesh_file2;
//    if (out_mesh[0] != '/')
//        out_mesh = env_shift + out_mesh;

    World3d::Object3D sm(mesh_file1);
    World3d::Object3D tm(mesh_file2);

    for (auto& t: extra_point_tags) ObReadTag(sm, t.first, t.second);

    auto x1 = sm.m_x;
//    auto x2 = tm.m_x;
    if (!x1_tag.empty()){
        auto it = sm.m_mesh.property_map<World3d::V_ind, World3d::Point>(x1_tag);
        if (it.second) x1 = it.first;
        else throw std::runtime_error("Src mesh doesn't have tag = \"" + x1_tag + "\"");
    }
//    if (!x2_tag.empty()){
//        auto it = tm.m_mesh.property_map<World3d::V_ind, World3d::Point>(x2_tag);
//        if (it.second) x2 = it.first;
//        else throw std::runtime_error("Target mesh doesn't have tag = \"" + x2_tag + "\"");
//    }
    auto dist_tag = sm.m_mesh.add_property_map<World3d::V_ind, double>(tag_name).first;
    Tree tree(faces(tm.m_mesh).first, faces(tm.m_mesh).second, tm.m_mesh);

    for (auto v: sm.m_mesh.vertices()){
        auto query = x1[v];
        auto pp = tree.closest_point_and_primitive(query);
        auto dir = query - pp.first;
        double dist = sqrt(dir.squared_length());
        if (dir * tm.m_normal[pp.second] < 0) dist *= -1;
        dist_tag[v] = dist;
    }

    std::set<std::string> save_tags;
    save_tags.insert(tag_name);
    save_tags.insert("v:point");
    if (!x1_tag.empty()) save_tags.insert(x1_tag);
    for (auto& i: extra_point_tags) save_tags.insert(i.second);
    save_tags.insert("v:boundary_lbl");
    save_tags.insert("f:normal");
    save_tags.erase("v:x");

    World3d::VTKSaver vs;
    vs.setObj(sm);
    vs.setX(sm.m_x);
    auto nodetg = sm.m_mesh.properties<World3d::V_ind>();
    for (auto& n: nodetg) if (save_tags.count(n)) vs.setTagByName(n, n);
    auto facetg = sm.m_mesh.properties<World3d::F_ind>();
    for (auto& n: facetg) if (save_tags.count(n)) vs.setTagByName(n, n);
    vs.save(out_mesh);

    if (!quiet) {
        double max_dist = 0;
        for (auto v: sm.m_mesh.vertices()) max_dist = std::max(max_dist, abs(static_cast<double>(dist_tag[v])));
        std::cout << "Write res into " << out_mesh << "\n\thausdorff_dist = " << max_dist << std::endl;
    }

    return 0;
}

