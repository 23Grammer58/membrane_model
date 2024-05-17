#include "MetaTriFormat.h"
#include <limits>

using namespace World3d;
using TMH = TMH_File; 

TMH::Section TMH::get_section(const char* str){
    if (!strcmp(str, "POINTS")) return POINTS_SC;
    if (!strcmp(str, "EDGES")) return EDGES_SC;
    if (!strcmp(str, "CELLS")) return CELLS_SC;
    if (!strcmp(str, "POINT_DATA")) return POINT_DATA_SC;
    if (!strcmp(str, "EDGE_DATA")) return EDGE_DATA_SC;
    if (!strcmp(str, "CELL_DATA")) return CELL_DATA_SC;
    if (!strcmp(str, "HALFEDGE_DATA")) return HALFEDGE_DATA_SC;
    if (!strcmp(str, "LOAD")) return LOAD_SC;
    throw std::runtime_error("Faced unknown section \"" + std::string(str) + "\"");
}

const char* TMH::get_section_name(Section sec){
    switch(sec){
        case POINTS_SC: return "POINTS";
        case EDGES_SC: return "EDGES";
        case CELLS_SC: return "CELLS";
        case POINT_DATA_SC: return "POINT_DATA";
        case EDGE_DATA_SC: return "EDGE_DATA";
        case CELL_DATA_SC: return "CELL_DATA";
        case HALFEDGE_DATA_SC: return "HALFEDGE_DATA";
        case LOAD_SC: return "LOAD";
        default:
             throw std::runtime_error("Faced unknown section code \"" + std::to_string(int(sec)) + "\"");
    }
}

TMH::Section TMH::get_data_subsection_type(const char* str){
    if (!strcmp(str, "SCALARS")) return SCALARS_SSC;
    if (!strcmp(str, "CONSTANTS")) return CONSTANTS_SSC;
    if (!strcmp(str, "SPARSE_SCALARS")) return SPARSE_SCALARS_SSC;
    return NSECTIONS_SC;
}

const char* TMH::get_subsection_name(Section ssec){
    switch(ssec){
        case SCALARS_SSC: return "SCALARS";
        case CONSTANTS_SSC: return "CONSTANTS";
        case SPARSE_SCALARS_SSC: return "SPARSE_SCALARS";
        default: 
            throw std::runtime_error("Faced unknown subsection code \"" + std::to_string(int(ssec)) + "\"");
    }
}

bool TMH::is_binary_file(std::string file_name){
    char val[5] = "\0";
    std::ifstream in(file_name);
    if (in.eof()) return false;
    in.read(val, 4); val[4] = '\0';
    in.close();
    return strcmp(val, "btmh") == 0;
}

bool TMH::read(std::string file_name, MetaTriMesh& m, unsigned char status){ 
    m_is_binary_file = is_binary_file(file_name); m_file_name = file_name; m_title_text.resize(0); is_loaded_state = false;
    return m_is_binary_file ? read_bin(file_name, m, status) : read_txt(file_name, m, status);
}

static bool skip_empty(std::istream& in){
    in >> std::ws;
    while (!in.eof()){
        auto p = in.tellg();
        char c = '\0';
        in >> c;
        if (c == '#'){ //ignore comment line
            in.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
            in >> std::ws;  
        } else{ 
            in.seekg(p, std::ios_base::beg);
            return true;
        }
    }
    return false;
}

template<int N>
static bool read_mesh_graph_section_ascii(std::ifstream& in, MetaTriMesh::TStorage<N>& s, MetaTriMesh::DataType default_type, bool ignore_section = false){
    MetaTriMesh::TStorage<N> to_ignore;
    MetaTriMesh::TStorage<N>* p = ignore_section ? &to_ignore : &s;
    p->m_type = default_type;
    in >> p->m_size; skip_empty(in);
    p->resize(p->m_size);
    std::string last_prt;
    std::getline(in, last_prt);
    if (!last_prt.empty()){
        std::stringstream ss(last_prt);
        ss >> std::ws;
        std::string dtype_str;
        if (!ss.eof() && ss >> dtype_str && !dtype_str.empty())
            p->m_type = MetaTriMesh::get_type_by_name(dtype_str.c_str());
    }
    skip_empty(in);

    p->read_ascii_data(in, [](std::istream& in) { skip_empty(in); });
    return true;
}

template<int N>
static bool read_mesh_graph_section_bin(std::ifstream& in, MetaTriMesh::TStorage<N>& s, bool ignore_section = false){
    MetaTriMesh::TStorage<N> to_ignore;
    MetaTriMesh::TStorage<N>* p = ignore_section ? &to_ignore : &s;
    unsigned long p_size = 0;
    in.read(reinterpret_cast<char*>(&p_size), sizeof(p_size));
    p->m_size = p_size;
    char dtype_str = MetaTriMesh::DataType::NONE_TP;
    in.read(reinterpret_cast<char*>(&dtype_str), sizeof(char));
    p->m_type = static_cast<MetaTriMesh::DataType>(dtype_str);
    p->read_binary_data(in);
    return true;
}

static std::string read_mesh_data_name_ascii(std::istream& in){
    std::string name;
    auto p = in.tellg();
    char c = '\0';
    in >> c;
    if (c == '\"')
        std::getline(in, name, '\"');
    else{
        in.seekg(p);
        in >> name;
    }
    return name;
}

static std::string read_path(std::istream& in){
    skip_empty(in);
    if (in.eof())
        throw std::runtime_error("Can't read filename");
    return read_mesh_data_name_ascii(in);
}

static bool read_mesh_data_ascii(std::ifstream& in, std::map<std::string, MetaTriMesh::MeshValue>& s, bool ignore_section = false){
    size_t size = 0;
    in >> size; skip_empty(in);
    bool continue_read = true;
    while (continue_read && !in.eof()){
        auto pin = in.tellg();
        std::string scal_str;
        in >> scal_str;
        TMH::Section tp = TMH::get_data_subsection_type(scal_str.c_str());
        switch(tp){
            case TMH::SCALARS_SSC:{
                MetaTriMesh::MeshValue val;
                val.m_type = MetaTriMesh::MeshValue::CONTIGUOUS;
                std::string name = read_mesh_data_name_ascii(in);
                std::string dtype_str; size_t dim = 0;
                in >> dtype_str >> dim; skip_empty(in);
                val.m_values.m_type = MetaTriMesh::get_type_by_name(dtype_str.c_str());
                val.m_values.m_dim = dim;
                val.m_values.m_size = size;
                val.m_values.read_ascii_data(in, [](std::istream& in) { skip_empty(in); });
                if (!ignore_section) s[name] = val;
                break;
            }
            case TMH::CONSTANTS_SSC:{
                MetaTriMesh::MeshValue val;
                val.m_type = MetaTriMesh::MeshValue::CONSTANT;
                std::string name = read_mesh_data_name_ascii(in);
                std::string dtype_str; size_t dim = 0;
                in >> dtype_str >> dim; skip_empty(in);
                val.m_values.m_type = MetaTriMesh::get_type_by_name(dtype_str.c_str());
                val.m_values.m_dim = dim;
                val.m_values.m_size = 1;
                val.m_values.read_ascii_data(in, [](std::istream& in) { skip_empty(in); });
                val.m_values.m_size = size;
                if (!ignore_section) s[name] = val;
                break;
            }
            case TMH::SPARSE_SCALARS_SSC:{
                MetaTriMesh::MeshValue val;
                val.m_type = MetaTriMesh::MeshValue::SPARSE;
                std::string name = read_mesh_data_name_ascii(in);
                std::string index_type_str, value_type_str;
                std::size_t nSparse = 0, dim = 0;
                in >> index_type_str >> nSparse >> value_type_str >> dim; skip_empty(in);
                val.m_sparse_indexes.m_type = MetaTriMesh::get_type_by_name(index_type_str.c_str());
                val.m_sparse_indexes.resize(nSparse);
                val.m_values.m_type = MetaTriMesh::get_type_by_name(value_type_str.c_str());
                val.m_values.m_dim = dim;
                val.m_values.m_size = 1;
                val.m_values.read_ascii_data(in, [](std::istream& in) { skip_empty(in); });
                val.m_values.resize(nSparse+1);
                for (std::size_t i = 0; i < nSparse; ++i){
                    MetaTriMesh::Value index = MetaTriMesh::Value::read_ascii_var(in, val.m_sparse_indexes.m_type);
                    val.m_sparse_indexes.setV(i, 0, index);
                    for (std::size_t k = 0; k < dim; ++k){
                        auto lcomp = MetaTriMesh::Value::read_ascii_var(in, val.m_values.m_type);
                        val.m_values.setV(i+1, k, lcomp);
                    }
                    skip_empty(in);
                }
                if (!ignore_section) s[name] = val;
                break;
            }
            default:{
                continue_read = false;
                in.seekg(pin);
            }
        }
    }
    return true;
}

static std::string read_string_bin(std::istream& in){
    unsigned size = 0;
    in.read(reinterpret_cast<char*>(&size), sizeof(unsigned));
    std::string res; res.resize(size);
    in.read(res.data(), size);
    return res;
}

static bool read_mesh_data_bin(std::ifstream& in, std::map<std::string, MetaTriMesh::MeshValue>& s, bool ignore_section = false){
    unsigned long size = 0;
    in.read(reinterpret_cast<char*>(&size), sizeof(size));
    bool continue_read = true;
    while (continue_read && !in.eof()){
        auto pin = in.tellg();
        char subsection = TMH::NSECTIONS_SC;
        in.read(reinterpret_cast<char*>(&subsection), sizeof(subsection));
        switch(subsection){
            case TMH::SCALARS_SSC:{
                auto name = read_string_bin(in);
                MetaTriMesh::MeshValue val;
                val.m_type = MetaTriMesh::MeshValue::CONTIGUOUS;
                char dtype_str = MetaTriMesh::DataType::NONE_TP;
                in.read(reinterpret_cast<char*>(&dtype_str), sizeof(char));
                val.m_values.m_type = static_cast<MetaTriMesh::DataType>(dtype_str);
                unsigned long dim = 0; 
                in.read(reinterpret_cast<char*>(&dim), sizeof(dim));
                val.m_values.m_dim = dim; val.m_values.resize(size); 
                val.m_values.read_binary_data(in);
                if (!ignore_section) s[name] = val;
                break;
            }
            case TMH::CONSTANTS_SSC:{
                auto name = read_string_bin(in);
                MetaTriMesh::MeshValue val;
                val.m_type = MetaTriMesh::MeshValue::CONSTANT;
                char dtype_str = MetaTriMesh::DataType::NONE_TP;
                in.read(reinterpret_cast<char*>(&dtype_str), sizeof(char));
                val.m_values.m_type = static_cast<MetaTriMesh::DataType>(dtype_str);
                unsigned long dim = 0; 
                in.read(reinterpret_cast<char*>(&dim), sizeof(dim));
                val.m_values.m_dim = dim; val.m_values.resize(1); 
                val.m_values.read_binary_data(in);
                if (!ignore_section) s[name] = val;
                break;
            }
            case TMH::SPARSE_SCALARS_SSC:{
                auto name = read_string_bin(in);
                MetaTriMesh::MeshValue val;
                val.m_type = MetaTriMesh::MeshValue::SPARSE;
                char index_type_str = MetaTriMesh::DataType::NONE_TP, 
                     value_type_str = MetaTriMesh::DataType::NONE_TP;
                unsigned long nSparse = 0, dim = 0;
                in.read(reinterpret_cast<char*>(&index_type_str), sizeof(char));
                val.m_sparse_indexes.m_type = static_cast<MetaTriMesh::DataType>(index_type_str);
                in.read(reinterpret_cast<char*>(&nSparse), sizeof(nSparse));
                val.m_sparse_indexes.resize(nSparse);
                in.read(reinterpret_cast<char*>(&value_type_str), sizeof(char));
                val.m_values.m_type = static_cast<MetaTriMesh::DataType>(value_type_str);
                in.read(reinterpret_cast<char*>(&dim), sizeof(dim));
                val.m_values.m_dim = dim; val.m_values.resize(nSparse + 1);
                val.m_values.read_binary_data(in);
                val.m_sparse_indexes.read_binary_data(in);
            }
            default: {
                continue_read = false;
                in.seekg(pin);
            }   
        }
    }
    return true;
}

bool TMH::read_txt(std::string file_name, MetaTriMesh& m, unsigned char status, unsigned char faced_sec){
    std::ifstream in(file_name);
    std::string welcome_str;
    std::getline(in, welcome_str);
    for (auto it = welcome_str.begin(); it != welcome_str.end(); ++it)
        if (!std::isspace(*it) && *it == '#') 
            break;
        else 
            throw std::runtime_error("This is not ascii \"*.tmh\" file");

    std::string title;
    std::getline(in, title);
    if (!is_loaded_state) m_title_text = title;
    skip_empty(in);

    std::string sec_str;
    while (!in.eof() && in >> sec_str){
        auto sec = get_section(sec_str.c_str());
        switch(sec){
            case POINTS_SC:{
                if ((faced_sec & IO_POINTS) && (status & IO_POINTS))
                    throw std::runtime_error("Allowed only one occurance of POINTS section");
                faced_sec |= (status & IO_POINTS);  
                read_mesh_graph_section_ascii<3>(in, m.points, MetaTriMesh::DOUBLE_TP, !(status & IO_POINTS)); 
                break;   
            }
            case EDGES_SC:{
                if ((faced_sec & IO_EDGES) && (status & IO_EDGES))
                    throw std::runtime_error("Allowed only one occurance of EDGES section");
                faced_sec |= (status & IO_EDGES);
                read_mesh_graph_section_ascii<2>(in, m.edges, MetaTriMesh::ULONG_TP, !(status & IO_EDGES));  
                break;  
            }
            case CELLS_SC:{
                if ((faced_sec & IO_CELLS) && (status & IO_CELLS))
                    throw std::runtime_error("Allowed only one occurance of CELLS section");
                faced_sec |= (status & IO_CELLS);
                read_mesh_graph_section_ascii<3>(in, m.cells, MetaTriMesh::ULONG_TP, !(status & IO_CELLS)); 
                break;   
            }
            case POINT_DATA_SC: {
                faced_sec |= (status & IO_POINT_DATA);
                read_mesh_data_ascii(in, m.point_data, !(status & IO_POINT_DATA)); 
                break;
            }
            case EDGE_DATA_SC: {
                faced_sec |= (status & IO_EDGE_DATA);
                read_mesh_data_ascii(in, m.edge_data, !(status & IO_EDGE_DATA)); 
                break;
            }
            case HALFEDGE_DATA_SC:{
                faced_sec |= (status & IO_HALFEDGE_DATA);
                read_mesh_data_ascii(in, m.halfedge_data, !(status & IO_HALFEDGE_DATA)); 
                break;
            }
            case CELL_DATA_SC: {
                faced_sec |= (status & IO_CELL_DATA);
                read_mesh_data_ascii(in, m.cell_data, !(status & IO_CELL_DATA)); 
                break;
            }
            case LOAD_SC:{
                auto load_file = read_path(in);
                if (!(status & IO_LOAD_FILES)) break;
                faced_sec |= IO_LOAD_FILES;
                if (load_file.empty()) throw std::runtime_error("Can't read module with empty module path");
                if (load_file[0] != '/')
                    load_file = file_name.substr(0, file_name.find_last_of('/')) + '/' + load_file;
                bool is_binary = is_binary_file(load_file);
                bool is_loaded_state_tmp = is_loaded_state; is_loaded_state = true;
                if (is_binary) read_bin(load_file, m, status, faced_sec);  
                else read_txt(load_file, m, status, faced_sec);
                is_loaded_state = is_loaded_state_tmp;
                break;
            }
            default:
                throw std::runtime_error("Faced unreachable code");
        }
        skip_empty(in);
    }
    return true;
}

void write_string_bin(std::ostream& out, const std::string& str){
    unsigned size = str.size();
    out.write(reinterpret_cast<char*>(&size), sizeof(unsigned));
    out.write(str.data(), size);;
}

bool TMH::read_bin(std::string file_name, MetaTriMesh& m, unsigned char status, unsigned char faced_sec){
    std::ifstream in(file_name, std::ios_base::binary);
    char welcome[5] = "\0";
    in.read(welcome, 4);
    unsigned major_version = 0, minor_version = 0;
    in.read(reinterpret_cast<char*>(&major_version), sizeof(unsigned));
    in.read(reinterpret_cast<char*>(&minor_version), sizeof(unsigned));
    auto title = read_string_bin(in);
    if (!is_loaded_state) m_title_text = title;
    char sec_str = static_cast<char>(NSECTIONS_SC);
    while(!in.eof() && in.read(&sec_str, sizeof(char))){
        switch (sec_str){
            case POINTS_SC:{
                if ((faced_sec & IO_POINTS) && (status & IO_POINTS))
                    throw std::runtime_error("Allowed only one occurance of POINTS section");
                faced_sec |= (status & IO_POINTS);
                read_mesh_graph_section_bin<3>(in, m.points, !(status & IO_POINTS));
                break;
            }
            case EDGES_SC:{
                if ((faced_sec & IO_EDGES) && (status & IO_EDGES))
                    throw std::runtime_error("Allowed only one occurance of EDGES section");
                faced_sec |= (status & IO_EDGES);
                read_mesh_graph_section_bin<2>(in, m.edges, !(status & IO_EDGES));
                break;
            }
            case CELLS_SC:{
                if ((faced_sec & IO_CELLS) && (status & IO_CELLS))
                    throw std::runtime_error("Allowed only one occurance of CELLS section");
                faced_sec |= (status & IO_CELLS);
                read_mesh_graph_section_bin<3>(in, m.cells, !(status & IO_CELLS));
                break;
            }
            case POINT_DATA_SC: {
                faced_sec |= (status & IO_POINT_DATA);
                read_mesh_data_bin(in, m.point_data, !(status & IO_POINT_DATA)); 
                break;
            }
            case EDGE_DATA_SC: {
                faced_sec |= (status & IO_EDGE_DATA);
                read_mesh_data_bin(in, m.edge_data, !(status & IO_EDGE_DATA)); 
                break;
            }
            case HALFEDGE_DATA_SC:{
                faced_sec |= (status & IO_HALFEDGE_DATA);
                read_mesh_data_bin(in, m.halfedge_data, !(status & IO_HALFEDGE_DATA)); 
                break;
            }
            case CELL_DATA_SC: {
                faced_sec |= (status & IO_CELL_DATA);
                read_mesh_data_bin(in, m.cell_data, !(status & IO_CELL_DATA)); 
                break;
            }
            case LOAD_SC:{
                auto load_file = read_string_bin(in);
                if (!(status & IO_LOAD_FILES)) break;
                faced_sec |= IO_LOAD_FILES;
                if (load_file.empty()) throw std::runtime_error("Can't read module with empty module path");
                if (load_file[0] != '/')
                    load_file = file_name.substr(0, file_name.find_last_of('/')) + '/' + load_file;
                bool is_binary = is_binary_file(load_file);
                bool is_loaded_state_tmp = is_loaded_state; is_loaded_state = true;
                if (is_binary) read_bin(load_file, m, status, faced_sec);  
                else read_txt(load_file, m, status, faced_sec);
                is_loaded_state = is_loaded_state_tmp;
                break;
            }
            default:
                throw std::runtime_error("Faced unknown section with code = \"" + std::to_string(int(sec_str)) + "\"");
        }
    }
    return true;
}

bool TMH::write_bin(std::string file_name, const MetaTriMesh& m, unsigned char status, const std::vector<std::string>& loaded_paths){
    m_file_name = file_name; m_is_binary_file = true;
    std::ofstream out(file_name, std::ios_base::binary);
    const char* welcome = "btmh";
    unsigned major_version = 1, minor_version = 0;
    out.write(welcome, 4);
    out.write(reinterpret_cast<const char*>(&major_version), sizeof(unsigned));
    out.write(reinterpret_cast<const char*>(&minor_version), sizeof(unsigned));
    write_string_bin(out, m_title_text);
    auto save_mesh_graph_section_bin = [status, &out](auto IO_SEC, auto SEC, auto& section){
        if (!(status & IO_SEC) || section.size() == 0) return;
        char sec_str = SEC;
        out.write(reinterpret_cast<const char*>(&sec_str), sizeof(sec_str));
        unsigned long section_size = section.m_size;
        out.write(reinterpret_cast<const char*>(&section_size), sizeof(section_size));
        char tp = section.dataType();
        out.write(reinterpret_cast<char*>(&tp), sizeof(tp));
        section.write_binary_data(out);
    };
    auto save_mesh_data_bin = [status, &out](unsigned char IO_SEC, Section SEC, std::size_t size, const std::map<std::string, MetaTriMesh::MeshValue>& sm){
        if (!(status & IO_SEC) || sm.empty() || size == 0) return;
        char sec_str = SEC;
        out.write(reinterpret_cast<const char*>(&sec_str), sizeof(sec_str));
        unsigned long _size = size;
        out.write(reinterpret_cast<const char*>(&_size), sizeof(_size));
        for (auto& s: sm){
            auto& v = s.second;
            char ssec_str = NSECTIONS_SC;
            switch(v.m_type){
                case MetaTriMesh::MeshValue::CONSTANT: ssec_str = CONSTANTS_SSC; break;
                case MetaTriMesh::MeshValue::CONTIGUOUS: ssec_str = SCALARS_SSC; break;
                case MetaTriMesh::MeshValue::SPARSE: ssec_str = SPARSE_SCALARS_SSC; break;
                default:{
                    std::cerr << "In data \"" << s.first << "\" faced unknown data type with code " << static_cast<int>(v.m_type) << ". It will not be saved\n";
                    continue;
                }
            }
            v.m_type;
            out.write(reinterpret_cast<const char*>(&ssec_str), sizeof(ssec_str));
            write_string_bin(out, s.first);
            if (v.m_type == MetaTriMesh::MeshValue::SPARSE){
                char index_type = v.m_sparse_indexes.dataType();
                out.write(reinterpret_cast<const char*>(&index_type), sizeof(index_type));
                unsigned long v_sparse_size = v.m_sparse_indexes.m_size;
                out.write(reinterpret_cast<const char*>(&v_sparse_size), sizeof(v_sparse_size));
            }
            char value_type = v.m_values.dataType();
            out.write(reinterpret_cast<const char*>(&value_type), sizeof(value_type));
            unsigned long v_values_dim = v.m_values.m_dim;
            out.write(reinterpret_cast<const char*>(&v_values_dim), sizeof(v_values_dim));
            MetaTriMesh::Storage::write_binary_data(out, v.m_values.dataType(), v.m_type == MetaTriMesh::MeshValue::CONSTANT ? std::size_t(1) : v.m_values.size(), v.m_values.dim(), v.m_values.m_data);
            if (v.m_type == MetaTriMesh::MeshValue::SPARSE)
                v.m_sparse_indexes.write_binary_data(out);
        }
    };

    save_mesh_graph_section_bin(IO_POINTS, POINTS_SC, m.points);
    save_mesh_graph_section_bin(IO_EDGES, EDGES_SC, m.edges);
    save_mesh_graph_section_bin(IO_CELLS, CELLS_SC, m.cells);
    save_mesh_data_bin(IO_POINT_DATA, POINT_DATA_SC, m.points.m_size, m.point_data);
    save_mesh_data_bin(IO_CELL_DATA, CELL_DATA_SC, m.cells.m_size, m.cell_data);
    save_mesh_data_bin(IO_EDGE_DATA, EDGE_DATA_SC, m.edges.m_size, m.edge_data);
    save_mesh_data_bin(IO_HALFEDGE_DATA, HALFEDGE_DATA_SC, 2*m.edges.m_size, m.halfedge_data);
    for (auto& p: loaded_paths){
        char sec_str = LOAD_SC;
        out.write(reinterpret_cast<char*>(&sec_str), sizeof(sec_str));
        write_string_bin(out, p);
    }
    return true;
}

static inline TMH_File::Section convert(MetaTriMesh::MeshValue::Type tp){
    switch (tp){
        case MetaTriMesh::MeshValue::CONSTANT: return TMH::CONSTANTS_SSC;
        case MetaTriMesh::MeshValue::CONTIGUOUS: return TMH::SCALARS_SSC;
        case MetaTriMesh::MeshValue::SPARSE: return TMH::SPARSE_SCALARS_SSC;
        default:
            throw std::runtime_error("Faced unknown MeshValue::Type with code = " + std::to_string(int(tp)));
    }
}

bool TMH::write_txt(std::string file_name, const MetaTriMesh& m, unsigned char status, const std::vector<std::string>& loaded_paths){
    m_file_name = file_name; m_is_binary_file = false;
    std::ofstream out(file_name);
    out << std::setprecision(std::numeric_limits<double>::digits10 + 1) << std::fixed;
    out << "# tmh Mesh file version 1.0\n";
    out << m_title_text << "\n\n";
    
    auto save_mesh_graph_section_txt = [status, &out](unsigned char IO_SEC, Section SEC, auto& section){
        if (!(status & IO_SEC) || section.size() == 0) return;
        out << get_section_name(SEC) << " " << section.m_size << " " << MetaTriMesh::get_type_name(section.m_type) << "\n ";
        section.write_ascii_data(out, " ", "\n ");
        out << "\n\n";
    };
    auto save_mesh_data_txt = [status, &out](unsigned char IO_SEC, Section SEC, std::size_t size, const std::map<std::string, MetaTriMesh::MeshValue>& sm){
        if (!(status & IO_SEC) || sm.empty() || size == 0) return;
        out << get_section_name(SEC) << " " << size << "\n ";
        for (auto& s: sm){
            auto ssec = convert(s.second.m_type);
            out << get_subsection_name(ssec) << " \"" << s.first << "\" ";
            auto& v = s.second;
            if (v.m_type == MetaTriMesh::MeshValue::SPARSE)
                out << MetaTriMesh::get_type_name(v.m_sparse_indexes.dataType()) << " " << v.m_sparse_indexes.size() << " ";
            out << MetaTriMesh::get_type_name(v.m_values.dataType()) << " " << v.m_values.dim() << "\n  ";
            switch (ssec){
                case SCALARS_SSC: v.m_values.write_ascii_data(out, " ", "\n  "); break;
                case CONSTANTS_SSC: MetaTriMesh::Storage::write_ascii_data(out, v.m_values.dataType(), 1, v.m_values.dim(), v.m_values.m_data, " ", "\n  "); break;
                case SPARSE_SCALARS_SSC:{
                    assert(v.m_values.size() >= v.m_values.dim());
                    for (std::size_t k = 0; k < v.m_values.dim(); ++k){
                        MetaTriMesh::Value::write_ascii_var(out, v.m_values.at(0, k)); 
                        out << " ";
                    }
                    out << "\n ";
                    for (std::size_t i = 0; i < v.m_sparse_indexes.size(); ++i){
                        MetaTriMesh::Value::write_ascii_var(out, v.m_sparse_indexes.at(i, 0)); 
                        out << " ";
                        for (std::size_t k = 0; k < v.m_values.dim(); ++k){
                            MetaTriMesh::Value::write_ascii_var(out, v.m_values.at(i+1, k)); 
                            out << " ";
                        }
                        out << "\n ";
                    }
                }
            }
            out << "\n";
        }
        out << "\n";
    };
    save_mesh_graph_section_txt(IO_POINTS, POINTS_SC, m.points);
    save_mesh_graph_section_txt(IO_EDGES, EDGES_SC, m.edges);
    save_mesh_graph_section_txt(IO_CELLS, CELLS_SC, m.cells);
    save_mesh_data_txt(IO_POINT_DATA, POINT_DATA_SC, m.points.m_size, m.point_data);
    save_mesh_data_txt(IO_CELL_DATA, CELL_DATA_SC, m.cells.m_size, m.cell_data);
    save_mesh_data_txt(IO_EDGE_DATA, EDGE_DATA_SC, m.edges.m_size, m.edge_data);
    save_mesh_data_txt(IO_HALFEDGE_DATA, HALFEDGE_DATA_SC, 2*m.edges.m_size, m.halfedge_data);
    for (auto& p: loaded_paths)
        out << get_section_name(LOAD_SC) << " \"" << p << "\"\n\n";

    return true;
}