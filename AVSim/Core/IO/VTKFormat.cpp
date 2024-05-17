#include "VTKFormat.h"

using namespace World3d;

VTK_File::Section VTK_File::get_section(const char* name){
    if (!strcmp(name, "DATASET") ) return DATASET_SC;
    else if (!strcmp(name, "POINT_DATA")) return POINT_DATA_SC;
    else if (!strcmp(name, "CELL_DATA")) return CELL_DATA_SC;
    else throw std::runtime_error("Faced unknown section = \"" + std::string(name) + "\"");
}
VTK_File::Section VTK_File::get_subsection(const char* name){
    static std::map<std::string, Section> ssmap{
        {"SCALARS", SCALARS_SSC},
        {"COLOR_SCALARS", COLOR_SCALARS_SSC},
        {"LOOKUP_TABLE", LOOKUP_TABLE_SSC},
        {"VECTORS", VECTORS_SSC},
        {"NORMALS", NORMALS_SSC},
        {"TEXTURE_COORDINATES", TEXTURE_COORDINATES_SSC},
        {"TENSORS", TENSORS_SSC},
        {"FIELD", FIELD_SSC}
    };
    auto it = ssmap.find(name);
    return it != ssmap.end() ? it->second : NSECTION_SC;
}

static const int VTK_TRIANGLE  = 5;

void VTK_File::read(const std::string& filename, MetaTriMesh& m){
    std::ifstream in(filename);
    if (!in.is_open())
        throw std::runtime_error("Can't open file \"" + filename + "\"");
    m_file_name = filename;
    std::string format_str;
    std::getline(in, format_str);
    {
        unsigned txt_start = 0;
        while(txt_start < format_str.size() && (std::isspace(format_str[txt_start]) || format_str[txt_start] == '#') ) ++txt_start;
        std::stringstream ss(format_str.substr(txt_start));
        std::string text_prt; text_prt.resize(21);
        ss.read(text_prt.data(), 21*sizeof(char));
        if (text_prt != "vtk DataFile Version ")
            throw std::runtime_error("\"" + filename + "\" is not \"vtk\" file but considered as vtk" );
        char pnt = '.';
        ss >> m_major_version >> pnt >> m_minor_version;  
    }
    std::getline(in, m_title_text);
    std::string storage_type;
    in >> storage_type;
    if (storage_type == "ASCII") m_is_binary = false;
    else if (storage_type == "BINARY") m_is_binary = true;
    else throw std::runtime_error("Expected ASCII|BINARY but found \"" + storage_type + "\"");
    in >> std::ws;
    std::string line;
    auto read_storage_data = [is_bin = m_is_binary](std::istream& in, auto& s) {
        if (is_bin) s.read_binary_data(in, MetaTriMesh::Endian::Big);
        else s.read_ascii_data(in);
    };
    while(!in.eof()){
        std::getline(in, line); in >> std::ws;
        std::stringstream ss(line);
        std::string section, section_type;
        ss >> section >> section_type;
        auto sec = get_section(section.c_str());
        switch (sec){
            case DATASET_SC:{
                if (section_type != "UNSTRUCTURED_GRID") 
                    throw std::runtime_error("Currently supported reading only vtk UNSTRUCTURED_GRID");
                for (int _ss = 0; _ss < 3; ++_ss){
                    std::getline(in, line); in >> std::ws;
                    std::stringstream oss(line);
                    std::string part; int size = 0;
                    oss >> part >> size;
                    if (part == "POINTS"){
                        std::string dtype_str;
                        oss >> dtype_str;
                        m.points.m_type = MetaTriMesh::get_type_by_name(dtype_str.c_str());
                        m.points.resize(size);
                        read_storage_data(in, m.points);
                    } else if (part == "CELLS"){
                        int msize = 0;
                        oss >> msize;
                        if (msize != 4*size) 
                            throw std::runtime_error("Supported only vtk files for exceptionaly meshes with triangular cells");
                        MetaTriMesh::TStorage<4> tcells;
                        tcells.m_type = MetaTriMesh::INT_TP;
                        tcells.resize(size);
                        read_storage_data(in, tcells); 
                        m.cells.m_type = MetaTriMesh::INT_TP;
                        m.cells.resize(size);
                        for (int j = 0; j < size; ++j){
                            auto dst = reinterpret_cast<int*>(m.cells.m_data.data());
                            auto src = reinterpret_cast<const int*>(tcells.m_data.data());
                            for (int k = 0; k < 3; ++k) dst[3*j + k] = src[4*j + k + 1];
                        }
                    } else if (part == "CELL_TYPES"){
                        MetaTriMesh::TStorage<1> tcells;
                        tcells.m_type = MetaTriMesh::INT_TP;
                        tcells.resize(size);
                        read_storage_data(in, tcells); 
                        for (int j = 0; j < size; ++j){
                            auto tp = reinterpret_cast<const int*>(tcells.m_data.data())[j];
                            if (tp != VTK_TRIANGLE)
                                throw std::runtime_error("Supported only triangles but faced cell with type code = " + std::to_string(tp));
                        }

                    } else 
                        throw std::runtime_error("Can't parse \"" + part + "\", expected POINTS|CELLS|CELL_TYPES");
                    in >> std::ws;    
                }
                break;
            }
            case POINT_DATA_SC:
            case CELL_DATA_SC:{
                auto* p = (sec == POINT_DATA_SC ? &m.point_data : &m.cell_data);
                int size = std::stoi(section_type);
                bool continue_read = true;
                while(continue_read && !in.eof()){
                    auto pin = in.tellg();
                    MetaTriMesh::MeshValue vv;
                    vv.m_type = MetaTriMesh::MeshValue::CONTIGUOUS;
                    auto& s = vv.m_values;
                    s.m_size = size;
                    std::getline(in, line);
                    std::stringstream oss(line);
                    std::string subsec, name;
                    oss >> subsec >> name;
                    auto ssec = get_subsection(subsec.c_str());
                    bool save_res = true;
                    switch (ssec){
                        case SCALARS_SSC:{
                            std::string dtype_str; 
                            oss >> dtype_str; oss >> std::ws;
                            int nComp = 1;
                            if (!oss.eof()) oss >> nComp;
                            s.m_dim = nComp;
                            s.m_type = MetaTriMesh::get_type_by_name(dtype_str.c_str());
                            auto ppin = in.tellg();
                            std::string check_line = "LOOKUP_TABLE";
                            in.read(check_line.data(), check_line.size());
                            if (check_line == "LOOKUP_TABLE") std::getline(in, line);
                            else in.seekg(ppin);
                            read_storage_data(in, s); 
                            break;
                        }
                        case COLOR_SCALARS_SSC:{
                            int nComp = 1;
                            oss >> nComp;
                            s.m_dim = nComp;
                            s.m_type = m_is_binary ? MetaTriMesh::UCHAR_TP : MetaTriMesh::FLOAT_TP;
                            read_storage_data(in, s);
                            break;
                        }
                        case LOOKUP_TABLE_SSC:{
                            save_res = false;
                            int nComp = 1;
                            oss >> nComp;
                            s.m_size = nComp;
                            s.m_dim = 4;
                            s.m_type = m_is_binary ? MetaTriMesh::UCHAR_TP : MetaTriMesh::FLOAT_TP;
                            read_storage_data(in, s);
                            std::cerr << "WARRNING: LOOKUP_TABLE dataset type is not supported and will be ignored\n";
                            break;
                        }
                        case VECTORS_SSC:
                        case NORMALS_SSC:{
                            std::string dtype_str; 
                            oss >> dtype_str;
                            s.m_type = MetaTriMesh::get_type_by_name(dtype_str.c_str());
                            s.m_dim = 3;
                            read_storage_data(in, s);
                            break;
                        }
                        case TEXTURE_COORDINATES_SSC:{
                            int nComp = 1;
                            oss >> nComp;
                            std::string dtype_str; 
                            oss >> dtype_str;
                            s.m_type = MetaTriMesh::get_type_by_name(dtype_str.c_str());
                            s.m_dim = nComp;
                            read_storage_data(in, s);
                            break;
                        }
                        case TENSORS_SSC:{
                            std::string dtype_str; 
                            oss >> dtype_str;
                            s.m_type = MetaTriMesh::get_type_by_name(dtype_str.c_str());
                            s.m_dim = 9;
                            read_storage_data(in, s);
                            break;
                        }
                        case FIELD_SSC:{
                            save_res = false;
                            std::cerr << "WARRNING: FIELD dataset type is not supported and will be ignored\n";
                            int narr = 0;
                            oss >> narr;
                            for (int l = 0; l < narr; ++l){
                                std::getline(in, line);
                                std::stringstream ox(line);
                                std::string arr_name, dtype_str;
                                int nComp = 0, nTupl = 0;
                                ox >> arr_name >> nComp >> nTupl >> dtype_str;
                                MetaTriMesh::Storage ss;
                                ss.m_dim = nComp, ss.m_size = nTupl;
                                ss.m_type = MetaTriMesh::get_type_by_name(dtype_str.c_str());
                                read_storage_data(in, ss);
                            }
                            break;
                        }
                        
                        default:{
                            continue_read = false;
                            in.seekg(pin);
                            break;
                        }
                    }
                    if (save_res && continue_read)
                        p->operator[](name) = std::move(vv);   
                    in >> std::ws;    
                }
                break;
            }
        }
        in >> std::ws;
    }
}

void VTK_File::write(const std::string& filename, const MetaTriMesh& m, bool is_binary){
    m_is_binary = is_binary;
    std::ofstream f(filename);
    if (!f.is_open())
        throw std::runtime_error("Can't open file \"" + filename + "\"");
    m_file_name = filename; 
    f << std::setprecision(std::numeric_limits<double>::digits10 + 1) << std::fixed;    
    f << "# vtk DataFile Version " << s_major_version << "." << s_minor_version << "\n";
    std::string title_text = m_title_text;
    if (title_text.size() > 255) title_text = m_title_text.substr(0, 254);
    f << title_text << "\n";
    f << (is_binary ? "BINARY\n" : "ASCII\n");
    auto write_data_set = [is_binary](std::ofstream& f, auto& s){
        if (is_binary) s.write_binary_data(f, MetaTriMesh::Endian::Big);
        else s.write_ascii_data(f);
    };
    if (m.points.size() > 0 && m.cells.size() > 0){
        f << "DATASET UNSTRUCTURED_GRID\n";
        f << "POINTS " << m.points.size() << " " << MetaTriMesh::get_type_name(m.points.m_type) << "\n";
        write_data_set(f, m.points); f << "\n";
        f << "CELLS " << m.cells.size() << " " << (4*m.cells.size()) << "\n";
        {
            MetaTriMesh::TStorage<4> tcell;
            tcell.m_type = MetaTriMesh::INT_TP;
            tcell.resize(m.cells.size());
            for (int i = 0; i < m.cells.size(); ++i){
                std::array<int, 3> elem;
                m.cells.get_element(i, elem);
                auto* dst = reinterpret_cast<int*>(tcell.m_data.data()) + 4*i;
                dst[0] = 3, dst[1] = elem[0], dst[2] = elem[1], dst[3] = elem[2];
            }
            write_data_set(f, tcell); f << "\n";
        }
        f << "CELL_TYPES " << m.cells.size() << "\n";
        {
            MetaTriMesh::TStorage<1> tcell;
            tcell.m_type = MetaTriMesh::INT_TP;
            tcell.resize(m.cells.size());
            auto* dst = reinterpret_cast<int*>(tcell.m_data.data());
            std::fill(dst, dst + m.cells.size(), VTK_TRIANGLE);
            write_data_set(f, tcell); f << "\n";
        }
    }
    std::size_t npoint = m.points.size();
    if (npoint == 0){
        for (auto& p: m.point_data)
            if (p.second.m_type == MetaTriMesh::MeshValue::CONTIGUOUS){
                npoint = p.second.m_values.m_size;
                break;
            }
    }
    if (!m.point_data.empty() && npoint > 0){
        f << "POINT_DATA " << npoint << "\n";
        for (auto& p: m.point_data){
            if (!p.second.isValid() || p.second.m_values.m_dim == 0) continue;
            auto& v = p.second;
            f << "SCALARS " << p.first << " " << MetaTriMesh::get_type_name(v.m_values.m_type) << " " << p.second.m_values.m_dim << "\n";
            f << "LOOKUP_TABLE default\n";
            if (v.m_type == MetaTriMesh::MeshValue::CONTIGUOUS)
                write_data_set(f, v.m_values);
            else {
                MetaTriMesh::Storage s = MetaTriMesh::MeshValue::create_contiguous(v, npoint);
                write_data_set(f, s);
            }     
            f << "\n";
        }
    }
    std::size_t ncell = m.cells.size();
    if (npoint == 0){
        for (auto& p: m.cell_data)
            if (p.second.m_type == MetaTriMesh::MeshValue::CONTIGUOUS){
                npoint = p.second.m_values.m_size;
                break;
            }
    }
    if (!m.cell_data.empty() && ncell > 0){
        f << "CELL_DATA " << ncell << "\n";
        for (auto& p: m.cell_data){
            if (!p.second.isValid() || p.second.m_values.m_dim == 0) continue;
            auto& v = p.second;
            f << "SCALARS " << p.first << " " << MetaTriMesh::get_type_name(v.m_values.m_type) << " " << p.second.m_values.m_dim << "\n";
            f << "LOOKUP_TABLE default\n";
            if (v.m_type == MetaTriMesh::MeshValue::CONTIGUOUS)
                write_data_set(f, v.m_values);
            else {
                MetaTriMesh::Storage s = MetaTriMesh::MeshValue::create_contiguous(v, ncell);
                write_data_set(f, s);
            }     
            f << "\n";
        }
    }
}