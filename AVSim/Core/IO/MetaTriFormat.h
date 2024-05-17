//
// Created by alex on 10.09.2023.
//

#ifndef AORTIC_VALVE_METATRIFORMAT_H
#define AORTIC_VALVE_METATRIFORMAT_H

#include "MetaTriMesh.h"
#include <fstream>

namespace World3d{

struct TMH_File{
    enum Section{
        POINTS_SC = 0,
        EDGES_SC,
        CELLS_SC,
        POINT_DATA_SC,
        EDGE_DATA_SC,
        CELL_DATA_SC,
        HALFEDGE_DATA_SC,
        LOAD_SC,

        NSECTIONS_SC,
        SCALARS_SSC = 9,
        CONSTANTS_SSC,
        SPARSE_SCALARS_SSC,
    };
    static Section get_data_subsection_type(const char* s);
    static const char* get_subsection_name(Section ssec);
    static const unsigned char IO_POINTS = 0x1;         ///< read/write mesh points
    static const unsigned char IO_EDGES = 0x2;          ///< read/write mesh edge connectivity
    static const unsigned char IO_CELLS = 0x4;          ///< read/write mesh cell connectivity
    static const unsigned char IO_POINT_DATA = 0x8;     ///< read/write mesh point data
    static const unsigned char IO_EDGE_DATA = 0x10;     ///< read/write mesh edge data
    static const unsigned char IO_CELL_DATA = 0x20;     ///< read/write mesh cell data
    static const unsigned char IO_HALFEDGE_DATA = 0x40; ///< read/write mesh halfedge data
    static const unsigned char IO_LOAD_FILES = 0x80;    ///< only in read mode, allow open load modules
    static const unsigned char IO_DEFAULT = 0xFF;
    
    std::string m_title_text = "";
    
    static const char* file_extension() { return "tmh"; }
    static bool is_binary_file(std::string file_name);

    bool read(std::string file_name, MetaTriMesh& m, unsigned char status = IO_DEFAULT);
    bool write_bin(std::string file_name, const MetaTriMesh& m, unsigned char status = IO_DEFAULT, const std::vector<std::string>& loaded_paths = std::vector<std::string>());
    bool write_txt(std::string file_name, const MetaTriMesh& m, unsigned char status = IO_DEFAULT, const std::vector<std::string>& loaded_paths = std::vector<std::string>());
    
    std::string get_last_file_name() const { return m_file_name; }
    std::string get_last_title_text() const { return m_title_text; }
    TMH_File&   set_title_text(const std::string& str) { return m_title_text = str, *this; }
protected:
    static Section get_section(const char* str);
    static const char* get_section_name(Section sec);
    bool read_txt(std::string file_name, MetaTriMesh& m, unsigned char status, unsigned char faced_sec = 0x0);
    bool read_bin(std::string file_name, MetaTriMesh& m, unsigned char status, unsigned char faced_sec = 0x0);

    bool is_loaded_state = false;
    std::string m_file_name;
    bool m_is_binary_file = false;
};

}


#endif //AORTIC_VALVE_METATRIFORMAT_H