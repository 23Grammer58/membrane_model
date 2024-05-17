//
// Created by alex on 25.09.2023.
//

#ifndef AORTIC_VALVE_VTKFORMAT_H
#define AORTIC_VALVE_VTKFORMAT_H

#include "MetaTriMesh.h"
#include <fstream>

namespace World3d{

struct VTK_File{
    enum Section{
        DATASET_SC,
        POINT_DATA_SC,
        CELL_DATA_SC,
        NSECTION_SC,

        SCALARS_SSC,
        COLOR_SCALARS_SSC,
        LOOKUP_TABLE_SSC,
        VECTORS_SSC,
        NORMALS_SSC,
        TEXTURE_COORDINATES_SSC,
        TENSORS_SSC,
        FIELD_SSC
    };
    static Section get_section(const char* name);
    static Section get_subsection(const char* name);
    void read(const std::string& filename, MetaTriMesh& m);
    void write(const std::string& filename, const MetaTriMesh& m, bool is_binary = false);

    std::string get_last_file_name() const { return m_file_name; }
    std::string get_last_title_text() const { return m_title_text; }
    std::pair<int, int> get_last_vtk_version() const { return {m_major_version, m_minor_version}; }
    VTK_File&   set_title_text(const std::string& str) { return m_title_text = str, *this; }
protected:
    std::string m_file_name;
    std::string m_title_text = "File written by membrane-model";
    bool m_is_binary = false; 
    int m_major_version = 3, m_minor_version = 0;
    static const int s_major_version = 3;
    static const int s_minor_version = 0;  
};

} 

#endif //AORTIC_VALVE_VTKFORMAT_H