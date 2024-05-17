//
// Created by alex on 10.09.2023.
//

#ifndef AORTIC_VALVE_METATRIMESH_H
#define AORTIC_VALVE_METATRIMESH_H

#include <vector>
#include <array>
#include <string>
#include <algorithm>
#include <locale>
#include <map>
#include <set>
#include <type_traits>
#include <typeindex>
#include <iterator>
#include "../AVMesh.h"
#include "../AVMeshConnectivityHelpers.h"
#include <cstdint>

namespace World3d{

struct MetaTriMesh{
    enum DataType{
        BIT_TP = 0,
        CHAR_TP = 1,
        SHORT_TP = 2,
        INT_TP = 3,
        LONG_TP = 4,
        UCHAR_TP = 5,
        USHORT_TP = 6,
        UINT_TP = 7,
        ULONG_TP = 8,
        FLOAT_TP = 9,
        DOUBLE_TP = 10,
        NONE_TP,
    };
    enum ElementSparsity{
        NODE = 1,
        EDGE = 2,
        FACE = 3,
        HALFEDGE = 4,
    };
    enum class Endian{
        Big,
        Little
    };
    static inline Endian getEndianOrder();
    static inline DataType    get_type_by_name(const char* str);
    static inline std::size_t get_type_sizeof(DataType dt);
    static inline std::string get_type_name(DataType dt);
    static inline std::string get_sparsity_name(ElementSparsity e);
    static inline ElementSparsity get_sparsity_by_name(std::string s);
    static inline bool type_is_indexable(DataType dt);
    static inline bool type_is_floating_point(DataType dt);
    struct Value{
        union Var{
            bool b;
            char c;
            short s;
            int i;
            long l;
            unsigned char uc;
            unsigned short us;
            unsigned int ui;
            unsigned long ul;
            float f;
            double d;
            Var(): d() {}
        };

        DataType m_type = NONE_TP;
        Var m_val;

        Value& operator=(bool v){ m_type = BIT_TP; m_val.b = v; return *this; }
        Value& operator=(char v){ m_type = CHAR_TP; m_val.c = v; return *this; }
        Value& operator=(short v){ m_type = SHORT_TP; m_val.s = v; return *this; }
        Value& operator=(int v){ m_type = INT_TP; m_val.i = v; return *this; }
        Value& operator=(long v){ m_type = LONG_TP; m_val.l = v; return *this; }
        Value& operator=(unsigned char v){ m_type = UCHAR_TP; m_val.uc = v; return *this; }
        Value& operator=(unsigned short v){ m_type = USHORT_TP; m_val.us = v; return *this; }
        Value& operator=(unsigned int v){ m_type = UINT_TP; m_val.ui = v; return *this; }
        Value& operator=(unsigned long v){ m_type = ULONG_TP; m_val.ul = v; return *this; }
        Value& operator=(float v){ m_type = FLOAT_TP; m_val.f = v; return *this; }
        Value& operator=(double v){ m_type = DOUBLE_TP; m_val.d = v; return *this; }
        template<typename T>
        Value& set(T v){ return (*this = v); }
        inline Value& set(DataType tp, const void* val);
        template<typename T>
        typename std::enable_if<!std::is_pointer<T>::value, Value&>::type set(DataType tp, T val);
        const void* get(void* dummy = nullptr) const { return &m_val; }
        void* get(void* dummy = nullptr) { return &m_val; }
        DataType type() const { return m_type; }
        template<typename T> T& as(){ return *reinterpret_cast<T*>(get(nullptr)); }
        template<typename T> T as() const { return *reinterpret_cast<T*>(get(nullptr)); }
        template<typename T> T get() const;
        void clear() { m_type = NONE_TP; m_val.d = 0; }
        static Value read_ascii_var(std::istream& in, DataType type);
        static Value read_binary_var(std::istream& in, DataType type, Endian order = getEndianOrder());
        static void write_ascii_var(std::ostream& out, Value val);
        static void write_binary_var(std::ostream& out, Value val, Endian order = getEndianOrder());
    };
    
    template<unsigned N> struct TStorage;
    struct Storage{
        DataType m_type = NONE_TP;
        std::size_t m_size = 0;
        std::size_t m_dim = 0;
        std::vector<char> m_data;

        void clear() { m_type = NONE_TP, m_size = 0, m_dim = 0, m_data.clear(); }
        void resize(std::size_t size){ m_size = size; std::size_t sz = m_dim*m_size*get_type_sizeof(m_type); if (m_type == BIT_TP) sz = (sz/8) + ((sz%8 > 0) ? 1 : 0); m_data.resize(sz); }
        Value at(std::size_t cont_index) const { return at(cont_index/m_dim, cont_index%m_dim); }
        Value at(std::size_t i, std::size_t d) const { return get_value(i, d, m_type, m_size, m_dim, m_data.data(), 0); }
        template<typename T> void set(std::size_t i, std::size_t d, T val){ set_value(i, d, val,  m_type, m_size, m_dim, m_data.data(), 0); }
        void setV(std::size_t i, std::size_t d, Value val){ setV(i, d, val, m_size, m_dim, m_data.data(), 0); }
        template<typename iterator> void set_values(iterator beg, iterator end, std::size_t start_cont_index) { copy(beg, end, start_cont_index, m_type, m_data.data()); }
        template<typename iterator> void get_values(std::size_t start_cont_index, std::size_t end_cont_index, iterator dst) const { copy(start_cont_index, end_cont_index, m_type, dst, m_data.data()); }
        
        std::size_t size() const { return m_size; }
        std::size_t dim() const { return m_dim; }
        DataType dataType() const { return m_type; }
        std::size_t get_cont_index(std::size_t i, std::size_t d) const { return i * m_dim + d; }
        template<typename T> void get_element(std::size_t i, T& dst) const;
        template<typename T> void set_element(std::size_t i, const T& val);

        template<typename T> static void set_value(std::size_t i, std::size_t d, T val, DataType type, std::size_t size, std::size_t dim, char* data, std::size_t off = 0);
        static void setV(std::size_t i, std::size_t d, Value val, std::size_t size, std::size_t dim, char* data, std::size_t off = 0);
        inline static Value get_value(std::size_t i, std::size_t d, DataType type, std::size_t size, std::size_t dim, const char* data, std::size_t off = 0);
        template<typename iterator> static iterator copy(std::size_t start, std::size_t end, DataType type, iterator dst, const char* src_data);
        template<typename iterator> static void* copy(iterator beg, iterator end, std::size_t start_dst, DataType type, char* dst_data);

        Storage() = default;
        template<unsigned N>
        Storage(TStorage<N>&& a);
        template<unsigned N>
        Storage(const TStorage<N>& a);

        struct ElementDataInfo{
            DataType m_type;
            std::size_t m_dim;
        };
        struct ElementView{
            Storage& m_base;
            std::size_t elem_index;
            
            ElementView(Storage& from, std::size_t index): m_base{from}, elem_index{index} {}

            template<typename T> void get_element(T& dst) const { m_base.get_element(elem_index, dst); }
            template<typename T> void set_element(const T& src) const { m_base.set_element(elem_index, src); }
        };
        ElementView get_element_view(std::size_t i) { return ElementView(*this, i); }

        void read_binary_data(std::istream& in, Endian order = getEndianOrder()){ read_binary_data(in, dataType(), size(), dim(), m_data, order); }
        void write_binary_data(std::ostream& out, Endian order = getEndianOrder()) const { write_binary_data(out, dataType(), size(), dim(), m_data, order); }
        void read_ascii_data(std::istream& in, std::function<void(std::istream& in)> skip_after_element = nullptr) { read_ascii_data(in, dataType(), size(), dim(), m_data, skip_after_element); }
        void write_ascii_data(std::ostream& out, std::string values_delimiter = " ", std::string elements_delimiter = "\n") const { write_ascii_data(out, dataType(), size(), dim(), m_data, values_delimiter, elements_delimiter); }

        static void read_binary_data(std::istream& in, DataType type, std::size_t size, std::size_t dim, std::vector<char>& out, Endian order = getEndianOrder());
        static void read_ascii_data(std::istream& in, DataType type, std::size_t size, std::size_t dim, std::vector<char>& out, std::function<void(std::istream& in)> skip_after_element = nullptr);
        static void write_binary_data(std::ostream& out, DataType type, std::size_t size, std::size_t dim, const std::vector<char>& data, Endian order = getEndianOrder());
        static void write_ascii_data(std::ostream& out, DataType type, std::size_t size, std::size_t dim, const std::vector<char>& data, std::string values_delimiter = " ", std::string elements_delimiter = "\n");
    };
    template<unsigned N>
    struct TStorage{
        DataType m_type = NONE_TP;
        std::size_t m_size = 0;
        std::vector<char> m_data;
        using m_dim = std::integral_constant<unsigned, N>;

        void clear() { m_type = NONE_TP, m_size = 0, m_data.clear(); }
        void resize(std::size_t size){ m_size = size; std::size_t sz = N*m_size*get_type_sizeof(m_type); if (m_type == BIT_TP) sz = (sz/8) + ((sz%8 > 0) ? 1 : 0); m_data.resize(sz); }
        Value at(std::size_t i, std::size_t d) const { return Storage::get_value(i, d, m_type, m_size, dim(), m_data.data(), 0); }
        Value at(std::size_t cont_index) const { return at(cont_index/dim(), cont_index%dim()); }
        template<typename T> void set(std::size_t i, std::size_t d, T val){ Storage::template set_value(i, d, val,  m_type, m_size, dim(), m_data.data(), 0); }
        void setV(std::size_t i, std::size_t d, Value val){ Storage::setV(i, d, val, size(), dim(), m_data.data(), 0); }
        template<typename iterator> void set_values(iterator beg, iterator end, std::size_t start_cont_index) { Storage::copy(beg, end, start_cont_index, m_type, m_data.data()); }
        template<typename iterator> void get_values(std::size_t start_cont_index, std::size_t end_cont_index, iterator dst) const { Storage::copy(start_cont_index, end_cont_index, m_type, dst, m_data.data()); }
        
        std::size_t size() const { return m_size; }
        std::size_t dim() const { return N; }
        DataType dataType() const { return m_type; }
        std::size_t get_cont_index(std::size_t i, std::size_t d) const { return i * dim() + d; }
        template<typename T> void get_element(std::size_t i, T& dst) const;
        template<typename T> void set_element(std::size_t i, const T& val);

        TStorage() = default;
        TStorage(Storage&& a);
        TStorage(const Storage& a);
        TStorage& operator=(Storage&& a){ return *this = std::move(TStorage<N>(a)); }
        TStorage& operator=(const Storage& a){ return *this = std::move(TStorage<N>(a)); }

        void read_binary_data(std::istream& in, Endian order = getEndianOrder()){ Storage::read_binary_data(in, dataType(), size(), dim(), m_data, order); }
        void write_binary_data(std::ostream& out, Endian order = getEndianOrder()) const { Storage::write_binary_data(out, dataType(), size(), dim(), m_data, order); }
        void read_ascii_data(std::istream& in, std::function<void(std::istream& in)> skip_after_element = nullptr) { Storage::read_ascii_data(in, dataType(), size(), dim(), m_data, skip_after_element); }
        void write_ascii_data(std::ostream& out, std::string values_delimiter = " ", std::string elements_delimiter = "\n") const { Storage::write_ascii_data(out, dataType(), size(), dim(), m_data, values_delimiter, elements_delimiter); }
    };
    
    struct MeshValue{
        enum Type{
            CONSTANT = 0,
            CONTIGUOUS = 1,
            SPARSE = 2,
            UNSET
        };
        Type m_type = UNSET;
        TStorage<1> m_sparse_indexes;
        Storage m_values;
        static const char* getTypeName(Type type) { static const char* stp[] = {"CONSTANT", "CONTIGUOUS", "SPARSE", "UNSET"}; return stp[type]; }

        bool isValid() const { return m_type != UNSET; }
        void clear() { m_type = UNSET, m_sparse_indexes.clear(), m_values.clear(); }
        std::size_t get_cont_index(std::size_t i, std::size_t d) const { return m_values.get_cont_index(i, d); }
        template<typename iterator> void get_values(std::size_t start_cont_index, std::size_t end_cont_index, iterator dst) const;
        template<typename T> void get_element(std::size_t i, T& val) const;
        static Storage create_contiguous(const MeshValue& v, std::size_t expect_size);
    };
    struct NamedValueView{
        const std::string& m_name;
        MeshValue& m_value;

        NamedValueView(const std::string& name, MeshValue& value): m_value{value}, m_name{name} {}

        const std::string& name() const { return m_name; }
        const MeshValue& value() const { return m_value; }
        MeshValue& value() { return m_value; }
    };

    template<typename T>
    struct DefaultConverter{
        Storage::ElementDataInfo dinfo;
        DefaultConverter() = default;
        DefaultConverter(DataType type, std::size_t dim): dinfo{type, dim} {}
        DefaultConverter(Storage::ElementDataInfo dinfo): dinfo{dinfo} {}

        Storage operator()(const Mesh& m, std::string tag_name, ElementSparsity sp);
    };
    
    using TagDataConverter = std::function<Storage(const Mesh& m, std::string tag_name, ElementSparsity sp)>;
    
    static std::map<std::type_index, MetaTriMesh::TagDataConverter>& get_converter();
    static void _init_tag_converters(bool enforce = false);
    static void register_tag_converter(std::type_index tindex, TagDataConverter conv);

    TStorage<3> points;
    TStorage<2> edges; 
    TStorage<3> cells; 
    std::map<std::string, MeshValue> point_data;
    std::map<std::string, MeshValue> halfedge_data;
    std::map<std::string, MeshValue> edge_data;
    std::map<std::string, MeshValue> cell_data;

    MetaTriMesh() { _init_tag_converters(); }
    MetaTriMesh(const Mesh& m, std::string point_tag = "v:point") { _init_tag_converters(); initFromMesh(m, point_tag); }
    void clear() { points.clear(); edges.clear(); cells.clear(); point_data.clear(); halfedge_data.clear(); edge_data.clear(); cell_data.clear(); }
    bool empty_mesh() const { return points.size() == 0 && edges.size() == 0 && cells.size() == 0; }
    bool empty_data() const { return point_data.size() == 0 && halfedge_data.size() == 0 && edge_data.size() == 0 && cell_data.size() == 0; }
    bool empty() const { return empty_mesh() && empty_data(); }
    std::ostream& print_content_info(std::ostream& out = std::cout, const std::string& prefix = "") const;

    void setup_edges();

    template<typename Geom_index_type, typename T>
    bool createTagOnMesh(Mesh* m, std::string meta_tag_name, std::string mesh_tag_name = "", bool overwrite = true) const;
    template<typename T = Point>
    bool createTagOnMeshFromSelfPoints(Mesh* m, std::string mesh_tag_name, bool overwrite = true) const;
    Mesh convertToMesh() const;

    template<typename Geom_index_type>
    static MeshValue constructTagDataFromMesh(const Mesh& m, std::string name);
    MetaTriMesh& initFromMesh(const Mesh& m, std::string point_tag = "v:point");
    template<typename Geom_index_type>
    MetaTriMesh& readTagDataFromMesh(const Mesh& m, const std::string& mesh_tag_name, const std::string& meta_tag_name);
    template<typename Geom_index_type>
    MetaTriMesh& readTagsDataFromMesh(const Mesh& m, const std::set<std::string>& exclude_names = std::set<std::string>());
    MetaTriMesh& readAllTagsFromMesh(const Mesh& m, const std::set<std::string>& exclude_names = std::set<std::string>());

    template<typename Geom_index_type>
    std::set<std::string> value_names() const; 
    inline std::map<std::string, MeshValue>& data_array(ElementSparsity e);
    inline const std::map<std::string, MeshValue>& data_array(ElementSparsity e) const;

private:
    TStorage<2> compute_edges() const;
};


}


#include "MetaTriMesh.inl"

#endif //METATRIMESH

