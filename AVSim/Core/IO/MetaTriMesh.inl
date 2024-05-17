//
// Created by alex on 29.01.2021.
//

#ifndef AORTIC_VALVE_METATRIMESH_INL
#define AORTIC_VALVE_METATRIMESH_INL

#include "MetaTriMesh.h"

namespace World3d{

inline bool MetaTriMesh::type_is_indexable(DataType dt){ return static_cast<unsigned>(dt) < static_cast<unsigned>(FLOAT_TP) && dt != BIT_TP; }
inline bool MetaTriMesh::type_is_floating_point(DataType dt) { return dt == FLOAT_TP || dt == DOUBLE_TP; }
inline std::size_t MetaTriMesh::get_type_sizeof(DataType dt){
    static size_t _data_type_names[]{sizeof(char), sizeof(char), sizeof(short), sizeof(int), sizeof(long), sizeof(unsigned char), sizeof(unsigned short), sizeof(unsigned int), sizeof(unsigned long), sizeof(float), sizeof(double), 0};
    return _data_type_names[static_cast<unsigned>(dt)];
}
inline std::string MetaTriMesh::get_type_name(DataType dt){
    static const char* _data_type_names[]{"bit", "char", "short", "int", "long", "uchar", "ushort", "uint", "ulong", "float", "double", "none"};
    return _data_type_names[static_cast<unsigned>(dt)];
}
inline std::string MetaTriMesh::get_sparsity_name(ElementSparsity e){
    static const char* res[]{"NODE", "EDGE", "FACE", "HALFEDGE"};
    int id = 0;
    switch (e){
        case NODE: id = 0; break;
        case EDGE: id = 1; break;
        case FACE: id = 2; break;
        case HALFEDGE: id = 3; break;
    } 
    return res[id];
}
inline MetaTriMesh::ElementSparsity MetaTriMesh::get_sparsity_by_name(std::string s){
    if (s == "NODE") return NODE;
    if (s == "EDGE") return EDGE;
    if (s == "FACE") return FACE;
    if (s == "HALFEDGE") return HALFEDGE;
    throw std::runtime_error("Faced unknown sparsity name = \"" + s + "\"");
}
inline MetaTriMesh::DataType MetaTriMesh::get_type_by_name(const char* str){
    std::string s = str;
    std::transform(s.begin(), s.end(), s.begin(), [](char c) {return std::tolower(c, std::locale()); });
    static std::map<std::string, MetaTriMesh::DataType> type_map{
        {"double", DOUBLE_TP}, {"int", INT_TP}, {"char", CHAR_TP}, {"long", LONG_TP}, {"float", FLOAT_TP}, {"bit", BIT_TP},
        {"short", SHORT_TP}, {"uchar", UCHAR_TP}, {"ushort", USHORT_TP}, {"uint", UINT_TP}, {"ulong", ULONG_TP}, {"none", NONE_TP}
    };
    auto it = type_map.find(s);
    if (it == type_map.end())
        throw std::runtime_error("Faced unknown data type \"" + std::string(str) + "\"");
    return it->second;
}
inline MetaTriMesh::Endian MetaTriMesh::getEndianOrder(){
    const union endian_tester {
        uint32_t   n;
        uint8_t    p[4];
    } test = {0x01020304};
    return test.p[0] == 0x04 ? Endian::Little : Endian::Big;
}

inline MetaTriMesh::Value& MetaTriMesh::Value::set(DataType tp, const void* val) {
    m_type = tp;
    switch (tp){
        case BIT_TP: set(*reinterpret_cast<const bool*>(val)); break;
        case CHAR_TP: set(*reinterpret_cast<const char*>(val)); break;
        case SHORT_TP: set(*reinterpret_cast<const short*>(val)); break;
        case INT_TP: set(*reinterpret_cast<const int*>(val)); break;
        case LONG_TP: set(*reinterpret_cast<const long*>(val)); break;
        case UCHAR_TP: set(*reinterpret_cast<const unsigned char*>(val)); break;
        case USHORT_TP: set(*reinterpret_cast<const unsigned short*>(val)); break;
        case UINT_TP: set(*reinterpret_cast<const unsigned int*>(val)); break;
        case ULONG_TP: set(*reinterpret_cast<const unsigned long*>(val)); break;
        case FLOAT_TP: set(*reinterpret_cast<const float*>(val)); break;
        case DOUBLE_TP: set(*reinterpret_cast<const double*>(val)); break;
        default:{}
    }
    return *this;
}

template<typename T> T MetaTriMesh::Value::get() const { 
    switch (m_type){
        case BIT_TP: return static_cast<T>(m_val.b);
        case CHAR_TP: return static_cast<T>(m_val.c);
        case SHORT_TP: return static_cast<T>(m_val.s);
        case INT_TP: return static_cast<T>(m_val.i);
        case LONG_TP: return static_cast<T>(m_val.l);
        case UCHAR_TP: return static_cast<T>(m_val.uc);
        case USHORT_TP: return static_cast<T>(m_val.us);
        case UINT_TP: return static_cast<T>(m_val.ui);
        case ULONG_TP: return static_cast<T>(m_val.ul);
        case FLOAT_TP: return static_cast<T>(m_val.f);
        case DOUBLE_TP: return static_cast<T>(m_val.d);
    }
    assert(0);
    return T(); 
}
template<typename T>
typename std::enable_if<!std::is_pointer<T>::value, MetaTriMesh::Value&>::type MetaTriMesh::Value::set(DataType tp, T val){
    switch (m_type){
        case BIT_TP: return set(static_cast<bool>(val));
        case CHAR_TP: return set(static_cast<char>(val));
        case SHORT_TP: return set(static_cast<short>(val));
        case INT_TP: return set(static_cast<int>(val));
        case LONG_TP: return set(static_cast<long>(val));
        case UCHAR_TP: return set(static_cast<unsigned char>(val));
        case USHORT_TP: return set(static_cast<unsigned short>(val));
        case UINT_TP: return set(static_cast<unsigned int>(val));
        case ULONG_TP: return set(static_cast<unsigned long>(val));
        case FLOAT_TP: return set(static_cast<float>(val));
        case DOUBLE_TP: return set(static_cast<double>(val));
    }
    assert(0);
    return *this;
}

template<typename T> 
void MetaTriMesh::Storage::set_value(std::size_t i, std::size_t d, T val, DataType type, std::size_t size, std::size_t dim, char* data, std::size_t off){
    auto dt_sz = get_type_sizeof(type);
    off = off + (i*dim + d) * dt_sz;
    switch (type){
        case CHAR_TP: *reinterpret_cast<char*>(data+off) = static_cast<char>(val); break;
        case SHORT_TP: *reinterpret_cast<short*>(data+off) = static_cast<short>(val); break;
        case INT_TP: *reinterpret_cast<int*>(data+off) = static_cast<int>(val); break;
        case LONG_TP: *reinterpret_cast<long*>(data+off) = static_cast<long>(val); break;
        case UCHAR_TP: *reinterpret_cast<unsigned char*>(data+off) = static_cast<unsigned char>(val); break;
        case USHORT_TP: *reinterpret_cast<unsigned short*>(data+off) = static_cast<unsigned short>(val); break;
        case UINT_TP: *reinterpret_cast<unsigned int*>(data+off) = static_cast<unsigned int>(val); break;
        case ULONG_TP: *reinterpret_cast<unsigned long*>(data+off) = static_cast<unsigned long>(val); break;
        case FLOAT_TP: *reinterpret_cast<float*>(data+off) = static_cast<float>(val); break;
        case DOUBLE_TP: *reinterpret_cast<double*>(data+off) = static_cast<double>(val); break;
        case BIT_TP: {
            if (val)
                data[off/8] |= (1 >> (off%8));
            else
                data[off/8] &= ~((char)(1) >> (off%8)); 
        }
        default:  assert(0);
    }
}  

namespace internals{
    template<typename iterator, typename Dummy = void>
    struct ReadBool{
        static iterator apply(std::size_t start, std::size_t end, iterator dst, const char* src_data){
            throw std::runtime_error(std::string("Can't convert types from bool to ") + typeid(decltype(*std::declval<iterator>())).name());
            return dst;
        }
    }; 
    template<typename iterator>
    struct ReadBool<iterator, std::void_t<decltype(std::declval<iterator>()->operator*() = std::declval<bool>())> >{
        static iterator apply(std::size_t start, std::size_t end, iterator dst, const char* src_data){
            for (std::size_t i = start; i < end; ++i)
                *(dst++) = ((src_data[i/8] & (1 >> (i%8))) > 0);
            return dst;
        }
    };
    template<typename iterator, typename Dummy = void>
    struct WriteBool{
        static char* apply(iterator beg, iterator end, std::size_t start_dst, char* dst_data){
            throw std::runtime_error(std::string("Can't convert type \"") + typeid(decltype(*std::declval<iterator>())).name() + "\" to \"bool\"");
            return dst_data;
        }
    };
    template<typename iterator>
    struct WriteBool<iterator, std::void_t<decltype( bool(*std::declval<iterator>()) )> >{
        static char* apply(iterator beg, iterator end, std::size_t start_dst, char* dst_data){
            std::size_t i = start_dst;
            for (auto it = beg; it != end; ++it){
                if (*it)
                    dst_data[i/8] |= (1 << (i%8));
                else
                    dst_data[i/8] &= ~(char(1) << (i%8));    
                i++;
            }
            return dst_data + (i/8);
        }
    };
    template<typename iterator>
    struct WriteBool<iterator, std::void_t<decltype( bool(*(*std::declval<iterator>())) )> >{
        static char* apply(iterator beg, iterator end, std::size_t start_dst, char* dst_data){
            std::size_t i = start_dst;
            for (auto it = beg; it != end; ++it){
                if (*(*it))
                    dst_data[i/8] |= (1 << (i%8));
                else
                    dst_data[i/8] &= ~(char(1) << (i%8));    
                i++;
            }
            return dst_data + (i/8);
        }
    };
}

template<typename iterator> iterator MetaTriMesh::Storage::copy(std::size_t start, std::size_t end, DataType type, iterator dst, const char* src_data){
    auto dt_sz = get_type_sizeof(type);
    auto src_beg = src_data + start*dt_sz;
    auto src_end = src_data + end*dt_sz;

    switch (type){
        case CHAR_TP: return std::copy(reinterpret_cast<const char*>(src_beg), reinterpret_cast<const char*>(src_end), dst);
        case SHORT_TP: return std::copy(reinterpret_cast<const short*>(src_beg), reinterpret_cast<const short*>(src_end), dst);
        case INT_TP: return std::copy(reinterpret_cast<const int*>(src_beg), reinterpret_cast<const int*>(src_end), dst);
        case LONG_TP: return std::copy(reinterpret_cast<const long*>(src_beg), reinterpret_cast<const long*>(src_end), dst);
        case UCHAR_TP: return std::copy(reinterpret_cast<const unsigned char*>(src_beg), reinterpret_cast<const unsigned char*>(src_end), dst);
        case USHORT_TP: return std::copy(reinterpret_cast<const unsigned short*>(src_beg), reinterpret_cast<const unsigned short*>(src_end), dst);
        case UINT_TP: return std::copy(reinterpret_cast<const unsigned int*>(src_beg), reinterpret_cast<const unsigned int*>(src_end), dst);
        case ULONG_TP: return std::copy(reinterpret_cast<const unsigned long*>(src_beg), reinterpret_cast<const unsigned long*>(src_end), dst);
        case FLOAT_TP: return std::copy(reinterpret_cast<const float*>(src_beg), reinterpret_cast<const float*>(src_end), dst);
        case DOUBLE_TP: return std::copy(reinterpret_cast<const double*>(src_beg), reinterpret_cast<const double*>(src_end), dst);
        case BIT_TP: {
            return internals::ReadBool<iterator>::apply(start, end, dst, src_data);    
        }
    }
    assert(0);
    return dst;
}
template<typename iterator> void* MetaTriMesh::Storage::copy(iterator beg, iterator end, std::size_t start_dst, DataType type, char* dst_data){
    auto dt_sz = get_type_sizeof(type);
    auto dst_beg = dst_data + start_dst*dt_sz;
    
    switch (type){
        case CHAR_TP: return reinterpret_cast<void*>(std::copy(beg, end, reinterpret_cast<char*>(dst_beg)));
        case SHORT_TP: return reinterpret_cast<void*>(std::copy(beg, end, reinterpret_cast<short*>(dst_beg)));
        case INT_TP: return reinterpret_cast<void*>(std::copy(beg, end, reinterpret_cast<int*>(dst_beg)));
        case LONG_TP: return reinterpret_cast<void*>(std::copy(beg, end, reinterpret_cast<long*>(dst_beg)));
        case UCHAR_TP: return reinterpret_cast<void*>(std::copy(beg, end, reinterpret_cast<unsigned char*>(dst_beg)));
        case USHORT_TP: return reinterpret_cast<void*>(std::copy(beg, end, reinterpret_cast<unsigned short*>(dst_beg)));
        case UINT_TP: return reinterpret_cast<void*>(std::copy(beg, end, reinterpret_cast<unsigned int*>(dst_beg)));
        case ULONG_TP: return reinterpret_cast<void*>(std::copy(beg, end, reinterpret_cast<unsigned long*>(dst_beg)));
        case FLOAT_TP: return reinterpret_cast<void*>(std::copy(beg, end, reinterpret_cast<float*>(dst_beg)));
        case DOUBLE_TP: return reinterpret_cast<void*>(std::copy(beg, end, reinterpret_cast<double*>(dst_beg)));
        case BIT_TP: {
            return internals::WriteBool<iterator>::apply(beg, end, start_dst, dst_data);  
        }
    }
    assert(0);
    return dst_data;
}

namespace internals{
    template<typename T, typename = void>
    struct is_container : std::false_type {};
    template<typename T>
    struct is_container<T, std::void_t<
        decltype(std::declval<T>().size()),
        decltype(std::declval<T>().begin()),
        decltype(std::declval<T>().end()),
        decltype(std::declval<T>().cbegin()),
        decltype(std::declval<T>().cend())
    >> : std::true_type {};

    template<typename T, typename = void>
    struct is_resizable_container : std::false_type {};
    template<typename T>
    struct is_resizable_container<T, std::void_t<decltype(std::declval<T>().resize())>>: is_container<T> {};

    template<typename T, unsigned isContainer = 0, typename Dummy = void>
    struct ProxyIterableImpl {
        T& m_val;

        ProxyIterableImpl(T& val): m_val{val} {};
        inline T* begin() { return &m_val; }
        inline T* end() { return (&m_val) + 1; }
        inline const T* cbegin() const { return &m_val; }
        inline const T* cend() const { return (&m_val) + 1; }
        inline const T* begin() const { return &m_val; }
        inline const T* end() const { return (&m_val) + 1; }
        inline std::size_t size() const { return 1; }
        inline void resize(std::size_t sz) { if (sz != size()) throw std::runtime_error("Can't resize to size = " + std::to_string(sz) + " object of size = " + std::to_string(size())); }
        inline void setup() {}
    };
    template<typename T>
    struct ProxyIterableImpl<T, 1, void>{
        T& m_val;
        ProxyIterableImpl(T& val): m_val{val} {};
        inline auto begin() { return m_val.begin(); }
        inline auto end() { return m_val.end(); }
        inline auto cbegin() const { return m_val.cbegin(); }
        inline auto cend() const { return m_val.cend(); }
        inline auto begin() const { return m_val.begin(); }
        inline auto end() const { return m_val.end(); }
        inline auto size() const { return m_val.size(); }
        inline void resize(std::size_t sz) { if (sz != size()) throw std::runtime_error("Can't resize to size = " + std::to_string(sz) + " object of size = " + std::to_string(size())); }
        inline void setup() {}
    };
    template<typename T>
    struct ProxyIterableImpl<T, 2, void>{
        T& m_val;
        ProxyIterableImpl(T& val): m_val{val} {};
        inline auto begin() { return m_val.begin(); }
        inline auto end() { return m_val.end(); }
        inline auto cbegin() const { return m_val.cbegin(); }
        inline auto cend() const { return m_val.cend(); }
        inline auto begin() const { return m_val.begin(); }
        inline auto end() const { return m_val.end(); }
        inline auto size() const { return m_val.size(); }
        inline void resize(std::size_t sz) { m_val.resize(sz); }
        inline void setup() {}
    };
    template<typename T>
    struct ProxyIterableImpl<T, 0, typename std::enable_if<std::is_same<T, Point>::value || std::is_same<T, Vector>::value>::type >{
        T& m_val;

        std::array<DReal, 3> m_proxy;
        ProxyIterableImpl(T& val): m_val{val}, m_proxy{val[0], val[1], val[2]} {};
        inline auto begin() { return m_proxy.begin(); }
        inline auto end() { return m_proxy.end(); }
        inline auto cbegin() const { return m_proxy.cbegin(); }
        inline auto cend() const { return m_proxy.cend(); }
        inline auto begin() const { return m_proxy.begin(); }
        inline auto end() const { return m_proxy.end(); }
        inline auto size() const { return m_proxy.size(); }
        inline void resize(std::size_t sz) { if (sz != size()) throw std::runtime_error("Can't resize to size = " + std::to_string(sz) + " object of size = " + std::to_string(size())); }
        inline void setup() { m_val = T(m_proxy[0], m_proxy[1], m_proxy[2]); }
    };
    template<typename T>
    struct ProxyIterableImpl<T, 0, typename std::enable_if<std::is_same<T, Point_2>::value || std::is_same<T, Vector_2>::value>::type >{
        T& m_val;

        std::array<DReal, 2> m_proxy;
        ProxyIterableImpl(T& val): m_val{val}, m_proxy{val[0], val[1]} {};
        inline auto begin() { return m_proxy.begin(); }
        inline auto end() { return m_proxy.end(); }
        inline auto cbegin() const { return m_proxy.cbegin(); }
        inline auto cend() const { return m_proxy.cend(); }
        inline auto begin() const { return m_proxy.begin(); }
        inline auto end() const { return m_proxy.end(); }
        inline auto size() const { return m_proxy.size(); }
        inline void resize(std::size_t sz) { if (sz != size()) throw std::runtime_error("Can't resize to size = " + std::to_string(sz) + " object of size = " + std::to_string(size())); }
        inline void setup() { m_val = T(m_proxy[0], m_proxy[1]); }
    };

    template<typename T>
    using ProxyIterable = ProxyIterableImpl<T, (is_container<T>::value ? 1 : 0) + (is_resizable_container<T>::value ? 1 : 0)>;

    template<typename Geom_index_type, bool Dummy = true>
    struct GetTagsContainer{
        MetaTriMesh* m;
        GetTagsContainer(MetaTriMesh* m): m{m} {}
        std::map<std::string, MetaTriMesh::MeshValue>& operator()() { static_assert("Unsupported CGAL::SurfaceMesh geometric element type"); return m->point_data; }
    };

    template<bool Dummy>
    struct GetTagsContainer<V_ind, Dummy>{
        MetaTriMesh* m;
        GetTagsContainer(MetaTriMesh* m): m{m} {}
        std::map<std::string, MetaTriMesh::MeshValue>& operator()() { return m->point_data; }
        const std::map<std::string, MetaTriMesh::MeshValue>& operator()() const { return m->point_data; }
        std::size_t number_of_element() const { return m->points.m_size; }
    };
    template<bool Dummy>
    struct GetTagsContainer<E_ind, Dummy>{
        MetaTriMesh* m;
        GetTagsContainer(MetaTriMesh* m): m{m} {}
        std::map<std::string, MetaTriMesh::MeshValue>& operator()() { return m->edge_data; }
        const std::map<std::string, MetaTriMesh::MeshValue>& operator()() const { return m->edge_data; }
        std::size_t number_of_element() const { return m->edges.m_size; }
    };
    template<bool Dummy>
    struct GetTagsContainer<HE_ind, Dummy>{
        MetaTriMesh* m;
        GetTagsContainer(MetaTriMesh* m): m{m} {}
        std::map<std::string, MetaTriMesh::MeshValue>& operator()() { return m->halfedge_data; }
        const std::map<std::string, MetaTriMesh::MeshValue>& operator()() const { return m->halfedge_data; }
        std::size_t number_of_element() const { return 2*m->edges.m_size; }
    };
    template<bool Dummy>
    struct GetTagsContainer<F_ind, Dummy>{
        MetaTriMesh* m;
        GetTagsContainer(MetaTriMesh* m): m{m} {}
        std::map<std::string, MetaTriMesh::MeshValue>& operator()() { return m->cell_data; }
        const std::map<std::string, MetaTriMesh::MeshValue>& operator()() const { return m->cell_data; }
        std::size_t number_of_element() const { return m->cells.m_size; }
    };

    template<typename Geom_index_type>
    struct GetMeshElements{
        const Mesh* m;
        GetMeshElements(const Mesh* m): m{m} {}
        std::size_t number() const { static_assert("Unsupported CGAL::SurfaceMesh geometric element type"); return 0; }
    };
    template<>
    struct GetMeshElements<V_ind>{
        const Mesh* m;
        GetMeshElements(const Mesh* m): m{m} {}
        std::size_t number() const { return m->number_of_vertices(); }
        auto begin() const { return m->vertices_begin(); }
        auto end() const { return m->vertices_end(); }
        auto range() const { return m->vertices(); }
        static MetaTriMesh::ElementSparsity inner_index() { return MetaTriMesh::NODE; }
    };
    template<>
    struct GetMeshElements<E_ind>{
        const Mesh* m;
        GetMeshElements(const Mesh* m): m{m} {}
        std::size_t number() const { return m->number_of_edges(); }
        auto begin() const { return m->edges_begin(); }
        auto end() const { return m->edges_end(); }
        auto range() const { return m->edges(); }
        static MetaTriMesh::ElementSparsity inner_index() { return MetaTriMesh::EDGE; }
    };
    template<>
    struct GetMeshElements<F_ind>{
        const Mesh* m;
        GetMeshElements(const Mesh* m): m{m} {}
        std::size_t number() const { return m->number_of_faces(); }
        auto begin() const { return m->faces_begin(); }
        auto end() const { return m->faces_end(); }
        auto range() const { return m->faces(); }
        static MetaTriMesh::ElementSparsity inner_index() { return MetaTriMesh::FACE; }
    };
    template<>
    struct GetMeshElements<HE_ind>{
        const Mesh* m;
        GetMeshElements(const Mesh* m): m{m} {}
        std::size_t number() const { return m->number_of_halfedges(); }
        auto begin() const { return m->halfedges_begin(); }
        auto end() const { return m->halfedges_end(); }
        auto range() const { return m->halfedges(); }
        static MetaTriMesh::ElementSparsity inner_index() { return MetaTriMesh::HALFEDGE; }
    };
}

template<typename T> void MetaTriMesh::MeshValue::get_element(std::size_t i, T& val) const {
    internals::ProxyIterable<T> p(val);
    p.resize(m_values.m_dim);
    get_values(get_cont_index(i, 0), get_cont_index(i+1, 0), p.begin());
    p.setup();
}

template<typename iterator> void MetaTriMesh::MeshValue::get_values(std::size_t start_cont_index, std::size_t end_cont_index, iterator dst) const {
    switch (m_type){
        case CONSTANT: {
            std::vector<typename std::iterator_traits<iterator>::value_type> vals; vals.resize(m_values.dim());
            m_values.get_values(0, m_values.m_dim, vals.begin());
            for (auto i = start_cont_index; i < end_cont_index; ++i)
                *(dst++) = vals[i%m_values.m_dim];
            return; 
        }
        case CONTIGUOUS: {
            m_values.get_values(start_cont_index, end_cont_index, dst);
            return;
        } 
        case SPARSE: {
            std::vector<typename std::iterator_traits<iterator>::value_type> def_vals(m_values.m_dim);
            std::vector<typename std::iterator_traits<iterator>::value_type> cpy_vals(m_values.m_dim);
            m_values.get_values(0, m_values.m_dim, def_vals.begin());
            auto sparse_id = end_cont_index / m_values.m_dim + 1;
            int index_id = 0;
            if (m_sparse_indexes.m_size)
                sparse_id = m_sparse_indexes.at(index_id).get<std::size_t>();
            auto val_id = start_cont_index/m_values.m_dim;
            
            auto move_to_next_sparse_id = [&](){
                for (; sparse_id < val_id && index_id < m_sparse_indexes.m_size; ++index_id)
                    sparse_id = m_sparse_indexes.at(index_id).get<std::size_t>();
                if (index_id == m_sparse_indexes.m_size) 
                    sparse_id = end_cont_index / m_values.m_dim + 1;
            };
            const typename std::iterator_traits<iterator>::value_type* from_data = def_vals.data();
            auto choose_data_array = [&](){
                move_to_next_sparse_id();
                if (val_id == sparse_id){
                    m_values.get_values(m_values.get_cont_index(index_id, 0), m_values.get_cont_index(index_id+1, 0), cpy_vals.begin());
                    from_data = cpy_vals.data();
                } else 
                    from_data = def_vals.data();
            };
            if (start_cont_index%m_values.m_dim != 0){
                choose_data_array();
                dst = std::copy(from_data + start_cont_index%m_values.m_dim, 
                                from_data + start_cont_index%m_values.m_dim + std::min(m_values.m_dim - start_cont_index%m_values.m_dim, end_cont_index - start_cont_index),
                                dst);    
                ++val_id;    
            }   
            for(auto end_val_id = end_cont_index/m_values.m_dim; val_id < end_val_id; ++val_id){
                choose_data_array();
                dst = std::copy(from_data, from_data + m_values.m_dim, dst);
            } 
            if (end_cont_index/m_values.m_dim != start_cont_index/m_values.m_dim && end_cont_index%m_values.m_dim != 0){
                choose_data_array();
                dst = std::copy(from_data + val_id*m_values.m_dim, from_data + end_cont_index, dst); 
                ++val_id;   
            }
        }
    }
    return;
}

inline MetaTriMesh::Value MetaTriMesh::Storage::get_value(std::size_t i, std::size_t d, DataType type, std::size_t size, std::size_t dim, const char* data, std::size_t off){
    auto dt_sz = get_type_sizeof(type);
    off = off + (i*dim + d) * dt_sz;
    const void* val = nullptr;
    bool bval = false;
    if (type == BIT_TP){
        bval = (data[off/8] & (1 >> (off%8)));
        val = &bval;
    } else 
        val = data + off;
    return Value().set(type, val);    
}

template<typename T> void MetaTriMesh::Storage::get_element(std::size_t i, T& dst) const {
    internals::ProxyIterable<T> p(dst);
    p.resize(dim());
    get_values(get_cont_index(i, 0), get_cont_index(i+1, 0), p.begin());
    p.setup();
}

template<unsigned N>
template<typename T> void MetaTriMesh::TStorage<N>::get_element(std::size_t i, T& dst) const {
    internals::ProxyIterable<T> p(dst);
    p.resize(dim());
    get_values(get_cont_index(i, 0), get_cont_index(i+1, 0), p.begin());
    p.setup();
}

template<typename T> void MetaTriMesh::Storage::set_element(std::size_t i, const T& val){
    internals::ProxyIterable<T> p(const_cast<T&>(val));
    if (p.size() != dim()) throw std::runtime_error("Expected " + std::to_string(m_dim) + " values, but get " + std::to_string(p.size()));
    set_values(p.cbegin(), p.cend(), get_cont_index(i, 0));
}

template<unsigned N>
template<typename T> void MetaTriMesh::TStorage<N>::set_element(std::size_t i, const T& val){
    internals::ProxyIterable<T> p(const_cast<T&>(val));
    if (p.size() != dim()) throw std::runtime_error("Expected " + std::to_string(dim()) + " values, but get " + std::to_string(p.size()));
    set_values(p.cbegin(), p.cend(), get_cont_index(i, 0));
}

template<typename Geom_index_type, typename T>
bool MetaTriMesh::createTagOnMesh(Mesh* m, std::string meta_tag_name, std::string mesh_tag_name, bool overwrite) const { 
    if (mesh_tag_name.empty()) mesh_tag_name = meta_tag_name;
    internals::GetTagsContainer<Geom_index_type> getter(const_cast<MetaTriMesh*>(this));
    auto& data = getter();
    std::size_t nelems = getter.number_of_element();
    auto it = data.find(meta_tag_name);
    if (it == data.end()) return false;
    auto& d = it->second;
    auto q = m->add_property_map<Geom_index_type, T>(mesh_tag_name);
    if (!q.second && !overwrite) return false;
    for (std::size_t i = 0; i < nelems; ++i){
        T p = T();
        d.get_element(i, p);
        q.first[Geom_index_type(i)] = p;
    }
    return true; 
}
template<typename T>
bool MetaTriMesh::createTagOnMeshFromSelfPoints(Mesh* m, std::string mesh_tag_name, bool overwrite) const{
    auto q = m->add_property_map<V_ind, T>(mesh_tag_name);
    if (!q.second && !overwrite) return false;
    for (std::size_t i = 0; i < points.m_size; ++i){
        T p = T();
        points.get_element(i, p);
        q.first[V_ind(i)] = p;
    }
    return true;
}

template<unsigned N>
MetaTriMesh::TStorage<N>::TStorage(MetaTriMesh::Storage&& a){
    if (N != a.dim()) throw std::runtime_error("Wrong type conversion");
    m_type = a.m_type; a.m_type = NONE_TP;
    m_size = a.m_size; a.m_size = 0;
    m_data = std::move(a.m_data);
    a.m_dim = 0;
}
template<unsigned N>
MetaTriMesh::TStorage<N>::TStorage(const MetaTriMesh::Storage& a){
    if (N != a.dim()) throw std::runtime_error("Wrong type conversion");
    m_type = a.m_type;
    m_size = a.m_size;
    m_data = a.m_data;
}

template<unsigned N>
MetaTriMesh::Storage::Storage(TStorage<N>&& a): m_type{a.m_type}, m_size{a.m_size}, m_dim{N}, m_data{std::move(a.m_data)}{ a.m_type = NONE_TP; a.m_size = 0; }
template<unsigned N>
MetaTriMesh::Storage::Storage(const TStorage<N>& a): m_type{a.m_type}, m_size{a.m_size}, m_dim{N}, m_data{a.m_data}{ a.m_type = NONE_TP; a.m_size = 0; }

template<typename T>
MetaTriMesh::Storage MetaTriMesh::DefaultConverter<T>::operator()(const Mesh& m, std::string tag_name, ElementSparsity sp){
    Storage val;
    val.m_dim = dinfo.m_dim;
    val.m_type = dinfo.m_type;
    switch (sp){
        case NODE:{
            auto tag_p = m.property_map<V_ind, T>(tag_name);
            if (!tag_p.second) return Storage();
            auto& tag = tag_p.first;
            val.resize(m.number_of_vertices());
            std::size_t i = 0;
            for (auto v: m.vertices())
                val.set_element(i++, tag[v]);
            return val;
        }
        case EDGE:{
            auto tag_p = m.property_map<E_ind, T>(tag_name);
            if (!tag_p.second) return Storage();
            auto& tag = tag_p.first;
            val.resize(m.number_of_edges());
            std::size_t i = 0;
            for (auto e: m.edges())
                val.set_element(i++, tag[e]);
            return val;
        }
        case FACE:{
            auto tag_p = m.property_map<F_ind, T>(tag_name);
            if (!tag_p.second) return Storage();
            auto& tag = tag_p.first;
            val.resize(m.number_of_faces());
            std::size_t i = 0;
            for (auto f: m.faces())
                val.set_element(i++, tag[f]);
            return val;
        }
        case HALFEDGE:{
            auto tag_p = m.property_map<HE_ind, T>(tag_name);
            if (!tag_p.second) return Storage();
            auto& tag = tag_p.first;
            val.resize(m.number_of_halfedges());
            std::size_t i = 0;
            for (auto h: m.halfedges())
                val.set_element(i++, tag[h]);
            return val;
        }
        default: 
            throw std::runtime_error("Achieved unreachable code");
    }
    val.clear();
    return val;
}

template<typename Geom_index_type>
MetaTriMesh::MeshValue MetaTriMesh::constructTagDataFromMesh(const Mesh& m, std::string name){
    MeshValue val;
    const std::type_info& info = const_cast<Mesh&>(m).property_type<Geom_index_type>(name);
    if (info == typeid(void)) return val;
    auto& conv = get_converter(); 
    auto it = conv.find(std::type_index(info));
    if (it == conv.end()) return val;
    val.m_type = MeshValue::CONTIGUOUS;
    val.m_values = it->second(m, name, internals::GetMeshElements<Geom_index_type>::inner_index());
    return val;
}
template<typename Geom_index_type>
MetaTriMesh& MetaTriMesh::readTagDataFromMesh(const Mesh& m, const std::string& name, const std::string& meta_tag_name){
    auto v = constructTagDataFromMesh<Geom_index_type>(m, name);
    if (v.isValid()) 
            internals::GetTagsContainer<Geom_index_type>(this)()[meta_tag_name] = std::move(v);
    else
        throw std::runtime_error("Can't construct MeshValue from tag \"" + name + "\"");
    return *this;            
}
template<typename Geom_index_type>
MetaTriMesh& MetaTriMesh::readTagsDataFromMesh(const Mesh& m, const std::set<std::string>& exclude_names){
    auto names = m.properties<Geom_index_type>();
    for (auto nm: names) if (exclude_names.find(nm) == exclude_names.end()){
        MeshValue v = constructTagDataFromMesh<Geom_index_type>(m, nm);
        if (v.isValid()) 
            internals::GetTagsContainer<Geom_index_type>(this)()[nm] = std::move(v);
    }
    return *this;
}

template<typename Geom_index_type>
std::set<std::string> MetaTriMesh::value_names() const {
    auto& p = internals::GetTagsContainer<Geom_index_type>(const_cast<MetaTriMesh*>(this))();
    std::set<std::string> res;
    for (auto v: p) res.insert(v.first);
    return res;
}

inline std::map<std::string, MetaTriMesh::MeshValue>& MetaTriMesh::data_array(ElementSparsity e){
    switch(e){
        case NODE: return point_data;
        case EDGE: return edge_data;
        case FACE: return cell_data;
        case HALFEDGE: return halfedge_data;
    }
    return point_data;
}
inline const std::map<std::string, MetaTriMesh::MeshValue>& MetaTriMesh::data_array(ElementSparsity e) const{
    switch(e){
        case NODE: return point_data;
        case EDGE: return edge_data;
        case FACE: return cell_data;
        case HALFEDGE: return halfedge_data;
    }
    return point_data;
}

}


#endif //AORTIC_VALVE_METATRIMESH_INL