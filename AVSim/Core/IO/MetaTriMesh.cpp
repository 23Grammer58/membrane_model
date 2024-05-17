#include "MetaTriMesh.h"

using namespace World3d;
using MTM = MetaTriMesh;
using namespace std;

static std::map<std::type_index, MetaTriMesh::TagDataConverter> s_converters;

std::map<std::type_index, MetaTriMesh::TagDataConverter>& MetaTriMesh::get_converter() { return s_converters; }
void MetaTriMesh::register_tag_converter(std::type_index tindex, TagDataConverter conv){ s_converters[tindex] = std::move(conv); }
void MetaTriMesh::_init_tag_converters(bool enforce) {
    if (!s_converters.empty() && !enforce) return;
#define CONV_VAL(FROM_T, TO_T, CNT) {std::type_index(typeid(FROM_T)), TagDataConverter(DefaultConverter<FROM_T>(TO_T, CNT))}   
#define CONV_ARR(FROM_T, TO_T, CNT) {std::type_index(typeid(std::array<FROM_T, CNT>)), TagDataConverter(DefaultConverter<std::array<FROM_T, CNT>>(TO_T, CNT))} 
    s_converters = std::map<std::type_index, TagDataConverter>{
        CONV_VAL(Point   , DOUBLE_TP, 3),
        CONV_VAL(Vector  , DOUBLE_TP, 3),
        CONV_VAL(double  , DOUBLE_TP, 1),
        CONV_VAL(float   , FLOAT_TP , 1),
        CONV_VAL(bool    , BIT_TP   , 1),
        CONV_VAL(char    , CHAR_TP  , 1),
        CONV_VAL(short   , SHORT_TP , 1),
        CONV_VAL(int     , INT_TP   , 1),
        CONV_VAL(long    , LONG_TP  , 1),
        CONV_VAL(Point_2 , DOUBLE_TP, 2),
        CONV_VAL(Vector_2, DOUBLE_TP, 2),
        CONV_ARR(double  , DOUBLE_TP, 1),
        CONV_ARR(double  , DOUBLE_TP, 2),
        CONV_ARR(double  , DOUBLE_TP, 3),
        CONV_ARR(double  , DOUBLE_TP, 4),
        CONV_ARR(double  , DOUBLE_TP, 5),
        CONV_ARR(double  , DOUBLE_TP, 6),
        CONV_ARR(double  , DOUBLE_TP, 7),
        CONV_ARR(double  , DOUBLE_TP, 8),
        CONV_ARR(double  , DOUBLE_TP, 9),
        CONV_ARR(double  , DOUBLE_TP, 12),
        CONV_ARR(double  , DOUBLE_TP, 16),
        CONV_ARR(double  , DOUBLE_TP, 24),
        CONV_ARR(double  , DOUBLE_TP, 27),
        CONV_ARR(double  , DOUBLE_TP, 32),
        CONV_ARR(double  , DOUBLE_TP, 36),
        CONV_ARR(double  , DOUBLE_TP, 45),
        CONV_ARR(double  , DOUBLE_TP, 64),
        CONV_ARR(double  , DOUBLE_TP, 81),
        CONV_ARR(float   , FLOAT_TP , 1),
        CONV_ARR(float   , FLOAT_TP , 2),
        CONV_ARR(float   , FLOAT_TP , 3),
        CONV_ARR(float   , FLOAT_TP , 4),
        CONV_ARR(float   , FLOAT_TP , 5),
        CONV_ARR(float   , FLOAT_TP , 6),
        CONV_ARR(float   , FLOAT_TP , 7),
        CONV_ARR(float   , FLOAT_TP , 8),
        CONV_ARR(float   , FLOAT_TP , 9),
        CONV_ARR(float   , FLOAT_TP , 12),
        CONV_ARR(float   , FLOAT_TP , 16),
        CONV_ARR(float   , FLOAT_TP , 24),
        CONV_ARR(float   , FLOAT_TP , 27),
        CONV_ARR(float   , FLOAT_TP , 32),
        CONV_ARR(float   , FLOAT_TP , 36),
        CONV_ARR(float   , FLOAT_TP , 45),
        CONV_ARR(float   , FLOAT_TP , 64),
        CONV_ARR(float   , FLOAT_TP , 81),
        CONV_ARR(int     , INT_TP   , 3),
        CONV_ARR(long    , LONG_TP  , 3),
        CONV_VAL(unsigned int , UINT_TP , 1),
        CONV_VAL(unsigned long, ULONG_TP, 1),
        CONV_ARR(unsigned int , UINT_TP , 3),
        CONV_ARR(unsigned long, ULONG_TP, 3),
    };
#undef CONV_ARR    
#undef CONV_VAL
}

MetaTriMesh::TStorage<2> MetaTriMesh::compute_edges() const {
    if (edges.m_size != 0) return edges;
    if (cells.m_size == 0) return TStorage<2>();
    TStorage<2> res;
    std::set<std::pair<std::size_t, std::size_t>> i_edges;
    for (std::size_t i = 0; i < cells.m_size; ++i){
        std::array<std::size_t, 3> tri = {
            cells.at(i, 0).get<std::size_t>(),
            cells.at(i, 1).get<std::size_t>(),
            cells.at(i, 2).get<std::size_t>()
        };
        for (int k = 0; k < 3; ++k)
            i_edges.insert({std::min(tri[k], tri[(k+1)%3]), std::max(tri[k], tri[(k+1)%3])});
    }
    res.m_type = cells.m_type;
    res.resize(i_edges.size());
    std::size_t i = 0;
    for (auto v: i_edges){
        res.set(i, 0, v.first);
        res.set(i, 1, v.second);
        ++i;
    }
    return res;
}

void MetaTriMesh::setup_edges(){
    if (cells.m_size == 0 || edges.m_size != 0) return;
    edges = std::move(compute_edges());
}

Mesh MetaTriMesh::convertToMesh() const {
    const TStorage<2>* pedges = &edges;
    TStorage<2> cedges;
    if (edges.m_size == 0){
        cedges = std::move(compute_edges());
        pedges = &cedges;
    }
    Mesh m;
    for (std::size_t i = 0; i < points.m_size; ++i){
        std::array<double, 3> p;
        points.get_values(points.get_cont_index(i, 0), points.get_cont_index(i, 3), p.begin());
        m.add_vertex(Point{p[0], p[1], p[2]});
    }

    std::map<std::pair<unsigned, unsigned>, HE_ind> hmap;
    for (std::size_t i = 0; i < pedges->m_size; ++i){ 
        std::array<std::size_t, 2> p;
        pedges->get_values(pedges->get_cont_index(i, 0), pedges->get_cont_index(i, 2), p.begin());
        std::array<V_ind, 2> v{V_ind(p[0]), V_ind(p[1])};
        auto hid  = m.add_edge(v[0], v[1]);
        hmap[{p[0], p[1]}] = hid;
        hmap[{p[1], p[0]}] = m.opposite(hid);;
    }
    for (std::size_t i = 0; i < cells.m_size; ++i){
        std::array<std::size_t, 3> p;
        cells.get_values(cells.get_cont_index(i, 0), cells.get_cont_index(i, 3), p.begin());
        using hmap_iterator = std::map<std::pair<unsigned, unsigned>, HE_ind>::iterator;
        std::array<hmap_iterator, 3> hit = {hmap.find({p[0], p[1]}), hmap.find({p[1], p[2]}), hmap.find({p[2], p[0]})};
        int k_ = -1;
        if (( (hit[0]==hmap.end()) && (k_ = 0)) || ( (hit[1]==hmap.end()) && (k_ = 1)) || ( (hit[2]==hmap.end()) && (k_ = 2)))
            throw std::runtime_error("Face {" + std::to_string(p[0]) + ", " + std::to_string(p[1]) + ", " + std::to_string(p[2]) + " contains undefined edge {"
                + std::to_string(p[k_]) + ", " + std::to_string(p[(k_+1)%3]) + "}");
        std::array<HE_ind, 3> h{hit[0]->second, hit[1]->second, hit[2]->second};        
        F_ind f = m.add_face();
        m.set_halfedge(f, h[0]);
        for (unsigned k = 0; k < h.size(); ++k)
            m.set_face(h[k], f);
        for (unsigned k = 0; k < h.size(); ++k)
            m.set_next(h[k], h[(k+1)%(h.size())]);
        for (unsigned k = 0; k < 3; ++k) 
            hmap.erase(hit[k]);  
    }
    std::map<unsigned, HE_ind> cycles;
    for (auto h: hmap)
        cycles.insert({h.first.first, h.second});
    hmap.clear();
    for (auto h: m.halfedges()){
        auto f = m.face(h);
        auto v = m.target(h);
        if (!f.is_valid()){
            m.set_halfedge(v, h);
            auto nh = m.next(h);
            if (!nh.is_valid()){
                auto cit = cycles.find(v.idx());
                if (cit == cycles.end())
                    throw std::runtime_error("Can't find next vertex for boundary halfedge");
                nh = cit->second;
                m.set_next(h, nh);
                cycles.erase(cit);
            }
        } else {
            auto h1 = m.halfedge(v);
            if (!h1.is_valid())
                m.set_halfedge(v, h);
        }
    }
    // for (auto& tt: point_data) {
    //     auto& t = tt.second;
    //     if (tt.first == "v:x"){
    //         auto m_x = m.add_property_map<V_ind, Point>("v:x").first;
    //         for (std::size_t i = 0; i < points.m_size; ++i){
    //             std::array<double, 3> p;
    //             t.get_values(t.get_cont_index(i, 0), t.get_cont_index(i, 3), p.begin());
    //             m_x[V_ind(i)] = Point(p[0], p[1], p[2]);
    //         }
    //     } else if (tt.first == "v:boundary_lbl"){
    //         auto m_lbl = m.add_property_map<V_ind, int>("v:boundary_lbl").first;
    //         for (std::size_t i = 0; i < points.m_size; ++i){
    //             int p = 0;
    //             t.get_values(t.get_cont_index(i, 0), t.get_cont_index(i, 1), &p);
    //             m_lbl[V_ind(i)] = p;
    //         }
    //     }
    // }
    return m;
}

MetaTriMesh& MetaTriMesh::initFromMesh(const Mesh& m, std::string point_tag){
    points = TStorage<3>(constructTagDataFromMesh<V_ind>(m, point_tag).m_values);
    edges.m_type = UINT_TP;
    cells.m_type = UINT_TP;
    edges.resize(m.number_of_edges());
    cells.resize(m.number_of_faces());
    auto vit = std::max_element(m.vertices().begin(), m.vertices().end());
    std::vector<int> loc_ind(*vit+1, -1);
    {
        int i = 0;
        for (auto v: m.vertices()) loc_ind[v] = i++;
    }
    std::size_t ei = 0;
    for (auto ec: m.edges()){
        auto v = vert_around(m, ec);
        std::array<int, 2> vi{loc_ind[v[0]], loc_ind[v[1]]};
        edges.set_element(ei++, vi);
    }
    std::size_t fi = 0;
    for (auto fc: m.faces()){
        auto v = vert_around(m, fc);
        std::array<int, 3> vi{loc_ind[v[0]], loc_ind[v[1]], loc_ind[v[2]]};
        cells.set_element(fi++, vi);
    }

    return *this;
}

MetaTriMesh& MetaTriMesh::readAllTagsFromMesh(const Mesh& m, const std::set<std::string>& exclude_names){
    readTagsDataFromMesh<V_ind>(m, exclude_names);
    readTagsDataFromMesh<E_ind>(m, exclude_names);
    readTagsDataFromMesh<F_ind>(m, exclude_names);
    readTagsDataFromMesh<HE_ind>(m, exclude_names);
    return *this;
}

void MetaTriMesh::Storage::read_binary_data(std::istream& in, DataType type, std::size_t size, std::size_t dim, std::vector<char>& out, Endian order){
    auto szof = get_type_sizeof(type);
    std::size_t sz = dim*size*szof; 
    if (type == BIT_TP) 
        sz = (sz/8) + ((sz%8 > 0) ? 1 : 0);
    out.resize(sz);
    in.read(out.data(), sz);
    if (szof > 1 && order != getEndianOrder()){
        for (std::size_t i = 0, cnt = dim*size; i < cnt; ++i){
            for (std::size_t k = 0, half_dim = szof/2; k < half_dim; ++k)
                std::swap(out[i*szof + k], out[(i+1)*szof - k - 1]);
        }
    }
}
void MetaTriMesh::Storage::read_ascii_data(std::istream& in, DataType type, std::size_t size, std::size_t dim, std::vector<char>& out, std::function<void(std::istream& in)> skip_after_element){
    std::size_t sz = dim*size*get_type_sizeof(type); 
    if (type == BIT_TP) 
        sz = (sz/8) + ((sz%8 > 0) ? 1 : 0);
    out.resize(sz);
    switch(type)
    {
#define CASE_TP(TYPE_ID, CPP_TP) \
        case TYPE_ID:{\
            auto* data = reinterpret_cast<CPP_TP*>(out.data());\
            for (size_t i = 0; i < size; ++i){\
                for (size_t k = 0; k < dim; ++k)\
                    in >> *(data++);\
                if (skip_after_element) skip_after_element(in);\
            }\
            break;\
        }
#define CASE_TP_R(TYPE_ID, CPP_TP, PROXY_TP) \
        case TYPE_ID:{\
            auto* data = reinterpret_cast<CPP_TP*>(out.data());\
            for (size_t i = 0; i < size; ++i){\
                for (size_t k = 0; k < dim; ++k){\
                    PROXY_TP t = 0;\
                    in >> t;\
                    *(data++) = static_cast<CPP_TP>(t);\
                }\
                if (skip_after_element) skip_after_element(in);\
            }\
            break;\
        }        

        CASE_TP(DOUBLE_TP, double)
        CASE_TP(INT_TP, int)
        CASE_TP_R(CHAR_TP, char, short)
        CASE_TP(LONG_TP, long)
        CASE_TP(FLOAT_TP, float)
        CASE_TP(SHORT_TP, short)
        CASE_TP_R(UCHAR_TP, unsigned char, unsigned short)
        CASE_TP(USHORT_TP, unsigned short)
        CASE_TP(UINT_TP, unsigned int)
        CASE_TP(ULONG_TP, unsigned long)
#undef CASE_TP_R        
#undef CASE_TP  
        case BIT_TP:{
            auto* data = reinterpret_cast<char*>(out.data());
            std::fill(out.begin(), out.end(), 0);
            for (size_t i = 0; i < size; ++i){
                for (size_t k = 0; k < dim; ++k){
                    std::string v;
                    in >> v;
                    std::transform(v.begin(), v.end(), v.begin(), [](char c) {return std::tolower(c, std::locale()); });
                    if (v != "0" && v != "false"){
                        auto ip = k + i * dim;
                        data[ip/8] |= (1 << (ip%8));
                    }
                }
                if (skip_after_element) skip_after_element(in);    
            }
            break;
        }
        default:
            throw std::runtime_error("Faced unknown data type");
    }
}
void MetaTriMesh::Storage::write_binary_data(std::ostream& out, DataType type, std::size_t size, std::size_t dim, const std::vector<char>& data, Endian order){
    auto szof = get_type_sizeof(type);
    std::size_t sz = dim*size*szof; 
    if (type == BIT_TP) 
        sz = (sz/8) + ((sz%8 > 0) ? 1 : 0);
    if (data.size() < sz)
        throw std::runtime_error("Wrong data size");
    if (szof <= 1 || order == getEndianOrder())    
        out.write(data.data(), sz); 
    else{
        char buf[sizeof(unsigned long) + sizeof(double)];
        for (std::size_t i = 0, cnt = dim*size; i < cnt; ++i){
            for (std::size_t k = 0; k < szof; ++k)
                buf[k] = data[(i+1)*szof - k - 1];
            out.write(buf, szof);    
        }
    }       
}
void MetaTriMesh::Storage::write_ascii_data(std::ostream& out, DataType type, std::size_t size, std::size_t dim, const std::vector<char>& data, std::string v_delim, std::string e_delim){
    std::size_t sz = dim*size*get_type_sizeof(type); 
    if (type == BIT_TP) 
        sz = (sz/8) + ((sz%8 > 0) ? 1 : 0);
    if (data.size() < sz)
        throw std::runtime_error("Wrong data size");
    if (dim == 0) return;    
    switch(type)
    {
#define CASE_TP_G(TYPE_ID, CPP_TP, CONV) \
        case TYPE_ID:{\
            auto* buf = reinterpret_cast<const CPP_TP*>(data.data());\
            for (size_t i = 0; i < size; ++i){\
                for (size_t k = 0; k < dim - 1; ++k)\
                    out << CONV(*(buf++)) << v_delim;\
                out << CONV(*(buf++)) << e_delim;   \
            }\
            break;\
        }
#define NO_OP         
#define CASE_TP(TYPE_ID, CPP_TP)  CASE_TP_G(TYPE_ID, CPP_TP, NO_OP)       
#define CASE_TP_R(TYPE_ID, CPP_TP, PROXY_TP) CASE_TP_G(TYPE_ID, CPP_TP, static_cast<PROXY_TP>)
        CASE_TP(DOUBLE_TP, double)
        CASE_TP(INT_TP, int)
        CASE_TP_R(CHAR_TP, char, short)
        CASE_TP(LONG_TP, long)
        CASE_TP(FLOAT_TP, float)
        CASE_TP(SHORT_TP, short)
        CASE_TP_R(UCHAR_TP, unsigned char, unsigned short)
        CASE_TP(USHORT_TP, unsigned short)
        CASE_TP(UINT_TP, unsigned int)
        CASE_TP(ULONG_TP, unsigned long)
#undef CASE_TP
#undef NO_OP
#undef CASE_TP_G        
        case BIT_TP:{
            auto* buf = reinterpret_cast<const char*>(data.data());
            for (size_t i = 0; i < size; ++i){
                for (size_t k = 0; k < dim - 1; ++k){
                    auto ip = k + i * dim;
                    out << ((buf[ip/8] & (1 << (ip%8))) ? '1' : '0') << v_delim;
                }
                auto ip = (dim - 1) + i * dim;
                out << ((buf[ip/8] & (1 << (ip%8))) ? '1' : '0') << e_delim;   
            }
            break;
        }
        default:
            throw std::runtime_error("Faced unknown data type");
    }
}

MetaTriMesh::Value MetaTriMesh::Value::read_ascii_var(std::istream& in, DataType type){
    Value res; res.m_type = type;
    switch (type)
    {
        case BIT_TP: in >> res.m_val.b; break;
        case CHAR_TP: { short t; in >> t; res.m_val.c = static_cast<char>(t); break; }
        case SHORT_TP: in >> res.m_val.s; break;
        case INT_TP: in >> res.m_val.i; break;
        case LONG_TP: in >> res.m_val.l; break;
        case UCHAR_TP: { unsigned short t; in >> t; res.m_val.uc = static_cast<unsigned char>(t); break; }
        case USHORT_TP: in >> res.m_val.us; break;
        case UINT_TP: in >> res.m_val.ui; break;
        case ULONG_TP: in >> res.m_val.ul; break;
        case FLOAT_TP: in >> res.m_val.f; break;
        case DOUBLE_TP: in >> res.m_val.d; break;
        default:
            throw std::runtime_error("Faced unknown datatype");
    }
    return res;
}
MetaTriMesh::Value MetaTriMesh::Value::read_binary_var(std::istream& in, DataType type, Endian order){
    Value res; res.m_type = type;
    auto szof = MetaTriMesh::get_type_sizeof(type);
    in.read(reinterpret_cast<char*>(&res.m_val), szof);
    if (order != getEndianOrder())
        for (std::size_t k = 0, half_szof = szof/2; k < half_szof; ++k){
            auto* data = reinterpret_cast<char*>(&res.m_val);
            std::swap(data[k], data[szof - k - 1]);
        }
    return res;
}
void MetaTriMesh::Value::write_ascii_var(std::ostream& out, Value val){
    switch (val.m_type)
    {
        case BIT_TP: out << val.m_val.b; break;
        case CHAR_TP: out << static_cast<short>(val.m_val.c); break;
        case SHORT_TP: out << val.m_val.s; break;
        case INT_TP: out << val.m_val.i; break;
        case LONG_TP: out << val.m_val.l; break;
        case UCHAR_TP: out << static_cast<unsigned short>(val.m_val.uc); break;
        case USHORT_TP: out << val.m_val.us; break;
        case UINT_TP: out << val.m_val.ui; break;
        case ULONG_TP: out << val.m_val.ul; break;
        case FLOAT_TP: out << val.m_val.f; break;
        case DOUBLE_TP: out << val.m_val.d; break;
        default:
            throw std::runtime_error("Faced unknown datatype");
    }
}
void MetaTriMesh::Value::write_binary_var(std::ostream& out, Value val, Endian order){
    auto szof = MetaTriMesh::get_type_sizeof(val.m_type);
    if (szof <= 1 || order == getEndianOrder())
        out.write(reinterpret_cast<const char*>(&val.m_val), szof);
    else {
        char buf[sizeof(unsigned long) + sizeof(double)];
        for (std::size_t k = 0; k < szof; ++k)
            buf[k] = reinterpret_cast<const char*>(&val.m_val)[szof - k - 1];
        out.write(buf, szof);    
    }    
}

void MetaTriMesh::Storage::setV(std::size_t i, std::size_t d, Value val, std::size_t size, std::size_t dim, char* data, std::size_t off){
    switch (val.m_type)
    {
        case BIT_TP: set_value(i, d, val.m_val.b, val.m_type, size, dim, data, off); break;
        case CHAR_TP: set_value(i, d, val.m_val.c, val.m_type, size, dim, data, off); break;
        case SHORT_TP: set_value(i, d, val.m_val.s, val.m_type, size, dim, data, off); break;
        case INT_TP: set_value(i, d, val.m_val.i, val.m_type, size, dim, data, off); break;
        case LONG_TP: set_value(i, d, val.m_val.l, val.m_type, size, dim, data, off); break;
        case UCHAR_TP: set_value(i, d, val.m_val.uc, val.m_type, size, dim, data, off); break;
        case USHORT_TP: set_value(i, d, val.m_val.us, val.m_type, size, dim, data, off); break;
        case UINT_TP: set_value(i, d, val.m_val.ui, val.m_type, size, dim, data, off); break;
        case ULONG_TP: set_value(i, d, val.m_val.ul, val.m_type, size, dim, data, off); break;
        case FLOAT_TP: set_value(i, d, val.m_val.f, val.m_type, size, dim, data, off); break;
        case DOUBLE_TP: set_value(i, d, val.m_val.d, val.m_type, size, dim, data, off); break;
        default:
            throw std::runtime_error("Faced unknown datatype");
    }
}

MetaTriMesh::Storage MetaTriMesh::MeshValue::create_contiguous(const MeshValue& v, std::size_t expect_size){
    if (v.m_type == CONTIGUOUS){
        assert(v.m_values.size() == expect_size); 
        return v.m_values;
    }
    if (v.m_type != CONSTANT && v.m_type != SPARSE)
        return Storage();
    Storage res;
    res.m_dim = v.m_values.m_dim; res.m_type = v.m_values.m_type;
    res.resize(expect_size); 
    if (res.m_type != MetaTriMesh::BIT_TP){
        auto szof = get_type_sizeof(res.m_type);
        auto sz_elem = szof * res.m_dim;
        for (std::size_t i = 0; i < expect_size; ++i)
            std::copy(v.m_values.m_data.data(), v.m_values.m_data.data() + sz_elem, res.m_data.data() + i*sz_elem);
    } else {
        for (std::size_t i = 0; i < expect_size; ++i)
            for (std::size_t k = 0; k < res.m_dim; ++k)
                res.set(i, k, v.m_values.at(0, k).as<bool>());
    }
    if (v.m_type == CONSTANT) return res; 
    for (std::size_t i = 0; i < v.m_sparse_indexes.size(); ++i){
        std::size_t id = 0;
        v.m_sparse_indexes.get_element(i, id);
        for (std::size_t k = 0; k < res.m_dim; ++k)
            res.setV(id, k, v.m_values.at(i+1, k));
    }
    return res;   
}

std::ostream& MetaTriMesh::print_content_info(std::ostream& out, const std::string& prefix) const {
    out << prefix << "Mesh #P = " << points.size() << " #E = " << edges.size() << " #C = " << cells.size() << "\n";
    auto print_data_info = [&](const char* data_name, const std::map<std::string, MeshValue>& data){
        if (data.empty()) return;
        out << prefix << data_name << ":\n";
        for (const auto& v: data){
            out << prefix << "  \"" << v.first << "\": " << MeshValue::getTypeName(v.second.m_type);
            switch(v.second.m_type){
                case MeshValue::CONSTANT: break;
                case MeshValue::CONTIGUOUS: out << " size = " << v.second.m_values.size(); break;
                case MeshValue::SPARSE: out << " sparse_count = " << v.second.m_sparse_indexes.size(); break;
            }
            out << " dim = " << v.second.m_values.dim() << "\n";
        }
    };
    print_data_info("POINTS_DATA", point_data);
    print_data_info("EDGE_DATA", edge_data);
    print_data_info("HALFEDGE_DATA", halfedge_data);
    print_data_info("CELLS_DATA", cell_data);
    return out;
}