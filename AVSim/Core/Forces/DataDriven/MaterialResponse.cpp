#include "MaterialResponse.h"
#include <numeric>
#include <fstream>

using namespace World3d;

void Response::toStressStrainTensor(const std::array<double, 3>& xi, const std::array<double, 3>& response,
                          std::array<double, 3>& E, std::array<double, 3>& S){
    double e_mk0 = exp(-2 * xi[0]), e_mk1 = exp(-2 * xi[1]), k2 = xi[2];
    double e_k0 = 1.0/e_mk0, e_k1 = 1.0 / e_mk1;
    E[0] = 0.5 * (expm1(2 * xi[0]));
    E[1] = 0.5 * (expm1(2 * xi[1]) + k2 * k2 * e_k0);
    E[2] = 0.5 * k2 * e_k0;
    
    S[0] = e_mk0 * (response[0] - 2 * k2 * response[2]) + e_mk1 * k2 * k2 * response[1];
    S[1] = e_mk1 * response[1];
    S[2] = -e_mk1 * k2 * response[1] + e_mk0 * response[2];
    
}

void Response::toStrain(const std::array<double, 3>& xi, std::array<double, 3>& E){
    double e_mk0 = exp(-2 * xi[0]), e_mk1 = exp(-2 * xi[1]), k2 = xi[2];
    double e_k0 = 1.0/e_mk0, e_k1 = 1.0 / e_mk1;
    E[0] = 0.5 * (expm1(2 * xi[0]));
    E[1] = 0.5 * (expm1(2 * xi[1]) + k2 * k2 * e_k0);
    E[2] = 0.5 * k2 * e_k0;
}

void Response::derivByStrain(const std::array<double, 3>& xi, std::array<std::array<double, 3>, 3>& dxi_dE){
    double e_mk0 = exp(-2 * xi[0]), e_mk1 = exp(-2 * xi[1]), k2 = xi[2];
    dxi_dE[0][0] = e_mk0, dxi_dE[0][1] = dxi_dE[0][2] = 0;
    dxi_dE[1][0] = k2*k2*e_mk1, dxi_dE[1][1] = e_mk1, dxi_dE[1][2] = -k2*e_mk1;
    dxi_dE[2][0] = -2*k2*e_mk0, dxi_dE[2][1] = 0, dxi_dE[2][2] = e_mk0;
}

void Response::secondDerivByStrain(const std::array<double, 3>& xi, std::array<std::array<double, 3>, 3>& dxi_dE, std::array<std::array<double, 6>, 3>& ddxi_dE){
    double e_mk0 = exp(-2 * xi[0]), e_mk1 = exp(-2 * xi[1]), k2 = xi[2];
    dxi_dE[0][0] = e_mk0, dxi_dE[0][1] = dxi_dE[0][2] = 0;
    dxi_dE[1][0] = k2*k2*e_mk1, dxi_dE[1][1] = e_mk1, dxi_dE[1][2] = -k2*e_mk1;
    dxi_dE[2][0] = -2*k2*e_mk0, dxi_dE[2][1] = 0, dxi_dE[2][2] = e_mk0;

    ddxi_dE[0] = std::array<double, 6>{-2*e_mk0*e_mk0, 0, 0, 0, 0, 0};
    ddxi_dE[1] = std::array<double, 6>{
        -2*k2*k2*e_mk1*(2*e_mk0 + k2*k2 * e_mk1),
        -2*e_mk1*e_mk1,
        -2*(k2*e_mk1)*(k2*e_mk1),
        -2*(k2*e_mk1)*(k2*e_mk1),
         2*k2*e_mk1*(e_mk0 + k2*k2*e_mk1),
         2*k2*e_mk1*e_mk1
    };
    ddxi_dE[2] = std::array<double, 6>{
        8*k2*e_mk0*e_mk0, 0, 0, 0, -2*e_mk0*e_mk0, 0
    };
}

void Response::fromStressStrainTensor(std::array<double, 3>& xi, std::array<double, 3>& response,
                            const std::array<double, 3>& E, const std::array<double, 3>& S){
    double a = E[2]*E[2] / (0.5 + E[0]);
    xi[0] = 0.5*log1p( 2*E[0] );
    xi[1] = 0.5*log1p( 2*(E[1] - a) );
    xi[2] = E[2] / (0.5 + E[0]);

    response[0] = S[0] + 2 * ( E[0]*S[0] + a*S[1] + 2*E[2]*S[2] );
    response[1] = S[1] * ( 1 + 2 * (E[1] - a) );
    response[2] = S[2] + 2*( S[2]*E[0] + S[1]*E[2] );
}

void Response::fromStrain(std::array<double, 3>& xi, const std::array<double, 3>& E){
    double a = E[2]*E[2] / (0.5 + E[0]);
    xi[0] = 0.5*log1p( 2*E[0] );
    xi[1] = 0.5*log1p( 2*(E[1] - a) );
    xi[2] = E[2] / (0.5 + E[0]);
}

void ResponseTable::unifyEqualElems(double atol, double rtol) {
    std::array<double, 3> tol = {atol, atol, atol};
    std::for_each(begin(), end(), [&tol, rtol](auto& r){ 
        for (int k = 0; k < 3; ++k) 
            tol[k] = std::max(tol[k], rtol * std::abs(r.xi[k])); 
    });
    auto resp_less = [tol](auto& a, auto& b){
        for (int k = 0; k < 3; ++k)
            if (std::abs(a.xi[k] - b.xi[k]) > tol[k])
                return a.xi[k] < b.xi[k];
        return false;
    };
    auto resp_eq = [tol](auto& a, auto& b){
        for (int k = 0; k < 3; ++k)
            if (std::abs(a.xi[k] - b.xi[k]) > tol[k])
                return false;
        return true;    
    };
    std::stable_sort(begin(), end(), resp_less);
    erase( std::unique( begin(), end(), resp_eq ), end() ); 
}

ResponseTable ResponseTable::concat(const ResponseTable* first, const ResponseTable* last){
    std::size_t sz = std::accumulate(first, last, std::size_t(0), [](std::size_t sum, const auto& r)->std::size_t{ return sum + r.size(); });
    ResponseTable r; r.reserve(sz);
    std::for_each(first, last, [&r](const auto& a) { r.append(a); });
    return r;
}

ResponseTable ResponseTable::concat(const std::shared_ptr<ResponseTable>* first, const std::shared_ptr<ResponseTable>* last){
    std::size_t sz = std::accumulate(first, last, std::size_t(0), [](std::size_t sum, const auto& r)->std::size_t{ return sum + r->size(); });
    ResponseTable r; r.reserve(sz);
    std::for_each(first, last, [&r](const auto& a) { r.append(*a); });
    return r;
}

ResponseTable::bbox ResponseTable::get_bbox() const {
    if (empty()) return {{0, 0, 0}, {0, 0, 0}};
    bbox b{begin()->xi, begin()->xi};
    std::for_each(begin(), end(), [&b](auto f){
        for (int k = 0; k < 3; ++k){
            b.first[k] = std::min(b.first[k], f.xi[k]);
            b.second[k] = std::max(b.second[k], f.xi[k]);
        }
    });
    return b;
}

void RegionalResponseTable::unifyEqualElems(double atol, double rtol){
    std::for_each(m_dat.begin(), m_dat.end(), [atol, rtol](auto& i) { i.unifyEqualElems(atol, rtol); });
}

void RegionalResponseTable::combine(int to) {
    ResponseTable superts = ResponseTable::concat(m_dat.data(), m_dat.data() + m_dat.size());
    if (to < 0) to = m_dat.size();
    if (to > m_dat.size()-1) m_dat.resize(to+1);
    m_common = to;
    m_dat[to] = std::move(superts);
}

void RegionalResponseTable::append(RegionalResponseTable rs){
    if (rs.m_dat.size() > 0) m_common = -1;
    m_dat.reserve(m_dat.size() + rs.m_dat.size());
    for (const auto& i: rs.m_dat) m_dat.emplace_back(std::move(i));
}

void RegionalResponseTable::merge_inplace(const RegionalResponseTable& r){
    m_common = -1;
    m_dat.resize(std::max(m_dat.size(), r.m_dat.size()));
    for (int i = 0, cnt = r.m_dat.size(); i < cnt; ++i) {
        int before = m_dat[i].size();
        m_dat[i].resize(m_dat[i].size() + r.m_dat[i].size());
        copy(r.m_dat[i].begin(), r.m_dat[i].end(), m_dat[i].begin() + before);
    }
}

RegionalResponseTable RegionalResponseTable::merge(const RegionalResponseTable* first, const RegionalResponseTable* last){
    RegionalResponseTable res;
    std::for_each(first, last, [&res](const auto& x) { res.merge_inplace(x); });
    return res;
}

RegionalResponseTable RegionalResponseTable::filter(const RegionalResponseTable& rs, std::function<bool(const std::array<double, 3>&)> f){
    RegionalResponseTable rs1 = rs;
    auto pred = [&f](const auto& e) { return f(e.xi); };
    for (int i = 0; i < rs.m_dat.size(); ++i){
        auto id = copy_if(rs.m_dat[i].begin(), rs.m_dat[i].end(), rs1.m_dat[i].begin(), pred);
        rs1.m_dat[i].resize(id - rs1.m_dat[i].begin());
    }
    return rs1;
}

RegionalResponseTable::bbox RegionalResponseTable::get_bbox() const{
    bool is_empty = true;
    bbox res;
    for (auto r: m_dat) if (!r.empty()){
        if (is_empty) {
            res = r.get_bbox();
            is_empty = false;
        } else {
            auto rb = r.get_bbox();
            for (int k = 0; k < 3; ++k){
                res.first[k] = std::min(res.first[k], rb.first[k]);
                res.second[k] = std::max(res.second[k], rb.second[k]);
            }
        }
    }
    if (is_empty) res = {{0, 0, 0}, {0, 0, 0}};
    return res;
} 

std::ostream& RegionalResponseTable::print_bbox(const bbox& b, std::ostream& out){
    out << "bbox{ " 
        << "P{ " << b.first[0] << ", " << b.first[1] << ", " << b.first[1] << " } -> P{ "
        << b.second[0] << ", " << b.second[1] << ", " << b.second[1] << " } }";
    return out;    
}

std::ostream& RegionalResponseTable::print_xi_cloud(std::ostream& out, int n){
     if (n < 0)
        for (const auto& i: m_dat)
            for (const auto& j: i)
                out << j.xi[0] << ", " << j.xi[1] << ", " << j.xi[2] << "\n";
    else
        for (const auto& j: (m_dat[n]))
            out << j.xi[0] << ", " << j.xi[1] << ", " << j.xi[2] << "\n";
    return out;        
}
std::ostream& RegionalResponseTable::print_xi_response_cloud(std::ostream& out, int n){
    if (n < 0)
        for (const auto& i: m_dat)
            for (const auto& j: i)
                out << j.xi[0] << ", " << j.xi[1] << ", " << j.xi[2] << ", \t"
                    << j.response[0] << ", " << j.response[1] << ", " << j.response[2] << "\n";
    else
        for (const auto& j: (m_dat[n]))
            out << j.xi[0] << ", " << j.xi[1] << ", " << j.xi[2] << ", \t"
                << j.response[0] << ", " << j.response[1] << ", " << j.response[2] << "\n";
    return out;            
}

void RegionalResponseTable::save_ascii(std::string fname){
    std::ofstream f(fname);
    f << "reg, xi_0, xi_1, xi_2, r_0, r_1, r_2\n";
    for (std::size_t i = 0; i < m_dat.size(); ++i)
        for (auto& j: (m_dat[i]))
            f << i << ", " << j.xi[0] << ", " << j.xi[1] << ", " << j.xi[2] << 
                ", " << j.response[0] << ", " << j.response[1] << ", " << j.response[2] << "\n";
    f.close();            
}
RegionalResponseTable RegionalResponseTable::read_ascii(std::string fname){
    RegionalResponseTable t;
    std::ifstream f(fname);
    f.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
    f >> std::ws;
    char pc = f.peek();
    while(!f.eof()){
        int r_id = 0;
        Response r;
        f >> r_id >> std::ws >> pc >> std::ws;
        for (int k = 0; k < 3; ++k)
            f >> r.xi[k] >> std::ws >> pc >> std::ws; 
        for (int k = 0; k < 2; ++k)
            f >> r.response[k] >> std::ws >> pc >> std::ws;
        f >> r.response[2] >> std::ws;
        pc = f.peek();   
        if (t.m_dat.size() <= r_id) t.m_dat.resize(r_id + 1);
        t.m_dat[r_id].emplace_back(std::move(r));
    }

    return t;
}

void RegionalResponseTable::save_binary(std::string fname){
    std::vector<int32_t> reg_vol(m_dat.size());
    for (std::size_t i = 0; i < m_dat.size(); ++i)
        reg_vol[i] = m_dat[i].size() + (i > 0 ? reg_vol[i-1] : 0);
    std::vector<double> raw_buf; raw_buf.reserve(6*reg_vol.back());
    for (const auto& i: m_dat)
        for (const auto& j: i){
            for (int k = 0; k < 3; ++k)
                raw_buf.push_back(j.xi[k]);
            for (int k = 0; k < 3; ++k)
                raw_buf.push_back(j.response[k]);
        }

    std::ofstream f(fname, std::ios::binary);
    int32_t st = 0, mone = -1;
    f.write(reinterpret_cast<const char*>(&st), sizeof(int32_t)*1);
    f.write(reinterpret_cast<const char*>(&mone), sizeof(int32_t)*1);
    f.write(reinterpret_cast<const char*>(raw_buf.data()), sizeof(double)*raw_buf.size());
    f.write(reinterpret_cast<const char*>(reg_vol.data()), sizeof(int32_t)*reg_vol.size());
    int32_t rsz = reg_vol.size();
    f.write(reinterpret_cast<const char*>(&rsz), sizeof(int32_t)*1);
    f.close();
}

RegionalResponseTable RegionalResponseTable::read_binary(std::string fname){
    std::ifstream f(fname, std::ios::binary);
    int32_t st = 0, mone = -1;
    f.read(reinterpret_cast<char*>(&st), sizeof(int32_t)*1);
    f.read(reinterpret_cast<char*>(&mone), sizeof(int32_t)*1);
    bool is_native_endian = (mone == -1);
    auto bflip = [](char* st, char* ed, std::size_t val_size){
        for (char* p = st; p < ed; p += val_size){
            for (std::size_t k = 0; k < val_size / 2; ++k)
                std::swap(p[k], p[val_size - 1 - k]);
        }
    };
    auto cp = f.tellg();
    f.seekg(-sizeof(int32_t), std::ios::end);
    int32_t rsz = 0;
    f.read(reinterpret_cast<char*>(&rsz), sizeof(int32_t)*1);
    if (!is_native_endian) bflip(reinterpret_cast<char*>(&rsz), reinterpret_cast<char*>(&rsz + 1), sizeof(int32_t));
    std::vector<int32_t> reg_vol(rsz);
    f.seekg(-(1 + rsz)*sizeof(int32_t), std::ios::end);
    f.read(reinterpret_cast<char*>(reg_vol.data()), sizeof(int32_t)*rsz);
    if (!is_native_endian) bflip(reinterpret_cast<char*>(reg_vol.data()), reinterpret_cast<char*>(reg_vol.data() + rsz), sizeof(int32_t));
    f.seekg(cp);
    std::vector<double> raw_buf(6*reg_vol.back());
    f.read(reinterpret_cast<char*>(raw_buf.data()), sizeof(double)*raw_buf.size());
    if (!is_native_endian) bflip(reinterpret_cast<char*>(raw_buf.data()), reinterpret_cast<char*>(raw_buf.data() + raw_buf.size()), sizeof(double));
    f.close();

    RegionalResponseTable t;
    t.m_dat.resize(rsz);
    for (int32_t i = 0; i < rsz; ++i){
        t.m_dat[i].reserve(reg_vol[i]);
        for (int32_t k = (i == 0 ? 0 : reg_vol[i-1]); k < reg_vol[i]; ++k){
            const double* p = raw_buf.data() + 6*k;
            Response r;
            for (int l = 0; l < 3; ++l) r.xi[l] = p[l], r.response[l] = p[l+3];
            t.m_dat[i].emplace_back(std::move(r));
        }
    }
    
    return t;
}

void RegionalResponseTable::read(std::string fname){
    std::ifstream f(fname, std::ios::binary);
    if (!f.is_open())
        throw std::runtime_error("Can't open file \"" + fname + "\"");
    int32_t st = 0;
    f.read(reinterpret_cast<char*>(&st), sizeof(int32_t)*1);
    f.close();
    RegionalResponseTable t;
    if (st == 0)
        t = read_binary(fname);
    else
        t = read_ascii(fname);
    merge_inplace(t);        
}

std::function<bool(const std::array<double, 3>&)> World3d::makeRespSphericalConstraint(double R){ return makeSphericalConstraint<std::array<double, 3>, 3>(R); }