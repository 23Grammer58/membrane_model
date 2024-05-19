//
// Created by alex on 29.01.2021.
//

#ifndef AORTIC_VALVE_RESPONSESTATISTICS_H
#define AORTIC_VALVE_RESPONSESTATISTICS_H

#include <array>
#include <vector>
#include <limits>
#include <algorithm>
#include <cmath>
#include <initializer_list>
#include <iostream>

#include "InterpolatorBase.h"
#include "KNearestInterpolant.h"

namespace World3d{

struct Response{
    std::array<double, 3> xi;      ///< meassure of deformation
    std::array<double, 3> response; ///< d(psi)/d(xi_k) where psi is energy potential (response-function) without thickness

    struct Traits{
        using PNT = std::array<double, 3>;
        using VAL = std::array<double, 3>;
        static PNT getPnt(const Response& d) { return d.xi; }
        static VAL getVal(const Response& d) { return d.response; }
        static void scalVAL(double a, VAL& x) { for (auto& i: x) i *= a; }
        static void saxpyVAL(VAL& y, double a, VAL x){ for (int i = 0; i < 3; ++i) y[i] += a*x[i]; }
    };

    /// Convert response to Lagrange Green strain E and Second Piola–Kirchhoff stress S
    /// xi = (xi_0, xi_1, xi_2), E = (E00, E11, E01), r = (d psi / dxi_0, d psi / dxi_1, d psi / dxi_2), S = (S00, S11, S01) 
    static void toStressStrainTensor(const std::array<double, 3>& xi, const std::array<double, 3>& response,
                                     std::array<double, 3>& E_strain, std::array<double, 3>& S_stress);
    /// Convert Lagrange Green strain E and Second Piola–Kirchhoff stress S to response
    static void fromStressStrainTensor(std::array<double, 3>& xi, std::array<double, 3>& response,
                                       const std::array<double, 3>& E_strain, const std::array<double, 3>& S_stress); 
    static void toStrain(const std::array<double, 3>& xi, std::array<double, 3>& E_strain);
    static void fromStrain(std::array<double, 3>& xi, const std::array<double, 3>& E_strain);  
    static void derivByStrain(const std::array<double, 3>& xi, std::array<std::array<double, 3>, 3>& dxi_dE);
    static void secondDerivByStrain(const std::array<double, 3>& xi, std::array<std::array<double, 3>, 3>& dxi_dE, std::array<std::array<double, 6>, 3>& ddxi_dE);                                                           
};

struct ResponseTable: public std::vector<Response>{
    using bbox = std::pair<std::array<double, 3>, std::array<double, 3>>;
    /// Remove duplicate points where performed approximate comparison based on absolute and relative tolerances
    /// rel_tol define tolerance relative to size of xi bounding box
    void unifyEqualElems(double abs_tol = std::numeric_limits<double>::epsilon(), double rel_tol = 0);
    /// Remove exactly same points
    void removeDuplicates() { unifyEqualElems(0, 0); }
    void append(const std::vector<Response>& add) { insert(end(), add.begin(), add.end()); }
    template<typename iterator>
    void append(iterator first, iterator last) { insert(end(), first, last); }
    /// Return AABB of xi point cloud
    bbox get_bbox() const ;

    static ResponseTable concat(const ResponseTable* first, const ResponseTable* last);
    static ResponseTable concat(const std::shared_ptr<ResponseTable>* first, const std::shared_ptr<ResponseTable>* last);
    static ResponseTable concat(std::initializer_list<ResponseTable> r) { return concat(r.begin(), r.end()); }
};

struct RegionalResponseTable{
    using bbox = std::pair<std::array<double, 3>, std::array<double, 3>>;

    std::vector<ResponseTable> m_dat;
    
    RegionalResponseTable() = default;

    ResponseTable& region(unsigned reg_id){ return m_dat[reg_id]; }
    const ResponseTable& region(unsigned reg_id) const { return m_dat[reg_id]; }
    const ResponseTable& cregion(unsigned reg_id) const { return m_dat[reg_id]; }
    /// Remove duplicate points based on approximate comparasion
    void unifyEqualElems(double abs_tol = std::numeric_limits<double>::epsilon(), double rel_tol = 0);
    /// Creates a common (including data from all others) table and saves it for the region with the "to" number
    /// if to < 0 will be created new region number, if to > current maximal region number then maximal region number will increased
    /// Last created common table region id is available by get_common_region_id() method
    void combine(int to = -1);
    /// Creates a common table and remove duplicates in it
    void qcombine(int to = -1, double atol = std::numeric_limits<double>::epsilon(), double rtol = 0){ combine(to); m_dat[get_common_region_id()].unifyEqualElems(atol, rtol); }
    /// Unite all table into one table
    void squashAll(double atol = std::numeric_limits<double>::epsilon(), double rtol = 0) { qcombine(0, atol, rtol); m_dat.resize(1); }
    /// Add new region, return region number
    unsigned append(ResponseTable r){ m_dat.emplace_back(std::move(r)); m_common = -1; return m_dat.size() - 1; }
    void append(RegionalResponseTable r);
    /// Merge tables region by region
    void merge_inplace(const RegionalResponseTable& r);

    static RegionalResponseTable merge(const RegionalResponseTable* first, const RegionalResponseTable* last);
    static RegionalResponseTable merge(std::initializer_list<RegionalResponseTable> tabs){ return merge(tabs.begin(), tabs.end()); }
    /// Create RegionalResponseTable where data in tables correspond to filter
    static RegionalResponseTable filter(const RegionalResponseTable& rs, std::function<bool(const std::array<double, 3>& xi)> f);
    static std::function<bool(const std::array<double, 3>& xi)> makeSphericalFilter(double R){ return makeSphericalConstraint<std::array<double, 3>, 3>(R); }


    int get_common_region_id() const { return m_common; }
    /// Return AABB of xi point cloud
    bbox get_bbox() const ;
    std::ostream& print_cloud_sparsity_info(std::ostream& out = std::cout) const { return print_bbox(get_bbox(), out); } 
    std::ostream& print_xi_cloud(std::ostream& out = std::cout, int reg_id = -1);
    std::ostream& print_xi_response_cloud(std::ostream& out = std::cout, int reg_id = -1);
    void save(std::string fname, bool binary = true){ if (binary) save_binary(fname); else save_ascii(fname); }
    void read(std::string fname);
    void clear() { m_dat.clear(); m_common = -1; }

    static std::ostream& print_bbox(const bbox& b, std::ostream& out = std::cout);
private:
    int m_common = -1;   

    void save_binary(std::string fname);
    void save_ascii(std::string fname);
    static RegionalResponseTable read_binary(std::string fname);
    static RegionalResponseTable read_ascii(std::string fname);
};

using ResponseInterpolant = Interpolant< std::vector<Response>, Response::Traits::PNT, Response::Traits::VAL, 3>;
using ResponseKNearest = KNearestSearcher< Response, Response::Traits::PNT, Response::Traits::VAL, Response::Traits >;
using ResponseLocKInvDist = LocalKInterpInverseDistance< Response, Response::Traits::PNT, Response::Traits::VAL, Response::Traits >;
using ResponseInterpFactory = MultiInterpolant< std::vector<Response>, Response::Traits::PNT, Response::Traits::VAL, 3>;

template <typename T, typename ...Ts>
ResponseInterpFactory makeMultiResponseInterpolant(T interp, Ts ... ts){
    return makeMultiInterpolant<Response::Traits::PNT, Response::Traits::VAL, std::vector<Response>>(interp, ts...);
}
template <typename T>
ResponseInterpFactory makeRespConstraintShell(T interp, std::function<bool(const std::array<double, 3>&)> constr){
    return makeConstraintShell<Response::Traits::PNT, Response::Traits::VAL, std::vector<Response>>(interp, constr);
}
std::function<bool(const std::array<double, 3>&)> makeRespSphericalConstraint(double R);

}

#endif //AORTIC_VALVE_RESPONSESTATISTICS_H