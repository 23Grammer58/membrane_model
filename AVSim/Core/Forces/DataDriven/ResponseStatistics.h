//
// Created by alex on 29.01.2021.
//

#ifndef AORTIC_VALVE_RESPONSESTATISTICS_H
#define AORTIC_VALVE_RESPONSESTATISTICS_H
#include <vector>
#include <array>
#include <cfloat>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/array.hpp>

#include "InterpolatorBase.h"
#include "KNearestSearcher.h"


struct ResponseStatistics{
    friend class boost::serialization::access;
    struct Elem{
        friend class boost::serialization::access;
        std::array<double, 3> ksi;      //meassure of deformation
        std::array<double, 3> response; //d(psi)/d(ksi_k) where psi is energy potential (response-function)

        struct Traits{
            using PNT = std::array<double, 3>;
            using VAL = std::array<double, 3>;
            static PNT getPnt(const Elem& d) { return d.ksi; }
            static VAL getVal(const Elem& d) { return d.response; }
            static void scalVAL(double a, VAL& x) { for (auto& i: x) i *= a; }
            static void saxpyVAL(VAL& y, double a, VAL x){ for (int i = 0; i < 3; ++i) y[i] += a*x[i]; }
        };
    private:
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version){ ar & ksi & response; }
    };
    struct ResponseData: public std::vector<Elem>{
        friend class boost::serialization::access;

        void unifyEqualElems(double tol = DBL_EPSILON);
        void concat(const std::vector<Elem>& add);
    private:
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version){ ar & *reinterpret_cast<std::vector<Elem>*>(this); }
    };
    using Rdat = ResponseData;
    using Stat = std::vector<ResponseData>;

    Stat stat;

    ResponseStatistics() = default;
    ResponseStatistics(const ResponseStatistics& ) = default;
    ResponseStatistics(ResponseStatistics&& ) noexcept = default;
    ResponseStatistics& operator=(const ResponseStatistics&) = default;
    ResponseStatistics& operator=(ResponseStatistics&&) noexcept = default;
    ResponseStatistics(std::string fname, int version = 0){ load_from_file(fname); }

    void unifyEqualsPoints(double tol = DBL_EPSILON){ std::for_each(stat.begin(), stat.end(), [tol](Rdat& i) { i.unifyEqualElems(tol); }); }
    void combine(int to = -1);
    void qcombine(int to = -1){ combine(to); stat[to < 0 ? stat.size() - 1 : to].unifyEqualElems(); }
    void combineAll() { qcombine(0); stat.resize(1); }
    void concat(const Rdat& dat){ common = -1; stat.push_back(dat); }
    void concat(const ResponseStatistics& rs);
    void merge(const ResponseStatistics& t);
    void merge(std::string fname){ merge(ResponseStatistics(fname)); }

    static ResponseStatistics filter(const ResponseStatistics rs, std::function<bool(std::array<double, 3> ksi)> f);
    Rdat get_common_cloud();

    void load_from_file(std::string fname, int version = 0);
    void save_to_file(std::string fname);
    void save_ksi_distribution(std::string fname, int n = -1);
    void save_distribution(std::string fname, int n = -1);
    void print_data_bbox();


private:
    int common = -1;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version){
        std::string mark = "TSS";
        ar & mark;
        if (mark != "TSS") throw std::runtime_error("Load error!!! File damaged\n");
        ar & stat;
    }
};

using InterpolantBase3D = InterpolantBase<
        ResponseStatistics::Elem,
        ResponseStatistics::Elem::Traits::PNT,
        ResponseStatistics::Elem::Traits::VAL
        >;
using KNearestSearcher3D = KNearestSearcher<
        ResponseStatistics::Elem,
        ResponseStatistics::Elem::Traits::PNT,
        ResponseStatistics::Elem::Traits::VAL,
        ResponseStatistics::Elem::Traits
        >;
using InterpolantFactory3D = InterpolantFactory<
        ResponseStatistics::Elem,
        ResponseStatistics::Elem::Traits::PNT,
        ResponseStatistics::Elem::Traits::VAL
        >;

template <typename T, typename ...Ts>
InterpolantFactory3D makeInterpolantFactory3D(T interp, Ts ... ts){
    using DAT = ResponseStatistics::Elem;
    using TR = ResponseStatistics::Elem::Traits;
    return makeInterpolantFactory<TR::PNT, TR::VAL, DAT>(interp, ts...);
}


#endif //AORTIC_VALVE_RESPONSESTATISTICS_H
