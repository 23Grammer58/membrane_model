//
// Created by alex on 01.02.2021.
//

#ifndef AORTIC_VALVE_KNEARESTSEARCHER_H
#define AORTIC_VALVE_KNEARESTSEARCHER_H

#include <type_traits>
#include "InterpolatorBase.h"

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/property_map.h>

///ENTRY should have common info to get PNT and VAL
///PNT should has operator[] with dimension 3
///struct ETraits should have
///  static PNT getPnt(ENTRY d);
///  static VAL getVal(ENTRY d);
///  static void scalVAL(double sc, VAL& vl)
///  static void saxpyVAL(VAL& y, double a, VAL x)
template<typename ENTRY, typename PNT, typename VAL, class ETraits, typename TraitCheck = void>
struct LocalKInterpolationBase{
    // static_assert(false, "Attempt to use not proper ETraits");
};

template<typename ENTRY, typename PNT, typename VAL, class ETraits>
struct LocalKInterpolationBase<ENTRY, PNT, VAL, ETraits,
    std::void_t<
        decltype( std::declval<PNT>() = ETraits::getPnt(std::declval<ENTRY>()) ),
        decltype( std::declval<VAL>() = ETraits::getVal(std::declval<ENTRY>()) ),
        decltype( ETraits::scalVAL(std::declval<double>(), std::declval<VAL&>()) ),
        decltype( ETraits::saxpyVAL(std::declval<VAL&>(), std::declval<double>(), std::declval<VAL>()) )
    >>
{
    typedef CGAL::Simple_cartesian<double> K;
    typedef K::Point_3                                     Point_3;
    typedef std::tuple<Point_3, int>                       Point_and_int;
    typedef CGAL::Search_traits_3<K>                       Traits_base;
    typedef CGAL::Search_traits_adapter<Point_and_int,
    CGAL::Nth_of_tuple_property_map<0, Point_and_int>, Traits_base> Traits;
    typedef CGAL::Orthogonal_k_neighbor_search<Traits>      K_neighbor_search;
    typedef ENTRY                                             iPnt;

    virtual VAL evaluate(const PNT& x, const K_neighbor_search& srch, const std::vector<iPnt>& cloud, int k_neighbour) = 0;
    virtual std::array<VAL,3> gradient(const PNT& x, const K_neighbor_search& srch, const std::vector<iPnt>& cloud, int k_neighbour) = 0;
    virtual std::unique_ptr<LocalKInterpolationBase<ENTRY, PNT, VAL, ETraits>> clone() const = 0;
};

/// Evaluate value by k points as p-inverse distance k-mean
template<typename ENTRY, typename PNT, typename VAL, class ETraits>
struct LocalKInterpInverseDistance: public LocalKInterpolationBase<ENTRY, PNT, VAL, ETraits>{
    using IBASE = LocalKInterpolationBase<ENTRY, PNT, VAL, ETraits>;
    using K_neighbor_search = typename IBASE::K_neighbor_search;
    using iPnt = ENTRY;

    double max_inverse_distance = 1e16;
    double metric_pow = 1.;

    LocalKInterpInverseDistance(double _metric_pow = 1, double max_inv_dist = 1e16): metric_pow{_metric_pow}, max_inverse_distance{max_inv_dist} {}
    VAL evaluate(const PNT& x, const K_neighbor_search& srch, const std::vector<iPnt>& cloud, int k_neighbour) override;
    std::array<VAL,3> gradient(const PNT& x, const K_neighbor_search& srch, const std::vector<iPnt>& cloud, int k_neighbour) override;
    std::unique_ptr<IBASE> clone() const override { return std::make_unique<LocalKInterpInverseDistance<ENTRY, PNT, VAL, ETraits>>(*this); }
};

/// Find K-nearest points (in euclidean norm of points) and call local interpolation
template<typename ENTRY, typename PNT, typename VAL, class DATTraits>
class KNearestSearcher: public Interpolant<std::vector<ENTRY>, PNT, VAL, 3>{
public:
    typedef CGAL::Simple_cartesian<double> K;
    typedef K::Point_3                                     Point_3;
    typedef std::tuple<Point_3, int>                       Point_and_int;
    typedef CGAL::Search_traits_3<K>                       Traits_base;
    typedef CGAL::Search_traits_adapter<Point_and_int,
    CGAL::Nth_of_tuple_property_map<0, Point_and_int>, Traits_base> Traits;
    typedef CGAL::Orthogonal_k_neighbor_search<Traits>      K_neighbor_search;
    typedef K_neighbor_search::Tree                         Tree;
    //typedef K_neighbor_search::Distance                         Distance;
    typedef ENTRY iPnt;
    using Derived = Interpolant<std::vector<ENTRY>, PNT, VAL, 3>;

    std::vector<ENTRY> m_cloud;
    Tree m_tree;
    unsigned k = 1;
    std::unique_ptr<K_neighbor_search> last_search;
    std::unique_ptr<LocalKInterpolationBase<ENTRY, PNT, VAL, DATTraits>> local_interp;

    KNearestSearcher(unsigned k = 1): k{k}, local_interp{std::make_unique<LocalKInterpInverseDistance<ENTRY, PNT, VAL, DATTraits>>()} { }
    KNearestSearcher(const KNearestSearcher<ENTRY, PNT, VAL, DATTraits>& a): k{a.k}{ fit(a.m_cloud); local_interp = a.local_interp->clone(); }
    KNearestSearcher(KNearestSearcher<ENTRY, PNT, VAL, DATTraits>&& a) = default;
    bool fit(const std::vector<ENTRY>& data) override;
    VAL evaluate(const PNT& x) override;
    std::array<VAL, 3> gradient(const PNT& x) override;
    /// Get last found nearest point
    PNT lastNearestX() const;
    /// Set interpolator in k-points
    KNearestSearcher& setLocalInterpolator(std::unique_ptr<LocalKInterpolationBase<ENTRY, PNT, VAL, DATTraits>> f){ return local_interp = std::move(f), *this; }
    /// Get last performed k-nearest query result
    const K_neighbor_search& lastSearch() const { return *last_search; }
    void info() { std::cout << "m_tree.size() = " << m_tree.size() << std::endl; }
    std::string getType() const{ return "KNearestSearcher"; }
    std::unique_ptr<Derived> clone() const override{ return std::make_unique<KNearestSearcher<ENTRY, PNT, VAL, DATTraits>>(*this); }
};

#include "KNearestInterpolant.inl"

#endif //AORTIC_VALVE_KNEARESTSEARCHER_H