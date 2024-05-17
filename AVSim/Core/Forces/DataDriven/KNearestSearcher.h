//
// Created by alex on 01.02.2021.
//

#ifndef AORTIC_VALVE_KNEARESTSEARCHER_H
#define AORTIC_VALVE_KNEARESTSEARCHER_H

#include "InterpolatorBase.h"

//DAT should have common info to get PNT and VAL
//PNT should has operator[] with dimension 3
//struct DATTraits should have
//  static PNT getPnt(DAT d);
//  static VAL getVal(DAT d);
//  static void scalVAL(double sc, VAL& vl)
//  static void saxpyVAL(VAL& y, double a, VAL x)
template<typename DAT, typename PNT, typename VAL, class DATTraits>
class KNearestSearcher: public InterpolantBase<DAT, PNT, VAL>{
public:
    typedef CGAL::Simple_cartesian<double> K;
    typedef K::Point_3                                     Point_3;
    typedef std::tuple<Point_3,int>                        Point_and_int;
    typedef CGAL::Search_traits_3<K>                       Traits_base;
    typedef CGAL::Search_traits_adapter<Point_and_int,
    CGAL::Nth_of_tuple_property_map<0, Point_and_int>,
    Traits_base>                                              Traits;
    typedef CGAL::Orthogonal_k_neighbor_search<Traits>          K_neighbor_search;
    typedef K_neighbor_search::Tree                             Tree;
    //typedef K_neighbor_search::Distance                         Distance;
    typedef DAT iPnt;
    using Derived = InterpolantBase<DAT, PNT, VAL>;

    Tree m_tree;
    unsigned k = 1;
    PNT last_x;

    std::function<VAL(const K_neighbor_search& srch, const std::vector<iPnt>& cloud, int k)> local_interp =
            [](const K_neighbor_search& srch, const std::vector<iPnt>& cloud, int k){
                auto metric = [](double len) -> double{ return std::sqrt(len); };
                const double drp_tol = 1.0e16;
                auto weight = [&drp_tol](double metr) -> double{ return (metr > 1.0/drp_tol) ? 1.0 / metr : drp_tol; };
                double sum = 0;
                std::vector<double> wws; wws.reserve(k);
                for (const auto &i: srch) {
                    wws.push_back(weight(metric(i.second)));
                    sum += wws.back();
                }
                VAL res = DATTraits::getVal(cloud[std::get<1>(srch.begin()->first)]);
                if (sum < drp_tol) {
                    int off = 0;
                    {
                        for (auto& i: wws) i /= sum;
                        DATTraits::scalVAL(wws[off], res);
                        ++off;
                    }
                    for (auto i = srch.begin() + 1; i != srch.end(); ++i) {
                        DATTraits::saxpyVAL(res, wws[off], DATTraits::getVal(cloud[std::get<1>(i->first)]));
                        ++off;
                    }
                }
                return res;
            };

    explicit KNearestSearcher(unsigned k = 1): k{k} { }
    KNearestSearcher(const KNearestSearcher& other): k{other.k} { if (other.m_cloud) init(other.m_cloud); };
    KNearestSearcher(KNearestSearcher&&) noexcept = default;
    KNearestSearcher& operator=(const KNearestSearcher& other) {
        k = other.k;
        if (other.m_cloud) init(other.m_cloud);
    };
    KNearestSearcher& operator=(KNearestSearcher&&) noexcept = default;

    void init(std::vector<DAT>* cloud) override {
        std::vector<std::tuple<Point_3, int>> pp;
        for (int i = 0; i < cloud->size(); ++i){
            auto p = DATTraits::getPnt(cloud->operator[](i));
            pp.emplace_back(Point_3{p[0], p[1], p[2]}, i);
        }
        m_tree.insert(pp.begin(), pp.end());
        Derived::init(cloud);
    }
    KNearestSearcher& setLocalInterpolator(std::function<VAL(const K_neighbor_search&, const std::vector<iPnt>&, int)> f){
        local_interp = f;
        return *this;
    }
    PNT lastX(){ return last_x; }
    VAL operator()(const PNT& x) final{
        K_neighbor_search srch(m_tree, {x[0], x[1], x[2]}, k);
        if (srch.begin() == srch.end()) {
            std::cout << "m_tree.size() = " << m_tree.size() << std::endl;
            throw std::runtime_error("CGAL broken");
        }
        auto pp = std::get<0>(srch.begin()->first);
        for (int i = 0; i < 3; ++i) last_x[i] = pp[i];
        if (k <= 1)
            return DATTraits::getVal((*Derived::m_cloud)[std::get<1>(srch.begin()->first)]);
        else return local_interp(srch, *Derived::m_cloud, k);
    }
    void info() {
        std::cout << "m_tree.size() = " << m_tree.size() << std::endl;
    }
    std::string getType() const override { return "KNearestSearcher"; }
    bool interpolated() final { return Derived::m_cloud->size() > 0; }
    Derived* clone() final{
        auto p = new KNearestSearcher(k);
        if (Derived::m_cloud) p->init(Derived::m_cloud);
        return p;
    }
};



#endif //AORTIC_VALVE_KNEARESTSEARCHER_H
