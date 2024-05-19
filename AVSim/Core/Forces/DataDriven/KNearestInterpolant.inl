
#include "KNearestInterpolant.h"

#ifndef AORTIC_VALVE_KNEARESTSEARCHER_INL
#define AORTIC_VALVE_KNEARESTSEARCHER_INL

template<typename ENTRY, typename PNT, typename VAL, class ETraits>
VAL LocalKInterpInverseDistance<ENTRY, PNT, VAL, ETraits>::evaluate(const PNT& x, const K_neighbor_search& srch, const std::vector<iPnt>& cloud, int k_neighbour){
    (void) x;
    VAL res = ETraits::getVal(cloud[std::get<1>(srch.begin()->first)]);
    if (std::next(srch.begin()) == srch.end()) 
        return res;
    auto metric = [metric_pow = this->metric_pow](double sqr_len) -> double{ return pow(sqr_len, metric_pow/2); };
    const double drp_tol = max_inverse_distance;
    auto weight = [drp_tol](double metr) -> double{ return (metr > 1.0/drp_tol) ? 1.0 / metr : drp_tol; };
    double sum = 0;
    std::vector<double> wws; wws.reserve(k_neighbour);
    for (const auto &i: srch) {
        wws.push_back(weight(metric(i.second)));
        sum += wws.back();
    }
    if (sum < drp_tol) {
        int off = 0;
        {
            for (auto& i: wws) i /= sum;
            ETraits::scalVAL(wws[off], res);
            ++off;
        }
        for (auto i = srch.begin() + 1; i != srch.end(); ++i) {
            ETraits::saxpyVAL(res, wws[off], ETraits::getVal(cloud[std::get<1>(i->first)]));
            ++off;
        }
    }
    return res;
}

template<typename ENTRY, typename PNT, typename VAL, class ETraits>
std::array<VAL,3> LocalKInterpInverseDistance<ENTRY, PNT, VAL, ETraits>::gradient(const PNT& x, const K_neighbor_search& srch, const std::vector<iPnt>& cloud, int k_neighbour) {
    VAL zero;
    ETraits::scalVAL(0, zero);
    std::array<VAL,3> res{zero, zero, zero};
    if (std::next(srch.begin()) == srch.end()) 
        return res;
    auto metric = [metric_pow = this->metric_pow](double sqr_len) -> double{ return metric_pow == 2 ? sqr_len : pow(sqr_len, metric_pow/2); };
    const double drp_tol = max_inverse_distance;
    auto weight = [drp_tol](double metr) -> double{ return (metr > 1.0/drp_tol) ? 1.0 / metr : drp_tol; };
    auto dweight = [drp_tol, metric, metric_pow = this->metric_pow](double sqr_len) -> std::pair<double, bool> {
        double nrm = metric(sqr_len);
        if (nrm <= 1.0/drp_tol) return std::pair<double, bool>{0.0, false};
        return std::pair<double, bool>{-metric_pow / (nrm * sqr_len), true}; //-metric_pow * pow(sqr_len, -(metric_pow+2)/2)
    };

    double sum = 0;
    std::vector<double> wws; wws.reserve(k_neighbour);
    for (const auto &i: srch) {
        wws.push_back(weight(metric(i.second)));
        sum += wws.back();
    }

    VAL evl = ETraits::getVal(cloud[std::get<1>(srch.begin()->first)]);
    if (sum < drp_tol) {
        int off = 0;
        {
            for (auto& i: wws) i /= sum;
            ETraits::scalVAL(wws[off], evl);
            ++off;
        }
        for (auto i = srch.begin() + 1; i != srch.end(); ++i) {
            ETraits::saxpyVAL(evl, wws[off], ETraits::getVal(cloud[std::get<1>(i->first)]));
            ++off;
        }
    }

    for (const auto &i: srch) {
        double sqr_len = i.second;
        auto dwww = dweight(sqr_len);
        if (!dwww.second) continue;
        double w = dwww.first / sum;
        PNT xi = ETraits::getPnt(cloud[std::get<1>(i.first)]);
        std::array<double, 3> dx = {x[0] - xi[0], x[1] - xi[1], x[2] - xi[2]};
        VAL dv = ETraits::getVal(cloud[std::get<1>(i.first)]);
        ETraits::saxpyVAL(dv, -1, evl); // v_i - v_mid
        for (int k = 0; k < 3; ++k)
            ETraits::saxpyVAL(res[k], w*dx[k], dv);
    }

    return res;
}

template<typename ENTRY, typename PNT, typename VAL, class DATTraits>
bool KNearestSearcher<ENTRY, PNT, VAL, DATTraits>::fit(const std::vector<ENTRY>& data){
    m_cloud = data;
    std::vector<std::tuple<Point_3, int>> pp;
    for (int i = 0; i < m_cloud.size(); ++i){
        auto p = DATTraits::getPnt(m_cloud[i]);
        pp.emplace_back(Point_3{p[0], p[1], p[2]}, i);
    }
    m_tree.insert(pp.begin(), pp.end());
    return true;
}
template<typename ENTRY, typename PNT, typename VAL, class DATTraits>
VAL KNearestSearcher<ENTRY, PNT, VAL, DATTraits>::evaluate(const PNT& x){
    last_search = std::make_unique<K_neighbor_search>(m_tree, Point_3{x[0], x[1], x[2]}, k);
    auto& srch = *last_search;
    if (srch.begin() == srch.end()) {
        std::cout << "m_tree.size() = " << m_tree.size() << std::endl;
        throw std::runtime_error("CGAL broken");
    }
    if (k <= 1)
        return DATTraits::getVal(m_cloud[std::get<1>(srch.begin()->first)]);
    else return local_interp->evaluate(x, srch, m_cloud, k);
}
template<typename ENTRY, typename PNT, typename VAL, class DATTraits>
std::array<VAL, 3> KNearestSearcher<ENTRY, PNT, VAL, DATTraits>::gradient(const PNT& x){
    if (k <= 1){
        VAL zero;
        DATTraits::scalVAL(0, zero);
        return std::array<VAL, 3>{zero, zero, zero};
    }

    last_search = std::make_unique<K_neighbor_search>(m_tree, Point_3{x[0], x[1], x[2]}, k);
    auto& srch = *last_search;
    if (srch.begin() == srch.end()) {
        std::cout << "m_tree.size() = " << m_tree.size() << std::endl;
        throw std::runtime_error("CGAL broken");
    } 
    return local_interp->gradient(x, srch, m_cloud, k);   
}
/// Get last found nearest point
template<typename ENTRY, typename PNT, typename VAL, class DATTraits>
PNT KNearestSearcher<ENTRY, PNT, VAL, DATTraits>::lastNearestX() const { 
    PNT last_x;
    auto pp = std::get<0>(last_search->begin()->first);
    for (int i = 0; i < 3; ++i) last_x[i] = pp[i];
    return last_x;
}

#endif //AORTIC_VALVE_KNEARESTSEARCHER_INL