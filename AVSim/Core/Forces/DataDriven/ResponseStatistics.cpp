//
// Created by alex on 29.01.2021.
//

#include "ResponseStatistics.h"
#include <gsl/gsl_multifit.h>

using namespace std;
inline std::ostream& operator<<(std::ostream& out, const std::array<double, 3>& ar){
    for (int i = 0; i < ar.size()-1; ++i)
        out << ar[i] << ", ";
    out << ar[ar.size()-1];
    return out;
}

void ResponseStatistics::ResponseData::unifyEqualElems(double tol) {
    std::sort(begin(), end(), [tol](auto a, auto b) {
        for (int i = 0; i < 3; ++i)
            if (fabs(a.ksi[i] - b.ksi[i]) > tol)
                return a.ksi[i] < b.ksi[i];
        return false;
    });
    int del = 0;
    for (int i = 1; i < size(); ++i){
        bool same = true;
        for (int j = 0; j < 3; ++j){
            if (fabs(operator[](i-del-1).ksi[j] - operator[](i).ksi[j]) > tol) same = false;
        }
        if (!same) operator[](i - del) = operator[](i);
        else del++;
    }
    resize(size() - del);
}

void ResponseStatistics::ResponseData::concat(const std::vector<Elem>& add) {
    reserve(size() + add.size());
    for (auto i: add) push_back(std::move(i));
}

void ResponseStatistics::load_from_file(std::string fname, int version) {
    std::ifstream ff(fname);
    boost::archive::binary_iarchive ib(ff);
    serialize(ib, version);
}

void ResponseStatistics::combine(int to) {
    if (to < 0) to = stat.size();
    int sz = 0;
    for (const auto& i: stat)
        sz += i.size();
    Rdat superts;
    superts.reserve(sz);
    for (const auto& i: stat)
        for (auto& j: i)
            superts.push_back(j);
    if (to > stat.size()-1) stat.resize(to+1);
    common = to;
    stat[to] = std::move(superts);
}

ResponseStatistics::Rdat ResponseStatistics::get_common_cloud() {
    if (common >= 0) return stat[common];
    int sz = 0;
    for (const auto& i: stat) sz += i.size();
    Rdat superts; superts.reserve(sz);
    for (const auto& i: stat)
        for (auto& j: i)
            superts.push_back(j);
    return superts;
}

void ResponseStatistics::save_to_file(std::string fname) {
    std::ofstream ff(fname);
    boost::archive::binary_oarchive ob(ff);
    serialize(ob, 0);
}

void ResponseStatistics::concat(const ResponseStatistics &rs) {
    if (rs.stat.size() > 0) common = -1;
    stat.reserve(stat.size() + rs.stat.size());
    for (const auto& i: rs.stat) stat.push_back(i);
}

void ResponseStatistics::merge(const ResponseStatistics &t) {
    common = -1;
    stat.resize(max(stat.size(), t.stat.size()));
    for (int i = 0, cnt = min(stat.size(), t.stat.size()); i < cnt; ++i) {
        int before = stat[i].size();
        stat[i].resize(stat[i].size() + t.stat[i].size());
        copy(t.stat[i].begin(), t.stat[i].end(), stat[i].begin() + before);
    }
}

ResponseStatistics ResponseStatistics::filter(const ResponseStatistics rs, std::function<bool(std::array<double, 3>)> f) {
    ResponseStatistics rs1;
    auto pred = [&f](const auto& e) { return f(e.ksi); };
    rs1.stat.resize(rs.stat.size());
    for (int i = 0; i < rs.stat.size(); ++i){
        rs1.stat[i].resize(rs.stat[i].size());
        auto id = copy_if(rs.stat[i].begin(), rs.stat[i].end(), rs1.stat[i].begin(), pred);
        rs1.stat[i].resize(id - rs1.stat[i].begin());
    }
    return rs1;
}

void ResponseStatistics::save_ksi_distribution(std::string fname, int n) {
    ofstream f(fname);
    if (n < 0)
        for (const auto& i: stat)
            for (const auto& j: i)
                f << j.ksi << "\n";
    else
        for (const auto& j: stat[n])
            f << j.ksi << "\n";
}

void ResponseStatistics::save_distribution(std::string fname, int n){
    ofstream f(fname);
    if (n < 0)
        for (const auto& i: stat)
            for (const auto& j: i)
                f << j.ksi << ", " << j.response << "\n";
    else
        for (const auto& j: stat[n])
            f << j.ksi << ", " << j.response << "\n";
}

void ResponseStatistics::print_data_bbox() {
    array<pair<double, double>, 3> bbox;
    bbox[0] = pair<double, double>{stat[0][0].ksi[0], stat[0][0].ksi[0]};
    bbox[1] = pair<double, double>{stat[0][0].ksi[1], stat[0][0].ksi[1]};
    bbox[2] = pair<double, double>{stat[0][0].ksi[2], stat[0][0].ksi[2]};
    for (const auto& i: stat)
        for (const auto& j: i)
            for (int k = 0; k < 3; ++k)
                bbox[k].first = std::min(bbox[k].first, j.ksi[k]), bbox[k].second = std::max(bbox[k].second, j.ksi[k]);
    cout << "range: " << endl;
    char crd[3] = {'x', 'y', 'z'};
    for (int k = 0; k < 3; ++k)
        cout << " " << crd[k] << ": "<< bbox[k].first << " " << bbox[k].second << endl;
}
