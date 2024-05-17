//
// Created by alex on 01.02.2021.
//

#ifndef AORTIC_VALVE_INTERPOLATORBASE_H
#define AORTIC_VALVE_INTERPOLATORBASE_H

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/property_map.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Delaunay_triangulation_cell_base_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>

#include <vector>
#include <tuple>
#include <functional>

template<typename DAT, typename PNT, typename VAL>
class InterpolantBase{
public:
    std::vector<DAT>* m_cloud = nullptr;

    virtual VAL operator()(const PNT& x) { return VAL(); }
    virtual bool interpolated() { return false; }
    virtual InterpolantBase* clone() = 0;
    virtual void init(std::vector<DAT>* cloud) { m_cloud = cloud; }
    virtual std::string getType() const { return "base"; }
};


template<typename DAT, typename PNT, typename VAL>
class InterpolantFactory: public InterpolantBase<DAT, PNT, VAL>{
public:
    using IBASE = InterpolantBase<DAT, PNT, VAL>;
    using Interpolant = std::shared_ptr<IBASE>;

    std::vector<Interpolant> factory;
    bool m_interpolated = false;
    std::function<VAL(InterpolantFactory&, const PNT& )> strategy =
            [](InterpolantFactory& fact, const PNT& x){
                fact.m_interpolated = false;
                VAL res{};
                for (int i = 0; i < fact.factory.size() && !fact.m_interpolated; ++i){
                    res = fact.factory[i].get()->operator()(x);
                    fact.m_interpolated = fact.factory[i].get()->interpolated();
                }
                return res;
            };

    void init(std::vector<DAT>* cloud) override {
        for (auto& i: factory) i.get()->init(cloud);
        InterpolantBase<DAT, PNT, VAL>::init(cloud);
    }
    InterpolantFactory& setStrategy(std::function<VAL(InterpolantFactory&, const PNT& )> newStrat) {
        strategy = std::move(newStrat);
        return *this;
    }
    InterpolantFactory& append(IBASE* inter) {
        factory.emplace_back(inter->clone());
        if (IBASE::m_cloud) factory.back().get()->init(IBASE::m_cloud);
        return *this;
    }
    InterpolantFactory& append(std::shared_ptr<IBASE> inter) {
        factory.push_back(inter);
        if (IBASE::m_cloud) factory.back().get()->init(IBASE::m_cloud);
        return *this;
    }
    template <typename INTERP>
    InterpolantFactory& append(INTERP& inter) {
        factory.emplace_back(inter.clone());
        if (IBASE::m_cloud) factory.back().get()->init(IBASE::m_cloud);
        return *this;
    }

    std::string getType() const override{
        std::string type = "Factory{ ";
        for (const auto& i: factory) type += i.get()->getType() + " -> ";
        type += "}";
        return type;
    }
    bool interpolated() override { return m_interpolated; }
    IBASE* clone() override { return new InterpolantFactory(*this); }
    VAL operator()(const PNT& x)  override  { return strategy(*this, x); }
};

template<typename DAT, typename PNT, typename VAL>
class InterpolantFactoryOwn: public InterpolantFactory<DAT, PNT, VAL>{
public:
    using IBASE = InterpolantBase<DAT, PNT, VAL>;
    using Interpolant = std::shared_ptr<IBASE>;

    std::shared_ptr<std::vector<DAT>> m_cloud;

    InterpolantFactoryOwn& setCloud(std::vector<DAT> dat){
        m_cloud = std::make_shared<std::vector<DAT>>(dat);
        IBASE::m_cloud = m_cloud.get();
        return *this;
    }
    IBASE* clone() override { return new InterpolantFactoryOwn(*this); }
};


template <typename DAT, typename PNT, typename VAL, typename T, typename ...Ts>
void InterpolantFactoryAppender(InterpolantFactoryOwn<DAT, PNT, VAL>& factory, T interp, Ts ... ts){
    factory.append(interp);
    if constexpr (sizeof...(ts) > 0) { InterpolantFactoryAppender(factory, ts...); }
}

template <typename DAT, typename PNT, typename VAL, typename T, typename ...Ts>
void InterpolantFactoryAppender(InterpolantFactory<DAT, PNT, VAL>& factory, T interp, Ts ... ts){
    factory.append(interp);
    if constexpr (sizeof...(ts) > 0) { InterpolantFactoryAppender(factory, ts...); }
}

template <typename PNT, typename VAL, typename DAT, typename T, typename ...Ts>
InterpolantFactory<DAT, PNT, VAL> makeInterpolantFactoryOwn(std::vector<DAT> dat, T interp, Ts ... ts){
    InterpolantFactoryOwn<DAT, PNT, VAL> factory;
    factory.setCloud(dat).append(interp);
    if constexpr (sizeof...(ts) > 0) { InterpolantFactoryAppender(factory, ts...); }
    return factory;
}

template <typename PNT, typename VAL, typename DAT, typename T, typename ...Ts>
InterpolantFactory<DAT, PNT, VAL> makeInterpolantFactory(T interp, Ts ... ts){
    InterpolantFactory<DAT, PNT, VAL> factory;
    factory.append(interp);
    if constexpr (sizeof...(ts) > 0) { InterpolantFactoryAppender(factory, ts...); }
    return factory;
}







#endif //AORTIC_VALVE_INTERPOLATORBASE_H
