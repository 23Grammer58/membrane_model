//
// Created by alex on 01.02.2021.
//

#ifndef AORTIC_VALVE_INTERPOLATORBASE_H
#define AORTIC_VALVE_INTERPOLATORBASE_H

#include <vector>
#include <tuple>
#include <functional>
#include <memory>
#include <cassert>

/// Base class for interpolation of values specified in 3D points
template<typename DAT, typename PNT, typename VAL, int N = 3>
class Interpolant{
public:    
    virtual bool fit(const DAT& data) = 0;
    virtual VAL evaluate(const PNT& x) = 0;
    virtual std::array<VAL, N> gradient(const PNT& x) = 0;
    virtual bool defined_on(const PNT& x) const { return true; }
    virtual std::unique_ptr<Interpolant<DAT, PNT, VAL, N>> clone() const = 0;
};

template<typename DAT, typename PNT, typename VAL, int N = 3>
class MultiInterpolant: public Interpolant<DAT, PNT, VAL, N>{
public:    
    using Interp = Interpolant<DAT, PNT, VAL, N>;
    using Strategy = std::function<int(const std::vector<std::unique_ptr<Interp>>& , const PNT&)>;
    using Constraint = std::function<bool(const MultiInterpolant<DAT, PNT, VAL, N>& , const PNT&)>;

    std::vector<std::unique_ptr<Interp>> m_factory;
    Strategy m_strategy = 
        [](const std::vector<std::unique_ptr<Interp>>& interps, const PNT& x)->int{
            for (int i = 0; i < interps.size(); ++i)
                if (interps[i]->defined_on(x))
                    return i;
            return -1;        
        };
    Constraint m_constraint = [](const MultiInterpolant<DAT, PNT, VAL, N>& back, const PNT& x)->bool { return back.m_strategy(back.m_factory, x) != -1; };    

    MultiInterpolant() = default;
    MultiInterpolant(const MultiInterpolant<DAT, PNT, VAL, N>& a): m_strategy(a.m_strategy), m_constraint{a.m_constraint} {
        m_factory.resize(a.m_factory.size());
        for (std::size_t i = 0; i < a.m_factory.size(); ++i) 
            m_factory[i] = a.m_factory[i]->clone();
    }
    MultiInterpolant(MultiInterpolant<DAT, PNT, VAL, N>&& a) = default;
    bool fit(const DAT& data) override { bool res = true; for(auto& i: m_factory) res &= i->fit(data); return res; }
    /// Set strategy of choose interpolant in current point
    MultiInterpolant& setStrategy(Strategy newStrat) { return m_strategy = std::move(newStrat), *this; }
    MultiInterpolant& setConstraint(Constraint newCstr) { return m_constraint = std::move(newCstr), *this; }
    /// Add new interpolation variant
    MultiInterpolant& append(Interp* inter) { return m_factory.emplace_back(inter->clone()), *this; }
    MultiInterpolant& append(std::unique_ptr<Interp> inter) { return m_factory.emplace_back(std::move(inter)), *this; }
    template <typename INTERP>
    typename std::enable_if<std::is_base_of<Interp, INTERP>::value, MultiInterpolant&>::type
      append(INTERP inter) { return m_factory.emplace_back(std::make_unique<INTERP>(std::move(inter))), *this; }
    VAL evaluate(const PNT& x) override { return m_factory[m_strategy(m_factory, x)]->evaluate(x); }
    std::array<VAL, N> gradient(const PNT& x) override { return m_factory[m_strategy(m_factory, x)]->gradient(x); }
    bool defined_on(const PNT& x) const override { return m_constraint(*this, x); }
    std::unique_ptr<Interp> clone() const override { return std::make_unique<MultiInterpolant<DAT, PNT, VAL, N>>(*this); }
};

namespace internals{
template <typename DAT, typename PNT, typename VAL, int N = 3, typename T, typename ...Ts>
void MultiInterpolantAppender(MultiInterpolant<DAT, PNT, VAL, N>& factory, T interp, Ts ... ts){
    factory.append(interp);
    if constexpr (sizeof...(ts) > 0) { InterpolantFactoryAppender(factory, ts...); }
}
}

template <typename PNT, typename VAL, typename DAT, int N = 3, typename T, typename ...Ts>
MultiInterpolant<DAT, PNT, VAL, N> makeMultiInterpolant(T interp, Ts ... ts){
    MultiInterpolant<DAT, PNT, VAL, N> factory;
    factory.append(interp);
    if constexpr (sizeof...(ts) > 0) { ::internals::MultiInterpolantAppender(factory, ts...); }
    return factory;
}

template <typename PNT, typename VAL, typename DAT, int N = 3, typename T>
MultiInterpolant<DAT, PNT, VAL, N> makeConstraintShell(T interp, std::function<bool(const PNT&)> constr){
    MultiInterpolant<DAT, PNT, VAL, N> factory;
    factory.append(interp);
    factory.setConstraint([constr](const MultiInterpolant<DAT, PNT, VAL, N>&, const PNT& x){ return constr(x); });
    return factory;
}

template <typename PNT, int N>
std::function<bool(const PNT&)> makeSphericalConstraint(double R){
    return [R](const PNT& x){
        double sqr_r = 0;
        for (int i = 0; i < N; ++i)
            sqr_r += x[i]*x[i];
        return sqr_r <= R;    
    };
}

#endif  //AORTIC_VALVE_INTERPOLATORBASE_H