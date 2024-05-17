//
// Created by alex on 30.09.2023.
//

#ifndef AORTIC_VALVE_HYPERELASTICFORMAT_INL
#define AORTIC_VALVE_HYPERELASTICFORMAT_INL

#include "HyperElasticFormat.h"
#include <casadi/casadi.hpp>
#include <cmath>

namespace World3d{

    template<typename T, bool Dummy> 
    struct Energy_Parser::DefaultMathFuncsApplier{
    #define REGISTER_FUNC1(FUNC)\
        template<typename TT, typename = void> struct have_global_##FUNC : std::false_type {};\
        template<typename TT> struct have_global_##FUNC<TT, std::void_t< decltype(TT(::FUNC(std::declval<TT>()))) >> : std::true_type {};\
        \
        template<typename TT, typename = void> struct have_nested_##FUNC : std::false_type {};\
        template<typename TT> struct have_nested_##FUNC<TT, std::void_t< decltype(TT(TT::FUNC(std::declval<TT>()))) >> : std::true_type {};\
        \
        template<typename TT, int num = 0> struct TryGetImpl_##FUNC{\
            static inline TT apply(TT a) { throw std::runtime_error(std::string("Function ") + #FUNC + " is not implemented by default for " + typeid(TT).name() + " type"); return a; }\
        };\
        template<typename TT> struct TryGetImpl_##FUNC<TT, 1>{ static inline TT apply(TT a) { return T::FUNC(a); } };\
        template<typename TT> struct TryGetImpl_##FUNC<TT, 2>{ static inline TT apply(TT a) { return ::FUNC(a); } };\
        template<typename TT> using Try_##FUNC = TryGetImpl_##FUNC<TT, have_nested_##FUNC<TT>::value ? 1 : (have_global_##FUNC<TT>::value ? 2 : 0)>;\
        \
        virtual T FUNC(T a) const { return Try_##FUNC<T>::apply(a); }

    #define REGISTER_FUNC2_0(FUNC, APPL0, APPL1)\
        template<typename TT, typename = void> struct have_global_##FUNC : std::false_type {};\
        template<typename TT> struct have_global_##FUNC<TT, std::void_t< decltype(TT(APPL0(std::declval<TT>(), std::declval<TT>()))) >> : std::true_type {};\
        \
        template<typename TT, typename = void> struct have_nested_##FUNC : std::false_type {};\
        template<typename TT> struct have_nested_##FUNC<TT, std::void_t< decltype( TT( TT:: APPL1(std::declval<TT>(), std::declval<TT>())) ) >> : std::true_type {};\
        \
        template<typename TT, int num = 0> struct TryGetImpl_##FUNC{\
            static inline T apply(TT a, TT b) { throw std::runtime_error(std::string("Function ") + #FUNC + " is not implemented by default for " + typeid(TT).name() + " type"); return a; }\
        };\
        template<typename TT> struct TryGetImpl_##FUNC<TT, 1>{ static inline TT apply(TT a, TT b) { return TT::APPL1(a, b); } };\
        template<typename TT> struct TryGetImpl_##FUNC<TT, 2>{ static inline TT apply(TT a, TT b) { return APPL0(a, b); } };\
        template<typename TT> using Try_##FUNC = TryGetImpl_##FUNC<TT, have_nested_##FUNC<TT>::value ? 1 : (have_global_##FUNC<TT>::value ? 2 : 0)>;\
        \
        virtual T FUNC(T a, T b) const { return Try_##FUNC<T>::apply(a, b); }    

    #define REGISTER_FUNC2_2(FUNC, APPL) REGISTER_FUNC2_0(FUNC, std::APPL<>(), std::template APPL<>())
    #define REGISTER_FUNC2(FUNC, APPL) REGISTER_FUNC2_0(FUNC, ::APPL, APPL)    

        REGISTER_FUNC1(abs);
        REGISTER_FUNC1(exp);
        REGISTER_FUNC1(expm1);
        REGISTER_FUNC1(log);
        REGISTER_FUNC1(log1p);
        REGISTER_FUNC1(sqrt);
        REGISTER_FUNC1(cbrt);
        REGISTER_FUNC1(sin);
        REGISTER_FUNC1(cos);
        REGISTER_FUNC1(tan);
        REGISTER_FUNC1(asin);
        REGISTER_FUNC1(acos);
        REGISTER_FUNC1(atan);
        REGISTER_FUNC1(sinh);
        REGISTER_FUNC1(cosh);
        REGISTER_FUNC1(tanh);
        REGISTER_FUNC1(erf);
        REGISTER_FUNC1(erfc);
        REGISTER_FUNC1(tgamma);
        REGISTER_FUNC1(lgamma);
        REGISTER_FUNC1(floor);
        REGISTER_FUNC1(ceil);       
        template<typename TT, typename Dummy2 = void> struct TryGetImpl_uminus{ static inline TT apply(TT a) { throw std::runtime_error(std::string("Function ") + "unary minus" + " is not implemented by default for " + typeid(TT).name() + " type"); return a; } };
        template<typename TT> struct TryGetImpl_uminus<TT, std::void_t<decltype(-std::declval<TT>())> >{ static inline TT apply(TT a) { return -a; } };
        virtual T uminus(T a) const { return TryGetImpl_uminus<T>::apply(a); }
        template<typename TT, typename Dummy2 = void> struct TryGetImpl_not_op{ static inline TT apply(TT a) { throw std::runtime_error(std::string("Function ") + "not_op" + " is not implemented by default for " + typeid(TT).name() + " type"); return a; } };
        template<typename TT> struct TryGetImpl_not_op<TT, std::void_t<decltype(!std::declval<TT>())> >{ static inline TT apply(TT a) { return !a; } };
        virtual T not_op(T a) const { return TryGetImpl_not_op<T>::apply(a); }

        REGISTER_FUNC2_2(plus, plus);
        REGISTER_FUNC2_2(minus, minus);
        REGISTER_FUNC2_2(mul, multiplies);
        virtual T sq(T a) const { return mul(a, a); }
        virtual T cube(T a) const { return mul(mul(a, a), a); }
        REGISTER_FUNC2_2(div, divides);
        REGISTER_FUNC2(pow, pow);
        REGISTER_FUNC2(atan2, atan2);
        REGISTER_FUNC2(fmod, fmod);
        REGISTER_FUNC2_2(lt, less);
        REGISTER_FUNC2_2(le, less_equal);
        REGISTER_FUNC2_2(eq, equal_to);
        REGISTER_FUNC2_2(ne, not_equal_to);
        REGISTER_FUNC2_2(gt, greater);
        REGISTER_FUNC2_2(ge, greater_equal);
        REGISTER_FUNC2_2(and_op, logical_and);
        REGISTER_FUNC2_2(or_op, logical_or);

        template<typename TT, typename = void> struct have_global_ifelse : std::false_type {};\
        template<typename TT> struct have_global_ifelse<TT, std::void_t< decltype(TT(ifelse(std::declval<TT>(), std::declval<TT>(), std::declval<TT>()))) >> : std::true_type {};
        template<typename TT, typename = void> struct have_nested_ifelse : std::false_type {};
        template<typename TT> struct have_nested_ifelse<TT, std::void_t< decltype(TT(TT::ifelse(std::declval<TT>(), std::declval<TT>(), std::declval<TT>()))) >> : std::true_type {};
        template<typename TT, int num = 0> struct TryGetImpl_ifelse{
            static T apply(T a, T b, T c) { throw std::runtime_error(std::string("Function ") + "ifelse" + " is not implemented by default for " + typeid(TT).name() + " type"); return a; }
        };
        template<typename TT> struct TryGetImpl_ifelse<TT, 1>{ static inline T apply(T a, T b, T c) { return T::ifelse(a, b, c); } };
        template<typename TT> struct TryGetImpl_ifelse<TT, 2>{ static inline T apply(T a, T b, T c) { return ifelse(a, b, c); } };
        template<typename TT> struct TryGetImpl_ifelse<TT, 3>{ static inline T apply(T a, T b, T c) { return a ? b : c; } };
        template<typename TT> using Try_ifelse = TryGetImpl_ifelse<TT, 
            std::is_arithmetic<TT>::value ? 3 : (have_nested_ifelse<TT>::value ? 1 : (have_global_ifelse<TT>::value ? 2 : 0))>;
        virtual T ifelse(T a, T b, T c) const { return Try_ifelse<T>::apply(a, b, c); } 
        virtual T sign(T a) const { return ifelse(gt(a, 0), T(1), ifelse(lt(a, 0), T(-1), T(0))); }   
        REGISTER_FUNC2_0(min, std::min, min);
        REGISTER_FUNC2_0(max, std::max, max);

        template<typename TT, int num = 0> struct TryGetImpl_norm{
            static inline T apply(T a, T b) { throw std::runtime_error(std::string("Function ") + "norm" + " is not implemented by default for " + typeid(TT).name() + " type"); return a; }
        };
        template<typename TT> struct TryGetImpl_norm<TT, 1>{ 
            static inline TT apply(std::vector<TT> a) { 
                if (a.empty()) throw std::runtime_error("Wrong argument in norm function"); 
                auto m = *std::max_element(a.begin(), a.end());
                TT v(0);
                for (auto i: a)
                    v += (i/m)*(i/m);
                return m*std::sqrt(v); 
            } 
        };
        template<typename TT> using Try_norm = TryGetImpl_norm<TT, std::is_floating_point<TT>::value ? 1 : 0>;
        virtual T norm(std::vector<T> a) const { return Try_norm<T>::apply(a); }

    #undef REGISTER_FUNC2_0
    #undef REGISTER_FUNC2_2
    #undef REGISTER_FUNC2
    #undef REGISTER_FUNC1 
    };

    template<bool Dummy>
    struct Energy_Parser::DefaultMathFuncsApplier<casadi::SX, Dummy>{
        using SX = casadi::SX;
        virtual SX abs(SX a) const { return SX::abs(a); }
        virtual SX exp(SX a) const { return SX::exp(a); }
        virtual SX expm1(SX a) const { return SX::exp(a) - 1; }
        virtual SX log(SX a) const { return SX::log(a); }
        virtual SX log1p(SX a) const { return SX::log(1+a); }
        virtual SX sqrt(SX a) const { return SX::sqrt(a); }
        virtual SX cbrt(SX a) const { return SX::pow(a, SX(1)/SX(3)); }
        virtual SX sq(SX a) const { return SX::sq(a); }
        virtual SX cube(SX a) const { return a*a*a; } //SX::constpow(a, 3); }//SX::sq(a)*a; }
        virtual SX sin(SX a) const { return SX::sin(a); }
        virtual SX cos(SX a) const { return SX::cos(a); }
        virtual SX tan(SX a) const { return SX::tan(a); }
        virtual SX asin(SX a) const { return SX::asin(a); }
        virtual SX acos(SX a) const { return SX::acos(a); }
        virtual SX atan(SX a) const { return SX::atan(a); }
        virtual SX sinh(SX a) const { return SX::sinh(a); }
        virtual SX cosh(SX a) const { return SX::cosh(a); }
        virtual SX tanh(SX a) const { return SX::tanh(a); }
        virtual SX erf(SX a) const { return SX::erf(a); }
        virtual SX erfc(SX a) const { return 1 - SX::erf(a); }
        virtual SX tgamma(SX a) const { throw std::runtime_error("tgamma is not implemented for SX type"); return a; }
        virtual SX lgamma(SX a) const { throw std::runtime_error("lgamma is not implemented for SX type"); return a; }
        virtual SX floor(SX a) const { return SX::floor(a); }
        virtual SX ceil(SX a) const { return SX::ceil(a); }
        virtual SX sign(SX a) const { return SX::sign(a); }
        virtual SX plus(SX a, SX b) const { return a + b; }
        virtual SX minus(SX a, SX b) const { return a - b; }
        virtual SX uminus(SX a) const { return -a; }
        virtual SX mul(SX a, SX b) const { return a*b; }
        virtual SX div(SX a, SX b) const { return a/b; }
        virtual SX pow(SX a, SX b) const { return SX::pow(a, b); }
        virtual SX atan2(SX a, SX b) const { return SX::atan2(a, b); }
        virtual SX fmod(SX a, SX b) const { return SX::mod(a, b); }
        virtual SX not_op(SX a) const { return !a; }
        virtual SX lt(SX a, SX b) const { return a < b; }
        virtual SX le(SX a, SX b) const { return a <= b; }
        virtual SX ne(SX a, SX b) const { return a != b; }
        virtual SX eq(SX a, SX b) const { return a == b; }
        virtual SX ge(SX a, SX b) const { return a >= b; }
        virtual SX gt(SX a, SX b) const { return a > b; }
        virtual SX and_op(SX a, SX b) const { return a && b; }
        virtual SX or_op(SX a, SX b) const { return a || b; }
        virtual SX min(SX a, SX b) const { return SX::fmin(a, b); }
        virtual SX max(SX a, SX b) const { return SX::fmax(a, b); }
        virtual SX ifelse(SX a, SX b, SX c) const { return SX::if_else(a, b, c); }
        virtual SX norm(std::vector<SX> a) const { 
            for (auto& s: a) s = SX::reshape(s, s.size1() * s.size2(), 1);
            return SX::norm_2(SX::vertcat(a)); 
        }
    };

    template<typename T, typename MathFuncsApplier>
    T Energy_Parser::evaluateSubtree(NodeID subtree, std::vector<std::pair<bool, T>>& evaluated, const std::vector<T>& params, const std::vector<T>& invs, const MathFuncsApplier& t){
        if (evaluated[subtree].first) return evaluated[subtree].second;
        auto evl = [&](NodeID id)->T{ return evaluateSubtree(id, evaluated, params, invs, t); };
        auto& res = evaluated[subtree].second;
        auto& n = m_expr[subtree];
        switch(n.m_val.m_op){
            case OP_ABS:  res = t.abs(evl(n.child(0))); break;
            case OP_EXP: res = t.exp(evl(n.child(0))); break;
            case OP_EXPM1: res = t.expm1(evl(n.child(0))); break;
            case OP_LOG: res = t.log(evl(n.child(0))); break;
            case OP_LOG1P: res = t.log1p(evl(n.child(0))); break;
            case OP_SQRT: res = t.sqrt(evl(n.child(0))); break;
            case OP_CBRT: res = t.cbrt(evl(n.child(0))); break;
            case OP_CUBE: res = t.cube(evl(n.child(0))); break;
            case OP_SQ: res = t.sq(evl(n.child(0))); break;
            case OP_SIN: res = t.sin(evl(n.child(0))); break;
            case OP_COS: res = t.cos(evl(n.child(0))); break;
            case OP_TAN: res = t.tan(evl(n.child(0))); break;
            case OP_ASIN: res = t.asin(evl(n.child(0))); break;
            case OP_ACOS: res = t.acos(evl(n.child(0))); break;
            case OP_ATAN: res = t.atan(evl(n.child(0))); break;
            case OP_SINH: res = t.sinh(evl(n.child(0))); break;
            case OP_COSH: res = t.cosh(evl(n.child(0))); break;
            case OP_TANH: res = t.tanh(evl(n.child(0))); break;
            case OP_ERF: res = t.erf(evl(n.child(0))); break;
            case OP_ERFC: res = t.erfc(evl(n.child(0))); break;
            case OP_TGAMMA: res = t.tgamma(evl(n.child(0))); break;
            case OP_LGAMMA: res = t.lgamma(evl(n.child(0))); break;
            case OP_FLOOR: res = t.floor(evl(n.child(0))); break;
            case OP_CEIL: res = t.ceil(evl(n.child(0))); break;
            case OP_SIGN: res = t.sign(evl(n.child(0))); break;
            case OP_PLUS: {
                for (auto ch: n.childs()) 
                    evl(ch);
                res = evaluated[n.child(0)].second;
                for (std::size_t i = 1; i < n.number_of_childs(); ++i)
                    evaluated[subtree].second = t.plus(evaluated[subtree].second, evaluated[n.child(i)].second);
                res = evaluated[subtree].second;  
                break;
            }
            case OP_MINUS: res = t.minus(evl(n.child(0)), evl(n.child(1))); break;
            case OP_UMINUS: res = t.uminus(evl(n.child(0))); break;
            case OP_MUL: {
                for (auto ch: n.childs()) 
                    evl(ch);
                res = evaluated[n.child(0)].second;
                for (std::size_t i = 1; i < n.number_of_childs(); ++i)
                    evaluated[subtree].second = t.mul(evaluated[subtree].second, evaluated[n.child(i)].second);
                res = evaluated[subtree].second;  
                break;
            }
            case OP_DIV: res = t.div(evl(n.child(0)), evl(n.child(1))); break;
            case OP_POW: res = t.pow(evl(n.child(0)), evl(n.child(1))); break;
            case OP_ATAN2: res = t.atan2(evl(n.child(0)), evl(n.child(1))); break;
            case OP_FMOD: res = t.fmod(evl(n.child(0)), evl(n.child(1))); break;
            case OP_NOT: res = t.not_op(evl(n.child(0))); break;
            case OP_AND: {
                for (auto ch: n.childs()) 
                    evl(ch);
                res = evaluated[n.child(0)].second;
                for (std::size_t i = 1; i < n.number_of_childs(); ++i)
                    evaluated[subtree].second = t.and_op(evaluated[subtree].second, evaluated[n.child(i)].second);
                res = evaluated[subtree].second;  
                break;
            }
            case OP_OR: {
                for (auto ch: n.childs()) 
                    evl(ch);
                res = evaluated[n.child(0)].second;
                for (std::size_t i = 1; i < n.number_of_childs(); ++i)
                    evaluated[subtree].second = t.or_op(evaluated[subtree].second, evaluated[n.child(i)].second);
                res = evaluated[subtree].second;  
                break;
            }
            case OP_LT: res = t.lt(evl(n.child(0)), evl(n.child(1))); break;
            case OP_LE: res = t.le(evl(n.child(0)), evl(n.child(1))); break;
            case OP_NE: res = t.ne(evl(n.child(0)), evl(n.child(1))); break;
            case OP_EQ: res = t.eq(evl(n.child(0)), evl(n.child(1))); break;
            case OP_GT: res = t.gt(evl(n.child(0)), evl(n.child(1))); break;
            case OP_GE: res = t.ge(evl(n.child(0)), evl(n.child(1))); break;
            case OP_MIN_N:{
                for (auto ch: n.childs()) 
                    evl(ch);
                res = evaluated[n.child(0)].second;
                for (std::size_t i = 1; i < n.number_of_childs(); ++i)
                    evaluated[subtree].second = t.min(evaluated[subtree].second, evaluated[n.child(i)].second);
                res = evaluated[subtree].second;  
                break;
            }
            case OP_MAX_N:{
                for (auto ch: n.childs()) 
                    evl(ch);
                res = evaluated[n.child(0)].second;
                for (std::size_t i = 1; i < n.number_of_childs(); ++i)
                    evaluated[subtree].second = t.max(evaluated[subtree].second, evaluated[n.child(i)].second);
                res = evaluated[subtree].second;  
                break;
            }
            case OP_HYPOT_N:{
                std::vector<T> to_evl; to_evl.reserve(n.number_of_childs());
                for (auto ch: n.childs()) 
                    to_evl.emplace_back(evl(ch));
                res = t.norm(to_evl); 
                break;
            }
            case OP_IFELSE: res = t.ifelse(evl(n.child(0)), evl(n.child(1)), evl(n.child(2))); break;
            case H_HEAD: res = evl(n.child(0)); break;
            case W_VAL: res = T(n.m_val.m_value.d); break;
            case W_IVAL: res = T(n.m_val.m_value.l); break;
            case W_VINIT: res = invs[n.m_val.m_value.l]; break;
            case W_PAR: res = params[n.m_val.m_value.l]; break;
            case W_VLET: res = evl(m_vlets[n.m_val.m_value.l].id); break;
            default: 
                throw std::runtime_error("Faced unknown operator type = " + std::to_string(n.m_val.m_op));
        }
        evaluated[subtree].first = true;
        return res;
    }

    template<typename T, typename MathFuncsApplier>
    T Energy_Parser::Evaluate(const std::vector<T>& params, const std::vector<T>& invs, const MathFuncsApplier& t){
        std::vector<std::pair<bool, T>> evaluated(m_expr.maxID(), {false, T()});
        return evaluateSubtree(0, evaluated, params, invs, t);
    }
}

#endif //AORTIC_VALVE_HYPERELASTICFORMAT_INL