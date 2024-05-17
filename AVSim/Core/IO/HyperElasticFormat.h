//
// Created by alex on 30.09.2023.
//

#ifndef AORTIC_VALVE_HYPERELASTICFORMAT_H
#define AORTIC_VALVE_HYPERELASTICFORMAT_H

#include <string>
#include <vector>
#include <stack>
#include <set>
#include <stdexcept>
#include <cmath>
#include <map>
#include <iostream>
#include <algorithm>
#include "MultiTree.h"

namespace World3d{

struct Energy_Parser{
    enum Lexems{
        KW_Var,     // "Var"
        KW_Let,     // "Let"
        KW_Fiber,   // "Fiber"
        KW_Potential, // "Potential"
        LX_SEMICOLON, // ';'
        LX_COMMA,   // ','
        LX_COLON,   // ':'
        LX_ASSIGN,  // '='
        LX_PLUS,    // '+'
        LX_MINUS,   // '-'
        LX_MUL,     // '*'
        LX_DIV,     // '/'
        LX_POW,     // '^'
        LX_NOT,     // "not"
        LX_LT,      // '<'
        LX_EQ,      // "=="
        LX_NE,      // "!="
        LX_AND,     // "and"
        LX_OR,      // "or"
        LX_LE,      // "<="
        LX_GE,      // ">=""
        LX_GT,      // '>'
        LX_NUM,     // (d+|d*.d+)((e|E)(+|-)?d+)?
        LX_INUM,    // d+
        LX_STR,     // "w*"
        LX_VAR,     // a+w*([(d+|)])?
        LX_OPBRACKET, // '('
        LX_CLBRACKET, // ')'
        LX_OPSQRBRACKET,//'['
        LX_CLSQRBRACKET,//']'
    };
    struct LexElem{
        Lexems lexem;
        std::size_t pos_st, pos_end;
        union {
            double d = 0;
            std::size_t l;
        } value; //< used only for LX_NUM
        LexElem() = default;
        LexElem(Lexems lexem, std::size_t pos_st, std::size_t pos_end, double value = 0): lexem{lexem}, pos_st{pos_st}, pos_end{pos_end}, value{value} {}
    };

    enum Operations{
        OP_ABS,
        OP_EXP,
        OP_EXPM1,
        OP_LOG,
        OP_LOG1P,
        OP_SQRT,
        OP_CBRT,
        OP_SQ,
        OP_CUBE,
        OP_SIN,
        OP_COS,
        OP_TAN,
        OP_ASIN,
        OP_ACOS,
        OP_ATAN,
        OP_SINH,
        OP_COSH,
        OP_TANH,
        OP_ERF,
        OP_ERFC,
        OP_TGAMMA,
        OP_LGAMMA,
        OP_FLOOR,
        OP_CEIL,
        OP_SIGN,

        OP_PLUS,
        OP_MINUS,
        OP_UMINUS, //unary minus
        OP_MUL,
        OP_DIV,
        OP_POW,
        // OP_ASSIGN,
        OP_ATAN2,
        OP_FMOD,

        OP_MIN_N,
        OP_MAX_N,
        OP_HYPOT_N,
        
        OP_IFELSE,
        OP_NOT,
        OP_LT,
        OP_EQ,
        OP_NE,
        OP_AND,
        OP_OR,
        OP_LE,
        OP_GE,
        OP_GT,

        S_ANY_VAR,   //temporary anonym symbol

        W_VAL,
        W_IVAL,
        W_VINIT, //predefined variable
        W_PAR,
        W_VLET,
        
        H_HEAD,

    };

    struct TreeValue{
        Operations m_op = H_HEAD;
        union {
            double d = 0;
            long l;
        } m_value;
        std::size_t ilex;

        TreeValue& set_value(double d, Operations op = W_VAL){ return m_value.d = d, m_op = op, *this; }
        TreeValue& set_value(long l, Operations op = W_IVAL){ return m_value.l = l, m_op = op, *this; }
        TreeValue& set_value(int l, Operations op = W_IVAL){ return m_value.l = l, m_op = op, *this; }
        TreeValue& set_value(bool l, Operations op = W_IVAL) { return m_value.l = l, m_op = op, *this; }
    };
    using NodeID = MultiTree<TreeValue>::NodeID;
    struct Parameter{
        double def_val = NAN;
        std::string def_tag = "";
        std::string name;

        Parameter() = default;
        Parameter(std::string name, double val): name{name}, def_val{val} {};
        Parameter(std::string name, std::string tag_name): name{name}, def_tag{tag_name} {};
        Parameter(std::string name, std::string tag_name, double val): name{name}, def_tag{tag_name}, def_val{val} {};
    };
    struct Fiber{
        std::string name;
        std::string def_tag = "";

        Fiber() = default;
        Fiber(std::string name, std::string tag_name = ""): name{name}, def_tag{tag_name} {};
    };
    struct Invariants{
        enum Type{
            I1_T,
            I2_T,
            J_T,
            If_T,
            Ifs_T,
        };
        Type type = I1_T;
        std::size_t v1 = std::numeric_limits<std::size_t>::max();
        std::size_t v2 = std::numeric_limits<std::size_t>::max();
    };
    struct SubVar{
        std::string name;
        NodeID id = 0;
    };
    struct MathFuncMeta{
        int supposed_nargs = -1;
        Operations operation = H_HEAD;

        bool isValid() const { return operation != H_HEAD;}
        operator bool() const { return isValid(); }
        MathFuncMeta() = default;
        MathFuncMeta(Operations op, int nargs): supposed_nargs{nargs}, operation{op} {}
    };
    struct TrivialValue{
        union {
            double d = 0;
            long l;
        } m_value;
        enum Type{
            IVAL,
            RVAL
        };
        Type m_type = RVAL;
        TrivialValue() = default;
        TrivialValue(double d): m_type{RVAL} { m_value.d = d; }
        TrivialValue(long l): m_type{IVAL} { m_value.l = l; }
        TrivialValue& set(double d) { m_type = RVAL, m_value.d = d; return *this; }
        TrivialValue& set(long l) { m_type = IVAL, m_value.l = l; return *this; }
        double Real() const { return m_value.d; }
        long Integer() const { return m_value.l; }
        bool Bool() const { 
            switch (m_type) {
                case IVAL: return m_value.l != 0;
                case RVAL: return m_value.d != 0.0;
            }
            return false;
        }
        template<typename T>
        typename std::enable_if<std::is_integral<T>::value || std::is_floating_point<T>::value, T>::type get() const { 
            switch (m_type) {
                case IVAL: return static_cast<T>(m_value.l);
                case RVAL: return static_cast<T>(m_value.d);
            }
            return T(); 
        }
    };

    Energy_Parser(const std::string& code = ""): m_code(code) {}
    Energy_Parser& setCode(const std::string& code) { return m_code = code, *this; }
    Energy_Parser& setCodeFromFile(const std::string& filename);
    const std::string& getCode() const { return m_code; }
    Energy_Parser& LexicalAnalysis() { if (!performLexicalAnalysis()) throw std::runtime_error(m_err_msg); return *this; }
    Energy_Parser& SyntaxAnalysis() { if (!performSyntaxAnalysis()) throw std::runtime_error(m_err_msg); return *this; } 
    Energy_Parser& Simplify();
    
    template<typename T, bool Dummy = true> 
    struct DefaultMathFuncsApplier;
    template<typename T, typename MathFuncsApplier = DefaultMathFuncsApplier<T>>
    T Evaluate(const std::vector<T>& params, const std::vector<T>& invs, const MathFuncsApplier& t = MathFuncsApplier());

    bool performLexicalAnalysis() noexcept;
    bool performSyntaxAnalysis() noexcept;
    std::ostream& debud_print_lexical_buffer(std::ostream& out = std::cout) const;
    std::ostream& print_parsed_state(std::ostream& out = std::cout, const std::string& prefix = "") const;
    
    const std::vector<Parameter>& getParameters() const { return m_params; }
    std::vector<Parameter>& getParameters() { return m_params; }
    const std::vector<Fiber>& getFibers() const { return m_fibers; }
    std::vector<Fiber>& getFibers() { return m_fibers; }
    const std::vector<Invariants>& getInvariants() const { return m_invs; }
protected:
    std::string m_code;
    std::string m_err_msg;
    std::vector<LexElem> m_lex;
    int m_ErrorStatus = 0;

    std::vector<Parameter> m_params;
    std::vector<Fiber> m_fibers;
    std::vector<Invariants> m_invs;
    std::vector<SubVar> m_vlets;
    struct VarMeta{
        enum Type{
            PAR_T,
            LET_T,
            FIB_T,
            FUNC_T,
        };
        Type m_type = PAR_T;
        std::size_t id = std::numeric_limits<std::size_t>::max();
        VarMeta() = default;
        VarMeta(Type type, std::size_t id = std::numeric_limits<std::size_t>::max()) : m_type{type}, id{id} {}
    };
    std::map<std::string, VarMeta> m_symbol_map;

    void register_param(Parameter p);
    void register_fiber(Fiber p);
    std::size_t find_fiber(const std::string& name) const;
    void register_let_variable(SubVar p);
    void register_invariant(Invariants p);
    std::size_t find_or_register_invariant(Invariants p);
    void register_math_functions();

    MultiTree<TreeValue> m_expr;

    bool isSkipableOrEnd(std::size_t pos) const;
    bool parseSpace(std::size_t& pos) const;
    bool parseComment(std::size_t& pos) const;
    bool parseEmpty(std::size_t& pos) const;
    bool parseKeyWord(std::size_t& pos);
    bool parseString(std::size_t& pos);
    bool parseVar(std::size_t& pos);
    bool parseUInt(std::size_t& pos);
    bool parseNum(std::size_t& pos);
    void create_lexical_error_message(std::size_t pos);
    bool parseLexem(std::size_t& pos);

    bool testNoErrors() { return m_ErrorStatus == 0; }
    static const std::map<std::string, MathFuncMeta>& getMathFunctionList();
    static const std::map<Operations, std::string>& getMathBackFunctionList();
    static MathFuncMeta getMathFunctionMeta(std::string name);
    static const std::map<Operations, std::string>& getOperationNameList();
    static std::string getOperationName(Operations op) { auto& s = getOperationNameList(); auto it = s.find(op); return it != s.end() ? it->second : ""; }

    NodeID interp_VAL(std::size_t& ilex);
    NodeID interp_ANY_VAR(std::size_t& ilex);
    NodeID interp_MATH_FUNC(std::size_t& ilex);
    NodeID interp_PREDEF_VAR(std::size_t& ilex);
    NodeID interp_ExprInBrackets(std::size_t& ilex);
    NodeID interp_Simple(std::size_t& ilex);
    NodeID interp_UMINUS(std::size_t& ilex);
    NodeID interp_BracketOP(std::size_t& ilex);
    NodeID interp_POW(std::size_t& ilex);
    NodeID interp_MUL_DIV(std::size_t& ilex);
    NodeID interp_PLUS_MINUS(std::size_t& ilex);
    NodeID interp_COMPARE(std::size_t& ilex);
    NodeID interp_NOT(std::size_t& ilex);
    NodeID interp_AND_OR(std::size_t& ilex);
    NodeID interp_Expression(std::size_t& ilex);
    
    NodeID interp_LetAssign(std::size_t& ilex);
    NodeID interp_FiberAssign(std::size_t& ilex);
    NodeID interp_VarAssign(std::size_t& ilex);
    NodeID interp_Let(std::size_t& ilex);
    NodeID interp_Fiber(std::size_t& ilex);
    NodeID interp_Var(std::size_t& ilex);
    NodeID interp_Potential(std::size_t& ilex);

    NodeID GeneralSyntaxError(std::size_t operand_ilex, std::string prenum_text, std::string afternum_text);
    NodeID NotFoundAfterOperandError(std::size_t operand_ilex, std::string expected_object);
    NodeID NotFoundOperandError(std::size_t operand_ilex, std::string expected_operand);
    NodeID FacedUnknownSymbolError(std::size_t operand_ilex, const std::string& object_type = "undefined symbol");
    NodeID RedefinitionSymbolError(std::size_t operand_ilex);
    NodeID NotFindPotentialSection();
    NodeID SomeCodeAfterSection(std::size_t operand_ilex);
    
    double evaluate_trivial_expression(std::size_t operand_ilex, NodeID subtree);
    std::pair<TrivialValue, bool> getTrivialValue(NodeID subtree) const;
    bool simplify_subexpr_iter(NodeID subtree);
    void simplify_subexpr(NodeID subtree) { while(simplify_subexpr_iter(subtree)); }
    void unmap_anon_variables(NodeID subtree);
    
    template<typename T, typename MathFuncsApplier>
    T evaluateSubtree(NodeID subtree, std::vector<std::pair<bool, T>>& evaluated, const std::vector<T>& params, const std::vector<T>& invs, const MathFuncsApplier& t);

    // op E
    NodeID interp_UNARY_OP(std::size_t& ilex, std::initializer_list<std::pair<Lexems, Operations>> op_map, NodeID(Energy_Parser::*next)(std::size_t& ilex), const std::string& next_name = "expression");
    // (op)? E
    NodeID interp_UNARY_OP_or_OP(std::size_t& ilex, std::initializer_list<std::pair<Lexems, Operations>> op_map, NodeID(Energy_Parser::*next)(std::size_t& ilex), const std::string& next_name = "expression");
    // E (op E)*
    NodeID interp_FUSED_BIN_OP(std::size_t& ilex, std::initializer_list<std::pair<Lexems, Operations>> op_map, NodeID(Energy_Parser::*next)(std::size_t& ilex), const std::string& next_name = "expression");
    // E (op E)?
    NodeID interp_BIN_OP(std::size_t& ilex, std::initializer_list<std::pair<Lexems, Operations>> op_map, NodeID(Energy_Parser::*next)(std::size_t& ilex), const std::string& next_name = "expression");
    // left( E right) 
    NodeID interp_BRACKET(std::size_t& ilex, Lexems left, Lexems right, NodeID(Energy_Parser::*next)(std::size_t& ilex), const std::string& closed_bracket, const std::string& next_name = "expression");
    // E1 | E2 | E3 | ... | En
    NodeID interp_OR_EXPR(std::size_t& ilex, std::initializer_list<NodeID(Energy_Parser::*)(std::size_t& ilex)> nexts);
    // st (OP (, OP)+)? end
    NodeID interp_LISTED_OP(std::size_t& ilex, Lexems lx_st, Lexems lx_sep, Lexems lx_end, NodeID(Energy_Parser::*next)(std::size_t& ilex), const std::string& end_seq, const std::string& next_name = "expression");

    struct CodePositionInfo{
        std::size_t iline = 0;
        std::size_t line_st = 0, line_end = 0;
        std::size_t part_st = 0, part_end = 0;
        std::size_t pos = 0;
    };
    CodePositionInfo computePositionInfo(std::size_t pos) const;

    std::ostream& print_invariant(std::ostream& out, std::size_t inv_order_num) const;
    std::ostream& print_expression(std::ostream& out, NodeID id, bool not_substitute_let_exprs) const;
};

}

#include "HyperElasticFormat.inl"

#endif //AORTIC_VALVE_HYPERELASTICFORMAT_H