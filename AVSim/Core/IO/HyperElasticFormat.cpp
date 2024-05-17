#include "HyperElasticFormat.h"
#include <cassert>
#include <fstream>

using namespace World3d;

bool Energy_Parser::parseComment(std::size_t& pos) const  {
    if (pos < m_code.size() && m_code[pos] == '#') {
        do{
            ++pos;
        } while (m_code[pos] != '\n' && pos < m_code.size());
        if (m_code[pos] == '\n' && pos < m_code.size()) 
            ++pos;
        return true;    
    }
    return false;
}

bool Energy_Parser::parseSpace(std::size_t& pos) const  {
    std::size_t lpos = pos;
    for (; pos < m_code.size() && std::isspace(m_code[lpos]); ++lpos);
    bool change = (lpos > pos);
    pos = lpos;
    return change;
}

bool Energy_Parser::isSkipableOrEnd(std::size_t pos) const{
    if (pos == m_code.size()) return true;
    if (pos > m_code.size()) return false;
    if (m_code[pos] == '#' || std::isspace(m_code[pos])) return true;
    return false;
}

bool Energy_Parser::parseEmpty(std::size_t& pos) const{
    bool lchange = false, gchange = false;
    do{
        lchange = parseSpace(pos);
        if (!lchange) lchange = parseComment(pos);
        gchange |= lchange;
    } while(lchange);
    return gchange;
}

bool Energy_Parser::parseKeyWord(std::size_t& pos){
    static const std::vector<std::pair<std::string, Lexems>> kws{
        {"Var", KW_Var},
        {"Let", KW_Let},
        {"Fiber", KW_Fiber},
        {"Potential", KW_Potential},
        {"not", LX_NOT},
        {"and", LX_AND},
        {"or", LX_OR},
    };
    for (auto& v: kws){
        if (pos + v.first.size() > m_code.size()) continue;
        auto lpos = pos + v.first.size();
        if (!(lpos == m_code.size() || isSkipableOrEnd(lpos))) continue;
        if (v.first == std::string(m_code.data() + pos, m_code.data() + pos + v.first.size())){
            m_lex.emplace_back(v.second, pos, pos + v.first.size(), 0.0); return pos += v.first.size(), true;
        }
    }
    return false;
}
static std::string replaceEscapeSequences(const char* st, const char* end){
    std::string res; res.reserve(end - st);
    while(st < end){
        if (*st != '\\' || st + 1 == end){
            res.push_back(*st);
            ++st;
            continue;
        }
        ++st;
        // \a, \b, \f, \n, \r, \t, \v, \", \\, \[0-7]
        switch (*st){
            case 'a': res.push_back('\a'); break;
            case 'b': res.push_back('\b'); break;
            case 'f': res.push_back('\f'); break;
            case 'n': res.push_back('\n'); break;
            case 'r': res.push_back('\r'); break;
            case 't': res.push_back('\t'); break;
            case 'v': res.push_back('\v'); break;
            case '\"': res.push_back('\"'); break;
            case '\\': res.push_back('\\'); break;
            case '0': res.push_back('\0'); break;
            case '1': res.push_back('\1'); break;
            case '2': res.push_back('\2'); break;
            case '3': res.push_back('\3'); break;
            case '4': res.push_back('\4'); break;
            case '5': res.push_back('\5'); break;
            case '6': res.push_back('\6'); break;
            case '7': res.push_back('\7'); break;
            default:{
                res.push_back('\\');
                res.push_back(*st);
            }
        }
        ++st;
    }
    return res;
}
bool Energy_Parser::parseString(std::size_t& pos){
    if (pos >= m_code.size()) return false;
    if (m_code[pos] != '\"') return false;
    auto lpos = pos + 1;
    for (; lpos < m_code.size(); ++lpos){
        if (m_code[lpos] == '\n') break;
        if (m_code[lpos] != '\"') continue;
        std::size_t k = 0;
        for ( ; m_code[lpos-1-k] == '\\'; ++k);
        if (k%2 != 0) continue;
        m_lex.emplace_back(LX_STR, pos, lpos + 1, 0.0);
        return pos = lpos + 1, true;
    }
    m_ErrorStatus = 1;
    return false;
}

bool Energy_Parser::parseUInt(std::size_t& pos){
    if (pos >= m_code.size()) return false;
    if (!std::isdigit(m_code[pos])) return false;
    auto lpos = pos;
    std::size_t val = 0;
    while (pos < m_code.size() && std::isdigit(m_code[pos])){
        val = (m_code[pos] - '0') + val*10;
        ++pos; 
    }
    m_lex.emplace_back(LX_INUM, lpos, pos, 0);
    m_lex.back().value.l = val;
    return true;
}

bool Energy_Parser::parseNum(std::size_t& pos){
    if (pos >= m_code.size()) return false;
    if (!std::isdigit(m_code[pos]) && m_code[pos] != '.') return false;
    if (m_code[pos] == '.' && (pos + 1 >= m_code.size() || !std::isdigit(m_code[pos+1]))) return false;
    auto lpos = pos;
    double val = 0;
    while (pos < m_code.size() && std::isdigit(m_code[pos])){
        val = (m_code[pos] - '0') + val*10;
        ++pos; 
    }
    int deg_exp = 0; bool is_int = true;
    if (pos < m_code.size() && m_code[pos] == '.'){
        is_int = false;
        ++pos;
        while (pos < m_code.size() && std::isdigit(m_code[pos])){
            val = (m_code[pos] - '0') + val*10;
            ++pos; 
            --deg_exp;
        }
    }
    if (pos+1 >= m_code.size() || (m_code[pos] != 'e' && m_code[pos] != 'E')){
        if (!is_int)
            m_lex.emplace_back(LX_NUM, lpos, pos, val*std::pow(10, deg_exp));
        else  {
            m_lex.emplace_back(LX_INUM, lpos, pos, val);
            m_lex.back().value.l = val;
        }   
        return true;
    }
    //m_code[pos] == 'e' || m_code[pos] == 'E'
    if (std::isdigit(m_code[pos + 1])){
        ++pos;
        int add_exp = 0;
        while (pos < m_code.size() && std::isdigit(m_code[pos])){
            add_exp = (m_code[pos] - '0') + add_exp*10;
            ++pos; 
        }
        m_lex.emplace_back(LX_NUM, lpos, pos, val*std::pow(10, deg_exp + add_exp));
        return true;
    }
    if (pos+2 < m_code.size() && (m_code[pos + 1] == '+' || m_code[pos + 1] == '-') &&  std::isdigit(m_code[pos + 2]) ){
        ++pos;
        bool sign = (m_code[pos] == '+');
        ++pos;
        int add_exp = 0;
        while (pos < m_code.size() && std::isdigit(m_code[pos])){
            add_exp = (m_code[pos] - '0') + add_exp*10;
            ++pos; 
        }
        m_lex.emplace_back(LX_NUM, lpos, pos, val*std::pow(10, deg_exp + (sign ? 1 : -1)*add_exp));
        return true;
    }
    m_lex.emplace_back(LX_NUM, lpos, pos, val*std::pow(10, deg_exp));
    return true;
}

bool Energy_Parser::parseVar(std::size_t& pos){
    if (pos >= m_code.size()) return false;
    auto myisalpha = [](char c) { return c == '_' || std::isalnum(c); };
    if (!myisalpha(m_code[pos])) return false;
    auto lpos = pos;
    ++pos;
    while(pos < m_code.size() && myisalpha(m_code[pos])) ++pos;
    m_lex.emplace_back(LX_VAR, lpos, pos, 0);
    return true;
}

bool Energy_Parser::parseLexem(std::size_t& pos){
    if (pos >= m_code.size()) return false;
    switch (m_code[pos]){
        case ';': m_lex.emplace_back(LX_SEMICOLON, pos, pos + 1, 0.0); return ++pos, true;
        case ',': m_lex.emplace_back(LX_COMMA, pos, pos + 1, 0.0); return ++pos, true;
        case ':': m_lex.emplace_back(LX_COLON, pos, pos + 1, 0.0); return ++pos, true;
        case '=': {
            if (pos+1 == m_code.size() || m_code[pos+1] != '=') {
                m_lex.emplace_back(LX_ASSIGN, pos, pos + 1, 0.0);
                ++pos;
            } else {
                m_lex.emplace_back(LX_EQ, pos, pos + 2, 0.0); 
                pos += 2;
            }
            return true;   
        }
        case '+': m_lex.emplace_back(LX_PLUS, pos, pos + 1, 0.0); return ++pos, true;
        case '-': m_lex.emplace_back(LX_MINUS, pos, pos + 1, 0.0); return ++pos, true;
        case '*': m_lex.emplace_back(LX_MUL, pos, pos + 1, 0.0); return ++pos, true;
        case '/': m_lex.emplace_back(LX_DIV, pos, pos + 1, 0.0); return ++pos, true;
        case '^': m_lex.emplace_back(LX_POW, pos, pos + 1, 0.0); return ++pos, true;
        case '(': m_lex.emplace_back(LX_OPBRACKET, pos, pos + 1, 0.0); return ++pos, true;
        case ')': m_lex.emplace_back(LX_CLBRACKET, pos, pos + 1, 0.0); return ++pos, true;
        case '[': m_lex.emplace_back(LX_OPSQRBRACKET, pos, pos + 1, 0.0); return ++pos, true;
        case ']': m_lex.emplace_back(LX_CLSQRBRACKET, pos, pos + 1, 0.0); return ++pos, true;
        case '<': {
            if (pos+1 == m_code.size() || m_code[pos+1] != '=') {
                m_lex.emplace_back(LX_LT, pos, pos + 1, 0.0);
                ++pos;
            } else {
                m_lex.emplace_back(LX_LE, pos, pos + 2, 0.0); 
                pos += 2;
            }
            return true;   
        }
        case '>': {
            if (pos+1 == m_code.size() || m_code[pos+1] != '=') {
                m_lex.emplace_back(LX_GT, pos, pos + 1, 0.0);
                ++pos;
            } else {
                m_lex.emplace_back(LX_GE, pos, pos + 2, 0.0); 
                pos += 2;
            }
            return true;   
        }
    }
    if (parseKeyWord(pos)) return true;
    if (parseString(pos)) return true;
    if (parseNum(pos)) return true;
    if (parseVar(pos)) return true;
    return false;
}

Energy_Parser::CodePositionInfo Energy_Parser::computePositionInfo(std::size_t pos) const{
    assert(pos < m_code.size() && "Wrong pos");
    std::size_t line_st = 0, line_end = 0, iline = 0, lpos = 0;
    while(line_end < pos && lpos < m_code.size()){
        if (m_code[lpos] == '\n'){
            line_st = (line_end > 0 ? line_end+1 : 0);
            line_end = lpos;
            iline++; 
        }
        ++lpos;
    }
    if (line_end < pos) line_end = m_code.size();
    std::size_t part_st = line_st, part_end = line_end;
    if (line_end - line_st > 256){
        part_st = std::max(pos - 128, line_st);
        part_end = std::min(pos + 128, line_end);
    }
    CodePositionInfo info;
    info.iline = iline > 0 ? iline - 1 : 0;
    info.line_st = line_st, info.line_end = line_end;
    info.part_st = part_st, info.part_end = part_end;
    info.pos = pos;
    return info;
}

void Energy_Parser::create_lexical_error_message(std::size_t pos){
    auto info = computePositionInfo(pos);
    if (m_ErrorStatus == 1 && m_code[pos] == '\"') {
        m_err_msg = "Lexical Error: Start of line character \" encountered at " + std::to_string(info.iline+1) + ":" + std::to_string(pos - info.line_st + 1) + ", but end of line not found :\n";
    } else 
        m_err_msg = "Lexical Error: Can't parse literal at " + std::to_string(info.iline+1) + ":" + std::to_string(pos - info.line_st + 1) + " : \n";
    m_err_msg += std::string(m_code.data() + info.part_st, m_code.data() + info.part_end);
    if (m_err_msg.back() != '\n') m_err_msg += '\n';
    m_err_msg += std::string(pos - info.part_st, ' ') + "^\n";
}

bool Energy_Parser::performLexicalAnalysis() noexcept{
    m_lex.resize(0); m_ErrorStatus = 0;
    for (std::size_t pos = 0, pos_end = m_code.size(); pos < pos_end; ){
        parseEmpty(pos);
        bool read = parseLexem(pos);
        if (!read && pos != pos_end){
            create_lexical_error_message(pos);
            return false;
        }
    }
    return true;
}

std::ostream& Energy_Parser::debud_print_lexical_buffer(std::ostream& out) const{
    static const std::map<Lexems, std::string> lex_map{
        {KW_Var, "KW_Var"},    
        {KW_Let, "KW_Let"},     
        {KW_Fiber, "KW_Fiber"},   
        {KW_Potential, "KW_Potential"}, 
        {LX_SEMICOLON, "LX_SEMICOLON"}, 
        {LX_COMMA, "LX_COMMA"},   
        {LX_COLON, "LX_COLON"},   
        {LX_ASSIGN, "LX_ASSIGN"},  
        {LX_PLUS, "LX_PLUS"},    
        {LX_MINUS, "LX_MINUS"},   
        {LX_MUL, "LX_MUL"},     
        {LX_DIV, "LX_DIV"},     
        {LX_POW, "LX_POW"},     
        {LX_NOT, "LX_NOT"},     
        {LX_LT, "LX_LT"},      
        {LX_EQ, "LX_EQ"},      
        {LX_NE, "LX_NE"},      
        {LX_AND, "LX_AND"},     
        {LX_OR, "LX_OR"},      
        {LX_LE, "LX_LE"},      
        {LX_GE, "LX_GE"},      
        {LX_GT, "LX_GT"},      
        {LX_NUM, "LX_NUM"},     
        {LX_STR, "LX_STR"},     
        {LX_VAR, "LX_VAR"},     
        {LX_OPBRACKET, "LX_OPBRACKET"}, 
        {LX_CLBRACKET, "LX_CLBRACKET"},
        {LX_OPSQRBRACKET, "LX_OPSQRBRACKET"}, 
        {LX_CLSQRBRACKET, "LX_CLSQRBRACKET"},
        {LX_INUM, "LX_INUM"}
    };
    int ii = 0;
    for (auto& v: m_lex){
        auto str_num = std::to_string(ii++);
        if (str_num.size() < 4) str_num.resize(4, ' ');
        auto it = lex_map.find(v.lexem);
        unsigned shift = (v.lexem == LX_STR) ? 1 : 0;
        auto str = it->second; str.resize(16, ' ');
        out << str_num << " : " << str << " : \"" << std::string(m_code.data() + v.pos_st + shift, m_code.data() + v.pos_end - shift) + "\"";
        if (v.lexem == LX_NUM) out << " : " << v.value.d;
        else if (v.lexem == LX_INUM) out << " : " << v.value.l;
        else if (v.lexem == LX_STR) out << " : \"" << replaceEscapeSequences(m_code.data() + v.pos_st + shift, m_code.data() + v.pos_end - shift) << "\"";
        out << "\n";
    }
    return out;
}

Energy_Parser& Energy_Parser::setCodeFromFile(const std::string& filename){
    std::ifstream f(filename, std::ios::in | std::ios::binary | std::ios::ate);
    if (!f.is_open())
        throw std::runtime_error("Can't open file = \"" + filename + "\"");
    std::ifstream::pos_type fileSize = f.tellg();
    m_code.resize(fileSize);
    f.seekg(0, std::ios::beg);
    f.read(m_code.data(), fileSize);

    return *this;
}

const std::map<std::string, Energy_Parser::MathFuncMeta>& Energy_Parser::getMathFunctionList(){
    static std::map<std::string, MathFuncMeta> s_funcs_meta{
        {"abs", {OP_ABS, 1}},
        {"exp", {OP_EXP, 1}},
        {"expm1", {OP_EXPM1, 1}},
        {"log", {OP_LOG, 1}},
        {"log1p", {OP_LOG1P, 1}},
        {"sqrt", {OP_SQRT, 1}},
        {"cbrt", {OP_CBRT, 1}},
        {"sq", {OP_SQ, 1}},
        {"cube", {OP_CUBE, 1}},
        {"sin", {OP_SIN, 1}},
        {"cos", {OP_COS, 1}},
        {"tan", {OP_TAN, 1}},
        {"asin", {OP_ASIN, 1}},
        {"acos", {OP_ACOS, 1}},
        {"atan", {OP_ATAN, 1}},
        {"sinh", {OP_SINH, 1}},
        {"cosh", {OP_COSH, 1}},
        {"tanh", {OP_TANH, 1}},
        {"erf", {OP_ERF, 1}},
        {"erfc", {OP_ERFC, 1}},
        {"tgamma", {OP_TGAMMA, 1}},
        {"lgamma", {OP_LGAMMA, 1}},
        {"floor", {OP_FLOOR, 1}},
        {"ceil", {OP_CEIL, 1}},
        {"sign", {OP_SIGN, 1}},
        {"atan2", {OP_ATAN2, 2}},
        {"fmod", {OP_FMOD, 2}},
        {"min", {OP_MIN_N, -1}},
        {"max", {OP_MAX_N, -1}},
        {"hypot", {OP_HYPOT_N, -1}},
        {"ifelse", {OP_IFELSE, 3}},
    };
    return s_funcs_meta;
}
const std::map<Energy_Parser::Operations, std::string>& Energy_Parser::getMathBackFunctionList(){
    static std::map<Operations, std::string> s_back_fmap {
        {OP_ABS, "abs"},
        {OP_EXP, "exp"},
        {OP_EXPM1, "expm1"},
        {OP_LOG, "log"},
        {OP_LOG1P, "log1p"},
        {OP_SQRT, "sqrt"},
        {OP_CBRT, "cbrt"},
        {OP_SQ, "sq"},
        {OP_CUBE, "cube"},
        {OP_SIN, "sin"},
        {OP_COS, "cos"},
        {OP_TAN, "tan"},
        {OP_ASIN, "asin"},
        {OP_ACOS, "acos"},
        {OP_ATAN, "atan"},
        {OP_SINH, "sinh"},
        {OP_COSH, "cosh"},
        {OP_TANH, "tanh"},
        {OP_ERF, "erf"},
        {OP_ERFC, "erfc"},
        {OP_TGAMMA, "tgamma"},
        {OP_LGAMMA, "lgamma"},
        {OP_FLOOR, "floor"},
        {OP_CEIL, "ceil"},
        {OP_SIGN, "sign"},
        {OP_ATAN2, "atan2"},
        {OP_FMOD, "fmod"},
        {OP_MIN_N, "min"},
        {OP_MAX_N, "max"},
        {OP_HYPOT_N, "hypot"},
        {OP_IFELSE, "ifelse"},
    };
    return s_back_fmap;
}
const std::map<Energy_Parser::Operations, std::string>& Energy_Parser::getOperationNameList(){
    static std::map<Operations, std::string> s_op_map {
        {OP_ABS, "abs"},
        {OP_EXP, "exp"},
        {OP_EXPM1, "expm1"},
        {OP_LOG, "log"},
        {OP_LOG1P, "log1p"},
        {OP_SQRT, "sqrt"},
        {OP_CBRT, "cbrt"},
        {OP_SQ, "sq"},
        {OP_CUBE, "cube"},
        {OP_SIN, "sin"},
        {OP_COS, "cos"},
        {OP_TAN, "tan"},
        {OP_ASIN, "asin"},
        {OP_ACOS, "acos"},
        {OP_ATAN, "atan"},
        {OP_SINH, "sinh"},
        {OP_COSH, "cosh"},
        {OP_TANH, "tanh"},
        {OP_ERF, "erf"},
        {OP_ERFC, "erfc"},
        {OP_TGAMMA, "tgamma"},
        {OP_LGAMMA, "lgamma"},
        {OP_FLOOR, "floor"},
        {OP_CEIL, "ceil"},
        {OP_SIGN, "sign"},
        {OP_ATAN2, "atan2"},
        {OP_FMOD, "fmod"},
        {OP_MIN_N, "min"},
        {OP_MAX_N, "max"},
        {OP_HYPOT_N, "hypot"},
        {OP_IFELSE, "ifelse"},
        {OP_PLUS, "+"},
        {OP_MINUS, "-"},
        {OP_UMINUS, "-"},
        {OP_MUL, "*"},
        {OP_DIV, "/"},
        {OP_POW, "^"},
        {OP_NOT, "not"},
        {OP_LT, "<"},
        {OP_EQ, "=="},
        {OP_NE, "!="},
        {OP_AND, "and"},
        {OP_OR, "or"},
        {OP_LE, "<="},
        {OP_GE, ">="},
        {OP_GT, ">"}
    };
    return s_op_map;
}

Energy_Parser::MathFuncMeta Energy_Parser::getMathFunctionMeta(std::string name){
    auto& m = getMathFunctionList();
    auto it = m.find(name);
    if (it != m.end()) return it->second;
    return MathFuncMeta();
}

void Energy_Parser::register_param(Parameter p){
    m_params.push_back(p);
    m_symbol_map.insert({p.name, VarMeta(VarMeta::PAR_T, m_params.size() - 1)});
}
void Energy_Parser::register_fiber(Fiber p){
    m_fibers.push_back(p);
    m_symbol_map.insert({p.name, VarMeta(VarMeta::FIB_T, m_fibers.size() - 1)});
}
void Energy_Parser::register_let_variable(SubVar p){
    m_vlets.push_back(p);
    m_symbol_map.insert({p.name, VarMeta(VarMeta::LET_T, m_vlets.size() - 1)});
}
void Energy_Parser::register_invariant(Invariants p){ m_invs.push_back(p); }
void Energy_Parser::register_math_functions(){
    for(const auto& v: getMathFunctionList())
        m_symbol_map.insert({v.first, VarMeta(VarMeta::FUNC_T, static_cast<std::size_t>(v.second.operation))}); 
}
std::size_t Energy_Parser::find_or_register_invariant(Invariants p){
    if (p.type == Invariants::Type::Ifs_T && p.v1 > p.v2) std::swap(p.v1, p.v2);
    for (std::size_t i = 0; i < m_invs.size(); ++i){
        if (p.type == m_invs[i].type){
            if (p.type == Invariants::Type::I1_T || p.type == Invariants::Type::I2_T || p.type == Invariants::Type::J_T)
                return i;
            if (p.type == Invariants::Type::If_T && p.v1 == m_invs[i].v1) return i;
            if (p.type == Invariants::Type::Ifs_T && p.v1 == m_invs[i].v1 && p.v2 == m_invs[i].v2) return i;    
        }
    }
    m_invs.push_back(p);
    return m_invs.size()-1;
}
std::size_t Energy_Parser::find_fiber(const std::string& name) const { 
    auto it = std::find_if(m_fibers.begin(), m_fibers.end(), [&name](const Fiber& f){ return f.name == name; }); 
    return it != m_fibers.end() ? static_cast<std::size_t>(it - m_fibers.begin()) : std::numeric_limits<std::size_t>::max();
}

Energy_Parser::NodeID Energy_Parser::GeneralSyntaxError(std::size_t operand_ilex, std::string prenum_text, std::string afternum_text){
    if (operand_ilex >= m_lex.size()) operand_ilex = m_lex.size()-1;
    m_ErrorStatus = -1;
    auto info = computePositionInfo(m_lex[operand_ilex].pos_st);
    m_err_msg = std::string("Syntax Error: ") + prenum_text + " at " + std::to_string(info.iline+1) + ":" + std::to_string(info.pos - info.line_st + 1) + " " + afternum_text + "\n"; 
    m_err_msg += std::string(m_code.data() + info.part_st, m_code.data() + info.part_end); if (m_err_msg.back() != '\n') m_err_msg += '\n';
    m_err_msg += std::string(info.pos - info.part_st, ' ') + "^\n";

    return NodeID();
}

Energy_Parser::NodeID Energy_Parser::NotFoundAfterOperandError(std::size_t operand_ilex, std::string expected_object){
    if (operand_ilex >= m_lex.size()) operand_ilex = m_lex.size()-1;
    m_ErrorStatus = -2;
    auto info = computePositionInfo(m_lex[operand_ilex].pos_st);
    std::string op_name = std::string(m_code.data() + m_lex[operand_ilex].pos_st, m_code.data() + m_lex[operand_ilex].pos_end);
    m_err_msg = std::string("Syntax Error: After \"") + op_name + "\" operand at " + std::to_string(info.iline+1) + ":" + std::to_string(info.pos - info.line_st + 1) + " expected " + expected_object + ", but not found\n"; 
    m_err_msg += std::string(m_code.data() + info.part_st, m_code.data() + info.part_end); if (m_err_msg.back() != '\n') m_err_msg += '\n';
    m_err_msg += std::string(info.pos - info.part_st, ' ') + "^\n";

    return NodeID();
}
Energy_Parser::NodeID Energy_Parser::NotFoundOperandError(std::size_t operand_ilex, std::string expected_operand){
    if (operand_ilex >= m_lex.size()) operand_ilex = m_lex.size()-1;
    m_ErrorStatus = -3;
    auto info = computePositionInfo(m_lex[operand_ilex].pos_st);
    m_err_msg = std::string("Syntax Error: Expected " + expected_operand + " at ") + std::to_string(info.iline+1) + ":" + std::to_string(info.pos - info.line_st + 1) + " , but not found\n"; 
    m_err_msg += std::string(m_code.data() + info.part_st, m_code.data() + info.part_end); if (m_err_msg.back() != '\n') m_err_msg += '\n';
    m_err_msg += std::string(info.pos - info.part_st, ' ') + "^\n";
    return NodeID();
}

Energy_Parser::NodeID Energy_Parser::FacedUnknownSymbolError(std::size_t ilex, const std::string& object_type){
    if (ilex >= m_lex.size()) ilex = m_lex.size()-1;
    m_ErrorStatus = -4;
    std::string symbol = std::string(m_code.data() + m_lex[ilex].pos_st, m_code.data() + m_lex[ilex].pos_end);
    auto info = computePositionInfo(m_lex[ilex].pos_st);
    m_err_msg = "Syntax Error: Faced " + object_type + " \"" + symbol + "\" at " + std::to_string(info.iline+1) + ":" + std::to_string(info.pos - info.line_st + 1) + ": \n";
    m_err_msg += std::string(m_code.data() + info.part_st, m_code.data() + info.part_end);
    if (m_err_msg.back() != '\n') m_err_msg += '\n';
    m_err_msg += std::string(info.pos - info.part_st, ' ') + "^\n";

    return NodeID();
}

Energy_Parser::NodeID Energy_Parser::RedefinitionSymbolError(std::size_t operand_ilex){
    if (operand_ilex >= m_lex.size()) operand_ilex = m_lex.size()-1;
    m_ErrorStatus = -5;
    auto info = computePositionInfo(m_lex[operand_ilex].pos_st); 
    auto symbol = std::string(m_code.data() + m_lex[operand_ilex].pos_st, m_code.data() + m_lex[operand_ilex].pos_end);
    m_err_msg = std::string("Syntax Error: Redefinition at ") + std::to_string(info.iline+1) + ":" + std::to_string(info.pos - info.line_st + 1) + " of symbol \"" + symbol + "\" is not allowed \n"; 
    m_err_msg += std::string(m_code.data() + info.part_st, m_code.data() + info.part_end); if (m_err_msg.back() != '\n') m_err_msg += '\n';
    m_err_msg += std::string(info.pos - info.part_st, ' ') + "^\n";
    return NodeID();
}

Energy_Parser::NodeID Energy_Parser::NotFindPotentialSection(){
    m_ErrorStatus = -6;
    m_err_msg = "Reached end of code but not find Potential section\n";
    return NodeID();
}

Energy_Parser::NodeID Energy_Parser::SomeCodeAfterSection(std::size_t operand_ilex){
    m_ErrorStatus = -6;
    auto info = computePositionInfo(m_lex[operand_ilex].pos_st); 
    m_err_msg = "Potential section should be last section but found lexems after at " + std::to_string(info.iline+1) + ":" + std::to_string(info.pos - info.line_st + 1) + "\n";
    m_err_msg += std::string(m_code.data() + info.part_st, m_code.data() + info.part_end); if (m_err_msg.back() != '\n') m_err_msg += '\n';
    m_err_msg += std::string(info.pos - info.part_st, ' ') + "^\n";
    return NodeID();
}

Energy_Parser::NodeID Energy_Parser::interp_UNARY_OP(std::size_t& ilex, std::initializer_list<std::pair<Lexems, Operations>> op_map, NodeID(Energy_Parser::*next)(std::size_t& ilex), const std::string& next_name){
    if (!testNoErrors() || ilex >= m_lex.size()) return NodeID(); 
    for (auto& pr: op_map){
        auto lx_op = pr.first; 
        auto op = pr.second;
        if (m_lex[ilex].lexem != lx_op) continue;
        auto lilex = ilex;
        auto id = m_expr.createFreeNode();
        m_expr[id].m_val.ilex = ilex;
        m_expr[id].m_val.m_op = op;
        ++ilex;
        NodeID id1 = NodeID();
        if (ilex < m_lex.size())
            id1 = (this->*next)(ilex);
        if (id1 <= 0)
            return NotFoundAfterOperandError(lilex, next_name);   
        m_expr.addLink(id, id1);
        return id;
    }
    return NodeID();
}

Energy_Parser::NodeID Energy_Parser::interp_UNARY_OP_or_OP(std::size_t& ilex, std::initializer_list<std::pair<Lexems, Operations>> op_map, NodeID(Energy_Parser::*next)(std::size_t& ilex), const std::string& next_name){
    auto lilex = ilex;
    auto id = interp_UNARY_OP(ilex, op_map, next, next_name);
    if (!testNoErrors() || id > 0) return id;
    ilex = lilex;
    return (this->*next)(ilex);
}

Energy_Parser::NodeID Energy_Parser::interp_FUSED_BIN_OP(std::size_t& ilex, std::initializer_list<std::pair<Lexems, Operations>> op_map, NodeID(Energy_Parser::*next)(std::size_t& ilex), const std::string& next_name){
    NodeID id1 = (this->*next)(ilex);
    if (!testNoErrors() || id1 <= 0) return NodeID();
    bool finded = true;
    while (ilex < m_lex.size()){
        bool lfind = false;
        for (auto it = op_map.begin(); it != op_map.end() && !lfind; ++it){
            const auto& pr = *it;
            auto lx_op = pr.first; 
            auto op = pr.second;
            if (m_lex[ilex].lexem != lx_op) continue;
            lfind = true;
            auto id = m_expr.createFreeNode();
            m_expr[id].m_val.ilex = ilex;
            m_expr[id].m_val.m_op = op;
            m_expr.addLink(id, id1);
            auto lilex = ilex;
            ++ilex;
            NodeID id2 = (this->*next)(ilex);
            if (!testNoErrors()) return id2;
            if (id2 <= 0) 
                return NotFoundAfterOperandError(lilex, next_name);
            m_expr.addLink(id, id2);
            id1 = id;
        }
        if (!lfind) break;
    }
    return id1;
}

Energy_Parser::NodeID Energy_Parser::interp_BIN_OP(std::size_t& ilex, std::initializer_list<std::pair<Lexems, Operations>> op_map, NodeID(Energy_Parser::*next)(std::size_t& ilex), const std::string& next_name){
    NodeID id1 = (this->*next)(ilex);
    if (!testNoErrors() || id1 <= 0) return NodeID();
    bool finded = true;
    if (ilex < m_lex.size()){
        bool lfind = false;
        for (auto it = op_map.begin(); it != op_map.end() && !lfind; ++it){
            const auto& pr = *it;
            auto lx_op = pr.first; 
            auto op = pr.second;
            if (m_lex[ilex].lexem != lx_op) continue;
            lfind = true;
            auto id = m_expr.createFreeNode();
            m_expr[id].m_val.ilex = ilex;
            m_expr[id].m_val.m_op = op;
            m_expr.addLink(id, id1);
            auto lilex = ilex;
            ++ilex;
            NodeID id2 = (this->*next)(ilex);
            if (!testNoErrors()) return id2;
            if (id2 <= 0) 
                return NotFoundAfterOperandError(lilex, next_name);
            m_expr.addLink(id, id2);
            id1 = id;
        }
    }
    return id1;
}

Energy_Parser::NodeID Energy_Parser::interp_BRACKET(std::size_t& ilex, Lexems left, Lexems right, NodeID(Energy_Parser::*next)(std::size_t& ilex), const std::string& closed_bracket, const std::string& next_name){
    if (!testNoErrors() || ilex >= m_lex.size()) return NodeID();
    if ( m_lex[ilex].lexem == left){
        auto lilex = ilex;
        ++ilex;
        auto id = (this->*next)(ilex);
        if (id <= 0)
            return NotFoundAfterOperandError(lilex, next_name);
        if (m_lex[ilex].lexem != right)
            return NotFoundOperandError(lilex, "closing bracket \"" + closed_bracket + "\"");
        ilex++; 
        return id;
    }
    return NodeID();
}

Energy_Parser::NodeID Energy_Parser::interp_OR_EXPR(std::size_t& ilex, std::initializer_list<NodeID(Energy_Parser::*)(std::size_t& ilex)> nexts){
    if (!testNoErrors() || ilex >= m_lex.size()) return NodeID();
    NodeID id = NodeID();
    auto lilex = ilex;
    for (auto& n: nexts){
        id = (this->*n)(ilex);
        if (!testNoErrors()) return id;
        if (id > 0) return id;
        ilex = lilex;
    }
    return id;
}

Energy_Parser::NodeID Energy_Parser::interp_VAL(std::size_t& ilex){ 
    if (!testNoErrors() || ilex >= m_lex.size() || (m_lex[ilex].lexem != LX_NUM && m_lex[ilex].lexem != LX_INUM))
        return NodeID(); 
    auto id = m_expr.createFreeNode();   
    if (m_lex[ilex].lexem == LX_NUM){
        m_expr[id].value().m_op = W_VAL; 
        m_expr[id].value().m_value.d = m_lex[ilex].value.d;
    }else { 
        m_expr[id].value().m_op = W_IVAL; 
        m_expr[id].value().m_value.l = m_lex[ilex].value.l;
    }
    m_expr[id].value().ilex = ilex++;
    return id;
}

Energy_Parser::NodeID Energy_Parser::interp_ANY_VAR(std::size_t& ilex){
    if (!testNoErrors() || ilex >= m_lex.size() || (m_lex[ilex].lexem != LX_VAR)
        || (ilex+1 < m_lex.size() && m_lex[ilex+1].lexem == LX_OPBRACKET || m_lex[ilex+1].lexem == LX_OPSQRBRACKET)) //FUNCTION | PREDEFINED VAR
        return NodeID(); 
    auto id = m_expr.createFreeNode();    
    m_expr[id].value().m_op = S_ANY_VAR;   
    m_expr[id].value().ilex = ilex++; 

    return id;
}

Energy_Parser::NodeID Energy_Parser::interp_MATH_FUNC(std::size_t& ilex){
    if (!testNoErrors() || ilex+2 >= m_lex.size() || m_lex[ilex].lexem != LX_VAR || m_lex[ilex+1].lexem != LX_OPBRACKET)
        return NodeID();
    std::string symbol = std::string(m_code.data() + m_lex[ilex].pos_st, m_code.data() + m_lex[ilex].pos_end);    
    auto info = getMathFunctionMeta(symbol);
    if (!info)
        return FacedUnknownSymbolError(ilex, "unknown function name");  
    auto id = m_expr.createFreeNode();
    m_expr[id].m_val.ilex = ilex;
    m_expr[id].m_val.m_op = info.operation;
    ilex+=2;
    if (info.supposed_nargs >= 0){
        for (std::size_t i = 0; i < info.supposed_nargs; ++i){
            auto lilex = ilex;
            auto id1 = interp_Expression(ilex);
            if (!testNoErrors()) return id1;
            if (id1 <= 0) 
                return NotFoundOperandError(lilex, "function argument");
            m_expr.addLink(id, id1);
            if (i != info.supposed_nargs - 1){
                if (ilex >= m_lex.size() || m_lex[ilex].lexem != LX_COMMA)
                    return NotFoundOperandError(ilex, "comma \",\" spliting arguments of " + std::to_string(info.supposed_nargs) + "-args function \"" + symbol + "\"");
                ilex++;
            }
        }
        if (ilex >= m_lex.size() || m_lex[ilex].lexem != LX_CLBRACKET)
            return NotFoundOperandError(ilex, "closing bracket \")\" after function arguments");
        ilex++;
    } else {
        auto _ilex= ilex; 
        int cnt = 0;
        while(ilex < m_lex.size() && m_lex[ilex].lexem != LX_CLBRACKET){
            ++cnt;
            auto lilex = ilex;
            auto id1 = interp_Expression(ilex);
            if (!testNoErrors()) return id1;
            if (id1 <= 0) 
                return NotFoundOperandError(ilex, "argument of variadic argument function \"" + symbol + "\"");
            m_expr.addLink(id, id1);
            if (ilex < m_lex.size() && m_lex[ilex].lexem == LX_CLBRACKET) continue;
            if (ilex >= m_lex.size() || m_lex[ilex].lexem != LX_COMMA )
                return NotFoundOperandError(ilex, "comma \",\" spliting function arguments or close bracket \")\"");
            ++ilex;
        }
        if (ilex < m_lex.size() && m_lex[ilex].lexem == LX_CLBRACKET) 
            ++ilex;
    }
    return id;
}

Energy_Parser::NodeID Energy_Parser::interp_PREDEF_VAR(std::size_t& ilex){
    if (!testNoErrors() || ilex+2 >= m_lex.size() || m_lex[ilex].lexem != LX_VAR || m_lex[ilex+1].lexem != LX_OPSQRBRACKET)
        return NodeID();
    auto symbol_ilex = ilex;    
    std::string symbol = std::string(m_code.data() + m_lex[ilex].pos_st, m_code.data() + m_lex[ilex].pos_end);
    struct Index{
        enum Type {
            FIB_NAME,
            INT_VAL,
        };
        Type m_type = INT_VAL;
        std::string name;
        int val = 0;
        std::size_t ilex;
    };
    std::vector<Index> ids;
    ilex += 2;
    bool find_end = false;
    while (ilex < m_lex.size() && !find_end){
        if (m_lex[ilex].lexem == LX_CLSQRBRACKET) {find_end = true; break;}
        auto faced = std::string(m_code.data() + m_lex[ilex].pos_st, m_code.data() + m_lex[ilex].pos_end);
        if (m_lex[ilex].lexem != LX_INUM && m_lex[ilex].lexem != LX_VAR)
            return GeneralSyntaxError(ilex, "For predefined values inside brackets [] allowed only integers and fiber names", ", but faced \"" + faced + "\"");
        Index ix;
        if (m_lex[ilex].lexem == LX_INUM) { 
            ix.val = m_lex[ilex].value.l;
            ix.m_type = Index::INT_VAL;
        }
        else {
            ix.name = faced;
            ix.m_type = Index::FIB_NAME;
        }
        ix.ilex = ilex;
        ids.emplace_back(ix);
        ++ilex;
        if (ilex < m_lex.size() && m_lex[ilex].lexem != LX_CLSQRBRACKET){
            if (m_lex[ilex].lexem == LX_COMMA)
                ilex++; 
            else 
                return NotFoundOperandError(ilex, "comma \",\"");   
        }   
    }
    if (!find_end)
        return NotFoundOperandError(ilex, "closing bracket \"]\"");
    else 
        ++ilex;    
    if (symbol == "I"){
        if (ids.size() == 1 || (ids.size() == 2 && ids[0].m_type == Index::FIB_NAME && ids[1].m_type == Index::FIB_NAME && ids[0].name == ids[1].name)){
            if (ids[0].m_type == Index::INT_VAL && (ids[0].val == 1 || ids[0].val == 2)){
                Invariants p; p.type = (ids[0].val == 1) ? Invariants::I1_T : Invariants::I2_T;
                auto id = m_expr.createFreeNode();
                m_expr[id].m_val.ilex = symbol_ilex;
                m_expr[id].m_val.m_op = W_VINIT;
                m_expr[id].m_val.m_value.l = find_or_register_invariant(p);
                return id;
            }
            if (ids[0].m_type == Index::FIB_NAME){
                Invariants p; p.type = Invariants::If_T;
                auto fib_id = find_fiber(ids[0].name);
                if (fib_id == std::numeric_limits<std::size_t>::max())
                    return FacedUnknownSymbolError(symbol_ilex+2, "undefined fiber with name = \"" + ids[0].name + "\"");
                p.v1 = fib_id;
                auto id = m_expr.createFreeNode();
                m_expr[id].m_val.ilex = symbol_ilex;
                m_expr[id].m_val.m_op = W_VINIT;
                m_expr[id].m_val.m_value.l = find_or_register_invariant(p);
                return id;
            }
        } else if (ids.size() == 2 && ids[0].m_type == Index::FIB_NAME && ids[1].m_type == Index::FIB_NAME){
            auto fib_id1 = find_fiber(ids[0].name);
            auto fib_id2 = find_fiber(ids[1].name);
            std::size_t iilex = 0;
            if (fib_id1 == std::numeric_limits<std::size_t>::max()) iilex = symbol_ilex + 2;
            else if (fib_id2 == std::numeric_limits<std::size_t>::max()) iilex = symbol_ilex + 4;
            if (iilex > 0)
                return FacedUnknownSymbolError(iilex, "undefined fiber with name = \"" + ids[0].name + "\"");
            Invariants p; p.type = Invariants::Ifs_T;
            p.v1 = std::min(fib_id1, fib_id2), p.v2 = std::max(fib_id1, fib_id2);
            auto id = m_expr.createFreeNode();
            m_expr[id].m_val.ilex = symbol_ilex;
            m_expr[id].m_val.m_op = W_VINIT;
            m_expr[id].m_val.m_value.l = find_or_register_invariant(p);
            return id;
        }
        return GeneralSyntaxError(symbol_ilex, "For I predefined variable", "expected syntax I[1], I[2], I[f] or I[f,s]");
    } else if (symbol == "J"){
        if (ids.empty()){
            Invariants p; p.type = Invariants::J_T;
            auto id = m_expr.createFreeNode();
            m_expr[id].m_val.ilex = symbol_ilex;
            m_expr[id].m_val.m_op = W_VINIT;
            m_expr[id].m_val.m_value.l = find_or_register_invariant(p);
            return id;
        } else 
            return GeneralSyntaxError(symbol_ilex, "For J predefined variable", "expected empty braces []");
    } else 
        return FacedUnknownSymbolError(ilex, "unknown name \"" + symbol + "\" of predefined variable");

    return NodeID();
}
Energy_Parser::NodeID Energy_Parser::interp_ExprInBrackets(std::size_t& ilex){
    return interp_BRACKET(ilex, LX_OPBRACKET, LX_CLBRACKET, &Energy_Parser::interp_Expression, ")");
}
Energy_Parser::NodeID Energy_Parser::interp_Simple(std::size_t& ilex){
    return interp_OR_EXPR(ilex, {&Energy_Parser::interp_MATH_FUNC, &Energy_Parser::interp_ANY_VAR, &Energy_Parser::interp_VAL, &Energy_Parser::interp_PREDEF_VAR});
}
Energy_Parser::NodeID Energy_Parser::interp_UMINUS(std::size_t& ilex){
    return interp_UNARY_OP_or_OP(ilex, {{LX_MINUS, OP_UMINUS}}, &Energy_Parser::interp_Simple, "");
}
Energy_Parser::NodeID Energy_Parser::interp_BracketOP(std::size_t& ilex){
    return interp_OR_EXPR(ilex, {&Energy_Parser::interp_UMINUS, &Energy_Parser::interp_ExprInBrackets});
}
Energy_Parser::NodeID Energy_Parser::interp_POW(std::size_t& ilex){
    return interp_BIN_OP(ilex, {{LX_POW, OP_POW}}, &Energy_Parser::interp_BracketOP);
}
Energy_Parser::NodeID Energy_Parser::interp_MUL_DIV(std::size_t& ilex){
    return interp_FUSED_BIN_OP(ilex, {{LX_MUL, OP_MUL}, {LX_DIV, OP_DIV}}, &Energy_Parser::interp_POW);
}
// Energy_Parser::NodeID Energy_Parser::interp_PLUS_MINUS(std::size_t& ilex){
//     return interp_FUSED_BIN_OP(ilex, {{LX_PLUS, OP_PLUS}, {LX_MINUS, OP_MINUS}}, &Energy_Parser::interp_MUL_DIV);
// }

Energy_Parser::NodeID Energy_Parser::interp_PLUS_MINUS(std::size_t& ilex){
    NodeID id1 = interp_MUL_DIV(ilex); 
    if (!testNoErrors() || id1 <= 0) return id1;
    if (ilex < m_lex.size() && (m_lex[ilex].lexem == LX_PLUS || m_lex[ilex].lexem == LX_MINUS)){
        auto id = m_expr.createFreeNode();
        m_expr[id].m_val.ilex = ilex;
        m_expr[id].m_val.m_op = OP_PLUS;
        m_expr.addLink(id, id1);

        while(ilex < m_lex.size() && (m_lex[ilex].lexem == LX_PLUS || m_lex[ilex].lexem == LX_MINUS)){
            bool is_plus = m_lex[ilex].lexem == LX_PLUS;
            auto lilex = ilex;
            ++ilex;
            auto id2 = interp_MUL_DIV(ilex); 
            if (!testNoErrors()) return id2;
            if (id2 <= 0) 
                return NotFoundAfterOperandError(lilex, "expression");
            if (is_plus) 
                m_expr.addLink(id, id2);
            else {
                auto id3 = m_expr.createFreeNode();
                m_expr[id3].m_val.m_op = OP_UMINUS;
                m_expr[id3].m_val.ilex = lilex;
                m_expr.addLink(id3, id2).addLink(id, id3);
            }    
        }
        id1 = id;
    } 
    return id1;
}

Energy_Parser::NodeID Energy_Parser::interp_COMPARE(std::size_t& ilex){
    return interp_FUSED_BIN_OP(ilex, {{LX_LT, OP_LT}, {LX_LE, OP_LE}, {LX_NE, OP_NE}, {LX_EQ, OP_EQ}, {LX_GT, OP_GT}, {LX_GE, OP_GE}}, &Energy_Parser::interp_PLUS_MINUS);
}
Energy_Parser::NodeID Energy_Parser::interp_NOT(std::size_t& ilex){
    return interp_UNARY_OP_or_OP(ilex, {{LX_NOT, OP_NOT}}, &Energy_Parser::interp_COMPARE);
}
Energy_Parser::NodeID Energy_Parser::interp_AND_OR(std::size_t& ilex){
    return interp_FUSED_BIN_OP(ilex, {{LX_AND, OP_AND}, {LX_OR, OP_OR}}, &Energy_Parser::interp_NOT, "");
}
Energy_Parser::NodeID Energy_Parser::interp_Expression(std::size_t& ilex){
    return interp_AND_OR(ilex);
}

Energy_Parser::NodeID Energy_Parser::interp_LetAssign(std::size_t& ilex){
    if (!testNoErrors() || ilex >= m_lex.size() || m_lex[ilex].lexem != LX_VAR) return NodeID();
    auto symbol = std::string(m_code.data() + m_lex[ilex].pos_st, m_code.data() + m_lex[ilex].pos_end);
    if (m_symbol_map.find(symbol) != m_symbol_map.end())
        return RedefinitionSymbolError(ilex);
    ++ilex;
    if (ilex >= m_lex.size() || m_lex[ilex].lexem != LX_ASSIGN)
        return NotFoundOperandError(ilex, "assign \"=\"");
    ++ilex;
    auto id = interp_Expression(ilex);
    if (!testNoErrors() || id <= 0) return id;
    unmap_anon_variables(id);
    if (!testNoErrors()) return NodeID();
    SubVar sv;
    sv.name = symbol;
    sv.id = id;
    register_let_variable(sv);
    return id;
}
Energy_Parser::NodeID Energy_Parser::interp_VarAssign(std::size_t& ilex){
    if (!testNoErrors() || ilex >= m_lex.size() || m_lex[ilex].lexem != LX_VAR) return NodeID();
    auto symbol = std::string(m_code.data() + m_lex[ilex].pos_st, m_code.data() + m_lex[ilex].pos_end);
    if (m_symbol_map.find(symbol) != m_symbol_map.end())
        return RedefinitionSymbolError(ilex);
    ++ilex;
    Parameter p; p.name = symbol;
    if (ilex < m_lex.size() && m_lex[ilex].lexem == LX_COLON){
        ++ilex;
        if (ilex >= m_lex.size() || m_lex[ilex].lexem != LX_STR)
            return NotFoundOperandError(ilex, "string");
        p.def_tag = std::string(m_code.data() + m_lex[ilex].pos_st + 1, m_code.data() + m_lex[ilex].pos_end - 1);
        ++ilex;    
    } 
    if (ilex < m_lex.size() && m_lex[ilex].lexem == LX_ASSIGN){
        auto lilex = ilex;
        ++ilex;
        auto id = interp_Expression(ilex);
        if (!testNoErrors()) return NodeID();
        if (id <= 0)
            return NotFoundAfterOperandError(lilex, "trivial expression");
        p.def_val = evaluate_trivial_expression(lilex, id);
        if (!testNoErrors()) return NodeID();
        m_expr.removeBranch(id);    
    }
    register_param(p);
    return NodeID(1); 
}
Energy_Parser::NodeID Energy_Parser::interp_FiberAssign(std::size_t& ilex){
    if (!testNoErrors() || ilex >= m_lex.size() || m_lex[ilex].lexem != LX_VAR) return NodeID();
    auto symbol = std::string(m_code.data() + m_lex[ilex].pos_st, m_code.data() + m_lex[ilex].pos_end);
    if (m_symbol_map.find(symbol) != m_symbol_map.end())
        return RedefinitionSymbolError(ilex);
    ++ilex;
    Fiber fib; fib.name = symbol;
    if (ilex < m_lex.size() && m_lex[ilex].lexem == LX_COLON){
        ++ilex;
        if (ilex >= m_lex.size() || m_lex[ilex].lexem != LX_STR)
            return NotFoundOperandError(ilex, "string");
        fib.def_tag = std::string(m_code.data() + m_lex[ilex].pos_st + 1, m_code.data() + m_lex[ilex].pos_end - 1);
        ++ilex;    
    }
    register_fiber(fib);
    return NodeID(1);   
}

Energy_Parser::NodeID Energy_Parser::interp_LISTED_OP(std::size_t& ilex, Lexems lx_st, Lexems lx_sep, Lexems lx_end, NodeID(Energy_Parser::*next)(std::size_t& ilex), const std::string& end_seq, const std::string& next_name){
    if (!testNoErrors() || ilex >= m_lex.size() || m_lex[ilex].lexem != lx_st) return NodeID();
    bool interpret = true;
    auto lilex = ilex;
    ++ilex;
    auto id = (this->*next)(ilex);
    if (!testNoErrors()) return NodeID();
    if (id > 0){
        while(ilex < m_lex.size() && m_lex[ilex].lexem == lx_sep){
            ilex++;
            auto plex = std::min(ilex, m_lex.size()-1);
            if (ilex < m_lex.size())
                id = (this->*next)(ilex);
            if (!testNoErrors()) return NodeID();
            if (id <= 0)
                return NotFoundOperandError(ilex, next_name);
        }
    }
    if (ilex >= m_lex.size() || m_lex[ilex].lexem != lx_end)
        return NotFoundOperandError(std::min(ilex, m_lex.size()-1), end_seq);
    ++ilex;    
    return NodeID(1);
}

Energy_Parser::NodeID Energy_Parser::interp_Let(std::size_t& ilex){
    return interp_LISTED_OP(ilex, KW_Let, LX_COMMA, LX_SEMICOLON, &Energy_Parser::interp_LetAssign, "semicolon \";\"" );
}
Energy_Parser::NodeID Energy_Parser::interp_Fiber(std::size_t& ilex){
    return interp_LISTED_OP(ilex, KW_Fiber, LX_COMMA, LX_SEMICOLON, &Energy_Parser::interp_FiberAssign, "semicolon \";\"" );
}
Energy_Parser::NodeID Energy_Parser::interp_Var(std::size_t& ilex){
    return interp_LISTED_OP(ilex, KW_Var, LX_COMMA, LX_SEMICOLON, &Energy_Parser::interp_VarAssign, "semicolon \";\"" );
}
Energy_Parser::NodeID Energy_Parser::interp_Potential(std::size_t& ilex){
    Lexems lx_st = KW_Potential, lx_end = LX_SEMICOLON;
    if (!testNoErrors() || ilex >= m_lex.size() || m_lex[ilex].lexem != lx_st) return NodeID();
    auto lilex = ilex;
    ++ilex;
    if (ilex >= m_lex.size() || m_lex[ilex].lexem != LX_ASSIGN)
        return NotFoundOperandError(ilex, "assign \"=\"");
    ++ilex;    
    auto id = interp_Expression(ilex);
    if (!testNoErrors()) return NodeID();
    if (id <= 0)
        return NotFoundAfterOperandError(lilex, "expression");
    if (ilex >= m_lex.size() || m_lex[ilex].lexem != lx_end)
        return NotFoundOperandError(ilex, "semicolon \";\"");
    ++ilex;
    unmap_anon_variables(id); 
    if (!testNoErrors()) return NodeID();    
    return id;       
}

std::pair<Energy_Parser::TrivialValue, bool> Energy_Parser::getTrivialValue(NodeID subtree) const{
    auto& n = m_expr[subtree];
    auto op = n.m_val.m_op;
    switch (op){
        case W_VAL: return {TrivialValue(n.m_val.m_value.d), true};
        case W_IVAL: return {TrivialValue(n.m_val.m_value.l), true};
        case W_VLET: return getTrivialValue(m_vlets[n.m_val.m_value.l].id);
    }
    return {TrivialValue(), false};
}
bool Energy_Parser::simplify_subexpr_iter(NodeID subtree){
    auto& n = m_expr[subtree];
    auto op = n.m_val.m_op;
    auto apply_unary_func = [&](auto func_int, auto func_real)->bool{
        if (!n.have_child()) return false;
        bool changed = simplify_subexpr_iter(n.child(0));
        auto v = getTrivialValue(n.child(0));
        if (v.second){
            if (v.first.m_type == TrivialValue::IVAL){
                n.m_val.set_value( func_int(v.first.Integer()) ); 
            } else {
                n.m_val.set_value( func_real(v.first.get<double>()) );
            }
            auto ch = n.child(0);
            m_expr.removeLink(subtree, ch);
            return true;
        } else 
            return changed;
    };
    auto apply_bin_func = [&](auto func_int, auto func_real)->bool{
        if (n.number_of_childs() != 2) return false;
        bool changed = simplify_subexpr_iter(n.child(0));
        changed |= simplify_subexpr_iter(n.child(1));
        auto v1 = getTrivialValue(n.child(0)), v2 = getTrivialValue(n.child(1));
        if (v1.second && v2.second){
            if (v1.first.m_type == TrivialValue::IVAL && v2.first.m_type == TrivialValue::IVAL){
                n.m_val.set_value( func_int(v1.first.Integer(), v2.first.Integer()) );
            } else {
                n.m_val.set_value( func_real(v1.first.get<double>(), v2.first.get<double>()) );
            }
            auto ch1 = n.child(0), ch2 = n.child(1);
            m_expr.removeLink(subtree, ch1).removeLink(subtree, ch2);
            return true;
        } else 
            return changed;
    };
    auto apply_narg_sym = [&](auto func_int, auto func_real)->bool{
        if (n.number_of_childs() == 0) return false;
        if (n.number_of_childs() == 1) {
            auto ch = n.child(0);
            m_expr.removeLink(subtree, ch);
            n.m_val = m_expr[ch].m_val;
            return true;    
        }
        bool changed = false;
        for (std::size_t i = 0; i < n.number_of_childs(); ++i)
            changed |= simplify_subexpr_iter(n.child(i));
        std::vector<NodeID> not_eval;
        NodeID some_eval = NodeID();
        std::vector<TrivialValue> eval;
        for (std::size_t i = 0; i < n.number_of_childs(); ++i){
            auto ch = n.child(i);
            auto v = getTrivialValue(ch);
            if (v.second){
                if (eval.empty()) some_eval = ch;
                eval.push_back(v.first);
            } else 
                not_eval.push_back(ch);
        }
        if (eval.size() <= 1) return changed;
        changed = true;
        bool is_int = true;
        for (auto& v: eval)
            is_int &= (v.m_type == TrivialValue::IVAL);
        if (is_int){
            std::vector<long> to_eval(eval.size());
            for (std::size_t i = 0; i < eval.size(); ++i)
                to_eval[i] = eval[i].Integer();
            m_expr[some_eval].m_val.set_value(func_int(to_eval));  
        } else {
            std::vector<double> to_eval(eval.size());
            for (std::size_t i = 0; i < eval.size(); ++i)
                to_eval[i] = eval[i].get<double>();
            m_expr[some_eval].m_val.set_value(func_real(to_eval));   
        }
        auto _childs = n.childs();
        for (auto ch: _childs)
            m_expr.removeLink(subtree, ch);
        for (auto& ch: not_eval) m_expr.addLink(subtree, ch); 
        m_expr.addLink(subtree, some_eval);
        return changed;   
    };
    auto perform_aggregation = [&](){
        if (!n.have_child()) return false;
        std::vector<NodeID> plus_ch;
        for (auto ch: n.childs())
            if (m_expr[ch].m_val.m_op == n.m_val.m_op)
                plus_ch.push_back(ch);
        if (!plus_ch.empty()){        
            for (auto ch: plus_ch) m_expr.removeLink(subtree, ch);
            for (auto v: plus_ch) for (auto ch: m_expr[v].childs()) m_expr.addLink(subtree, ch);
            return true; 
        }
        return false;
    };
    #define APPLY_UNARY(FUNC) apply_unary_func([](long v){ return FUNC(v); }, [](double v){ return FUNC(v); });
    #define APPLY_BINARY(FUNC) apply_bin_func([](long v1, long v2){ return FUNC(v1, v2); }, [](double v1, double v2){ return FUNC(v1, v2); });
    switch (op){
        case S_ANY_VAR: return false;
        case W_VAL: return false;
        case W_IVAL: return false;
        case W_VINIT: return false;
        case W_PAR: return false;
        case W_VLET: return false;
        case H_HEAD: return simplify_subexpr_iter(n.child(0));
        case OP_ABS: return APPLY_UNARY(abs);
        case OP_EXP: return APPLY_UNARY(exp);
        case OP_EXPM1: return APPLY_UNARY(expm1);
        case OP_LOG: return APPLY_UNARY(log);
        case OP_LOG1P: return APPLY_UNARY(log1p);
        case OP_SQRT: return APPLY_UNARY(sqrt);
        case OP_CBRT: return APPLY_UNARY(cbrt);
        case OP_SQ: return apply_unary_func([](long v){ return v*v; }, [](double v){ return v*v; });
        case OP_CUBE: return apply_unary_func([](long v){ return v*v*v; }, [](double v){ return v*v*v; });
        case OP_SIN: return APPLY_UNARY(sin);
        case OP_COS: return APPLY_UNARY(cos);
        case OP_TAN: return APPLY_UNARY(tan);
        case OP_ASIN: return APPLY_UNARY(asin);
        case OP_ACOS: return APPLY_UNARY(acos);
        case OP_ATAN: return APPLY_UNARY(atan);
        case OP_SINH: return APPLY_UNARY(sinh);
        case OP_COSH: return APPLY_UNARY(cosh);
        case OP_TANH: return APPLY_UNARY(tanh);
        case OP_ERF: return APPLY_UNARY(erf);
        case OP_ERFC: return APPLY_UNARY(erfc);
        case OP_TGAMMA: return APPLY_UNARY(tgamma);
        case OP_LGAMMA: return APPLY_UNARY(lgamma);
        case OP_FLOOR: return apply_unary_func([](long v){ return v; }, [](double v){ return floor(v); });
        case OP_CEIL: return apply_unary_func([](long v){ return v; }, [](double v){ return ceil(v); });
        case OP_SIGN: return apply_unary_func([](long v){ return (v > 0) - (v < 0); }, [](double v){ return (v > 0) - (v < 0); });
        case OP_NOT: return apply_unary_func([](long v){ return v == 0; }, [](double v){ return v == 0; });
        case OP_UMINUS: return apply_unary_func([](long v){ return -v; }, [](double v){ return -v; });
        case OP_POW: {
            bool changed = false;
            if (n.number_of_childs() != 2) return false;
            {
                changed = simplify_subexpr_iter(n.child(0));
                changed |= simplify_subexpr_iter(n.child(1));
                auto v1 = getTrivialValue(n.child(0)), v2 = getTrivialValue(n.child(1));
                if (v1.second && v2.second){
                    if (v2.first.m_type == TrivialValue::IVAL){
                        if (v1.first.m_type == TrivialValue::IVAL)
                            n.m_val.set_value( pow(v1.first.Integer(), v2.first.Integer()) );
                        else    
                            n.m_val.set_value( pow(v1.first.Real(), v2.first.Integer()) );    
                    } else
                        n.m_val.set_value( pow(v1.first.get<double>(), v2.first.get<double>()) );
                    auto ch1 = n.child(0), ch2 = n.child(1);
                    m_expr.removeLink(subtree, ch1).removeLink(subtree, ch2);
                    changed = true;
                }
            }
            // auto changed = APPLY_BINARY(pow);
            if (!changed && n.number_of_childs() == 2){
                auto v1 = getTrivialValue(n.child(0)), v2 = getTrivialValue(n.child(1));
                if (v1.first.m_type == TrivialValue::IVAL && v1.first.Integer() == 0){
                    auto _childs = n.childs();
                    for (auto ch: _childs)
                        m_expr.removeLink(subtree, ch);
                    n.m_val.set_value(0L); 
                    changed = true;   
                } else if (v2.first.m_type == TrivialValue::IVAL){
                    if (v2.first.Integer() == 0){
                        auto _childs = n.childs();
                        for (auto ch: _childs)
                            m_expr.removeLink(subtree, ch);
                        n.m_val.set_value(1L);
                        changed = true;
                    } else if (v2.first.Integer() == 2) {
                        m_expr.removeLink(subtree, n.child(1));
                        n.m_val.m_op = OP_SQ;
                        changed = true;
                    } else if (v2.first.Integer() == 3) {
                        m_expr.removeLink(subtree, n.child(1));
                        n.m_val.m_op = OP_CUBE;
                        changed = true;
                    }
                }
            }
            return changed;
        }
        case OP_ATAN2: return APPLY_BINARY(atan2);
        case OP_FMOD: return apply_bin_func([](long v1, long v2){ return v1 % v2; }, [](double v1, double v2){ return fmod(v1, v2); });
        case OP_LT: return apply_bin_func([](long v1, long v2){ return v1 < v2; }, [](double v1, double v2){ return v1 < v2; });
        case OP_LE: return apply_bin_func([](long v1, long v2){ return v1 <= v2; }, [](double v1, double v2){ return v1 <= v2; });
        case OP_EQ: return apply_bin_func([](long v1, long v2){ return v1 == v2; }, [](double v1, double v2){ return v1 == v2; });
        case OP_NE: return apply_bin_func([](long v1, long v2){ return v1 != v2; }, [](double v1, double v2){ return v1 != v2; });
        case OP_GT: return apply_bin_func([](long v1, long v2){ return v1 > v2; }, [](double v1, double v2){ return v1 > v2; });
        case OP_GE: return apply_bin_func([](long v1, long v2){ return v1 >= v2; }, [](double v1, double v2){ return v1 >= v2; });
        case OP_IFELSE:{
            if (n.number_of_childs() != 3) return false;
            auto r0 = simplify_subexpr_iter(n.child(0));
            auto r1 = simplify_subexpr_iter(n.child(1));
            auto r2 = simplify_subexpr_iter(n.child(2));
            bool changed = r0 || r1 || r2;
            auto v_cond = getTrivialValue(n.child(0));
            if (!v_cond.second) 
                return changed;
            auto ch = n.child(1);
            if (v_cond.first.get<long>() == 0)
                ch = n.child(2);
            auto _childs = n.childs();
            for (auto c: _childs) m_expr.removeLink(subtree, c);
            n.m_val = m_expr[ch].m_val;
            for (auto c: m_expr[ch].childs()) m_expr.addLink(subtree, c);
            return true;
        }
        case OP_AND:{
            if (n.number_of_childs() == 0) return false;
            bool changed = perform_aggregation();
            for (auto ch: n.childs())
                changed |= simplify_subexpr_iter(ch);  
            std::vector<NodeID> true_ch;
            for (auto ch: n.childs()){
                auto v = getTrivialValue(ch);
                if (!v.second) continue;
                if (v.first.get<bool>() == false){
                    n.m_val.set_value(0L);
                    return true;
                } else {
                    true_ch.push_back(ch);
                }
            }
            if (true_ch.empty()) return changed;
            if (true_ch.size() == n.number_of_childs()){
                n.m_val.set_value(1L);
                return true;
            }
            for (auto ch: true_ch)
                m_expr.removeLink(subtree, ch);  
            if (n.number_of_childs() == 1){
                n.m_val.m_op = OP_NE;
                auto id = m_expr.createFreeNode();
                m_expr[id].m_val.set_value(0L);
                m_expr[id].m_val.ilex = m_expr[n.child(0)].m_val.ilex;
                m_expr.addLink(subtree, id);
            }  
            return true;  
        }
        case OP_OR:{
            if (n.number_of_childs() == 0) return false;
            bool changed = perform_aggregation();
            for (auto ch: n.childs())
                changed |= simplify_subexpr_iter(ch);  
            std::vector<NodeID> false_ch;
            for (auto ch: n.childs()){
                auto v = getTrivialValue(ch);
                if (!v.second) continue;
                if (v.first.get<bool>() == true){
                    n.m_val.set_value(1L);
                    return true;
                } else {
                    false_ch.push_back(ch);
                }
            }
            if (false_ch.empty()) return changed;
            if (false_ch.size() == n.number_of_childs()){
                n.m_val.set_value(0L);
                return true;
            }
            for (auto ch: false_ch)
                m_expr.removeLink(subtree, ch);  
            if (n.number_of_childs() == 1){
                n.m_val.m_op = OP_NE;
                auto id = m_expr.createFreeNode();
                m_expr[id].m_val.set_value(0L);
                m_expr[id].m_val.ilex = m_expr[n.child(0)].m_val.ilex;
                m_expr.addLink(subtree, id);
            }  
            return true;  
        }
        case OP_PLUS: {
            bool changed = perform_aggregation();      
            auto res = apply_narg_sym([](const std::vector<long>& v){ return std::accumulate(v.begin(), v.end(), 0L); }, [](const std::vector<double>& v){ return std::accumulate(v.begin(), v.end(), 0.0); });
            return changed || res;
        }
        case OP_MUL: {
            bool changed = perform_aggregation();
            auto res = apply_narg_sym([](const std::vector<long>& v){ return std::accumulate(v.begin(), v.end(), 1L, std::multiplies<>()); }, [](const std::vector<double>& v){ return std::accumulate(v.begin(), v.end(), 1.0, std::multiplies<>()); });
            return changed || res;
        }
        case OP_MIN_N: {
            bool changed = perform_aggregation();
            auto res = apply_narg_sym([](const std::vector<long>& v){ return *std::min_element(v.begin(), v.end()); }, [](const std::vector<double>& v){ return *std::min_element(v.begin(), v.end()); });
            return changed || res;
        }
        case OP_MAX_N: {
            bool changed = perform_aggregation();
            auto res =  apply_narg_sym([](const std::vector<long>& v){ return *std::max_element(v.begin(), v.end()); }, [](const std::vector<double>& v){ return *std::max_element(v.begin(), v.end()); });
            return changed || res;
        }
        case OP_HYPOT_N:{
            bool changed = perform_aggregation();
            auto hypot_x = [](const auto& arr)->double{
                if (arr.size() == 0) return 0.0;
                if (arr.size() == 1) return static_cast<double>(arr[0]);
                if (arr.size() == 2) return std::hypot(arr[0], arr[1]);
                if (arr.size() == 3) return std::hypot(arr[0], arr[1], arr[2]);
                auto m = *std::max_element(arr.begin(), arr.end());
                auto v = 0.0;
                for (auto i: arr)
                    v += (i/m)*(i/m);
                return m*sqrt(v);
            };
            auto res =  apply_narg_sym([&hypot_x](const std::vector<long>& v){ return hypot_x(v); }, [&hypot_x](const std::vector<double>& v){ return hypot_x(v); });
            return changed || res; 
        }
        case OP_MINUS:{
            if (n.number_of_childs() == 0) return false;
            if (n.number_of_childs() == 1) {
                auto ch = n.child(0);
                m_expr.removeLink(subtree, ch);
                n.m_val = m_expr[ch].m_val;
                for (auto c: m_expr[ch].childs()) m_expr.addLink(subtree, c);
                return true;
            }
            bool changed = false;
            for (auto ch: n.childs())
                changed |= simplify_subexpr_iter(ch);
            if (n.number_of_childs() > 2) return changed;
            auto res =  apply_bin_func([](long v1, long v2){ return v1 - v2; }, [](double v1, double v2){ return v1 - v2; });
            return changed || res;
        }
        case OP_DIV:{
            if (n.number_of_childs() != 2) return false;
            bool changed = false;
            for (auto ch: n.childs())
                changed |= simplify_subexpr_iter(ch);
            auto v1 = getTrivialValue(n.child(0)), v2 = getTrivialValue(n.child(1));
            if (v1.second && v2.second){
                if (v1.first.m_type == TrivialValue::IVAL && v2.first.m_type == TrivialValue::IVAL && v1.first.Integer() % v2.first.Integer() == 0){
                    n.m_val.set_value( v1.first.Integer() / v2.first.Integer() );  
                } else {
                    n.m_val.set_value( v1.first.get<double>() / v2.first.get<double>() );
                }
                auto ch1 = n.child(0), ch2 = n.child(1);
                m_expr.removeLink(subtree, ch1).removeLink(subtree, ch2);
                return true;
            } else if (v2.second){
                auto id = m_expr.createFreeNode();
                m_expr[id].m_val = m_expr[n.child(1)].m_val;
                double val = 1.0;
                if (m_expr[n.child(1)].m_val.m_op == W_VAL)
                    val = 1. / m_expr[n.child(1)].m_val.m_value.d;
                else
                    val = 1. / m_expr[n.child(1)].m_val.m_value.l;    
                m_expr[id].m_val.set_value(val);
                m_expr.removeLink(subtree, n.child(1));
                m_expr.addLink(subtree, id);
                n.m_val.m_op = OP_MUL;
                return true;
            } else 
                return changed;               
        }
        default:
            return false;
    }
    #undef APPLY_BINARY
    #undef APPLY_UNARY
}

double Energy_Parser::evaluate_trivial_expression(std::size_t ilex, NodeID subtree){
    simplify_subexpr(subtree);
    auto& v = m_expr[subtree].m_val;
    double res = 0;
    if (v.m_op == W_VAL){
        res = v.m_value.d;
    } else if (v.m_op == W_IVAL){
        res = v.m_value.l;
    } else
        GeneralSyntaxError(ilex, "Expression ", "is not trivial");
    return res;        
}

void Energy_Parser::unmap_anon_variables(NodeID subtree){
    for (auto it = m_expr.begin(subtree); it != m_expr.end(subtree); ++it){
        auto& v = it->m_val;
        if (v.m_op == S_ANY_VAR){
            auto name = std::string(m_code.data() + m_lex[v.ilex].pos_st, m_code.data() + m_lex[v.ilex].pos_end);
            auto it = m_symbol_map.find(name);
            if (it == m_symbol_map.end() || (it->second.m_type != VarMeta::PAR_T && it->second.m_type != VarMeta::LET_T)){
                if (it == m_symbol_map.end()) FacedUnknownSymbolError(v.ilex, "undefined symbol \"" + name + "\"");
                else RedefinitionSymbolError(v.ilex);
                return;
            }
            v.m_op = (it->second.m_type == VarMeta::PAR_T) ? W_PAR : W_VLET;
            v.m_value.l = it->second.id;
        }
    }
}

bool Energy_Parser::performSyntaxAnalysis() noexcept{
    if (m_lex.empty()) return false;
    auto id = m_expr.createFreeNode();
    m_expr[id].m_val.m_op = H_HEAD;
    std::size_t ilex = 0;
    while(interp_OR_EXPR(ilex, {&Energy_Parser::interp_Var, &Energy_Parser::interp_Let, &Energy_Parser::interp_Fiber}) > 0);
    if (!testNoErrors()) return false;
    auto lilex = ilex;
    auto id_p = interp_Potential(ilex);
    if (!testNoErrors()) return false;
    if (id_p <= 0) {
        NotFoundOperandError(lilex, "Potential");
        return false;
    }
    if (ilex != m_lex.size()){
        SomeCodeAfterSection(ilex);
        return false;
    }
    m_expr.addLink(id, id_p);
    return true;
}

Energy_Parser& Energy_Parser::Simplify(){
    for (auto se: m_vlets)
        simplify_subexpr(se.id);
    auto pid = m_expr[0].child(0);    
    simplify_subexpr(pid);
    std::map<NodeID, bool> remote_let_id;
    std::set<NodeID> loc_let_id;
    std::set<NodeID> par_id, pre_id;
    bool add_let_id = false;
    auto analyze_expr = [&](NodeID id){
        for (auto it = m_expr.begin(id); it != m_expr.end(id); ++it){
            if (it->m_val.m_op == W_VINIT)  pre_id.insert(it.getID());
            else if (it->m_val.m_op == W_PAR) par_id.insert(it.getID());
            else if (it->m_val.m_op == W_VLET) {
                loc_let_id.insert(it.getID());
                auto remote_id = m_vlets[it->m_val.m_value.l].id;
                auto lt = remote_let_id.find(remote_id);
                if (lt == remote_let_id.end()) {
                    add_let_id = true;
                    remote_let_id.insert({remote_id, false});
                }
            }
        }
    };
    analyze_expr(0);
    while(add_let_id){
        add_let_id = false;
        for (auto& li: remote_let_id)
            if (!li.second){
                li.second = true;
                analyze_expr(li.first);
            }
    }
    std::vector<NodeID> head_nodes; head_nodes.reserve(remote_let_id.size() + 1);
    head_nodes.push_back(0);
    for (auto& v: remote_let_id) head_nodes.push_back(v.first);
    m_expr.removeUnreacheableNodes(head_nodes.begin(), head_nodes.end());

    //TODO: below uneffective algo
    std::vector<bool> params_present(m_params.size(), false);
    std::vector<bool> invs_present(m_invs.size(), false);
    std::vector<bool> vlets_present(m_vlets.size(), false);
    for(auto n: par_id) params_present[m_expr[n].m_val.m_value.l] = true;
    for(auto n: pre_id) invs_present[m_expr[n].m_val.m_value.l] = true;
    for(auto n: loc_let_id) vlets_present[m_expr[n].m_val.m_value.l] = true;
    auto remove_not_present = [&](const std::vector<bool>& present, auto& repl){
        std::size_t cnt = 0;
        for (std::size_t i = 0; i < present.size(); ++i) cnt += present[i];
        typename std::remove_reference<decltype(repl)>::type res; res.reserve(cnt);
        for (std::size_t i = 0; i < present.size(); ++i) if (present[i]) res.push_back(repl[i]);
        repl = std::move(res);
    };
    auto remove_repetitive = [&](const std::vector<bool>& present, auto arr, auto& repl, auto conv = [](auto i) { return i; })->void{
        for(auto cn: arr) {
            auto n = conv(cn);
            auto l = m_expr[n].m_val.m_value.l;
            long cnt = 0;
            for (long i = 0; i < l; ++i) cnt += present[i];
            m_expr[n].m_val.m_value.l = cnt;
        }
        remove_not_present(present, repl);
    };
    remove_repetitive(params_present, par_id, m_params, [](auto i) { return i; });
    remove_repetitive(invs_present, pre_id, m_invs, [](auto i) { return i; });
    remove_repetitive(vlets_present, loc_let_id, m_vlets, [](auto i){ return i; }); 

    std::vector<bool> fibers_present(m_fibers.size(), false);
    for (auto v: m_invs)
        if (v.type == Invariants::If_T) fibers_present[v.v1] = true;
        else if (v.type == Invariants::Ifs_T) fibers_present[v.v1] = fibers_present[v.v2] = true;
    remove_not_present(fibers_present, m_fibers);
    std::map<std::size_t, std::size_t> fiber_num_map;
    for (std::size_t i = 0, cnt = 0; i < fibers_present.size(); ++i) 
        if (fibers_present[i]) 
            fiber_num_map[i] = cnt++;
    for (auto v: m_invs)
        if (v.type == Invariants::If_T) v.v1 = fiber_num_map[v.v1];
        else if (v.type == Invariants::Ifs_T) v.v1 = fiber_num_map[v.v1], v.v2 = fiber_num_map[v.v2];  

    return *this;          
}

std::ostream& Energy_Parser::print_invariant(std::ostream& out, std::size_t inv_order_num) const {
    auto& p = m_invs[inv_order_num];
    switch(p.type){
        case Invariants::I1_T: return out << "I[1]";
        case Invariants::I2_T: return out << "I[2]";
        case Invariants::J_T: return out << "J[]";
        case Invariants::If_T: return out << "I[" << m_fibers[p.v1].name << "]";
        case Invariants::Ifs_T: return out << "I[" << m_fibers[p.v1].name << "," << m_fibers[p.v2].name << "]";
        default: return out << "I_UNKNOWN";
    }
}

std::ostream& Energy_Parser::print_expression(std::ostream& out, NodeID id, bool not_substitute_let_exprs) const {
    auto& n = m_expr[id];
    auto get_name = [&](const MultiTree<TreeValue>::Node& n) {
        return std::string(m_code.data() + m_lex[n.m_val.ilex].pos_st, m_code.data() + m_lex[n.m_val.ilex].pos_end);
    };
    static const std::map<Operations, int> priority = {
        {S_ANY_VAR, 10},
        {W_VAL, 10},
        {W_VINIT, 10},
        {W_PAR, 10},
        {W_VLET, 10},
        {OP_UMINUS, 0},
        {OP_PLUS, 1},
        {OP_MINUS, 1},
        {OP_MUL, 2},
        {OP_DIV, 2},
        {OP_POW, 3},
        {OP_LT, 4},
        {OP_EQ, 4},
        {OP_NE, 4},
        {OP_LE, 4},
        {OP_GE, 4},
        {OP_GT, 4},
        {OP_NE, 5},
        {OP_AND, 6},
        {OP_OR, 6},
    };
    auto get_priority = [](Operations op){
        auto it = priority.find(op);
        return (it != priority.end()) ? it->second : 10;
    };
    switch (n.m_val.m_op){
        case S_ANY_VAR: return out << "ANY_VAR(" << get_name(n) << ")";
        case W_VLET: {
            if (not_substitute_let_exprs)
                return out << m_vlets[n.m_val.m_value.l].name;
            else {
                out << "( "; print_expression(out, m_vlets[n.m_val.m_value.l].id, not_substitute_let_exprs) << " )";
                return out;
            }    
        }
        case W_VAL: {
            if (n.m_val.m_value.d < 0 && n.have_parent()) out << '(';
            out << n.m_val.m_value.d;
            if (n.m_val.m_value.d < 0 && n.have_parent()) out << ')';
            return out;
        }
        case W_IVAL: {
            if (n.m_val.m_value.l < 0 && n.have_parent()) out << '(';
            out << n.m_val.m_value.l;
            if (n.m_val.m_value.l < 0 && n.have_parent()) out << ')';
            return out;
        }
        case W_VINIT: return print_invariant(out, n.m_val.m_value.l);
        case W_PAR: return out << m_params[n.m_val.m_value.l].name;
        case OP_PLUS: 
        case OP_MINUS:
        {
            auto op_p1 = get_priority(n.m_val.m_op);
            for (std::size_t i = 0; i < n.number_of_childs(); ++i){
                auto ie = n.child(i);
                bool is_plus = (n.m_val.m_op == OP_PLUS);
                if (i != 0 && m_expr[ie].m_val.m_op == OP_UMINUS){
                    is_plus = (is_plus) ? false : true;
                    ie = m_expr[ie].child(0);
                }
                auto& n2 = m_expr[ie];
                if (n2.m_val.m_op == W_VAL){
                    if (n2.m_val.m_value.d < 0) is_plus = !is_plus;
                    out << ' ' << (is_plus ? '+' : '-') << ' ' << abs(n2.m_val.m_value.d);
                } else if (n2.m_val.m_op == W_IVAL){
                    if (n2.m_val.m_value.l < 0) is_plus = !is_plus;
                    out << ' ' << (is_plus ? '+' : '-') << ' ' << abs(n2.m_val.m_value.l);
                } else {
                    auto op_p2 = get_priority(n2.m_val.m_op);
                    if (i != 0)
                        out << ' ' << (is_plus ? '+' : '-') << ' ';
                    if (op_p1 > op_p2) out << "( ";
                    print_expression(out, ie, not_substitute_let_exprs);
                    if (op_p1 > op_p2) out << " )";
                }
            }
            return out;
        }
        case OP_DIV:
        case OP_MUL:
        case OP_AND:
        case OP_OR:
        {
            auto op_p1 = get_priority(n.m_val.m_op);
            for (std::size_t i = 0; i < n.number_of_childs(); ++i){
                auto ie = n.child(i);
                auto& n2 = m_expr[ie];
                auto op_p2 = get_priority(n2.m_val.m_op);
                if (i != 0)
                    out << ' ' << getOperationName(n.m_val.m_op) << ' ';
                bool in_bracket = (op_p1 > op_p2 || (op_p1 == op_p2 && (n.m_val.m_op != n2.m_val.m_op || n.m_val.m_op == OP_DIV)));     
                //bool in_bracket = (op_p1 > op_p2 || (op_p1 == op_p2 && n.m_val.m_op != n2.m_val.m_op && n.m_val.m_op != OP_MUL));    
                if (in_bracket) out << "( ";
                print_expression(out, ie, not_substitute_let_exprs);
                if (in_bracket) out << " )";
            }
            return out;
        }
        case OP_UMINUS:{
            auto ie = n.child(0);
            auto& n2 = m_expr[ie];
            auto op_p2 = get_priority(n2.m_val.m_op);
            out << '-';
            if (op_p2 != 10) out << '(';
            print_expression(out, ie, not_substitute_let_exprs);
            if (op_p2 != 10) out << ')';
            return out;
        }
        case OP_LT:
        case OP_EQ:
        case OP_NE:
        case OP_LE:
        case OP_GE:
        case OP_GT:
        case OP_POW:{
            auto op_p0 = get_priority(n.m_val.m_op);
            auto ie1 = n.child(0), ie2 = n.child(1);
            auto &n1 = m_expr[ie1], &n2 = m_expr[ie2];
            auto op_p1 = get_priority(n1.m_val.m_op), op_p2 = get_priority(n2.m_val.m_op);
            if (op_p1 <= op_p0) out << "( ";
            print_expression(out, ie1, not_substitute_let_exprs);
            if (op_p1 <= op_p0) out << " )";
            out << ' ' << getOperationName(n.m_val.m_op) << ' ';
            if (op_p2 <= op_p0) out << "( ";
            print_expression(out, ie2, not_substitute_let_exprs);
            if (op_p2 <= op_p0) out << " )";
            return out;
        }
        case OP_NOT:{
            auto op_p0 = get_priority(n.m_val.m_op);
            auto ie1 = n.child(0);
            auto &n1 = m_expr[ie1];
            auto op_p1 = get_priority(n1.m_val.m_op);
            out << "not ";
            if (op_p1 <= op_p0) out << "( ";
            print_expression(out, ie1, not_substitute_let_exprs);
            if (op_p1 <= op_p0) out << " )";
            return out;
        }
        default: { //functions
            out << ((n.m_val.m_op == H_HEAD) ? std::string("HEAD") : getOperationName(n.m_val.m_op));
            out << "(";
            for (std::size_t i = 0; i < n.number_of_childs(); ++i){
                auto ie = n.child(i);
                print_expression(out, ie, not_substitute_let_exprs);
                if (i + 1 != n.number_of_childs())
                    out << ", ";
            }
            out << ")";
            return out;
        }
    }
}

std::ostream& Energy_Parser::print_parsed_state(std::ostream& out, const std::string& prefix) const{
    out << std::scientific;
    if (!m_params.empty()){
        auto print_param = [&](std::size_t i){ 
            auto& p = m_params[i]; 
            out << p.name;
            if (!p.def_tag.empty())
                out << ":\"" << p.def_tag << "\"";
            if (!std::isnan(p.def_val))
                out << " = " << p.def_val;
        };
        out << prefix << "Var ";
        print_param(0);
        for (std::size_t i = 1; i < m_params.size(); ++i){
            out << ",\n" << prefix << "    ";
            print_param(i);
        }
        out << ";\n";
    }
    if (!m_fibers.empty()){
        auto print_fib = [&](std::size_t i){
            auto& p = m_fibers[i];
            out << p.name;
            if (!p.def_tag.empty())
                out << ":\"" << p.def_tag << "\"";
        };
        out << prefix << "Fiber ";
        print_fib(0);
        for (std::size_t i = 1; i < m_fibers.size(); ++i){
            out << ",\n" << prefix << "      ";
            print_fib(i);
        }
        out << ";\n";
    }
    if (!m_vlets.empty()){
        auto print_let = [&](std::size_t i){
            auto& p = m_vlets[i];
            out << p.name << " = ";
            print_expression(out, p.id, false);
        };
        out << prefix << "Let ";
        print_let(0);
        for (std::size_t i = 1; i < m_vlets.size(); ++i){
            out << ",\n" << prefix << "    ";
            print_let(i);
        }
        out << ";\n";
    }
    if (m_expr[0].have_child()){
        out << prefix << "Potential = ";
        print_expression(out, m_expr[0].child(0), true);
        out << ";\n";
    }
    return out;
}