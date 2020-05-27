// This file is a part of gEcon.
//
// (c) Chancellery of the Prime Minister of the Republic of Poland 2012-2015
// (c) Grzegorz Klima, Karol Podemski, Kaja Retkiewicz-Wijtiwiak 2015-2018
// License terms can be found in the file 'LICENCE'
//
// Author: Grzegorz Klima

grammar gEcon;

options
{
    language = Cpp;
}


@lexer::includes
{
#include <iostream>
#include <stdexcept>
#include <vector>
#include <string>
#include <model.h>

extern Model model_obj;
extern std::vector<std::string> errors;

}

@parser::includes
{
#include "gEconLexer.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <climits>
#include <stdexcept>
#include <model.h>


extern Model model_obj;
using symbolic::triplet;
using symbolic::ex;
using symbolic::idx_ex;
using symbolic::internal::num2str;


}

@lexer::namespace { parser }
@parser::namespace{ parser }

@lexer::traits
{
    class gEconLexer;
    class gEconParser;

    template<class ImplTraits>
    class UserTraits : public antlr3::CustomTraitsBase<ImplTraits>
    {
      public:
        static const bool TOKENS_ACCESSED_FROM_OWNING_RULE = true;
        static void displayRecognitionError(const std::string &s) { errors.push_back(s); };
    };

    typedef antlr3::Traits< gEconLexer, gEconParser, UserTraits > gEconLexerTraits;
    typedef gEconLexerTraits gEconParserTraits;
}


model
    : opts? sets? tryreduce? (block)+ EOF
    ;


opts
    : OPTS LBRACE (opt)+ RBRACE SEMI?
    ;

opt
    : OUTPUT opt_output
    | SILENT EQ b = atom_bool { model_obj.set_option(Model::silent, b); } SEMI
    | VERBOSE EQ b = atom_bool { model_obj.set_option(Model::verbose, b); } SEMI
    | WARNINGS EQ b = atom_bool { model_obj.set_option(Model::warnings_opt, b); } SEMI
    ;

opt_output
    : LOGF EQ b = atom_bool { model_obj.set_option(Model::output_logf, b); } SEMI
    | R opt_output_R
    | LATEX opt_output_latex
    ;

opt_output_R
    : EQ b = atom_bool { model_obj.set_option(Model::output_r, b); } SEMI
    | LONG EQ b = atom_bool { model_obj.set_option(Model::output_r_long, b); } SEMI
    | JACOBIAN EQ b = atom_bool { model_obj.set_option(Model::output_r_jacobian, b); } SEMI
    | RCPP EQ b = atom_bool { model_obj.set_option(Model::output_r_rcpp, b); } SEMI
    ;

opt_output_latex
    : EQ b = atom_bool { model_obj.set_option(Model::output_latex, b); } SEMI
    | LONG EQ b = atom_bool { model_obj.set_option(Model::output_latex_long, b); } SEMI
    | LANDSCAPE EQ b = atom_bool { model_obj.set_option(Model::output_latex_landscape, b); } SEMI
    ;

sets
    : SETS LBRACE (seteq | setvalid)+ RBRACE SEMI?
    ;

seteq
    : ids = id_str EQ s = setex SEMI {
        if (!model_obj.add_set(idx_set(s, ids.first)))
            model_obj.error("set \"" + ids.first + "\" already declared", ids.second);
      }
    ;

setvalid
@init {
    bool ok;
    int test;
}
    : a = setex (LEQ { test = 1; } | DEQ { test = 2; } | NEQ { test = 3; }) b = setex QUESTION {
        switch (test) {
            case 1: ok = (a <= b); break;
            case 2: ok = (a == b); break;
            case 3: ok = (a != b); break;
        }
        if (!ok) {
            std::string mes = "index set test failed (";
            switch (test) {
                case 1: mes += "lhs is not a subset of rhs"; break;
                case 2: mes += "sets are not equal"; break;
                case 3: mes += "sets are equal"; break;
            }
            mes += ")";
            model_obj.error(mes, $QUESTION.line);
        }
      }
   ;

setex returns [idx_set is]
    : s = setex_add { $is = s; }
    ;

setex_add returns [idx_set is]
    : a = setex_intersect { $is = a; }
        ((OR b = setex_intersect { $is = $is.set_sum(b); })
         |(BACKSLASH b = setex_intersect { $is = $is.set_diff(b); }))*
    ;

setex_intersect returns [idx_set is]
@init {
    a = idx_set("");
    b = idx_set("");
}
    : a = setex_cat { $is = a; } (AND b = setex_cat { $is = $is.set_intersect(b); })*
    ;

setex_cat returns [idx_set is]
@init {
    a = idx_set("");
    std::string pre, suf;
}
    : (QUOTE (IDU  { pre = $IDU.text; } | DIDU  { pre = $DIDU.text; } | idxp = idx_str { pre = idxp.first; }) QUOTE TILDE)? a = setex_atom
      (TILDE QUOTE (UID { suf = $UID.text; } | UDID { suf = $UDID.text; } | idxs = idx_str { suf = idxs.first; }) QUOTE)? { $is = a.prefix(pre).suffix(suf); }
    ;

setex_atom returns [idx_set is]
    : a = set { $is = a; }
    | LPAREN a = setex RPAREN { $is = a; }
    ;

set returns [idx_set iset]
@init {
    iset = idx_set("");
}
    : ZERO
    | is = id_str {
        if (!model_obj.is_set(is.first))
            model_obj.error("undefined set \"" + is.first + "\" used in expression", is.second);
        iset = model_obj.get_set(is.first);
      }
    | LBRACE ls = list_set {
        vec_strint::const_iterator it = ls.begin();
        for (; it != ls.end(); ++it) {
            if (!iset.add(it->first))
                model_obj.error("element \"" + it->first + "\" already present in the set",
                                it->second);
        }
      } RBRACE
    | LBRACE seq = seq_set {
        if (seq.first) {
            for (char i = seq.second; i <= (char) seq.third; ++i)
                iset.add(std::string() = i);
        } else {
            int d = 0;
            if (seq.third > 0) {
                d = 1 + (int) std::log10((double) seq.third);
            }
            for (unsigned i = seq.second; i <= seq.third; ++i)
                iset.add(num2str(i, d));
        }
      } RBRACE
    ;

list_set returns [vec_strint vsi]
    :  QUOTE ids = idx_str QUOTE { vsi.push_back(strint(ids.first, ids.second)); }
      (COMMA QUOTE ids = idx_str QUOTE { vsi.push_back(strint(ids.first, ids.second)); })*
    ;

seq_set returns [triplet<bool, unsigned, unsigned> seq]
    : QUOTE beg = SLETTER QUOTE DDOT QUOTE end = SLETTER QUOTE {
        seq.first = true; seq.second = $beg.text[0]; seq.third = $end.text[0];
        if (seq.second > seq.third)
            model_obj.error("decreasing sequence of elements in a set", $DDOT.line);
      }
    | QUOTE begc = capletter QUOTE DDOT QUOTE endc = capletter QUOTE {
        seq.first = true; seq.second = begc; seq.third = endc;
        if (seq.second > seq.third)
            model_obj.error("decreasing sequence of elements in a set", $DDOT.line);
      }
    | QUOTE begi = atom_int QUOTE DDOT QUOTE endi = atom_int QUOTE {
        seq.first = false; seq.second = begi; seq.third = endi;
        if (seq.second > seq.third)
            model_obj.error("decreasing sequence of elements in a set", $DDOT.line);
      }
    ;

capletter returns [unsigned c]
@init {
    $c = 0;
}
    : CLETTER { $c = $CLETTER.text[0]; }
    | E { $c = 'E'; }
    | R { $c = 'R'; }
    ;


tryreduce
    : TRYREDUCE LBRACE lv = list_livar { model_obj.add_red_vars(lv); } RBRACE SEMI?
    ;


block
    : BLOCK lie = list_indexing_ex
      ids = id_str {
        if (model_obj.block_declared(ids.first)) {
            model_obj.error("block \"" + ids.first + "\" already declared", ids.second);
        }
        model_obj.add_block(ids.first, ids.second, lie[0], lie[1]);
      }
      LBRACE
        block_definitions?
        (block_controls block_objective (block_constraints)? (block_identities)?
        | block_identities)
        block_shocks?
        block_calibr?
      RBRACE SEMI?
    ;

block_definitions
    : DEFS LBRACE definition+ RBRACE SEMI?
    ;

definition
    : lhs = atom_id_t EQ rhs = expr SEMI
        { model_obj.add_definition(lhs, rhs, $SEMI.line); }
    | lhs = atom_id_nt EQ rhs = expr SEMI
        { model_obj.add_definition(lhs, rhs, $SEMI.line); }
    ;

block_controls
    : CONTROLS LBRACE lv = list_ctr_var { model_obj.add_controls(lv); } RBRACE SEMI?
    ;

list_ctr_var returns [vec_exintstr listln]
@init {
    vec_ex liste;
    std::vector<std::string> lists;
    std::vector<int> listl;
}
    :  val = list_ctr_var_elem { liste.push_back(val.first); lists.push_back(val.second); }
        (COMMA { listl.push_back($COMMA.line); }
            val = list_ctr_var_elem { liste.push_back(val.first); lists.push_back(val.second); })*
        SEMI { listl.push_back($SEMI.line);
            unsigned i = 0, n = listl.size();
            listln.reserve(n);
            for (; i < n; ++i)
                listln.push_back(exintstr(liste[i], listl[i], lists[i]));
        }
    ;

list_ctr_var_elem returns [exstr val]
@init {
    std::string ref;
}
    : lie = list_indexing_ex
      v = atom_id_t (AT ids = id_str { ref = ids.first; })? {
        $val = exstr(ex(lie[0], ex(lie[1], v)), ref);
      }
    ;
block_objective
    : OBJ LBRACE objective RBRACE SEMI?
    ;

objective
    : obj = atom_id_t EQ obj_eq = expr (COLON lambda = atom_id_t)? SEMI
        { model_obj.add_objective(obj, obj_eq, lambda, $SEMI.line); }
    ;


block_constraints
    : CONSTRAINTS LBRACE constraint+ RBRACE SEMI?
    ;

constraint
    : lie = list_indexing_ex
      lhs = expr EQ rhs = expr (COLON lambda = atom_id_t)? SEMI
        { model_obj.add_constraint(ex(lie[0], ex(lie[1], lhs)), ex(lie[0], ex(lie[1], rhs)),
                                    ex(lie[0], ex(lie[1], lambda)), $SEMI.line); }
    | lr = list_ref AT ids = id_str SEMI {
            for (unsigned i = 0; i < lr.size(); ++i) {
                model_obj.add_constraint_ref(ids.first, lr[i].first, lr[i].second);
            }
      }
    ;

list_ref returns [vec_intint lr]
    : r = ref_sec { lr.push_back(r); } (COMMA r = ref_sec { lr.push_back(r); })*
    ;

ref_sec returns [std::pair<int, int> rs]
    : OBJ { return std::pair<int, int>(Model_block::objective, $OBJ.line); }
    | CONSTRAINTS { return std::pair<int, int>(Model_block::constraints, $CONSTRAINTS.line); }
    | FOCS { return std::pair<int, int>(Model_block::focs, $FOCS.line); }
    | IDS { return std::pair<int, int>(Model_block::identities, $IDS.line); }
    ;

block_identities
    : IDS LBRACE identity+ RBRACE SEMI?
    ;

identity
    : lie = list_indexing_ex
      lhs = expr EQ rhs = expr SEMI {
        model_obj.add_identity(ex(lie[0], ex(lie[1], lhs)), ex(lie[0], ex(lie[1], rhs)), $SEMI.line);
      }
    ;

block_shocks
    : SHOCKS LBRACE ls = list_var { model_obj.add_shocks(ls); } RBRACE SEMI?
    ;

shock
    : lie = list_indexing_ex
      s = atom_id_t SEMI { model_obj.add_shock(ex(lie[0], ex(lie[1], s)), $SEMI.line); }
    ;

block_calibr
    : CALIBR LBRACE calibr_eq+ RBRACE SEMI?
    ;

calibr_eq
    : lie = list_indexing_ex
      lhs = expr EQ rhs = expr ((RARROW lp = list_par) | SEMI) {
        vec_exint pln;
        pln.reserve(lp.size());
        for (unsigned i = 0; i < lp.size(); ++i) {
            pln.push_back(exint(ex(lie[0], ex(lie[1], lp[i].first)), lp[i].second));
        }
        model_obj.add_calibr(ex(lie[0], ex(lie[1], lhs)), ex(lie[0], ex(lie[1], rhs)), $EQ.line, pln);
      }
    ;



list_var returns [vec_exint listln]
@init {
    vec_ex liste;
    std::vector<int> listl;
}
    :  val = list_var_elem { liste.push_back(val); }
        (COMMA { listl.push_back($COMMA.line); }
            val = list_var_elem { liste.push_back(val); })*
        SEMI { listl.push_back($SEMI.line);
            unsigned i = 0, n = listl.size();
            listln.reserve(n);
            for (; i < n; ++i)
                listln.push_back(exint(liste[i], listl[i]));
        }
    ;



list_livar returns [vec_exint listln]
@init {
    vec_ex liste;
    std::vector<int> listl;
}
    :  val = list_livar_elem { liste.push_back(val); }
        (COMMA { listl.push_back($COMMA.line); }
            val = list_livar_elem { liste.push_back(val); })*
        SEMI { listl.push_back($SEMI.line);
            unsigned i = 0, n = listl.size();
            listln.reserve(n);
            for (; i < n; ++i)
                listln.push_back(exint(liste[i], listl[i]));
        }
    ;

list_var_elem returns [ex val]
    : lie = list_indexing_ex
        v = atom_id_t { $val = ex(lie[0], ex(lie[1], v)); }
    ;

list_livar_elem returns [ex val]
    : lie = list_lindexing_ex
        v = atom_id_t { $val = ex(lie[0], ex(lie[1], ex(lie[2], ex(lie[3], v)))); }
    ;

list_par returns [vec_exint listln]
@init {
    vec_ex liste;
    std::vector<int> listl;
}
    :  val = list_par_elem { liste.push_back(val); }
        (COMMA { listl.push_back($COMMA.line); }
            val = list_par_elem { liste.push_back(val); })*
        SEMI { listl.push_back($SEMI.line);
            unsigned i = 0, n = listl.size();
            listln.reserve(n);
            for (; i < n; ++i)
                listln.push_back(exint(liste[i], listl[i]));
        }
    ;

list_par_elem returns [ex val]
    : lie = list_indexing_ex
        p = atom_id_nt { $val = ex(lie[0], ex(lie[1], p)); }
    ;

expr returns [ex val]
    : a = expr_add { $val = a; }
    ;

expr_add returns [ex val]
    : a = expr_sum { $val = a; }
                ((PLUS b = expr_sum { $val = $val + b; })
                |(MINUS b = expr_sum { $val = $val - b; }))*
    ;

expr_sum returns [ex val]
@init {
    bool minus = false;
}
    : e = expr_mul { $val = e; }
    | (MINUS { minus = true;} )? SUM ie = indexing_ex e = expr_sum {
        $val = ie.first ? sum(ie.first, e) : e; $val = minus ? -$val : $val;
      }
    ;

expr_mul returns [ex val]
    : a = expr_prod { $val = a; }
                ((MUL b = expr_prod { $val = $val * b; })
                |(DIV b = expr_prod { $val = $val / b; }))*
    ;

expr_prod returns [ex val]
@init {
    bool minus = false;
}
    : e = expr_pow { $val = e; }
    | (MINUS { minus = true;} )? PROD ie = indexing_ex e = expr_prod {
        $val = ie.first ? prod(ie.first, e) : e; $val = minus ? -$val : $val;
      }
    ;


list_indexing_ex returns [vec_idx_ex lie]
@init {
    bool toomany = false;
}
    : (ie = indexing_ex {
        lie.push_back(ie.first);
        if (!toomany && (lie.size() > 2)) {
            toomany = true;
            errors.push_back("up to 2 indexing expressions in a template declaration are supported in this context; error near line "
                                + num2str(ie.second));
        }
      })*
      { while (lie.size() < 2) lie.push_back(idx_ex()); }
    ;


list_lindexing_ex returns [vec_idx_ex lie]
@init {
    bool toomany = false;
}
    : (ie = indexing_ex {
        lie.push_back(ie.first);
        if (!toomany && (lie.size() > 4)) {
            toomany = true;
            errors.push_back("up to 4 indexing expressions in a template declaration are supported in this context; error near line "
                                + num2str(ie.second));
        }
      })*
      { while (lie.size() < 4) lie.push_back(idx_ex()); }
    ;


indexing_ex returns [std::pair<idx_ex, int> val]
@init {
    idx_set iset;
}
    : LANGBR iv = id_str DBLCOLON is = id_str RANGBR {
        if (!model_obj.is_set(is.first))
            model_obj.error("undefined set \"" + is.first + "\" used in expression", $RANGBR.line);
        iset = model_obj.get_set(is.first);
        if (!iset.size())
            model_obj.warning("empty set \"" + is.first + "\" used in expression", $RANGBR.line);
        $val = std::pair<idx_ex, int>(idx_ex(iv.first, iset), $RANGBR.line);
      }
    | LANGBR iv = id_str DBLCOLON is = id_str BACKSLASH ei = id_str RANGBR {
        if (!model_obj.is_set(is.first))
            model_obj.error("undefined set \"" + is.first + "\" used in expression", $RANGBR.line);
        iset = model_obj.get_set(is.first);
        if (!iset.size())
            model_obj.warning("empty set \"" + is.first + "\" used in expression", $RANGBR.line);
        if (iv.first == ei.first)
            model_obj.error("excluded index (\"" + ei.first + "\") is the same as free index in indexing expression",
                            $RANGBR.line);
        $val = std::pair<idx_ex, int>(idx_ex(iv.first, iset, ei.first, false), $RANGBR.line);
      }
    | LANGBR iv = id_str DBLCOLON is = id_str BACKSLASH QUOTE ei = idx_str QUOTE RANGBR {
        if (!model_obj.is_set(is.first))
            model_obj.error("undefined set \"" + is.first + "\" used in expression",
                            $RANGBR.line);
        iset = model_obj.get_set(is.first);
        if (!iset.size())
            model_obj.warning("empty set \"" + is.first + "\" used in expression",
                              $RANGBR.line);
        $val = std::pair<idx_ex, int>(idx_ex(iv.first, iset, ei.first, true), $RANGBR.line);
      }
    ;

expr_pow returns [ex val]
@init {
    bool min = false;
    vec_ex ve;
    std::vector<bool> vs;
}
    : (MINUS { min = true; })?
        a = expr_atom { ve.push_back(a); vs.push_back(min); }
        (POW { min = false; } (MINUS { min = true; })? b = expr_atom
        { ve.push_back(b); vs.push_back(min); })*
        {
            int i = ve.size() - 1;
            val = ve[i];
            if (vs[i]) val = -val;
            for (--i; i >= 0; --i) {
                val = pow(ve[i], val);
                if (vs[i]) val = -val;
            }
        }
    ;


expr_atom returns [ex val]
    : a = atom_num { $val = a; }
    | a = atom_id { $val = a; }
    | a = atom_delta { $val = a; }
    | a = expr_func { $val = a; }
    | a = expr_e { $val = a; }
    | LPAREN (a = expr { $val = a; }) RPAREN
    ;

expr_func returns [ex val]
    : SQRT LPAREN a = expr RPAREN { $val = sqrt(a); }
    | EXP LPAREN a = expr RPAREN { $val = exp(a); }
    | LOG LPAREN a = expr RPAREN { $val = log(a); }
    | SIN LPAREN a = expr RPAREN { $val = sin(a); }
    | COS LPAREN a = expr RPAREN { $val = cos(a); }
    | TAN LPAREN a = expr RPAREN { $val = tan(a); }
    | ASIN LPAREN a = expr RPAREN { $val = asin(a); }
    | ACOS LPAREN a = expr RPAREN { $val = acos(a); }
    | ATAN LPAREN a = expr RPAREN { $val = atan(a); }
    | SINH LPAREN a = expr RPAREN { $val = sinh(a); }
    | COSH LPAREN a = expr RPAREN { $val = cosh(a); }
    | TANH LPAREN a = expr RPAREN { $val = tanh(a); }
//    | ERF LPAREN a = expr RPAREN { $val = erf(a); }
    | PNORM LPAREN a = expr RPAREN { $val = pnorm(a); }
    ;

expr_e returns [ex val]
    : E LBRACK l = time_idx RBRACK LBRACK a = expr RBRACK { $val = symbolic::E(a, l);}
    ;

atom_delta returns [ex val]
    : DELTA LANGBR i1 = id_str COMMA i2 = id_str RANGBR { $val = ex(false, i1.first, false, i2.first); }
    | DELTA LANGBR QUOTE i1 = idx_str QUOTE COMMA i2 = idx_str RANGBR { $val = ex(true, i1.first, false, i2.first); }
    | DELTA LANGBR i1 = id_str COMMA QUOTE i2 = idx_str QUOTE RANGBR { $val = ex(false, i1.first, true, i2.first); }
    | DELTA LANGBR QUOTE i1 = idx_str QUOTE COMMA QUOTE i2 = idx_str QUOTE RANGBR { $val = ex(true, i1.first, true, i2.first); }
    ;

atom_id returns [ex val]
    : v = atom_id_nt { $val = v; }
    | v = atom_id_t { $val = v; }
    ;

atom_id_t returns [ex val]
    : ids = id_str LBRACK l = time_idx RBRACK { $val = ex(ids.first, l); }
    | ids = id_str li = list_idx LBRACK l = time_idx RBRACK {
            switch(li.size()) {
                case 1:
                    $val = ex(ids.first, l, li[0].first, li[0].second);
                    break;
                case 2:
                    $val = ex(ids.first, l, li[0].first, li[0].second, li[1].first, li[1].second);
                    break;
                case 3:
                    $val = ex(ids.first, l, li[0].first, li[0].second, li[1].first, li[1].second,
                                            li[2].first, li[2].second);
                    break;
                case 4:
                    $val = ex(ids.first, l, li[0].first, li[0].second, li[1].first, li[1].second,
                                            li[2].first, li[2].second, li[3].first, li[3].second);
                    break;
            }
        }
    ;

time_idx returns [int val]
@init {
    $val = 0;
}
    :
    | l = atom_int { $val = l; }
    | PLUS l = atom_int { $val = l; }
    | MINUS l = atom_int { $val = -l; }
    | SS { $val = INT_MIN; }
    | MINUS INF { $val = INT_MIN; }
    ;

atom_id_nt returns [ex val]
    : ids = id_str { $val = ex(ids.first); }
    | ids = id_str li = list_idx {
            switch(li.size()) {
                case 1:
                    $val = ex(ids.first, li[0].first, li[0].second);
                    break;
                case 2:
                    $val = ex(ids.first, li[0].first, li[0].second, li[1].first, li[1].second);
                    break;
                case 3:
                    $val = ex(ids.first, li[0].first, li[0].second, li[1].first, li[1].second,
                                         li[2].first, li[2].second);
                    break;
                default:
                    $val = ex(ids.first, li[0].first, li[0].second, li[1].first, li[1].second,
                                         li[2].first, li[2].second, li[3].first, li[3].second);
                    break;
            }
        }
    ;

list_idx returns [std::vector<std::pair<bool, std::string> > listbs]
@init {
    typedef std::pair<bool, std::string> pbs;
    bool toomany = false;
}
    : LANGBR
        ((val = idx_str { listbs.push_back(pbs(false, val.first)); }) | (QUOTE val = idx_str QUOTE { listbs.push_back(pbs(true, val.first)); }))
        (COMMA ((val = idx_str { listbs.push_back(pbs(false, val.first)); }) | (QUOTE val = idx_str QUOTE { listbs.push_back(pbs(true, val.first)); })) {
            if (!toomany && (listbs.size() > 4)) {
                toomany = true;
                errors.push_back("up to 4 indices are supported; error near line "
                                 + num2str($COMMA.line) + ", pos: " + num2str($COMMA->get_charPositionInLine() + 1));
            }
          } )* RANGBR
    ;

idx_str returns [std::pair<std::string, int> val]
    : is = id_str { $val = is; }
    | E { $val = std::pair<std::string, int>($E.text, $E.line); }
    | INF { $val = std::pair<std::string, int>($INF.text, $INF.line); }
    | ZERO { $val = std::pair<std::string, int>($ZERO.text, $ZERO.line); }
    | INT { $val = std::pair<std::string, int>($INT.text, $INT.line); }
    | DID { $val = std::pair<std::string, int>($DID.text, $DID.line); }
    ;

id_str returns [std::pair<std::string, int> val]
    : ID { $val = std::pair<std::string, int>($ID.text, $ID.line); }
    | SLETTER { $val = std::pair<std::string, int>($SLETTER.text, $SLETTER.line); }
    | CLETTER { $val = std::pair<std::string, int>($CLETTER.text, $CLETTER.line); }
    | OUTPUT { $val = std::pair<std::string, int>($OUTPUT.text, $OUTPUT.line); }
    | R { $val = std::pair<std::string, int>($R.text, $R.line); }
    | LATEX { $val = std::pair<std::string, int>($LATEX.text, $LATEX.line); }
    | LANDSCAPE { $val = std::pair<std::string, int>($LANDSCAPE.text, $LANDSCAPE.line); }
    | LOGF { $val = std::pair<std::string, int>($LOGF.text, $LOGF.line); }
    | LONG { $val = std::pair<std::string, int>($LONG.text, $LONG.line); }
    | SHORT { $val = std::pair<std::string, int>($SHORT.text, $SHORT.line); }
    | BTRUE { $val = std::pair<std::string, int>($BTRUE.text, $BTRUE.line); }
    | BFALSE { $val = std::pair<std::string, int>($BFALSE.text, $BFALSE.line); }
    | SS { $val = std::pair<std::string, int>($SS.text, $SS.line); }
    | SILENT { $val = std::pair<std::string, int>($SILENT.text, $SILENT.line); }
    | VERBOSE { $val = std::pair<std::string, int>($VERBOSE.text, $VERBOSE.line); }
    | WARNINGS { $val = std::pair<std::string, int>($WARNINGS.text, $WARNINGS.line); }
    | JACOBIAN { $val = std::pair<std::string, int>($JACOBIAN.text, $JACOBIAN.line); }
    | RCPP { $val = std::pair<std::string, int>($RCPP.text, $RCPP.line); }
    ;

atom_num returns [ex val]
    : a = atom_int { $val = ex((double) a); }
    | a = atom_double { $val = ex(a); }
    ;

atom_int returns [int val]
@init {
    $val = 0;
}
    :  ZERO
    |  INT  { $val = atoi($INT.text.c_str()); }
    ;

atom_double returns [double val]
@init {
    $val = 0.;
}
    :  DOUBLE  { $val = atof($DOUBLE.text.c_str()); }
    ;

atom_bool returns [bool val]
@init {
    $val = false;
}
    : BTRUE { $val = true; }
    | BFALSE
    ;


// IDs that can be used as var names
SILENT      : 'silent';
VERBOSE     : 'verbose';
WARNINGS    : 'warnings';
OUTPUT      : 'output';
R           : 'R';
LOGF        : 'logfile';
LONG        : 'long';
SHORT       : 'short';
LATEX       : 'LaTeX'|'latex';
LANDSCAPE   : 'landscape';
JACOBIAN    : 'Jacobian'|'jacobian';
RCPP        : 'Rcpp'|'rcpp';

// Bool vals
BTRUE   : 'true'|'TRUE';
BFALSE  : 'false'|'FALSE';

// Keywords
OPTS        : 'options' ;
SETS        : 'indexsets' ;
TRYREDUCE   : 'tryreduce' ;
BLOCK       : 'block' ;
DEFS        : 'definitions';
CONTROLS    : 'controls';
OBJ         : 'objective';
CONSTRAINTS : 'constraints';
FOCS        : 'focs';
IDS         : 'identities';
SHOCKS      : 'shocks';
CALIBR      : 'calibration';

// Steady state
SS  : 'ss'|'SS';
// Infinity
INF : 'inf'|'Inf'|'INF';

// Expected val
E   : 'E';

// Sums, products and Kronecker deltas
SUM     : 'SUM';
PROD    : 'PROD';
DELTA   : 'KRONECKER_DELTA';

// Functions
SQRT    : 'sqrt';
EXP     : 'exp';
LOG     : 'log';
SIN     : 'sin';
COS     : 'cos';
TAN     : 'tan';
ASIN    : 'asin';
ACOS    : 'acos';
ATAN    : 'atan';
SINH    : 'sinh';
COSH    : 'cosh';
TANH    : 'tanh';
// ERF     : 'erf';
PNORM   : 'pnorm';

// Numbers
ZERO:   '0';

INT :   '0'
    |   '1'..'9'
    |   '1'..'9' '0'..'9'
    |   '1'..'9' '0'..'9' '0'..'9'
    |   '1'..'9' '0'..'9' '0'..'9' '0'..'9'
    |   '1'..'9' '0'..'9' '0'..'9' '0'..'9' '0'..'9'
    |   '1'..'9' '0'..'9' '0'..'9' '0'..'9' '0'..'9' '0'..'9'
    |   '1'..'9' '0'..'9' '0'..'9' '0'..'9' '0'..'9' '0'..'9' '0'..'9'
    |   '1'..'9' '0'..'9' '0'..'9' '0'..'9' '0'..'9' '0'..'9' '0'..'9' '0'..'9'
    |   '1'..'9' '0'..'9' '0'..'9' '0'..'9' '0'..'9' '0'..'9' '0'..'9' '0'..'9' '0'..'9'
    |   '1'..'9' '0'..'9' '0'..'9' '0'..'9' '0'..'9' '0'..'9' '0'..'9' '0'..'9' '0'..'9' '0'..'9'
    ;

DOUBLE
    :   '0'..'9'+ '.' '0'..'9'+ (('e' | 'E') ('+' | '-')? '0'..'9'+)?
    |   '0'..'9'+ '.' (('e' | 'E') ('+' | '-')? '0'..'9'+)?
    |   '.' '0'..'9'+ (('e' | 'E') ('+' | '-')? '0'..'9'+)?
    |   '0' '0'..'9'+
    |   '1'..'9' '0'..'9' '0'..'9' '0'..'9' '0'..'9' '0'..'9' '0'..'9' '0'..'9' '0'..'9' '0'..'9' '0'..'9'+
    ;

// Ids
SLETTER
    : 'a'..'z'
    ;

CLETTER
    : 'A'..'Z'
    ;

ID  : ('a'..'z'|'A'..'Z') (('_')?('a'..'z'|'A'..'Z'|'0'..'9'))*
    ;

UID
    : '_'? ('a'..'z'|'A'..'Z') (('_')?('a'..'z'|'A'..'Z'|'0'..'9'))*
    ;

DID
    : ('0'..'9')(('_')?('a'..'z'|'A'..'Z'|'0'..'9'))+
    ;

UDID
    : '_' ('0'..'9')(('_')?('a'..'z'|'A'..'Z'|'0'..'9'))*
    ;

IDU
    : ('a'..'z'|'A'..'Z') (('_')?('a'..'z'|'A'..'Z'|'0'..'9'))* '_'?
    ;

DIDU
    : ('0'..'9')(('_')?('a'..'z'|'A'..'Z'|'0'..'9'))* '_'?
    ;

// Whitespace
WS  :   (' '|'\t'|'\n'|'\r')+ { $channel = HIDDEN; } ;

// Symbols
TILDE    : '~' ;
QUESTION : '?' ;
EXCLAM   : '!' ;
DOLLAR   : '$' ;
AT       : '@' ;
AND      : '&' ;
OR       : '|' ;
DOR      : '||' ;
SEMI     : ';' ;
COLON    : ':' ;
DBLCOLON : '::' ;
DDOT     : '..' ;
COMMA    : ',' ;
RARROW   : '->' ;
PLUS     : '+' ;
MINUS    : '-' ;
MUL      : '*' ;
DIV      : '/' ;
POW      : '^' ;
EQ       : '=' ;
DEQ      : '==' ;
NEQ      : '!=' ;
LEQ      : '<=' ;
QUOTE    : '\'' ;
DQUOTE   : '"' ;
BACKSLASH: '\\' ;


// Brackets
LBRACE  : '{' ;
RBRACE  : '}' ;
LPAREN  : '(' ;
RPAREN  : ')' ;
LBRACK  : '[' ;
RBRACK  : ']' ;
LANGBR  : '<' ;
RANGBR  : '>' ;

// Comment
COMMENT
    :   '#' ~('\n'|'\r')* { $channel = HIDDEN; }
    |   '%' ~('\n'|'\r')* { $channel = HIDDEN; }
    |   '//' ~('\n'|'\r')* { $channel = HIDDEN; }
    ;

