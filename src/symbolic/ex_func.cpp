/***********************************************************
 * (c) Kancelaria Prezesa Rady Ministrów 2012-2015         *
 * Treść licencji w pliku 'LICENCE'                        *
 *                                                         *
 * (c) Chancellery of the Prime Minister 2012-2015         *
 * License terms can be found in the file 'LICENCE'        *
 *                                                         *
 * Author: Grzegorz Klima                                  *
 ***********************************************************/

/** \file ex_func.cpp
 * \brief Functions.
 */


#include <ex_func.h>
#include <ex_num.h>
#include <ops.h>
#include <cmp.h>
#include <error.h>
#include <cmath>
#include <utils.h>
#include <climits>



namespace {

const char *fname[] =
    {
    "exp",     //    EXP
    "log",     //    LOG
    "sin",     //    SIN
    "cos",     //    COS
    "tan",     //    TAN
    "asin",    //    ASIN
    "acos",    //    ACOS
    "atan",    //    ATAN
    "sinh",    //    SINH
    "cosh",    //    COSH
    "tanh"     //    TANH
    // "erf"      //    ERF
    };

symbolic::Number (*feval[])(const symbolic::Number&) =
    {
    symbolic::exp,   //    EXP
    symbolic::log,   //    LOG
    symbolic::sin,   //    SIN
    symbolic::cos,   //    COS
    symbolic::tan,   //    TAN
    symbolic::asin,  //    ASIN
    symbolic::acos,  //    ACOS
    symbolic::atan,  //    ATAN
    symbolic::sinh,  //    SIN
    symbolic::cosh,  //    COS
    symbolic::tanh   //    TAN
    // symbolic::erf    //    ERF
    };

} /* namespace */


using namespace symbolic;
using namespace symbolic::internal;


ex_func::ex_func(func_code c, const ptr_base &arg)
    : ex_base(FUN | (arg->flag() & HAST) | (arg->flag() & HASIDX) | SINGLE),
      m_code(c), m_arg(arg)
{
}



ptr_base
ex_func::create(func_code c, const ptr_base &arg)
{
    return ptr_base(new ex_func(c, arg));
}


void
ex_func::destroy(ex_base *ptr)
{
    delete ptr;
}


int
ex_func::compare(const ex_func &b) const
{
    int c = compareT(m_code, b.m_code);
    if (c) return c;
    c = symbolic::internal::compare(m_arg, b.m_arg);
    return c;
}


std::string
ex_func::str(int pflag) const
{
    if (pflag & INDEXING_ONLY) return std::string();
    std::string name(fname[m_code]);
    return name + '(' + m_arg->str(pflag) + ')';
}


std::string
ex_func::strmap(const map_str_str &mss) const
{
    std::string name(fname[m_code]);
    return name + '(' + m_arg->strmap(mss) + ')';
}


std::string
ex_func::tex(int pflag) const
{
    if (pflag & INDEXING_ONLY) return std::string();
    if (m_code == EXP) return "e^{" + m_arg->tex(pflag) + "}";
    std::string name = func2tex(fname[m_code]);
    if (m_arg->flag() & SINGLE) return name + '{' + m_arg->tex(pflag) + '}';
    return name + "\\left(" + m_arg->tex(pflag) + "\\right)";
}


symbolic::Number
ex_func::eval(func_code c, const symbolic::Number &x)
{
    return feval[c](x);
}


int
ex_func::get_lag_max(bool stop_on_E) const
{
    return m_arg->get_lag_max(stop_on_E);
}


int
ex_func::get_lag_min(bool stop_on_E) const
{
    return m_arg->get_lag_min(stop_on_E);
}


#ifndef M_2_SQRTPI
#define M_2_SQRTPI 1.12837916709551257389615890312154517
#endif /* M_2_SQRTPI */


ptr_base
ex_func::diff(const ptr_base &p) const
{
    ptr_base d = m_arg->diff(p);
    if (d->is0()) return d;

    ptr_base x = m_arg;

    switch (m_code) {
        case EXP:
            return mk_mul(d, mk_func(EXP, x));
        case LOG:
            return mk_div(d, x);
        case SIN:
            return mk_mul(d, mk_func(COS, x));
        case COS:
            return mk_mul(d, mk_neg(mk_func(SIN, x)));
        case TAN:
            return mk_mul(d, mk_pow(mk_func(COS, x), ex_num::create(-2)));
        case ASIN:
            return mk_mul(d, mk_pow(mk_sub(ex_num::one(), mk_pow(x, ex_num::create(2))),
                                    ex_num::create(-0.5)));
        case ACOS:
            return mk_mul(d, mk_neg(mk_pow(mk_sub(ex_num::one(), mk_pow(x, ex_num::create(2))),
                                            ex_num::create(-0.5))));
        case ATAN:
            return mk_div(d, mk_add(ex_num::one(), mk_pow(x, ex_num::create(2))));
        case SINH:
            return mk_mul(d, mk_func(COSH, x));
        case COSH:
            return mk_mul(d, mk_func(SINH, x));
        case TANH:
            return mk_mul(d, mk_sub(ex_num::one(), mk_pow(mk_func(TANH, x), ex_num::create(2))));
        // case ERF:
            // return mk_mul(mk_mul(d, ex_num::create(M_2_SQRTPI)),
                            // mk_func(EXP, mk_neg(mk_pow(x, ex_num::create(2)))));
        default:
            INTERNAL_ERROR
    }
    return ex_num::zero();
}


bool
ex_func::has(const ptr_base &what, search_flag f, bool exact_idx) const
{
    return m_arg->has(what, f, exact_idx);
}


bool
ex_func::hasidx(int idx) const
{
    return m_arg->hasidx(idx);
}



