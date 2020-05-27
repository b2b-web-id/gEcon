/*****************************************************************************
 * This file is a part of gEcon.                                             *
 *                                                                           *
 * (c) Chancellery of the Prime Minister of the Republic of Poland 2012-2015 *
 * (c) Grzegorz Klima, Karol Podemski, Kaja Retkiewicz-Wijtiwiak 2015-2018   *
 * License terms can be found in the file 'LICENCE'                          *
 *                                                                           *
 * Author: Grzegorz Klima                                                    *
 *****************************************************************************/

/** \file ex_pow.cpp
 * \brief Powers.
 */

#include <ex_pow.h>
#include <ex_num.h>
#include <ops.h>
#include <cmp.h>
#include <cmath>
#include <climits>
#include <algorithm>


using namespace symbolic::internal;


ex_pow::ex_pow(const ptr_base &a, const ptr_base &b)
    : ex_base(POW | (a->flag() & HAST) | (b->flag() & HAST)
                  | (a->flag() & HASIDX) | (b->flag() & HASIDX) | SINGLE),
      m_base(a), m_exp(b)
{
}


ptr_base
ex_pow::create(const ptr_base &a, const ptr_base &b)
{
    return ptr_base(new ex_pow(a, b));
}

void
ex_pow::destroy(ex_base *ptr)
{
    delete ptr;
}


int
ex_pow::compare(const ex_pow &b) const
{
    int c = symbolic::internal::compare(m_base, b.m_base);
    if (c) return c;
    return symbolic::internal::compare(m_exp, b.m_exp);
}


std::string
ex_pow::str(int pflag, bool c_style) const
{
    if (pflag & INDEXING_ONLY) return std::string();
    bool braceb = true, braceex = true;
    std::string res;
    unsigned bt = m_base->type(), et = m_exp->type();

    // braces
    if (m_base->val() > 0.) braceb = false;
    else if ((bt == DELTA) || (bt == SYMB) || (bt == SYMBIDX)
             || (bt == VART) || (bt == VARTIDX)) braceb = false;
    if ((et == NUM) || (et == DELTA) || (et == SYMB) || (et == SYMBIDX)
        || (et == VART) || (et == VARTIDX)) braceex = false;

    if (c_style) {
        res += "pow(";
        res += m_base->str(pflag, c_style);
        res += ", ";
        res += m_exp->str(pflag, c_style);
        res += ")";
    } else {
        if (braceb) res += "(";
        res += m_base->str(pflag, c_style);
        if (braceb) res += ")";
        res += '^';
        if (braceex) res += "(";
        res += m_exp->str(pflag, c_style);
        if (braceex) res += ")";
    }
    return res;
}


std::string
ex_pow::strmap(const map_str_str &mss, bool c_style) const
{
    bool braceb = true, braceex = true;
    std::string res;
    unsigned bt = m_base->type(), et = m_exp->type();

    // braces
    if (m_base->val() > 0.) braceb = false;
    else if ((bt == DELTA) || (bt == SYMB) || (bt == SYMBIDX)
             || (bt == VART) || (bt == VARTIDX)) braceb = false;
    if ((et == NUM) || (et == DELTA) || (et == SYMB) || (et == SYMBIDX)
        || (et == VART) || (et == VARTIDX)) braceex = false;

    if (c_style) {
        res += "pow(";
        if (braceb) res += "(";
        res += m_base->strmap(mss, c_style);
        if (braceb) res += ")";
        res += ", ";
        if (braceex) res += "(";
        res += m_exp->strmap(mss, c_style);
        if (braceex) res += ")";
        res += ")";
    } else {
        if (braceb) res += "(";
        res += m_base->strmap(mss, c_style);
        if (braceb) res += ")";
        res += '^';
        if (braceex) res += "(";
        res += m_exp->strmap(mss, c_style);
        if (braceex) res += ")";
    }

    return res;
}


std::string
ex_pow::tex(int pflag) const
{
    if (pflag & INDEXING_ONLY) return std::string();
    bool braceb = true;
    std::string res;
    unsigned bt = m_base->type();

    // braces
    if (m_base->val() > 0.) braceb = false;
    else if ((bt == DELTA) || (bt == SYMB) || (bt == SYMBIDX)
             || (bt == VART) || (bt == VARTIDX)) braceb = false;

    if (braceb) res += "\\left("; else res += '{';
    res += m_base->tex(pflag);
    if (braceb) res += "\\right)"; else res += '}';
    res += "^{";
    res += m_exp->tex(pflag);
    res += "}";

    return res;
}


int
ex_pow::get_lag_max(bool stop_on_E) const
{
    return std::max(m_base->get_lag_max(stop_on_E), m_exp->get_lag_max(stop_on_E));
}


int
ex_pow::get_lag_min(bool stop_on_E) const
{
    return std::min(m_base->get_lag_min(stop_on_E), m_exp->get_lag_min(stop_on_E));
}


ptr_base
ex_pow::diff(const ptr_base &p) const
{
    ptr_base f = m_base, g = m_exp, fp = m_base->diff(p), gp = m_exp->diff(p);
    bool f0 = fp->is0(), g0 = gp->is0();
    if (f0 && g0) return fp;
    if (f0) return mk_mul(mk_pow(f, g), mk_mul(mk_func(LOG, f), gp));
    if (g0) return mk_mul(mk_pow(f, mk_sub(g, ex_num::one())), mk_mul(g, fp));
    return mk_mul(mk_add(mk_mul(gp, mk_func(LOG, f)), mk_mul(mk_div(g, f), fp)), mk_pow(f, g));
}



bool
ex_pow::has(const ptr_base &what, search_flag f, bool exact_idx) const
{
    return m_base->has(what, f, exact_idx) || m_exp->has(what, f, exact_idx);
}



bool
ex_pow::hasidx(int idx) const
{
    return m_base->hasidx(idx) || m_exp->hasidx(idx);
}












