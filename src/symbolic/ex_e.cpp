/*****************************************************************************
 * This file is a part of gEcon.                                             *
 *                                                                           *
 * (c) Chancellery of the Prime Minister of the Republic of Poland 2012-2015 *
 * (c) Grzegorz Klima, Karol Podemski, Kaja Retkiewicz-Wijtiwiak 2015-2018   *
 * License terms can be found in the file 'LICENCE'                          *
 *                                                                           *
 * Author: Grzegorz Klima                                                    *
 *****************************************************************************/

/** \file ex_e.cpp
 * \brief Expected value (conditional - time indexed).
 */


#include <error.h>
#include <ex_e.h>
#include <cmp.h>
#include <ops.h>
#include <utils.h>
#include <climits>


using namespace symbolic::internal;


ex_e::ex_e(const ptr_base &e, int l)
    : ex_base(EX | HAST | (e->flag() & HASIDX)), m_arg(e), m_lag(l)
{
}

ptr_base
ex_e::create(const ptr_base &e, int l)
{
    return ptr_base(new ex_e(e, l));
}


void
ex_e::destroy(ex_base *ptr)
{
    delete ptr;
}


int
ex_e::compare(const ex_e &b) const
{
    if (m_lag < b.m_lag) return -1;
    if (m_lag > b.m_lag) return 1;
    return symbolic::internal::compare(m_arg, b.m_arg);
}


std::string
ex_e::str(int pflag, bool c_style) const
{
    if (pflag & INDEXING_ONLY) return std::string();
    if (pflag & CONVERT_T) {
        if (m_lag < 0) return "E__tm" + num2str(-m_lag) + '[' + m_arg->str(pflag, c_style) + ']';
        if (m_lag == 0) return "E__t[" + m_arg->str(pflag, c_style) + ']';
        return "E__tp" + num2str(m_lag) + '[' + m_arg->str(pflag, c_style) + ']';
    } else if (pflag & DROP_T) {
        USER_ERROR("invalid print flag in ex_e::str()")
    } else {
        return (m_lag) ? "E[" + num2str(m_lag) + "][" + m_arg->str(pflag, c_style) + ']':
                            "E[][" + m_arg->str(pflag, c_style) + ']';
    }
}


std::string
ex_e::strmap(const map_str_str &mss, bool c_style) const
{
    return m_arg->strmap(mss);
}


std::string
ex_e::tex(int pflag) const
{
    if (pflag & INDEXING_ONLY) return std::string();
    std::string res;
    if ((pflag & CONVERT_T) || (pflag & DROP_T)) {
        USER_ERROR("invalid print flag in ex_e::tex()")
    }
    if (m_lag == 0) res = "\\mathrm{E}_{t}";
    else if (m_lag < 0) res = "\\mathrm{E}_{t" + num2str(m_lag) + '}';
    else res = "\\mathrm{E}_{t+" + num2str(m_lag) + '}';
    res += "\\left[";
    res += m_arg->tex(pflag);
    res += "\\right]";
    return res;
}


int
ex_e::get_lag_max(bool stop_on_E) const
{
    if (stop_on_E) return m_lag;
    return std::max(m_lag, m_arg->get_lag_max());
}


int
ex_e::get_lag_min(bool stop_on_E) const
{
    if (stop_on_E) return m_lag;
    return std::min(m_lag, m_arg->get_lag_min());
}



ptr_base
ex_e::diff(const ptr_base &p) const
{
    return mk_E(m_arg->diff(p), m_lag);
}


bool
ex_e::has(const ptr_base &what, search_flag f, bool exact_idx) const
{
    if (what->type() == EX) {
        const ex_e *w = what.get<ex_e>();
        int ld = m_lag - w->m_lag;
        bool ok = !symbolic::internal::compare(symbolic::internal::lag(m_arg, ld), w->m_arg);
        if (!ok) return false;
        switch (f) {
            case EXACT_T:
                return (ld == 0) ? true : false;
            case ANY_T:
                return true;
            case DIFF_T:
                return (ld != 0) ? true : false;
            case LEAD_T:
                return (ld > 0) ? true : false;
            case LAG_T:
                return (ld < 0) ? true : false;
            default:
                INTERNAL_ERROR
        }
        return false;
    }
    return m_arg->has(what, f, exact_idx);
}


bool
ex_e::hasidx(int idx) const
{
    return m_arg->hasidx(idx);
}


ptr_base
ex_e::lag(int l) const
{
    if (l == INT_MIN) return symbolic::internal::lag(m_arg, l);
    return create(symbolic::internal::lag(m_arg, l), m_lag + l);
}








