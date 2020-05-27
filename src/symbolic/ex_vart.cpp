/*****************************************************************************
 * This file is a part of gEcon.                                             *
 *                                                                           *
 * (c) Chancellery of the Prime Minister of the Republic of Poland 2012-2015 *
 * (c) Grzegorz Klima, Karol Podemski, Kaja Retkiewicz-Wijtiwiak 2015-2018   *
 * License terms can be found in the file 'LICENCE'                          *
 *                                                                           *
 * Author: Grzegorz Klima                                                    *
 *****************************************************************************/

/** \file ex_vart.cpp
 * \brief Time indexed variables.
 */

#include <error.h>
#include <ex_vart.h>
#include <ex_num.h>
#include <ops.h>
#include <cmp.h>
#include <stringhash.h>
#include <utils.h>
#include <climits>


using namespace symbolic;
using namespace symbolic::internal;

ex_vart::ex_vart(const std::string &n, int l)
    : ex_base(VART | ((l == INT_MIN) ? 0 : HAST) | SINGLE), m_lag(l)
{
    m_hash = stringhash::get_instance().get_hash(n);
}

ex_vart::ex_vart(unsigned i, int l)
    : ex_base(VART | ((l == INT_MIN) ? 0 : HAST) | SINGLE), m_lag(l)
{
    m_hash = i;
}



ptr_base
ex_vart::create(const std::string &n, int l)
{
    return ptr_base(new ex_vart(n, l));
}


ptr_base
ex_vart::create(unsigned i, int l)
{
    return ptr_base(new ex_vart(i, l));
}

void
ex_vart::destroy(ex_base *ptr)
{
    delete ptr;
}


ptr_base
ex_vart::copy0() const
{
    return ptr_base(new ex_vart(m_hash, 0));
}


int
ex_vart::compare(const ex_vart &b) const
{
    if (m_lag < b.m_lag) return -1;
    if (m_lag > b.m_lag) return 1;
    return compareT(m_hash, b.m_hash);
}


int
ex_vart::compare_name(const ex_vart &b) const
{
    return compareT(m_hash, b.m_hash);
}



std::string
ex_vart::get_name() const
{
    return stringhash::get_instance().get_str(m_hash);
}



std::string
ex_vart::str(int pflag, bool c_style) const
{
    if (pflag & INDEXING_ONLY) return std::string();
    std::string res = get_name();
    if (pflag & DROP_T) {
        return res;
    } else if (pflag & CONVERT_T) {
        if (m_lag == INT_MIN) return res + "__ss";
        if (m_lag < 0) return res + "__tm" + num2str(-m_lag);
        if (m_lag == 0) return res + "__t";
        return res + "__tp" + num2str(m_lag);
    } else {
        if (m_lag == INT_MIN) return res + "[ss]";
        else if (m_lag == 0) return res + "[]";
        else return res + '[' + num2str(m_lag) + ']';
    }
}


std::string
ex_vart::strmap(const map_str_str &mss, bool c_style) const
{
    std::string name = get_name();
    map_str_str::const_iterator it;
    it = mss.find(name);
    if (it == mss.end()) {
        return name;
    }
    return it->second;
}


std::string
ex_vart::tex(int pflag) const
{
    if (pflag & INDEXING_ONLY) return std::string();
    std::string name = get_name();
    if (pflag & DROP_T) {
        return str2tex(name);
    } else {
        if (m_lag == INT_MIN) return str2tex(name) + "_\\mathrm{ss}";
        if (m_lag == 0) return str2tex(name) + "_{t}";
        if (m_lag < 0) return str2tex(name) + "_{t" + num2str(m_lag) + '}';
        return str2tex(name) + "_{t+" + num2str(m_lag) + '}';
    }
}


int
ex_vart::get_lag_max(bool) const
{
    return m_lag;
}


int
ex_vart::get_lag_min(bool) const
{
    if (m_lag == INT_MIN) return INT_MAX;
    return m_lag;
}


ptr_base
ex_vart::diff(const ptr_base &p) const
{
    if (p->type() == VART) {
        const ex_vart *pp = p.get<ex_vart>();
        if (!compare(*pp)) return ex_num::one();
    }
    return ex_num::zero();
}



bool
ex_vart::has(const ptr_base &what, search_flag f, bool) const
{
    if (what->type() != VART) return false;
    const ex_vart *w = what.get<ex_vart>();
    if (compareT(m_hash, w->m_hash)) return false;
    switch (f) {
        case EXACT_T:
            return (m_lag == w->m_lag) ? true : false;
        case ANY_T:
            return true;
        case DIFF_T:
            return (m_lag != w->m_lag) ? true : false;
        case LEAD_T:
            return (m_lag > w->m_lag) ? true : false;
        case LAG_T:
            return (m_lag < w->m_lag) ? true : false;
        default:
            INTERNAL_ERROR
    }
    return false;
}


bool
ex_vart::hasidx(int) const
{
    return false;
}



ptr_base
ex_vart::lag(int l) const
{
    if (l == INT_MIN) return create(m_hash, l);
    return create(m_hash, m_lag + l);
}








