/*****************************************************************************
 * This file is a part of gEcon.                                             *
 *                                                                           *
 * (c) Chancellery of the Prime Minister of the Republic of Poland 2012-2015 *
 * (c) Grzegorz Klima, Karol Podemski, Kaja Retkiewicz-Wijtiwiak 2015-2018   *
 * License terms can be found in the file 'LICENCE'                          *
 *                                                                           *
 * Author: Grzegorz Klima                                                    *
 *****************************************************************************/

/** \file ex_idx.cpp
 * \brief Indexed expression.
 */


#include <ex_idx.h>
#include <ops.h>
#include <cmp.h>
#include <error.h>
#include <cmath>
#include <utils.h>
#include <climits>


using namespace symbolic;
using namespace symbolic::internal;



ptr_base
ex_idx::create(const idx_ex &ie, const ptr_base &e)
{
    if (!ie) return e;
    if (!e->hasidx(ie.m_id)) return e;
    if (e->type() == IDX) {
        if (e.get<ex_idx>()->has_ie(ie)) return e;
    }
    return ptr_base(new ex_idx(ie, e));
}

void
ex_idx::destroy(ex_base *ptr)
{
    delete ptr;
}


int
ex_idx::compare(const ex_idx &b) const
{
    int c = m_ie.compare(b.m_ie, false);
    if (c) return c;
    return symbolic::internal::compare(m_e, substidx(b.m_e, b.m_ie.get_id(), m_ie.get_id()));
}


std::string
ex_idx::str(int pflag, bool c_style) const
{
    if (pflag & INDEXING_ONLY) {
        return m_ie.str() + m_e->str(pflag, c_style);
    } else if (pflag & DROP_INDEXING) {
        return get_ptr()->str(pflag, c_style);
    }

    std::string res = m_ie.str();
    res += " ";
    return res + m_e->str(pflag, c_style);
}


std::string
ex_idx::strmap(const map_str_str &mss, bool c_style) const
{
    USER_ERROR("ex_idx::strmap() called")
    return m_e->strmap(mss, c_style);
}


std::string
ex_idx::tex(int pflag) const
{
    if ((pflag & DROP_IDX) && (pflag & CONVERT_IDX))
        USER_ERROR("invalid print flag in ex_idx::tex()")

    if (pflag & INDEXING_ONLY) {
        return m_ie.tex() + "\\colon\\quad " + m_e->tex(pflag);
    } else if (pflag & DROP_INDEXING) {
        return get_ptr()->tex(pflag);
    }
    int t = get_ptr()->type();
    if ((t == SYMBIDX) || (t == VARTIDX)) {
        return "\\left(" + m_e->tex(pflag) + "\\right)_{" + m_ie.tex() + '}';
    }
    return m_ie.tex() + "\\colon\\quad " + m_e->tex(pflag);
}


int
ex_idx::get_lag_max(bool stop_on_E) const
{
    return m_e->get_lag_max(stop_on_E);
}


int
ex_idx::get_lag_min(bool stop_on_E) const
{
    return m_e->get_lag_min(stop_on_E);
}



ptr_base
ex_idx::diff(const ptr_base &p) const
{
    return create(m_ie, m_e->diff(p));
}


bool
ex_idx::has(const ptr_base &what, search_flag f, bool exact_idx) const
{
    return m_e->has(what, f, exact_idx);
}


bool
ex_idx::hasidx(int idx) const
{
    return ((int) m_ie.m_id == idx) || (m_ie.m_excl_id == idx) || m_e->hasidx(idx);
}

bool
ex_idx::has_ie(const idx_ex &ie) const
{
    if (!m_ie.compare(ie)) return true;
    if (m_e->type() == IDX) return m_e.get<ex_idx>()->has_ie(ie);
    return false;
}


const ptr_base&
ex_idx::get_ptr() const
{
    if (m_e->type() == IDX) return m_e.get<ex_idx>()->get_ptr();
    return m_e;
}



