/***********************************************************
 * (c) Kancelaria Prezesa Rady Ministrów 2012-2015         *
 * Treść licencji w pliku 'LICENCE'                        *
 *                                                         *
 * (c) Chancellery of the Prime Minister 2012-2015         *
 * License terms can be found in the file 'LICENCE'        *
 *                                                         *
 * Author: Grzegorz Klima                                  *
 ***********************************************************/

/** \file ex_add.cpp
 * \brief Addition.
 */

#include <ops.h>
#include <ex_add.h>
#include <ex_num.h>
#include <utils.h>
#include <error.h>
#include <iostream>



using namespace symbolic::internal;



ex_add::ex_add(const Number &s, const ptr_base &p)
        : ex_base(ADD)
{
    unsigned t = p->type();

    if (t == ADD) {
        m_ops = num_ex_pair_vec(s, p.get<ex_add>()->m_ops, ADD);
        m_ops.reduce(ADD);
        update_flags();
    } else {
        m_ops = num_ex_pair_vec(s, p);
        m_type |= SINGLE;
        if (p->flag() & HAST) m_type |= HAST;
        if (p->flag() & HASIDX) m_type |= HASIDX;
        if (s == 1.) m_type |= SCAL1;
    }
}


ex_add::ex_add(const ptr_base &a, const ptr_base &b)
       : ex_base(ADD)
{
    num_ex_pair_vec aa, bb;
    if (a->type() == ADD) {
        const ex_add *p = a.get<ex_add>();
        aa = p->m_ops;
    } else {
        aa = num_ex_pair_vec(a);
    }
    if (b->type() == ADD) {
        const ex_add *p = b.get<ex_add>();
        bb = p->m_ops;
    } else {
        bb = num_ex_pair_vec(b);
    }

    m_ops = num_ex_pair_vec(aa, bb, ADD);
    update_flags();
}



ex_add::ex_add(const num_ex_pair_vec &ops, bool try_reduce)
    : ex_base(ADD), m_ops(ops)
{
    if (try_reduce) m_ops.reduce(ADD);
    update_flags();
}


void
ex_add::update_flags()
{
    if (m_ops.size() == 0) {
        m_type |= EMPTY;
        return;
    }

    if (m_ops.size() == 1) {
        m_type |= SINGLE;
        if (m_ops[0].first == 1.) m_type |= SCAL1;
    }

    if (m_ops.hast()) m_type |= HAST;
    if (m_ops.hasidx()) m_type |= HASIDX;
}



ptr_base
ex_add::create(const Number &s, const ptr_base &p)
{
    return ptr_base(new ex_add(s, p));
}


ptr_base
ex_add::create(const ptr_base &a, const ptr_base &b)
{
    return ptr_base(new ex_add(a, b));
}


ptr_base
ex_add::create(const num_ex_pair_vec &ops, bool try_reduce)
{
    return ptr_base(new ex_add(ops, try_reduce));
}


void
ex_add::destroy(ex_base *ptr)
{
    delete ptr;
}



int
ex_add::compare(const ex_add &b) const
{
    return m_ops.compare(b.m_ops);
}



std::string
ex_add::str(int pflag) const
{
    if (pflag & INDEXING_ONLY) return std::string();
    size_t i = 0, n = m_ops.size();
    Number f;
    std::string res;

    f = m_ops[i].first;
    if (f == 1.) res += m_ops[i].second->str(pflag);
    else if (f == -1.) res += '-' + m_ops[i].second->str(pflag);
    else {
        res += f.str();
        res += " * ";
        res += m_ops[i].second->str(pflag);
    }

    for (++i; i < n; ++i) {
        f = m_ops[i].first;
        if (f == 1.) {
            res += " + ";
            res += m_ops[i].second->str(pflag);
        } else if (f == -1.) {
            res += " - ";
            res += m_ops[i].second->str(pflag);
        } else if (f > 0.) {
            res += " + ";
            res += f.str();
            res += " * ";
            res += m_ops[i].second->str(pflag);
        } else {
            res += " - ";
            res += (-f).str();
            res += " * ";
            res += m_ops[i].second->str(pflag);
        }
    }

    return res;
}



std::string
ex_add::strmap(const map_str_str &mss) const
{
    size_t i = 0, n = m_ops.size();
    Number f;
    std::string res;

    f = m_ops[i].first;
    if (f == 1.) res += m_ops[i].second->strmap(mss);
    else if (f == -1.) res += '-' + m_ops[i].second->strmap(mss);
    else {
        res += f.str();
        res += " * ";
        res += m_ops[i].second->strmap(mss);
    }

    for (++i; i < n; ++i) {
        f = m_ops[i].first;
        if (f == 1.) {
            res += " + ";
            res += m_ops[i].second->strmap(mss);
        } else if (f == -1.) {
            res += " - ";
            res += m_ops[i].second->strmap(mss);
        } else if (f > 0.) {
            res += " + ";
            res += f.str();
            res += " * ";
            res += m_ops[i].second->strmap(mss);
        } else {
            res += " - ";
            res += (-f).str();
            res += " * ";
            res += m_ops[i].second->strmap(mss);
        }
    }

    return res;
}



std::string
ex_add::tex(int pflag) const
{
    if (pflag & INDEXING_ONLY) return std::string();
    size_t i = 0, n = m_ops.size();
    Number f;
    std::string res;

    f = m_ops[i].first;
    if (f == 1.) res += m_ops[i].second->tex(pflag);
    else if (f == -1.) res += '-' + m_ops[i].second->tex(pflag);
    else {
        res += f.tex();
        res += m_ops[i].second->tex(pflag);
    }

    for (++i; i < n; ++i) {
        f = m_ops[i].first;
        if (f == 1.) {
            res += " + ";
            res += m_ops[i].second->tex(pflag);
        } else if (f == -1.) {
            res += " - ";
            res += m_ops[i].second->tex(pflag);
        } else if (f > 0.) {
            res += " + ";
            res += f.tex();
            res += m_ops[i].second->tex(pflag);
        } else {
            res += " - ";
            res += (-f).tex();
            res += m_ops[i].second->tex(pflag);
        }
    }

    return res;
}


int
ex_add::get_lag_max(bool stop_on_E) const
{
    return m_ops.get_lag_max(stop_on_E);
}


int
ex_add::get_lag_min(bool stop_on_E) const
{
    return m_ops.get_lag_min(stop_on_E);
}



ptr_base
ex_add::diff(const ptr_base &p) const
{
    unsigned i, n = m_ops.size();
    if (!n) return ex_num::zero();

    num_ex_pair_vec ders;
    ders.reserve(n);

    for (i = 0; i < n; ++i) {
        ders.push_back(num_ex_pair(m_ops[i].first, m_ops[i].second->diff(p)));
    }

    return mk_add(ders);
}


bool
ex_add::has(const ptr_base &what, search_flag f, bool exact_idx) const
{
    return m_ops.has(what, f, exact_idx);
}


bool
ex_add::hasidx(int idx) const
{
    return m_ops.hasidx(idx);
}




















































