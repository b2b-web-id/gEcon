/*****************************************************************************
 * This file is a part of gEcon.                                             *
 *                                                                           *
 * (c) Chancellery of the Prime Minister of the Republic of Poland 2012-2015 *
 * (c) Grzegorz Klima, Karol Podemski, Kaja Retkiewicz-Wijtiwiak 2015-2018   *
 * License terms can be found in the file 'LICENCE'                          *
 *                                                                           *
 * Author: Grzegorz Klima                                                    *
 *****************************************************************************/

/** \file ex_mul.cpp
 * \brief Multiplication.
 */

#include <ex_add.h>
#include <ex_mul.h>
#include <ex_pow.h>
#include <ex_num.h>
#include <ops.h>
#include <utils.h>
#include <cmp.h>


using namespace symbolic::internal;




ex_mul::ex_mul(const Number &s, const ptr_base &p)
        : ex_base(MUL)
{
    if (p->type() != MUL) {
        m_ops = num_ex_pair_vec(s, p);
        m_type |= SINGLE;
        if (s == 1.) m_type |= SCAL1;
        if (p->flag() & HAST) m_type |= HAST;
        if (p->flag() & HASIDX) m_type |= HASIDX;
        return;
    }

    m_ops = num_ex_pair_vec(s, p.get<ex_mul>()->m_ops, MUL);
    update_flags();
}



ex_mul::ex_mul(const ptr_base &a, const ptr_base &b)
       : ex_base(MUL)
{
    num_ex_pair_vec aa, bb;
    if (a->type() == MUL) {
        const ex_mul *p = a.get<ex_mul>();
        aa = p->m_ops;
    } else {
        aa = num_ex_pair_vec(a);
    }
    if (b->type() == MUL) {
        const ex_mul *p = b.get<ex_mul>();
        bb = p->m_ops;
    } else {
        bb = num_ex_pair_vec(b);
    }

    m_ops = num_ex_pair_vec(aa, bb, MUL);

    update_flags();
}



ex_mul::ex_mul(const num_ex_pair_vec &ops, bool try_reduce) : ex_base(MUL), m_ops(ops)
{
    if (try_reduce) m_ops.reduce(MUL);
    update_flags();
}



void
ex_mul::update_flags()
{
    m_type &= 0xffffff00;
    if (m_ops.size() == 0) m_type |= EMPTY;
    else if (m_ops.size() == 1) {
        m_type |= SINGLE;
        if (m_ops[0].first == 1.) m_type |= SCAL1;
    } else if ((m_ops.size() == 2) && (m_ops[0].second->type() == NUM)
        && (m_ops[1].first == 1.)) {
        m_type |= SPOW1;
    }

    if (m_ops.hast()) m_type |= HAST;
    if (m_ops.hasidx()) m_type |= HASIDX;
}




ptr_base
ex_mul::create(const Number &s, const ptr_base &p)
{
    return ptr_base(new ex_mul(s, p));
}



ptr_base
ex_mul::create(const ptr_base &a, const ptr_base &b)
{
    return ptr_base(new ex_mul(a, b));
}



ptr_base
ex_mul::create(const num_ex_pair_vec &ops, bool try_reduce)
{
    return ptr_base(new ex_mul(ops, try_reduce));
}



void
ex_mul::destroy(ex_base *ptr)
{
    delete ptr;
}



int
ex_mul::compare(const ex_mul &b) const
{
    return m_ops.compare(b.m_ops);
}



std::string
ex_mul::str(int pflag, bool c_style) const
{
    if (pflag & INDEXING_ONLY) return std::string();
    size_t i = 0, n = m_ops.size();
    Number f;
    bool brace;
    std::string res;

    if ((n) && (m_ops[i].second->ism1())) {
        res = '-';
        ++i;
    }

    brace = (m_ops[i].second->flag() & SINGLE) ? false : true;
    if ((m_ops[i].second->type() == POW) && (m_ops[i].first != 1.)) brace = true;
    if (brace) {
        if ((m_ops[i].second->type() == EX) && (m_ops[i].first == 1.))
            brace = false;
    }
    f = m_ops[i].first;

    if (c_style) {
        if (f != 1.) res += "pow(";
        if (brace) res += '(';
        res += m_ops[i].second->str(pflag);
        if (brace) res += ')';
        if (f != 1.) {
            res += ", ";
            res += f.str();
            res += ")";
        }
    } else {
        if (brace) res += '(';
        res += m_ops[i].second->str(pflag, c_style);
        if (brace) res += ')';
        if (f != 1.) {
            res += "^";
            res += f.str();
        }
    }

    for (++i; i < n; ++i) {
        res += " * ";
        brace = (m_ops[i].second->flag() & SINGLE) ? false : true;
        if ((m_ops[i].second->type() == POW) && (m_ops[i].first != 1.)) brace = true;
        if (brace) {
            if ((m_ops[i].second->type() == EX) && (m_ops[i].first == 1.))
                brace = false;
        }
        f = m_ops[i].first;
        if (c_style) {
            if (f != 1.) res += "pow(";
            if (brace) res += '(';
            res += m_ops[i].second->str(pflag, c_style);
            if (brace) res += ')';
            if (f != 1.) {
                res += ", ";
                res += f.str();
                res += ")";
            }
        } else {
            if (brace) res += '(';
            res += m_ops[i].second->str(pflag, c_style);
            if (brace) res += ')';
            if (f != 1.) {
                res += "^";
                res += f.str();
            }
        }
    }

    return res;
}



std::string
ex_mul::strmap(const map_str_str &mss, bool c_style) const
{
    size_t i = 0, n = m_ops.size();
    Number f;
    bool brace;
    std::string res;

    if ((n) && (m_ops[i].second->ism1())) {
        res = '-';
        ++i;
    }

    brace = (m_ops[i].second->flag() & SINGLE) ? false : true;
    if ((m_ops[i].second->type() == POW) && (m_ops[i].first != 1.)) brace = true;
    if (brace) {
        if ((m_ops[i].second->type() == EX) && (m_ops[i].first == 1.))
            brace = false;
    }
    f = m_ops[i].first;

    if (c_style) {
        if (f != 1.) res += "pow(";
        if (brace) res += '(';
        res += m_ops[i].second->strmap(mss, c_style);
        if (brace) res += ')';
        if (f != 1.) {
            res += ", ";
            res += f.str();
            res += ")";
        }
    } else {
        if (brace) res += '(';
        res += m_ops[i].second->strmap(mss, c_style);
        if (brace) res += ')';
        if (f != 1.) {
            res += "^";
            res += f.str();
        }
    }

    for (++i; i < n; ++i) {
        res += " * ";
        brace = (m_ops[i].second->flag() & SINGLE) ? false : true;
        if ((m_ops[i].second->type() == POW) && (m_ops[i].first != 1.)) brace = true;
        if (brace) {
            if ((m_ops[i].second->type() == EX) && (m_ops[i].first == 1.))
                brace = false;
        }
        f = m_ops[i].first;
        if (c_style) {
            if (f != 1.) res += "pow(";
            if (brace) res += '(';
            res += m_ops[i].second->strmap(mss, c_style);
            if (brace) res += ')';
            if (f != 1.) {
                res += ", ";
                res += f.str();
                res += ")";
            }
        } else {
            if (brace) res += '(';
            res += m_ops[i].second->strmap(mss, c_style);
            if (brace) res += ')';
            if (f != 1.) {
                res += "^";
                res += f.str();
            }
        }
    }

    return res;
}



std::string
ex_mul::tex(int pflag) const
{
    if (pflag & INDEXING_ONLY) return std::string();
    size_t i = 0, n = m_ops.size();
    Number f;
    bool brace;
    std::string res;

    if ((n) && (m_ops[i].second->ism1())) {
        res = '-';
        ++i;
    }
    --n;
    for (; i <= n; ++i) {
        brace = (m_ops[i].second->flag() & SINGLE) ? false : true;
        if ((m_ops[i].second->type() == POW) && (m_ops[i].first != 1.)) brace = true;
        if (brace) {
            if ((m_ops[i].second->type() == EX) && (m_ops[i].first == 1.))
                brace = false;
        }
        f = m_ops[i].first;
        if (brace) res += "\\left("; else res += '{';
        res += m_ops[i].second->tex(pflag);
        if (brace) res += "\\right)"; else res += '}';
        if (f != 1.) {
            res += "^{";
            res += f.tex();
            res += "}";
        }
        if (i != n) res += ' ';
    }

    return res;
}


int
ex_mul::get_lag_max(bool stop_on_E) const
{
    return m_ops.get_lag_max(stop_on_E);
}


int
ex_mul::get_lag_min(bool stop_on_E) const
{
    return m_ops.get_lag_min(stop_on_E);
}



ptr_base
ex_mul::diff(const ptr_base &p) const
{
    unsigned i, n = m_ops.size();
    if (!n) return ex_num::zero();

    if (n == 1) {
        return mk_mul(mk_mul(m_ops[0].second->diff(p), ex_num::create(m_ops[0].first)),
                      mk_pow(m_ops[0].second, ex_num::create(m_ops[0].first - 1)));
    }

    num_ex_pair_vec prods, ders;
    ders.reserve(n);

    for (i = 0; i < n; ++i) {
        prods = m_ops;
        prods.erase(prods.begin() + i);
        ptr_base der = mk_mul(mk_mul(m_ops[i].second->diff(p),
                                        ex_num::create(m_ops[i].first)),
                              mk_mul(mk_pow(m_ops[i].second, ex_num::create(m_ops[i].first - 1)),
                                        mk_mul(prods)));
        ders.push_back(num_ex_pair(1, der));
    }

    return mk_add(ders);
}


bool
ex_mul::has(const ptr_base &what, search_flag f, bool exact_idx) const
{
    return m_ops.has(what, f, exact_idx);
}


bool
ex_mul::hasidx(int idx) const
{
    return m_ops.hasidx(idx);
}











