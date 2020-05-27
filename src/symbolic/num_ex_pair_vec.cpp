/*****************************************************************************
 * This file is a part of gEcon.                                             *
 *                                                                           *
 * (c) Chancellery of the Prime Minister of the Republic of Poland 2012-2015 *
 * (c) Grzegorz Klima, Karol Podemski, Kaja Retkiewicz-Wijtiwiak 2015-2018   *
 * License terms can be found in the file 'LICENCE'                          *
 *                                                                           *
 * Author: Grzegorz Klima                                                    *
 *****************************************************************************/

/** \file num_ex_pair_vec.cpp
 * \brief Vector of scalar-expression pairs used in construction
 * of classes representing addition and multiplication.
 */

#include <num_ex_pair_vec.h>
#include <ex_num.h>
#include <ex_pow.h>
#include <ops.h>
#include <cmp.h>
#include <cmath>
#include <climits>
#include <algorithm>


using namespace symbolic;
using namespace symbolic::internal;


num_ex_pair_vec::num_ex_pair_vec(const ptr_base &p)
{
    push_back(std::pair<Number, ptr_base>(1., p));
}


num_ex_pair_vec::num_ex_pair_vec(const Number &s, const ptr_base &p)
{
    push_back(std::pair<Number, ptr_base>(s, p));
}


num_ex_pair_vec::num_ex_pair_vec(const Number &s, const num_ex_pair_vec &p, ex_type t)
{
    int i, n = p.size();
    *this = p;
    for (i = 0; i < n; ++i) (*this)[i].first *= s;
}


num_ex_pair_vec::num_ex_pair_vec(const num_ex_pair_vec &a,
                                 const num_ex_pair_vec &b, ex_type t)
{
    num_ex_pair_vec_base::const_iterator ia, ib;

    insert(end(), a.begin(), a.end());
    insert(end(), b.begin(), b.end());
    reduce(t);
}


int
num_ex_pair_vec::compare(const num_ex_pair_vec &b) const
{
    int c;
    if ((c = compareT(size(), b.size()))) return c;

    for (num_ex_pair_vec_base::const_iterator i = begin(), ib = b.begin();
         i != end(); ++i, ++ib) {
        c = symbolic::internal::compare((*i).second, (*ib).second);
        if (c) return c;
        c = compareT((*i).first, (*ib).first);
        if (c) return c;
    }

    return 0;
}


namespace {

struct cmp_pair {
    bool operator()(const num_ex_pair &a, const num_ex_pair &b) {
        return (compare(a.second, b.second) < 0) ? true : false;
    }
}; /* struct cpm_pair */

} /* namespace */


void
num_ex_pair_vec::reduce(ex_type t)
{
    // sort
    std::sort(begin(), end(), cmp_pair());

    // reduce numbers
    num_ex_pair_vec_base::iterator i = begin();
    if (i->second->type() == NUM) {
        Number x = 0.;
        if (t == MUL) x = 1.;
        if (i->first != 1.) {
            if (t == ADD) x = i->first * i->second->val();
            if (t == MUL) x = pow(i->second->val(), i->first);
        } else x = i->second->val();
        ++i;
        while ((size() > 1) && (i->second->type() == NUM)) {
            if (t == ADD) x += i->first * i->second->val();
            if (t == MUL) x *= pow(i->second->val(), i->first);
            erase(i);
        }
        if ((t == ADD) && (x == 0.)) {
            erase(begin());
        } else if ((t == MUL) && (x == 0.)) {
            clear();
            return;
        } else if ((t == MUL) && (x == 1.)) {
            erase(begin()); return;
        } else {
            i = begin();
            i->first = 1.;
            i->second = ex_num::create(x);
        }
    }

    // if add handle case 2 * (0.5 * x* y) and a - b * c
    if (t == ADD) {
        for (num_ex_pair_vec_base::iterator it = begin(); it < end(); ++it) {
            if (it->second->type() == MUL) {
                ptr_base fst = it->second.get<ex_mul>()->get_ops()[0].second;
                if (fst->type() == NUM) {
                    it->first *= fst->val();
                    const ex_mul *pmul = it->second.get<ex_mul>();
                    num_ex_pair_vec newops = pmul->get_ops();
                    newops.erase(newops.begin());
                    it->second = ex_mul::create(newops);
                }
            }
        }
    }

    // collect similar terms
    for (num_ex_pair_vec_base::iterator ia = begin(), ib = begin() + 1; ib < end();) {
        int c = symbolic::internal::compare(ia->second, ib->second);
        if (c) {
            ++ia;
            ++ib;
        } else {
            ia->first += ib->first;
            erase(ib);
        }
    }

    // remove terms with coeff 0
    for (i = begin(); i != end();) {
        if (i->first == 0.) {
            erase(i);
        } else ++i;
    }

    // if mul collect remove powers of Kronecker deltas
    if (t == MUL) {
        for (i = begin(); i != end(); ++i) {
            if ((i->second->type() == DELTA) && (i->first != 1.)) {
                i->first = 1.;
            }
        }
    }

    // if mul collect powers with same base
    if (t == MUL) {
        if (size() < 2) return;
        for (num_ex_pair_vec_base::iterator ia = begin(), ib = begin() + 1; ib < end(); ++ib, ++ia) {
            if ((ia->second->type() == POW) && (ib->second->type() == POW)) {
                const ex_pow *p1 = ia->second.get<ex_pow>();
                const ex_pow *p2 = ib->second.get<ex_pow>();
                if (!symbolic::internal::compare(p1->get_base(), p2->get_base())) {
                    ptr_base tmp = mk_pow(p1->get_base(),
                                    mk_add(ex_add::create(ia->first, p1->get_exp()),
                                        ex_add::create(ib->first, p2->get_exp())));
                    ia->second = tmp;
                    ia->first = 1.;
                    --ia;
                    erase(ib);
                }
            }
        }
    }

    // if mul that is add reduce to add
    if ((t == MUL) && (size() == 2)) {
        num_ex_pair_vec_base::iterator it = begin();
        if ((it->second->type() == NUM) && ((++it)->second->type() == ADD)
            && (it->first == 1.)) {
            ptr_base pb = ex_add::create(begin()->second->val(), it->second);
            erase(end() - 1);
            begin()->first = 1.;
            begin()->second = pb;
        }
    }
}


int
num_ex_pair_vec::get_lag_max(bool stop_on_E) const
{
    int l = INT_MIN;
    if (!size()) return l;
    num_ex_pair_vec_base::const_iterator i;
    for (i = begin(); i < end(); ++i) l = std::max(l, i->second->get_lag_max(stop_on_E));
    return l;
}


int
num_ex_pair_vec::get_lag_min(bool stop_on_E) const
{
    int l = INT_MAX;
    if (!size()) return l;
    num_ex_pair_vec_base::const_iterator i;
    for (i = begin(); i < end(); ++i) l = std::min(l, i->second->get_lag_min(stop_on_E));
    return l;
}


bool
num_ex_pair_vec::has(const ptr_base &what, search_flag f, bool exact_idx) const
{
    if (!size()) return false;
    num_ex_pair_vec_base::const_iterator i;
    for (i = begin(); i != end(); ++i) {
        if (i->second->has(what, f, exact_idx)) return true;
    }
    return false;
}


bool
num_ex_pair_vec::hast() const
{
    if (!size()) return false;
    num_ex_pair_vec_base::const_iterator i;
    for (i = begin(); i != end(); ++i) {
        if (i->second->flag() & HAST) return true;
    }
    return false;
}


bool
num_ex_pair_vec::hasidx(int idx) const
{
    if (!size()) return false;
    num_ex_pair_vec_base::const_iterator i;
    if (idx) {
        for (i = begin(); i != end(); ++i) {
            if (i->second->hasidx(idx)) return true;
        }
    } else {
        for (i = begin(); i != end(); ++i) {
            if (i->second->flag() & HASIDX) return true;
        }
    }
    return false;
}




