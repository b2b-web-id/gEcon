/***********************************************************
 * (c) Kancelaria Prezesa Rady Ministrów 2012-2015         *
 * Treść licencji w pliku 'LICENCE'                        *
 *                                                         *
 * (c) Chancellery of the Prime Minister 2012-2015         *
 * License terms can be found in the file 'LICENCE'        *
 *                                                         *
 * Author: Grzegorz Klima                                  *
 ***********************************************************/

/** \file ex_sum.cpp
 * \brief Sum over indices.
 */


#include <ex_sum.h>
#include <cmp.h>
#include <ops.h>
#include <error.h>
#include <utils.h>


using namespace symbolic;
using namespace symbolic::internal;


ex_sum::ex_sum(const idx_ex &ie, const ptr_base &e)
    : ex_base(SUM | (e->flag() & HAST)), m_ie(ie), m_e(e)
{
    if (m_ie.m_excl_id > 0) m_type |= HASIDX;
    if (!(m_type & HASIDX)) {
        std::set<unsigned> is;
        collect_idx(m_e, is);
        is.erase(m_ie.m_id);
        if (is.size()) m_type |= HASIDX;
    }
}


ptr_base
ex_sum::create(const idx_ex &ie, const ptr_base &e)
{
    if (e->is0()) return e;
    if (!e->hasidx(ie.get_id())) {
        return mk_mul(ptr_base(new ex_sum(ie, ex_num::one())), e);
    }
    if (e->type() == ADD) {
        const num_ex_pair_vec &in = e.get<ex_add>()->get_ops();
        num_ex_pair_vec out;
        unsigned i, n = in.size();
        out.reserve(n);
        for (i = 0; i < n; ++i) {
            out.push_back(num_ex_pair(in[i].first, ex_sum::create(ie, in[i].second)));
        }
        return mk_add(out);
    } else if (e->type() == MUL) {
        const num_ex_pair_vec &in = e.get<ex_mul>()->get_ops();
        num_ex_pair_vec out1, out2;
        unsigned i, n = in.size();
        out1.reserve(n);
        out2.reserve(n);
        for (i = 0; i < n; ++i) {
            if (in[i].second->hasidx(ie.get_id())) {
                out2.push_back(in[i]);
            } else {
                out1.push_back(in[i]);
            }
        }
        out1.push_back(num_ex_pair(1., ptr_base(new ex_sum(ie, mk_mul(out2)))));
        return mk_mul(out1);
    }
    return ptr_base(new ex_sum(ie, e));
}



void
ex_sum::destroy(ex_base *ptr)
{
    delete ptr;
}



ptr_base
ex_sum::reindex(bool force) const
{
    unsigned id, idn;
    id = m_ie.m_id;
    stringhash &ref = stringhash::get_instance();
    idn = id;
    if (force) {
        idn = ref.append_underscore(idn);
    }
    while (m_e->hasidx(idn)) { idn = ref.append_underscore(idn); }
    return create(idx_ex(idn, m_ie.m_set, m_ie.m_excl_id), substidx(m_e, id, (int) idn));
}



int
ex_sum::compare(const ex_sum &b) const
{
    int c = m_ie.compare(b.m_ie);
    if (c) return c;
    return symbolic::internal::compare(m_e, b.m_e);
}


std::string
ex_sum::str(int pflag) const
{
    if (pflag & INDEXING_ONLY) return std::string();
    if ((pflag & DROP_IDX) || (pflag & CONVERT_IDX))
        USER_ERROR("invalid print flag in ex_sum::str()")

    bool brace = true;
    if ((m_e->type() == SUM) || (m_e->type() == PROD) || (m_e->type() == MUL)
        || (m_e->flag() & SINGLE)) brace = false;
    std::string res = "SUM" + m_ie.str() + ' ';
    if (brace) res += "(";
    res += m_e->str(pflag);
    if (brace) res += ")";
    return res;
}


std::string
ex_sum::strmap(const map_str_str &mss) const
{
    USER_ERROR("ex_sum::strmap() called")
}


std::string
ex_sum::tex(int pflag) const
{
    if (pflag & INDEXING_ONLY) return std::string();
    if ((pflag & DROP_IDX) || (pflag & CONVERT_IDX))
        USER_ERROR("invalid print flag in ex_sum::str()")

    bool brace = true;
    if ((m_e->type() == SUM) || (m_e->type() == PROD) || (m_e->type() == MUL)
        || (m_e->flag() & SINGLE)) brace = false;
    std::string res = "\\sum_{" + m_ie.tex() + "} ";
    if (brace) res += "\\left(";
    res += m_e->tex(pflag);
    if (brace) res += "\\right)";
    return res;
}


int
ex_sum::get_lag_max(bool stop_on_E) const
{
    return m_e->get_lag_max(stop_on_E);
}


int
ex_sum::get_lag_min(bool stop_on_E) const
{
    return m_e->get_lag_min(stop_on_E);
}



ptr_base
ex_sum::diff(const ptr_base &p) const
{
    if (p->hasidx(m_ie.m_id)) return reindex(true)->diff(p);
    return create(m_ie, m_e->diff(p));
}


bool
ex_sum::has(const ptr_base &what, search_flag f, bool exact_idx) const
{
    return m_e->has(what, f, exact_idx);
}


bool
ex_sum::hasidx(int idx) const
{
    return (m_ie.m_excl_id == idx) || m_e->hasidx(idx);
}


ptr_base
ex_sum::expand() const
{
    num_ex_pair_vec out;
    out.reserve(m_ie.m_set.size());
    set<unsigned>::const_iterator it;
    for (it = m_ie.m_set->begin(); it != m_ie.m_set->end(); ++it) {
        if (-m_ie.m_excl_id != (int) *it) {
            out.push_back(std::pair<Number, ptr_base>(1,
                    expandsp(substidx(m_e, m_ie.m_id, -*it))));
        }
    }
    return mk_add(out);
}


ptr_base
ex_sum::reduceidx(const idx_ex &ie) const
{
    unsigned t = m_e->type(), f = m_e->flag();
    if (!m_ie.compare(ie, false)) {
        int i1 = ie.m_id, i2 = m_ie.m_id;
        bool red = false;
        if (t == MUL) {
            const num_ex_pair_vec &in = m_e.get<ex_mul>()->get_ops();
            num_ex_pair_vec out;
            unsigned i, n = in.size();
            out.reserve(n);
            for (i = 0; i < n; ++i) {
                if (in[i].second->type() == DELTA) {
                    if (in[i].second->hasidx(i1) && in[i].second->hasidx(i2)) {
                        red = true;
                    } else {
                        out.push_back(in[i]);
                    }
                } else {
                    out.push_back(in[i]);
                }
            }
            if (red) {
                return symbolic::internal::reduceidx(substidx(mk_mul(out), i2, i1), ie);
            }
        } else if ((t == ADD) && (f & SINGLE)) {
            const num_ex_pair_vec &in = m_e.get<ex_add>()->get_ops();
                if (in[0].second->type() == DELTA) {
                    if (in[0].second->hasidx(i1) && in[0].second->hasidx(i2)) {
                        return ex_num::create(in[0].first);
                    }
                }
        } else if (t == DELTA) {
            if (m_e->hasidx(i1) && m_e->hasidx(i2)) return ex_num::one();
        }
    }
    return create(m_ie, symbolic::internal::reduceidx(m_e, ie));
}













