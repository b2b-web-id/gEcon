/***********************************************************
 * (c) Kancelaria Prezesa Rady Ministrów 2012-2015         *
 * Treść licencji w pliku 'LICENCE'                        *
 *                                                         *
 * (c) Chancellery of the Prime Minister 2012-2015         *
 * License terms can be found in the file 'LICENCE'        *
 *                                                         *
 * Author: Grzegorz Klima                                  *
 ***********************************************************/

/** \file ex_prod.cpp
 * \brief PRoduct over indices.
 */


#include <ex_prod.h>
#include <cmp.h>
#include <ops.h>
#include <error.h>
#include <utils.h>


using namespace symbolic;
using namespace symbolic::internal;


ex_prod::ex_prod(const idx_ex &ie, const ptr_base &e)
    : ex_base(PROD | (e->flag() & HAST)), m_ie(ie), m_e(e)
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
ex_prod::create(const idx_ex &ie, const ptr_base &e)
{
    if (e->is0()) return e;
    if (!e->hasidx(ie.get_id())) {
        return mk_pow(e, ptr_base(new ex_sum(ie, ex_num::one())));
    }
    return ptr_base(new ex_prod(ie, e));
}


void
ex_prod::destroy(ex_base *ptr)
{
    delete ptr;
}



ptr_base
ex_prod::reindex(bool force) const
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
ex_prod::compare(const ex_prod &b) const
{
    int c = m_ie.compare(b.m_ie);
    if (c) return c;
    return symbolic::internal::compare(m_e, b.m_e);
}


std::string
ex_prod::str(int pflag) const
{
    if (pflag & INDEXING_ONLY) return std::string();
    if ((pflag & DROP_IDX) || (pflag & CONVERT_IDX))
        USER_ERROR("invalid print flag in ex_prod::str()")

    bool brace = true;
    if ((m_e->type() == PROD) ||( m_e->type() == POW) ||
        ((m_e->type() != ADD) && (m_e->flag() & SINGLE))) brace = false;
    std::string res = "PROD" + m_ie.str() + ' ';
    if (brace) res += "(";
    res += m_e->str(pflag);
    if (brace) res += ")";
    return res;
}


std::string
ex_prod::strmap(const map_str_str &mss) const
{
    USER_ERROR("ex_prod::strmap() called")
}


std::string
ex_prod::tex(int pflag) const
{
    if (pflag & INDEXING_ONLY) return std::string();
    if ((pflag & DROP_IDX) || (pflag & CONVERT_IDX))
        USER_ERROR("invalid print flag in ex_prod::tex()")

    bool brace = true;
    if ((m_e->type() == POW) ||
        ((m_e->type() != ADD) && (m_e->flag() & SINGLE))) brace = false;
    std::string res = "\\prod_{" + m_ie.tex() + "} ";
    if (brace) res += "\\left(";
    res += m_e->tex(pflag);
    if (brace) res += "\\right)";
    return res;
}


int
ex_prod::get_lag_max(bool stop_on_E) const
{
    return m_e->get_lag_max(stop_on_E);
}


int
ex_prod::get_lag_min(bool stop_on_E) const
{
    return m_e->get_lag_min(stop_on_E);
}



ptr_base
ex_prod::diff(const ptr_base &p) const
{
    if (p->hasidx(m_ie.m_id)) return reindex(true)->diff(p);
    ptr_base cpy = create(m_ie, m_e);
    ptr_base pder = ex_sum::create(m_ie, mk_div(m_e->diff(p), m_e));
    return mk_mul(cpy, pder);
}


bool
ex_prod::has(const ptr_base &what, search_flag f, bool exact_idx) const
{
    return m_e->has(what, f, exact_idx);
}


bool
ex_prod::hasidx(int idx) const
{
    return (m_ie.m_excl_id == idx) || m_e->hasidx(idx);
}


ptr_base
ex_prod::expand() const
{
    num_ex_pair_vec out;
    out.reserve(m_ie.m_set.size());
    set<unsigned>::const_iterator it;
    for (it = m_ie.m_set->begin(); it != m_ie.m_set->end(); ++it) {
        if (-m_ie.m_excl_id != (int) *it) {
            out.push_back(std::pair<Number, ptr_base>(1,
                    symbolic::internal::expandsp(substidx(m_e, m_ie.m_id, -*it))));
        }
    }
    return mk_mul(out);
}



