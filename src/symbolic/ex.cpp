/***********************************************************
 * (c) Kancelaria Prezesa Rady Ministrów 2012-2015         *
 * Treść licencji w pliku 'LICENCE'                        *
 *                                                         *
 * (c) Chancellery of the Prime Minister 2012-2015         *
 * License terms can be found in the file 'LICENCE'        *
 *                                                         *
 * Author: Grzegorz Klima                                  *
 ***********************************************************/

/** \file ex.cpp
 * \brief Expression.
 */

#include <ex.h>
#include <ex_num.h>
#include <ex_delta.h>
#include <ex_symb.h>
#include <ex_symbidx.h>
#include <ex_vart.h>
#include <ex_vartidx.h>
#include <ex_sum.h>
#include <ex_idx.h>
#include <ops.h>
#include <cmp.h>
#include <error.h>
#ifdef R_DLL
#include <R.h>
#else /* R_DLL */
#include <cmath>
#endif /* R_DLL */

using namespace symbolic;
using namespace symbolic::internal;


ex::ex() : m_ptr(ex_num::zero())
{
}


ex::ex(const Number &v) : m_ptr(ex_num::create(v))
{
}

ex::ex(bool fixed1, const std::string &i1, bool fixed2, const std::string &i2)
        : m_ptr(ex_delta::create(fixed1, i1, fixed2, i2))
{
}


ex::ex(const std::string &n) : m_ptr(ex_symb::create(n))
{
}


ex::ex(const std::string &n, bool fixed1, const std::string &i1)
        : m_ptr(ex_symbidx::create(n, fixed1, i1))
{
}


ex::ex(const std::string &n, bool fixed1, const std::string &i1,
       bool fixed2, const std::string &i2)
        : m_ptr(ex_symbidx::create(n, fixed1, i1, fixed2, i2))
{
}


ex::ex(const std::string &n, bool fixed1, const std::string &i1,
       bool fixed2, const std::string &i2, bool fixed3, const std::string &i3)
        : m_ptr(ex_symbidx::create(n, fixed1, i1, fixed2, i2, fixed3, i3))
{
}


ex::ex(const std::string &n, bool fixed1, const std::string &i1,
       bool fixed2, const std::string &i2, bool fixed3, const std::string &i3,
       bool fixed4, const std::string &i4)
        : m_ptr(ex_symbidx::create(n, fixed1, i1, fixed2, i2, fixed3, i3,
                                   fixed4, i4))
{
}


ex::ex(const std::string &n, int l) : m_ptr(ex_vart::create(n, l))
{
}


ex::ex(const std::string &n, int l, bool fixed1, const std::string &i1)
        : m_ptr(ex_vartidx::create(n, l, fixed1, i1))
{
}


ex::ex(const std::string &n, int l, bool fixed1, const std::string &i1,
       bool fixed2, const std::string &i2)
        : m_ptr(ex_vartidx::create(n, l, fixed1, i1, fixed2, i2))
{
}


ex::ex(const std::string &n, int l, bool fixed1, const std::string &i1,
       bool fixed2, const std::string &i2, bool fixed3, const std::string &i3)
        : m_ptr(ex_vartidx::create(n, l, fixed1, i1, fixed2, i2, fixed3, i3))
{
}


ex::ex(const std::string &n, int l, bool fixed1, const std::string &i1,
       bool fixed2, const std::string &i2, bool fixed3, const std::string &i3,
       bool fixed4, const std::string &i4)
        : m_ptr(ex_vartidx::create(n, l, fixed1, i1, fixed2, i2, fixed3, i3,
                                   fixed4, i4))
{
}


ex::ex(const idx_ex &ie, const ex &e) : m_ptr(ex_idx::create(ie, e.m_ptr))
{
}


ex&
ex::operator=(const ex &rhs)
{
    m_ptr = rhs.m_ptr;
    return *this;
}


bool
ex::operator==(const ex &rhs)
{
    return compare(m_ptr, rhs.m_ptr) ? false : true;
}


bool
ex::operator!=(const ex &rhs)
{
    return compare(m_ptr, rhs.m_ptr) ? true : false;
}


bool
ex::has(const ex &what, search_flag f, bool exact_idx) const
{
    unsigned t = what.m_ptr->type();
    if (t == IDX) {
        return m_ptr->has(what.m_ptr.get<ex_idx>()->get_ptr(), f, exact_idx);
    }
    if ((t != SYMB) && (t != SYMBIDX) && (t != VART) && (t != VARTIDX) && (t != EX)) {
        USER_ERROR("trying to search for something that is not a constant, nor a variable, nor an expected value")
    }
    return m_ptr->has(what.m_ptr, f, exact_idx);
}


bool
ex::has_Es() const
{
    return symbolic::internal::has_Es(m_ptr);
}


ex
ex::subst(const ex &what, const ex &with, bool all_leads_lags) const
{
    return ex(symbolic::internal::subst(m_ptr, what.m_ptr, with.m_ptr, all_leads_lags));
}


ex
ex::substidx(const std::string &what, const std::string &with, bool fix) const
{
    stringhash &ref = stringhash::get_instance();
    unsigned iwhat = ref.get_hash(what);
    int iwith = ref.get_hash(with);
    if (fix) iwith = -iwith;
    return ex(symbolic::internal::substidx(m_ptr, iwhat, iwith));
}


std::string
ex::str(int pflag) const
{
    return m_ptr->str(pflag);
}

std::string
ex::strmap(const map_str_str &mss) const
{
    return m_ptr->strmap(mss);
}

std::string
ex::tex(int pflag) const
{
    return m_ptr->tex(pflag);
}


int
ex::get_lag_max(bool stop_on_E) const
{
    return m_ptr->get_lag_max(stop_on_E);
}


int
ex::get_lag_min(bool stop_on_E) const
{
    return m_ptr->get_lag_min(stop_on_E);
}



bool
ex::validnum() const
{
    Number v = m_ptr->val();
#ifdef R_DLL
    return R_FINITE(v.val());
#else /* R_DLL */
    return std::isfinite(v.val());
#endif /* R_DLL */
}



std::ostream&
symbolic::operator<<(std::ostream &s, const ex &e)
{
    s << e.str() << '\n';
    return s;
}


ex
ex::operator-() const
{
    return ex(mk_neg(m_ptr));
}

ex
symbolic::operator+(const ex &lhs, const ex &rhs)
{
    return ex(mk_add(lhs.m_ptr, rhs.m_ptr));
}

ex
symbolic::operator-(const ex &lhs, const ex &rhs)
{
    return ex(mk_sub(lhs.m_ptr, rhs.m_ptr));
}

ex
symbolic::operator*(const ex &lhs, const ex &rhs)
{
    return ex(mk_mul(lhs.m_ptr, rhs.m_ptr));
}

ex
symbolic::operator/(const ex &lhs, const ex &rhs)
{
    return ex(mk_div(lhs.m_ptr, rhs.m_ptr));
}

ex
symbolic::pow(const ex &lhs, const ex &rhs)
{
    return ex(mk_pow(lhs.m_ptr, rhs.m_ptr));
}


ex
symbolic::func(symbolic::internal::func_code fc, const ex &e)
{
    return ex(mk_func(fc, e.m_ptr));
}


ex
symbolic::sum(const idx_ex &ie, const ex &e)
{
    return ex(ex_sum::create(ie, e.m_ptr));
}


ex
symbolic::sum(const ex &e)
{
    return ex(mk_sum(e.m_ptr));
}


ex
symbolic::prod(const idx_ex &ie, const ex &e)
{
    return ex(ex_prod::create(ie, e.m_ptr));
}


ex
symbolic::prod(const ex &e)
{
    return ex(mk_prod(e.m_ptr));
}



ex
symbolic::diff(const ex &exp, const ex &var)
{
    return ex(mk_diff(exp.m_ptr, var.m_ptr));
}


ex
symbolic::append_name(const ex &e, const std::string &s)
{
    return ex(internal::append_name(e.m_ptr, s));
}


ex
symbolic::add_idx(const ex &e, const idx_ex &ie)
{
    return ex(internal::add_idx(e.m_ptr, ie));
}


ex
symbolic::add_idx(const ex &e, const ex &ie)
{
    return ex(internal::add_idx(e.m_ptr, ie.m_ptr));
}


ex
symbolic::apply_idx(const ex &e, const ex &of)
{
    return ex(internal::apply_idx(e.m_ptr, of.m_ptr));
}


ex
symbolic::lag(const ex &e, int l)
{
    return ex(lag(e.m_ptr, l));
}


ex
symbolic::ss(const ex &e)
{
    return ex(ss(e.m_ptr));
}


ex
symbolic::E(const ex &e, int l)
{
    return ex(mk_E(e.m_ptr, l));
}


vec_ex
symbolic::expand(const ex &e)
{
    vec_ex res;
    unsigned t = e.m_ptr->type();
    if (t != IDX) {
        res.push_back(ex(symbolic::internal::expandsp(e.m_ptr)));
        return res;
    }
    const ex_idx *pt = e.m_ptr.get<ex_idx>();
    idx_ex ie = pt->get_ie();
    ptr_base ee = pt->get_e();
    if (ie.m_excl_id > 0)
        USER_ERROR("excluded index not evalued (substituted for) in idx_ex::evalidx")
    res.reserve(ie.m_set->size());
    idx_set_impl::const_iterator it;
    for (it = ie.m_set->begin(); it != ie.m_set->end(); ++it) {
        if (-ie.m_excl_id != (int) *it) {
            vec_ex res0 = expand(ex(substidx(ee, ie.m_id, -(int)*it)));
            for (unsigned j = 0; j < res0.size(); ++j) {
                res.push_back(res0[j]);
            }
        }
    }
    return res;
}

// vec_ex
// symbolic::expand(const vec_ex &ve)
// {
// }


ex
symbolic::drop_Es(const ex &e)
{
    return ex(internal::drop_Es(e.m_ptr));
}



void
symbolic::find_Es(const ex &e, set_ex &sex)
{
    find_Es(e.m_ptr, sex);
}


void
symbolic::collect(const ex &e, set_ex &vars, set_ex &parms)
{
    collect(e.m_ptr, vars, parms);
}

void
symbolic::collect(const vec_ex &ve, set_ex &vars, set_ex &parms)
{
    for (vec_ex::const_iterator it = ve.begin(); it != ve.end(); ++it)
        collect(*it, vars, parms);
}


void
symbolic::collect_lags(const ex &e, map_ex_int &map)
{
    if (e.get_lag_min() >= -1) return;
    collect_lags(e.m_ptr, map);
}


void
symbolic::collect_idx(const ex &e, std::set<unsigned> &iset)
{
    collect_idx(e.m_ptr, iset);
}


void
symbolic::collect_fidx(const ex &e, std::map<unsigned, unsigned> &imap)
{
    collect_fidx(e.m_ptr, imap);
}

void
symbolic::collect_fidx(const idx_ex &ie, std::map<unsigned, unsigned> &imap)
{
    symbolic::internal::collect_fidx(ie, imap);
}





triplet<bool, ex, ex>
symbolic::find_subst(const ex &e, const set_ex &se)
{
    ptr_base pb = e.m_ptr;

    if (pb->type() == VART) {
        ex v(pb.get<ex_vart>()->copy0());
        if (se.find(v) != se.end()) return triplet<bool, ex, ex>(true, e, ex());
        return triplet<bool, ex, ex>(false, ex(), ex());
    } else if (pb->type() == VARTIDX) {
        ex v(pb.get<ex_vartidx>()->copy0());
        if (se.find(v) != se.end()) return triplet<bool, ex, ex>(true, e, ex());
        return triplet<bool, ex, ex>(false, ex(), ex());
    }
    if (pb->type() != ADD) return triplet<bool, ex, ex>(false, ex(), ex());

    const ex_add *p = pb.get<ex_add>();
    unsigned n = p->no_ops(), i, c, m, l;

    for (i = 0, c = 0; i < n;) {
        if (p->get_ops()[i].second->type() == VART) {
            const ex_vart *pv = p->get_ops()[i].second.get<ex_vart>();
            ptr_base v = pv->copy0();
            if (se.find(ex(v)) != se.end()) {
                for (unsigned j = 0; j < n; ++j) {
                    if ((j != i) && (p->get_ops()[j].second->has(v, ANY_T, true))) goto next;
                }
                m = i;
                l = pv->get_lag();
                ++c;
                break;
            }
        } else if (p->get_ops()[i].second->type() == VARTIDX) {
            const ex_vartidx *pv = p->get_ops()[i].second.get<ex_vartidx>();
            ptr_base v = pv->copy0();
            if (se.find(ex(v)) != se.end()) {
                for (unsigned j = 0; j < n; ++j) {
                    if ((j != i) && (p->get_ops()[j].second->has(v, ANY_T, true))) goto next;
                }
                m = i;
                l = pv->get_lag();
                ++c;
                break;
            }
        }
next:
        ++i;
    }

    if (!c) return triplet<bool, ex, ex>(false, ex(), ex());

    ex what(p->get_ops()[m].second);
    ex coef(-p->get_ops()[m].first);
    ex with = e / coef + what;
    return triplet<bool, ex, ex>(true, lag(what, -l), lag(with, -l));
}



triplet<bool, ex, ex>
symbolic::find_par_eq_num(const ex &e)
{
    ptr_base pb = e.m_ptr;

    if (pb->type() == IDX) {
        return find_par_eq_num(ex(pb.get<ex_idx>()->get_ptr()));
    }
    if ((pb->type() == SYMB) || (pb->type() == SYMBIDX)) {
        return triplet<bool, ex, ex>(true, e, ex());
    }
    if (pb->type() != ADD) return triplet<bool, ex, ex>(false, ex(), ex());

    const ex_add *p = pb.get<ex_add>();
    unsigned n = p->no_ops();

    if (n != 2) return triplet<bool, ex, ex>(false, ex(), ex());

    if ((p->get_ops()[0].second->type() == NUM)
        && ((p->get_ops()[1].second->type() == SYMB) || (p->get_ops()[1].second->type() == SYMBIDX))
        && (p->get_ops()[1].first == 1.)) {
        ex par(p->get_ops()[1].second);
        ex val(p->get_ops()[0].second);
        return triplet<bool, ex, ex>(true, par, -val);
    }

    return triplet<bool, ex, ex>(false, ex(), ex());
}



bool
less_ex::operator()(const symbolic::ex &a, const symbolic::ex &b) const
{
    return (compare(a.m_ptr, b.m_ptr) < 0);
}







