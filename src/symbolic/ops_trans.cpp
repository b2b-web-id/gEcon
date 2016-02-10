/***********************************************************
 * (c) Kancelaria Prezesa Rady Ministrów 2012-2015         *
 * Treść licencji w pliku 'LICENCE'                        *
 *                                                         *
 * (c) Chancellery of the Prime Minister 2012-2015         *
 * License terms can be found in the file 'LICENCE'        *
 *                                                         *
 * Author: Grzegorz Klima                                  *
 ***********************************************************/

 /** \file ops.cpp
 * \brief Operations.
 */


#include <ex_num.h>
#include <ex_symb.h>
#include <ex_vart.h>
#include <ex_add.h>
#include <ex_mul.h>
#include <ex_pow.h>
#include <ex_func.h>
#include <ex_e.h>
#include <ex_sum.h>
#include <cmp.h>
#include <error.h>
#include <ops.h>
#include <iostream>
#include <cmath>
#include <climits>


using namespace symbolic;
using namespace symbolic::internal;




ptr_base
symbolic::internal::lag(const ptr_base &p, int l)
{
    if (!l) return p;
    if (!(p->flag() & HAST)) return p;
    unsigned t = p->type();
    if (t == VART) return p.get<ex_vart>()->lag(l);
    else if (t == VARTIDX) return p.get<ex_vartidx>()->lag(l);
    else if (t == EX) return p.get<ex_e>()->lag(l);
    else if (t == ADD) {
        const num_ex_pair_vec &in = p.get<ex_add>()->get_ops();
        num_ex_pair_vec out;
        unsigned i, n = in.size();
        out.reserve(n);
        for (i = 0; i < n; ++i)
            out.push_back(std::pair<Number, ptr_base>(in[i].first,
                                                        lag(in[i].second, l)));
        return mk_add(out);
    } else if (t == MUL) {
        const num_ex_pair_vec &in = p.get<ex_mul>()->get_ops();
        num_ex_pair_vec out;
        unsigned i, n = in.size();
        out.reserve(n);
        for (i = 0; i < n; ++i)
            out.push_back(std::pair<Number, ptr_base>(in[i].first,
                                                        lag(in[i].second, l)));
        return mk_mul(out);
    } else if (t == POW) {
        const ex_pow *pt = p.get<ex_pow>();
        return mk_pow(lag(pt->get_base(), l), lag(pt->get_exp(), l));
    } else if (t == FUN) {
        const ex_func *pt = p.get<ex_func>();
        return mk_func(pt->get_code(), lag(pt->get_arg(), l));
    } else if (t == IDX) {
        const ex_idx *pt = p.get<ex_idx>();
        return ex_idx::create(pt->get_ie(), lag(pt->get_e(), l));
    } else if (t == SUM) {
        const ex_sum *pt = p.get<ex_sum>();
        return ex_sum::create(pt->get_ie(), lag(pt->get_e(), l));
    } else if (t == PROD) {
        const ex_prod *pt = p.get<ex_prod>();
        return ex_prod::create(pt->get_ie(), lag(pt->get_e(), l));
    } else INTERNAL_ERROR
}



ptr_base
symbolic::internal::ss(const ptr_base &p)
{
    if (!(p->flag() & HAST)) return p;
    unsigned t = p->type();
    if (t == VART)
        return p.get<ex_vart>()->lag(INT_MIN);
    else if (t == VARTIDX) {
        return p.get<ex_vartidx>()->lag(INT_MIN);
    } else if (t == EX) return ss(p.get<ex_e>()->get_arg());
    else if (t == ADD) {
        const num_ex_pair_vec &in = p.get<ex_add>()->get_ops();
        num_ex_pair_vec out;
        unsigned i, n = in.size();
        out.reserve(n);
        for (i = 0; i < n; ++i)
            out.push_back(num_ex_pair(in[i].first, ss(in[i].second)));
        return mk_add(out);
    } else if (t == MUL) {
        const num_ex_pair_vec &in = p.get<ex_mul>()->get_ops();
        num_ex_pair_vec out;
        unsigned i, n = in.size();
        out.reserve(n);
        for (i = 0; i < n; ++i)
            out.push_back(std::pair<Number, ptr_base>(in[i].first, ss(in[i].second)));
        return mk_mul(out);
    } else if (t == POW) {
        const ex_pow *pt = p.get<ex_pow>();
        return mk_pow(ss(pt->get_base()), ss(pt->get_exp()));
    } else if (t == FUN) {
        const ex_func *pt = p.get<ex_func>();
        return mk_func(pt->get_code(), ss(pt->get_arg()));
    } else if (t == IDX) {
        const ex_idx *pt = p.get<ex_idx>();
        return ex_idx::create(pt->get_ie(), ss(pt->get_e()));
    } else if (t == SUM) {
        const ex_sum *pt = p.get<ex_sum>();
        return ex_sum::create(pt->get_ie(), ss(pt->get_e()));
    } else if (t == PROD) {
        const ex_prod *pt = p.get<ex_prod>();
        return ex_prod::create(pt->get_ie(), ss(pt->get_e()));
    } else INTERNAL_ERROR;
}






ptr_base
symbolic::internal::drop_Es(const ptr_base &p)
{
    if (!has_Es(p)) return p;
    unsigned t = p->type();
    if (t == EX) {
        return p.get<ex_e>()->get_arg();
    } else if (t == ADD) {
        const num_ex_pair_vec &in = p.get<ex_add>()->get_ops();
        num_ex_pair_vec out;
        unsigned i, n = in.size();
        out.reserve(n);
        for (i = 0; i < n; ++i)
            out.push_back(std::pair<Number, ptr_base>(in[i].first,
                                                      drop_Es(in[i].second)));
        return mk_add(out);
    } else if (t == MUL) {
        const num_ex_pair_vec &in = p.get<ex_mul>()->get_ops();
        num_ex_pair_vec out;
        unsigned i, n = in.size();
        out.reserve(n);
        for (i = 0; i < n; ++i)
            out.push_back(std::pair<Number, ptr_base>(in[i].first,
                                                      drop_Es(in[i].second)));
        return mk_mul(out);
    } else if (t == POW) {
        const ex_pow *pt = p.get<ex_pow>();
        return mk_pow(drop_Es(pt->get_base()), drop_Es(pt->get_exp()));
    } else if (t == FUN) {
        const ex_func *pt = p.get<ex_func>();
        return mk_func(pt->get_code(), drop_Es(pt->get_arg()));
    } else if (t == IDX) {
        const ex_idx *pt = p.get<ex_idx>();
        return ex_idx::create(pt->get_ie(), drop_Es(pt->get_e()));
    } else if (t == SUM) {
        const ex_sum *pt = p.get<ex_sum>();
        return ex_sum::create(pt->get_ie(), drop_Es(pt->get_e()));
    } else if (t == PROD) {
        const ex_prod *pt = p.get<ex_prod>();
        return ex_prod::create(pt->get_ie(), drop_Es(pt->get_e()));
    } else INTERNAL_ERROR;
}



ptr_base
symbolic::internal::append_name(const ptr_base &p, const std::string &s)
{
    unsigned t = p->type();
    if (t == IDX) {
        const ex_idx *pt = p.get<ex_idx>();
        return ex_idx::create(pt->get_ie(), append_name(pt->get_e(), s));
    } else if (t == SYMB) {
        return ex_symb::create(stringhash::get_instance().append(p.get<ex_symb>()->m_hash, s));
    } else if (t == SYMBIDX) {
        const ex_symbidx *pt = p.get<ex_symbidx>();
        return ex_symbidx::create(stringhash::get_instance().append(pt->m_hash, s),
                                  pt->m_noid, pt->m_idx1, pt->m_idx2, pt->m_idx3, pt->m_idx4);
    } else if (t == VART) {
        return ex_vart::create(stringhash::get_instance().append(p.get<ex_vart>()->m_hash, s), 0);
    } else if (t == VARTIDX) {
        const ex_vartidx *pt = p.get<ex_vartidx>();
        return ex_vartidx::create(stringhash::get_instance().append(pt->m_hash, s), 0,
                                  pt->m_noid, pt->m_idx1, pt->m_idx2, pt->m_idx3, pt->m_idx4);
    } else USER_ERROR("append_name expects parameter or variable")
}



ptr_base
symbolic::internal::add_idx(const ptr_base &p, const idx_ex &ie)
{
    if (!ie) return p;
    unsigned t = p->type();
    if (t == IDX) {
        const ex_idx *pt = p.get<ex_idx>();
        return ex_idx::create(pt->get_ie(), add_idx(pt->get_e(), ie));
//         return add_idx(pt->get_e(), ie);
    } else if (t == SYMB) {
        const ex_symb *pt = p.get<ex_symb>();
        return ex_symbidx::create(pt->m_hash, 1, ie.m_id, 0, 0, 0);
    } else if (t == SYMBIDX) {
        if (p->hasidx(ie.m_id)) return p;
        const ex_symbidx *pt = p.get<ex_symbidx>();
        unsigned no = pt->m_noid;
        int i1 = pt->m_idx1, i2 = pt->m_idx2, i3 = pt->m_idx3, i4 = pt->m_idx4;
        switch (no) {
            case 1: i2 = ie.m_id; break;
            case 2: i3 = ie.m_id; break;
            case 3: i4 = ie.m_id; break;
            default: INTERNAL_ERROR
        }
        return ex_symbidx::create(pt->m_hash, no + 1, i1, i2, i3, i4);
    } else if (t == VART) {
        const ex_vart *pt = p.get<ex_vart>();
        return ex_vartidx::create(pt->m_hash, pt->m_lag, 1, ie.m_id, 0, 0, 0);
    } else if (t == VARTIDX) {
        if (p->hasidx(ie.m_id)) return p;
        const ex_vartidx *pt = p.get<ex_vartidx>();
        unsigned no = pt->m_noid;
        int i1 = pt->m_idx1, i2 = pt->m_idx2, i3 = pt->m_idx3, i4 = pt->m_idx4;
        switch (no) {
            case 1: i2 = ie.m_id; break;
            case 2: i3 = ie.m_id; break;
            case 3: i4 = ie.m_id; break;
            default: INTERNAL_ERROR
        }
        return ex_vartidx::create(pt->m_hash, pt->m_lag, no + 1, i1, i2, i3, i4);
    } else USER_ERROR("add_idx expects parameter or variable")
}



ptr_base
symbolic::internal::add_idx(const ptr_base &p, const ptr_base &ie)
{
    if (ie->type() != IDX) return p;
    const ex_idx *pt = ie.get<ex_idx>();
    return add_idx(add_idx(p, pt->get_e()), pt->get_ie());
}


ptr_base
symbolic::internal::apply_idx(const ptr_base &p, const ptr_base &of)
{
    if (of->type() != IDX) return p;
    const ex_idx *pt = of.get<ex_idx>();
    return ex_idx::create(pt->get_ie(), apply_idx(p, pt->get_e()));
}


ptr_base
symbolic::internal::subst(const ptr_base &p,
                          const ptr_base &what, const ptr_base &with,
                          bool all_leads_lags)
{
    if (with->type() == IDX) return subst(p, what, with.get<ex_idx>()->get_ptr(), all_leads_lags);
    unsigned t = p->type(), tw = what->type();
    if (tw == IDX) {
        return subst(p, what.get<ex_idx>()->get_e(), with, all_leads_lags);
    }
    if ((tw != SYMB) && (tw != SYMBIDX) && (tw != VART) && (tw != VARTIDX) && (tw != EX)) {
        USER_ERROR("trying to substitute for something that is not a constant, nor variable, nor expected val expression")
    }
    if (!(with->flag() & HAST) && (with->flag() & HAST)) {
        USER_ERROR("trying to substitute for constant with an expression containing variable")
    }

    if (all_leads_lags) {
        if (!p->has(what, ANY_T, false)) return p;
    } else {
        if (!p->has(what, EXACT_T, false)) return p;
    }
    if (t == EX) {
        const ex_e *pt = p.get<ex_e>();
        const ex_e *ptw = what.get<ex_e>();
        if (tw == EX) {
            int l1 = pt->get_lag(), l2 = ptw->get_lag(), ld = l2 - l1;
            if (ld && (!all_leads_lags))
                return mk_E(subst(pt->get_arg(), what, with, all_leads_lags), pt->get_lag());
            if (!compare(lag(pt->get_arg(), ld), ptw->get_arg())) {
                return lag(with, -ld);
            }
        }
        return mk_E(subst(pt->get_arg(), what, with, all_leads_lags), pt->get_lag());
    } else if (t == tw) { // primitive type that contains what, so is equal
        if (t == SYMB) return with;
        else if (t == SYMBIDX) {
            ptr_base with2 = with;
            const ex_symbidx *pt = p.get<ex_symbidx>();
            const ex_symbidx *pwhat = what.get<ex_symbidx>();
            if (pwhat->m_idx4 > 0) with2 = substidx(with, (unsigned) pwhat->m_idx4, pt->m_idx4);
            if (pwhat->m_idx3 > 0) with2 = substidx(with, (unsigned) pwhat->m_idx3, pt->m_idx3);
            if (pwhat->m_idx2 > 0) with2 = substidx(with, (unsigned) pwhat->m_idx2, pt->m_idx2);
            if (pwhat->m_idx1 > 0) with2 = substidx(with, (unsigned) pwhat->m_idx1, pt->m_idx1);
            return with2;
        } else if (t == VART) {
            if (all_leads_lags) {
                int l1 = p->get_lag_max(), l2 = what->get_lag_max();
                if (l1 == INT_MIN) return ss(with);
                return lag(with, l1 - l2);
            } else return with;
        } else if (t == VARTIDX) {
            ptr_base with2 = with;
            const ex_vartidx *pt = p.get<ex_vartidx>();
            const ex_vartidx *pwhat = what.get<ex_vartidx>();
            if (pwhat->m_idx4 > 0) with2 = substidx(with, (unsigned) pwhat->m_idx4, pt->m_idx4);
            if (pwhat->m_idx3 > 0) with2 = substidx(with, (unsigned) pwhat->m_idx3, pt->m_idx3);
            if (pwhat->m_idx2 > 0) with2 = substidx(with, (unsigned) pwhat->m_idx2, pt->m_idx2);
            if (pwhat->m_idx1 > 0) with2 = substidx(with, (unsigned) pwhat->m_idx1, pt->m_idx1);
            if (all_leads_lags) {
                int l1 = p->get_lag_max(), l2 = what->get_lag_max();
                if (l1 == INT_MIN) return ss(with2);
                return lag(with2, l1 - l2);
            } else return with2;
        } else INTERNAL_ERROR
    } else if (t == SUM) {
        const ex_sum *pt = p.get<ex_sum>();
        return ex_sum::create(pt->get_ie(), subst(pt->get_e(), what, with, all_leads_lags));
    } else if (t == PROD) {
        const ex_prod *pt = p.get<ex_prod>();
        return ex_prod::create(pt->get_ie(), subst(pt->get_e(), what, with, all_leads_lags));
    } else if (t == ADD) {
        const num_ex_pair_vec &in = p.get<ex_add>()->get_ops();
        num_ex_pair_vec out;
        unsigned i, n = in.size();
        out.reserve(n);
        for (i = 0; i < n; ++i)
            out.push_back(std::pair<Number, ptr_base>(
                in[i].first, subst(in[i].second, what, with, all_leads_lags)));
        return mk_add(out);
    } else if (t == MUL) {
        const num_ex_pair_vec &in = p.get<ex_mul>()->get_ops();
        num_ex_pair_vec out;
        unsigned i, n = in.size();
        out.reserve(n);
        for (i = 0; i < n; ++i)
            out.push_back(std::pair<Number, ptr_base>(
                in[i].first, subst(in[i].second, what, with, all_leads_lags)));
        return mk_mul(out);
    } else if (t == POW) {
        const ex_pow *pt = p.get<ex_pow>();
        return mk_pow(subst(pt->get_base(), what, with, all_leads_lags),
                      subst(pt->get_exp(), what, with, all_leads_lags));
    } else if (t == FUN) {
        const ex_func *pt = p.get<ex_func>();
        return mk_func(pt->get_code(), subst(pt->get_arg(), what, with, all_leads_lags));
    } else if (t == IDX) {
        const ex_idx *pt = p.get<ex_idx>();
        return ex_idx::create(pt->get_ie(), subst(pt->get_e(), what, with, all_leads_lags));
    } else INTERNAL_ERROR
}




ptr_base
symbolic::internal::substidx(const ptr_base &p, unsigned what, int with)
{
    if (!p->hasidx(what)) return p;
    unsigned t = p->type();
    if (t == DELTA) {
        return p.get<ex_delta>()->substidx(what, with);
    } else if (t == SYMBIDX) {
        return p.get<ex_symbidx>()->substidx(what, with);
    } else if (t == VARTIDX) {
        return p.get<ex_vartidx>()->substidx(what, with);
    } else if (t == SUM) {
        const ex_sum *pt = p.get<ex_sum>();
        if (pt->get_ie().get_id() == what) return substidx(pt->reindex(), what, with);
        else return ex_sum::create(pt->get_ie().substidx(what, with), substidx(pt->get_e(), what, with));
    } else if (t == PROD) {
        const ex_prod *pt = p.get<ex_prod>();
        if (pt->get_ie().get_id() == what) return substidx(pt->reindex(), what, with);
        else return ex_prod::create(pt->get_ie().substidx(what, with), substidx(pt->get_e(), what, with));
    } else if (t == IDX) {
        const ex_idx *pt = p.get<ex_idx>();
        return ex_idx::create(pt->get_ie().substidx(what, with), substidx(pt->get_e(), what, with));
    } else if (t == EX) {
        const ex_e *pt = p.get<ex_e>();
        return mk_E(substidx(pt->get_arg(), what, with), pt->get_lag());
    } else if (t == ADD) {
        const num_ex_pair_vec &in = p.get<ex_add>()->get_ops();
        num_ex_pair_vec out;
        unsigned i, n = in.size();
        out.reserve(n);
        for (i = 0; i < n; ++i)
            out.push_back(std::pair<Number, ptr_base>(
                in[i].first, substidx(in[i].second, what, with)));
        return mk_add(out);
    } else if (t == MUL) {
        const num_ex_pair_vec &in = p.get<ex_mul>()->get_ops();
        num_ex_pair_vec out;
        unsigned i, n = in.size();
        out.reserve(n);
        for (i = 0; i < n; ++i)
            out.push_back(std::pair<Number, ptr_base>(
                in[i].first, substidx(in[i].second, what, with)));
        return mk_mul(out);
    } else if (t == POW) {
        const ex_pow *pt = p.get<ex_pow>();
        return mk_pow(substidx(pt->get_base(), what, with),
                      substidx(pt->get_exp(), what, with));
    } else if (t == FUN) {
        const ex_func *pt = p.get<ex_func>();
        return mk_func(pt->get_code(), substidx(pt->get_arg(), what, with));
    } else INTERNAL_ERROR;
}



ptr_base
symbolic::internal::expandsp(const ptr_base &p)
{
    unsigned t = p->type();
    if ((t == NUM) || (t == DELTA) || (t == SYMB) || (t == SYMBIDX)
        || (t == VART) || (t == VARTIDX)) {
        return p;
    } else if (t == SUM) {
        return p.get<ex_sum>()->expand();
    } else if (t == PROD) {
        return p.get<ex_prod>()->expand();
    } else if (t == EX) {
        const ex_e *pt = p.get<ex_e>();
        return mk_E(expandsp(pt->get_arg()), pt->get_lag());
    } else if (t == ADD) {
        const num_ex_pair_vec &in = p.get<ex_add>()->get_ops();
        num_ex_pair_vec out;
        unsigned i, n = in.size();
        out.reserve(n);
        for (i = 0; i < n; ++i)
            out.push_back(std::pair<Number, ptr_base>(in[i].first,
                                                      expandsp(in[i].second)));
        return mk_add(out);
    } else if (t == MUL) {
        const num_ex_pair_vec &in = p.get<ex_mul>()->get_ops();
        num_ex_pair_vec out;
        unsigned i, n = in.size();
        out.reserve(n);
        for (i = 0; i < n; ++i)
            out.push_back(std::pair<Number, ptr_base>(in[i].first,
                                                      expandsp(in[i].second)));
        return mk_mul(out);
    } else if (t == POW) {
        const ex_pow *pt = p.get<ex_pow>();
        return mk_pow(expandsp(pt->get_base()), expandsp(pt->get_exp()));
    } else if (t == FUN) {
        const ex_func *pt = p.get<ex_func>();
        return mk_func(pt->get_code(), expandsp(pt->get_arg()));
    } else INTERNAL_ERROR;
}






