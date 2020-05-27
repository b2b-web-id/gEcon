/*****************************************************************************
 * This file is a part of gEcon.                                             *
 *                                                                           *
 * (c) Chancellery of the Prime Minister of the Republic of Poland 2012-2015 *
 * (c) Grzegorz Klima, Karol Podemski, Kaja Retkiewicz-Wijtiwiak 2015-2018   *
 * License terms can be found in the file 'LICENCE'                          *
 *                                                                           *
 * Author: Grzegorz Klima                                                    *
 *****************************************************************************/

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
#include <cmp.h>
#include <error.h>
#include <ops.h>
#include <iostream>
#include <cmath>
#include <climits>


using namespace symbolic;
using namespace symbolic::internal;


ptr_base
symbolic::internal::reduce(const ptr_base &p)
{
    if (p->flag() == EMPTY) {
        if (p->type() == ADD) return ex_num::zero();
        if (p->type() == MUL) return ex_num::one();
        INTERNAL_ERROR;
    }
    if ((p->type() == MUL) && (p->flag() & SINGLE) && (p->flag() &  SCAL1)) {
        return p.get<ex_mul>()->get_ops()[0].second;
    }
    if ((p->type() == ADD) && (p->flag() & SINGLE) && (p->flag() &  SCAL1)) {
        return p.get<ex_add>()->get_ops()[0].second;
    }
    if ((p->type() == MUL) && (p->flag() & SPOW1)) {
        const ex_mul *pp = p.get<ex_mul>();
        return ex_add::create(pp->get_ops()[0].second->val(), pp->get_ops()[1].second);
    }

    return p;
}


ptr_base
symbolic::internal::mk_neg(const ptr_base &x)
{
    return mk_mul(ex_num::create(-1.), x);
}


ptr_base
symbolic::internal::mk_add(const ptr_base &lhs, const ptr_base &rhs)
{
    if (lhs->is0()) return rhs;
    if (rhs->is0()) return lhs;

    unsigned t1, t2;
    t1 = lhs->type();
    t2 = rhs->type();

    if ((t1 == IDX) && (t2 == IDX)) {
        const ex_idx *p1 = lhs.get<ex_idx>();
        const ex_idx *p2 = rhs.get<ex_idx>();
        if (p1->get_ie().compare(p2->get_ie())) {
            USER_ERROR("mismatched indexing in binary operation")
        }
        return ex_idx::create(p1->get_ie(), mk_add(p1->get_e(), p2->get_e()));
    } else if (t1 == IDX) {
        const ex_idx *pt = lhs.get<ex_idx>();
        return ex_idx::create(pt->get_ie(), mk_add(pt->get_e(), rhs));
    } else if (t2 == IDX) {
        const ex_idx *pt = rhs.get<ex_idx>();
        return ex_idx::create(pt->get_ie(), mk_add(lhs, pt->get_e()));
    }

    if ((t1 == NUM) && (t2 == NUM)) {
        return ex_num::create(lhs->val() + rhs->val());
    } else {
        return reduce(ex_add::create(lhs, rhs));
    }
}



namespace {

ptr_base
mk_rec_add(const num_ex_pair_vec &ops, unsigned b, unsigned e)
{
    if (b == e) return mk_mul(ex_num::create(ops[b].first), ops[b].second);
    if ((b + 1) == e) return mk_add(mk_rec_add(ops, b, b), mk_rec_add(ops, e, e));
    unsigned d = (e - b) / 2;
    return mk_add(mk_rec_add(ops, b, b + d), mk_rec_add(ops, b + d + 1, e));
}

} /* namespace */


ptr_base
symbolic::internal::mk_add(const num_ex_pair_vec &ops)
{
    if (!ops.size()) return ex_num::zero();
    return mk_rec_add(ops, 0, ops.size() - 1);
}




ptr_base
symbolic::internal::mk_sub(const ptr_base &lhs, const ptr_base &rhs)
{
    return mk_add(lhs, mk_neg(rhs));
}




ptr_base
symbolic::internal::mk_mul(const ptr_base &lhs, const ptr_base &rhs)
{
    unsigned t1, t2, f1, f2;
    t1 = lhs->type();
    t2 = rhs->type();
    f1 = lhs->flag();
    f2 = rhs->flag();

    if ((t1 == IDX) && (t2 == IDX)) {
        const ex_idx *p1 = lhs.get<ex_idx>();
        const ex_idx *p2 = rhs.get<ex_idx>();
        if (p1->get_ie().compare(p2->get_ie())) {
            USER_ERROR("mismatched indexing in binary operation")
        }
        return ex_idx::create(p1->get_ie(), mk_mul(p1->get_e(), p2->get_e()));
    } else if (t1 == IDX) {
        const ex_idx *pt = lhs.get<ex_idx>();
        return ex_idx::create(pt->get_ie(), mk_mul(pt->get_e(), rhs));
    } else if (t2 == IDX) {
        const ex_idx *pt = rhs.get<ex_idx>();
        return ex_idx::create(pt->get_ie(), mk_mul(lhs, pt->get_e()));
    }

    if ((t1 == NUM) && (t2 == NUM)) {
        return ex_num::create(lhs->val() * rhs->val());
    }

    if (lhs->is0()) return lhs;
    if (rhs->is0()) return rhs;
    if (lhs->is1()) return rhs;
    if (rhs->is1()) return lhs;

    if ((t1 == NUM) && (t2 == ADD)) {
        return ex_add::create(lhs->val(), rhs);
    } else if ((t1 == ADD) && (t2 == NUM)) {
        return ex_add::create(rhs->val(), lhs);
    } else if ((t1 == NUM) && ((t2 == SYMB) || (t2 == VART))) {
        return ex_add::create(lhs->val(), rhs);
    } else if (((t1 == SYMB) || (t1 == VART)) && (t2 == NUM)) {
        return ex_add::create(rhs->val(), lhs);
    } else if ((t1 == ADD) && (f1 & SINGLE) && (t2 == ADD) && (f2 & SINGLE)) {
        const ex_add *p1 = lhs.get<ex_add>();
        const ex_add *p2 = rhs.get<ex_add>();
        return reduce(ex_mul::create(
            ex_mul::create(
                ex_num::create(p1->get_ops()[0].first),
                p1->get_ops()[0].second),
            ex_mul::create(
                ex_num::create(p2->get_ops()[0].first),
                p2->get_ops()[0].second)));
    } else if ((t1 == ADD) && (f1 & SINGLE)) {
        const ex_add *p = lhs.get<ex_add>();
        return reduce(ex_mul::create(ex_num::create(p->get_ops()[0].first),
                    ex_mul::create(p->get_ops()[0].second, rhs)));
    } else if ((t2 == ADD) && (f2 & SINGLE)) {
        const ex_add *p = rhs.get<ex_add>();
        return reduce(ex_mul::create(ex_num::create(p->get_ops()[0].first),
                    ex_mul::create(p->get_ops()[0].second, lhs)));
//     } else if (t1 == ADD) {
//         const num_ex_pair_vec &in = lhs.get<ex_add>()->get_ops();
//         num_ex_pair_vec out;
//         unsigned i, n = in.size();
//         out.reserve(n);
//         for (i = 0; i < n; ++i) {
//             out.push_back(num_ex_pair(in[i].first, mk_mul(in[i].second, rhs)));
//         }
//         return mk_add(out);
//     } else if (t2 == ADD) {
//         const num_ex_pair_vec &in = rhs.get<ex_add>()->get_ops();
//         num_ex_pair_vec out;
//         unsigned i, n = in.size();
//         out.reserve(n);
//         for (i = 0; i < n; ++i) {
//             out.push_back(num_ex_pair(in[i].first, mk_mul(lhs, in[i].second)));
//         }
//         return mk_add(out);
    } else {
        return reduce(ex_mul::create(lhs, rhs));
    }
}



namespace {

ptr_base
mk_rec_mul(const num_ex_pair_vec &ops, unsigned b, unsigned e)
{
    if (b == e) return mk_pow(ops[b].second, ex_num::create(ops[b].first));
    if ((b + 1) == e) return mk_mul(mk_rec_mul(ops, b, b), mk_rec_mul(ops, e, e));
    unsigned d = (e - b) / 2;
    return mk_mul(mk_rec_mul(ops, b, b + d), mk_rec_mul(ops, b + d + 1, e));
}


} /* namespace */


ptr_base
symbolic::internal::mk_mul(const num_ex_pair_vec &ops)
{
    if (!ops.size()) return ex_num::one();
    return mk_rec_mul(ops, 0, ops.size() - 1);
}




ptr_base
symbolic::internal::mk_div(const ptr_base &lhs, const ptr_base &rhs)
{
    return mk_mul(lhs, mk_pow(rhs, ex_num::create(-1)));
}


ptr_base
symbolic::internal::mk_pow(const ptr_base &lhs, const ptr_base &rhs)
{
    unsigned t1, t2, f1;
    t1 = lhs->type();
    t2 = rhs->type();
    f1 = lhs->flag();

    if ((t1 == IDX) && (t2 == IDX)) {
        const ex_idx *p1 = lhs.get<ex_idx>();
        const ex_idx *p2 = rhs.get<ex_idx>();
        if (p1->get_ie().compare(p2->get_ie())) {
            USER_ERROR("mismatched indexing in binary operation")
        }
        return ex_idx::create(p1->get_ie(), mk_pow(p1->get_e(), p2->get_e()));
    } else if (t1 == IDX) {
        const ex_idx *pt = lhs.get<ex_idx>();
        return ex_idx::create(pt->get_ie(), mk_pow(pt->get_e(), rhs));
    } else if (t2 == IDX) {
        const ex_idx *pt = rhs.get<ex_idx>();
        return ex_idx::create(pt->get_ie(), mk_pow(lhs, pt->get_e()));
    }

    if ((t1 == NUM) && (t2 == NUM)) {
        return ex_num::create(pow(lhs->val(), rhs->val()));
    }

    if (lhs->is0()) return lhs;
    if (rhs->is0()) return ex_num::one();
    if (rhs->is1()) return lhs;

    if ((t1 == MUL) && (t2 == NUM)) {
        return ex_mul::create(rhs->val(), lhs);
    } else if (((t1 == SYMB) || (t1 == VART)) && (t2 == NUM)) {
        return ex_mul::create(rhs->val(), lhs);
    } else if ((t1 == ADD) && (f1 & SINGLE) && (t2 == NUM)) {
        const ex_add *p = lhs.get<ex_add>();
        if (p->get_ops()[0].first == 1.)
            return ex_pow::create(p->get_ops()[0].second, rhs);
        return ex_mul::create(
            ex_num::create(pow(p->get_ops()[0].first, rhs->val())),
            ex_mul::create(rhs->val(), p->get_ops()[0].second)
            );
    } else if (t1 == POW) {
        const ex_pow *p = lhs.get<ex_pow>();
        return mk_pow(p->get_base(), mk_mul(p->get_exp(), rhs));
    } else if (t2 == NUM) {
        return ex_mul::create(rhs->val(), lhs);
    } else {
        return ex_pow::create(lhs, rhs);
    }
}




ptr_base
symbolic::internal::mk_E(const ptr_base &p, int l)
{
    if (p->type() == IDX) {
        const ex_idx *pt = p.get<ex_idx>();
        return (ex_idx::create(pt->get_ie(), mk_E(pt->get_e(), l)));
    }
    if (!(p->flag() & HAST)) return p;
    if (l >= p->get_lag_max(true)) return p;
    unsigned t = p->type();
    if (t == VART) {
        return ex_e::create(p, l);
    } else if (t == EX) {
        const ex_e *pp = p.get<ex_e>();
        return ex_e::create(pp->get_arg(), std::min(l, pp->get_lag()));
    } else if (t == SUM) {
        const ex_sum *pt = p.get<ex_sum>();
        return (ex_sum::create(pt->get_ie(), mk_E(pt->get_e(), l)));
    } else if (t == ADD) {
        const num_ex_pair_vec &in = p.get<ex_add>()->get_ops();
        num_ex_pair_vec out;
        unsigned i, n = in.size();
        out.reserve(n);
        for (i = 0; i < n; ++i) {
            out.push_back(num_ex_pair(in[i].first, mk_E(in[i].second, l)));
        }
        return mk_add(out);
    } else if (t == MUL) {
        const num_ex_pair_vec &in = p.get<ex_mul>()->get_ops();
        num_ex_pair_vec out1, out2;
        unsigned i, n = in.size();
        out1.reserve(n);
        out2.reserve(n);
        for (i = 0; i < n; ++i) {
            if (in[i].second->get_lag_max(true) > l) {
                out2.push_back(in[i]);
            } else {
                out1.push_back(in[i]);
            }
        }
        out1.push_back(num_ex_pair(1., ex_e::create(mk_mul(out2), l)));
        return mk_mul(out1);
    } else return ex_e::create(p, l);
}




ptr_base
symbolic::internal::mk_func(func_code c, const ptr_base &arg)
{
    if (arg->type() == IDX) {
        const ex_idx *pt = arg.get<ex_idx>();
        return (ex_idx::create(pt->get_ie(), mk_func(c, pt->get_e())));
    }
    if (arg->type() == NUM) return ex_num::create(ex_func::eval(c, arg->val()));
    return ex_func::create(c, arg);
}


ptr_base
symbolic::internal::mk_sum(const ptr_base &arg)
{
    if (arg->type() != IDX) return arg;
    const ex_idx *pt = arg.get<ex_idx>();
    return ex_sum::create(pt->get_ie(), mk_sum(pt->get_e()));
}


ptr_base
symbolic::internal::mk_prod(const ptr_base &arg)
{
    if (arg->type() != IDX) return arg;
    const ex_idx *pt = arg.get<ex_idx>();
    return ex_prod::create(pt->get_ie(), mk_prod(pt->get_e()));
}




ptr_base
symbolic::internal::reduceidx(const ptr_base &p, const idx_ex &ie)
{
    unsigned t = p->type();
    if (!ie) {
        if (t == IDX) {
            const ex_idx *pt = p.get<ex_idx>();
            return ex_idx::create(pt->get_ie(), reduceidx(pt->get_e(), pt->get_ie()));
        } else return p;
    }

    if (t == IDX) {
        const ex_idx *pt = p.get<ex_idx>();
        ptr_base res = reduceidx(pt->get_e(), pt->get_ie());
        res = reduceidx(res, ie);
        res = ex_idx::create(pt->get_ie(), res);
        res = ex_idx::create(ie, res);
        return res;
    } else if ((t == NUM) || (t == DELTA) || (t == SYMB) || (t == SYMBIDX)
        || (t == VART) || (t == VARTIDX)) {
        return p;
    } else if (t == SUM) {
        return p.get<ex_sum>()->reduceidx(ie);
    } else if (t == PROD) {
        const ex_prod *pt = p.get<ex_prod>();
        return ex_prod::create(pt->get_ie(), reduceidx(pt->get_e(), ie));
    } else if (t == EX) {
        const ex_e *pt = p.get<ex_e>();
        return mk_E(reduceidx(pt->get_arg(), ie), pt->get_lag());
    } else if (t == ADD) {
        const num_ex_pair_vec &in = p.get<ex_add>()->get_ops();
        num_ex_pair_vec out;
        unsigned i, n = in.size();
        out.reserve(n);
        for (i = 0; i < n; ++i)
            out.push_back(std::pair<Number, ptr_base>(in[i].first,
                                                      reduceidx(in[i].second, ie)));
        return mk_add(out);
    } else if (t == MUL) {
        const num_ex_pair_vec &in = p.get<ex_mul>()->get_ops();
        num_ex_pair_vec out;
        unsigned i, n = in.size();
        out.reserve(n);
        for (i = 0; i < n; ++i)
            out.push_back(std::pair<Number, ptr_base>(in[i].first,
                                                      reduceidx(in[i].second, ie)));
        return mk_mul(out);
    } else if (t == POW) {
        const ex_pow *pt = p.get<ex_pow>();
        return mk_pow(reduceidx(pt->get_base(), ie), reduceidx(pt->get_exp(), ie));
    } else if (t == FUN) {
        const ex_func *pt = p.get<ex_func>();
        return ex_func::create(pt->get_code(), reduceidx(pt->get_arg(), ie));
    } else INTERNAL_ERROR;
}




ptr_base
symbolic::internal::mk_diff(const ptr_base &e, const ptr_base &v)
{
    unsigned t = v->type();
    if (t == IDX) {
        const ex_idx *pt = v.get<ex_idx>();
        return reduceidx(ex_idx::create(pt->get_ie(), mk_diff(e, pt->get_e())), idx_ex());
    } else if ((t != SYMB) && (t != SYMBIDX) && (t != VART) && (t != VARTIDX)) {
        USER_ERROR("trying to differentiate wrt something that is not a parameter nor a variable")
    }

    if (!e->has(v, EXACT_T, false)) return ex_num::zero();
    return e->diff(v);
}



