/***********************************************************
 * (c) Kancelaria Prezesa Rady Ministrów 2012-2015         *
 * Treść licencji w pliku 'LICENCE'                        *
 *                                                         *
 * (c) Chancellery of the Prime Minister 2012-2015         *
 * License terms can be found in the file 'LICENCE'        *
 *                                                         *
 * Author: Grzegorz Klima                                  *
 ***********************************************************/

/** \file cmp.h
 * \brief Comparing expressions (used for sorting and reduction).
 */

#ifndef SYMBOLIC_CMP_H

#define SYMBOLIC_CMP_H

#include <ptr_base.h>
#include <ex_num.h>
#include <ex_delta.h>
#include <ex_symb.h>
#include <ex_symbidx.h>
#include <ex_vart.h>
#include <ex_vartidx.h>
#include <ex_pow.h>
#include <ex_func.h>
#include <ex_e.h>
#include <ex_add.h>
#include <ex_mul.h>
#include <ex_sum.h>
#include <ex_prod.h>
#include <ex_idx.h>
#include <error.h>


namespace symbolic {
namespace internal {


template <typename T>
inline
int
compareT(T a, T b)
{
    if (a < b) return -1;
    if (a > b) return 1;
    return 0;
}


inline
int
compare(const ptr_base &a, const ptr_base &b)
{
    int c;
    unsigned t;
    if ((c = compareT(t = a->type(), b->type()))) return c;
#ifdef EXPAND_CASE
#undef EXPAND_CASE
#endif
#define EXPAND_CASE(TYPE, CLASS) \
        case TYPE: \
            return a.get<CLASS>()-> compare(*b.get<CLASS>());
    switch (t) {
        EXPAND_CASE(NUM, ex_num)
        EXPAND_CASE(DELTA, ex_delta)
        EXPAND_CASE(SYMB, ex_symb)
        EXPAND_CASE(SYMBIDX, ex_symbidx)
        EXPAND_CASE(VART, ex_vart)
        EXPAND_CASE(VARTIDX, ex_vartidx)
        EXPAND_CASE(EX, ex_e)
        EXPAND_CASE(POW, ex_pow)
        EXPAND_CASE(FUN, ex_func)
        EXPAND_CASE(ADD, ex_add)
        EXPAND_CASE(MUL, ex_mul)
        EXPAND_CASE(SUM, ex_sum)
        EXPAND_CASE(PROD, ex_prod)
        EXPAND_CASE(IDX, ex_idx)
        default:
            INTERNAL_ERROR
    }
#undef EXPAND_CASE
    return 0;
}


} /* namespace internal */
} /* namespace symbolic */

#endif /* SYMBOLIC_CMP_H */
