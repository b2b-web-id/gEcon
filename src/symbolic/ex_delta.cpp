/***********************************************************
 * (c) Kancelaria Prezesa Rady Ministrów 2012-2015         *
 * Treść licencji w pliku 'LICENCE'                        *
 *                                                         *
 * (c) Chancellery of the Prime Minister 2012-2015         *
 * License terms can be found in the file 'LICENCE'        *
 *                                                         *
 * Author: Grzegorz Klima                                  *
 ***********************************************************/

/** \file ex_delta.cpp
 * \brief Symbols.
 */

#include <ex_delta.h>
#include <ex_num.h>
#include <stringhash.h>
#include <utils.h>
#include <cmp.h>
#include <climits>


using namespace symbolic;
using namespace symbolic::internal;


ptr_base
ex_delta::create(bool fixed1, const std::string &i1,
                 bool fixed2, const std::string &i2)
{
    int idx1, idx2;
    stringhash &ref = stringhash::get_instance();
    idx1 = ref.get_hash(i1);
    if (fixed1) idx1 = -idx1;
    idx2 = ref.get_hash(i2);
    if (fixed2) idx2 = -idx2;
    return create(idx1, idx2);
}



ex_delta::ex_delta(int i1, int i2) : ex_base(DELTA | SINGLE | HASIDX)
{
    if (i1 < i2) {
        m_idx1 = i1;
        m_idx2 = i2;
    } else {
        m_idx1 = i2;
        m_idx2 = i1;
    }
}


ptr_base
ex_delta::create(int i1, int i2)
{
    if (i1 == i2) return ex_num::one();
    if ((i1 < 0) && (i2 < 0) && (i1 != i2)) return ex_num::zero();
    return ptr_base(new ex_delta(i1, i2));
}


void
ex_delta::destroy(ex_base *ptr)
{
    delete ptr;
}


int
ex_delta::compare(const ex_delta &b) const
{
    int t;
    if ((t = compareT(m_idx1, b.m_idx1))) return t;
    if ((t = compareT(m_idx2, b.m_idx2))) return t;
    return 0;
}



std::string
ex_delta::str(int pflag) const
{
    if (pflag & INDEXING_ONLY) return std::string();
    stringhash &ref = stringhash::get_instance();
    std::string res = "KRONECKER_DELTA";
    if (pflag & DROP_IDX) {
        return res;
    } else if (pflag & CONVERT_IDX) {
#ifdef EXPAND_PRINT_IDX
#undef EXPAND_PRINT_IDX
#endif
#define EXPAND_PRINT_IDX(id) \
if (id < 0) { \
    res += ref.get_str(-id); \
} else { \
    res += ref.get_str(id); \
}
            res += "__";
            EXPAND_PRINT_IDX(m_idx1)
            res += "__";
            EXPAND_PRINT_IDX(m_idx2)
            return res;
    } else {
#ifdef EXPAND_PRINT_IDX
#undef EXPAND_PRINT_IDX
#endif
#define EXPAND_PRINT_IDX(id) \
if (id < 0) { \
    res += '\''; \
    res += ref.get_str(-id); \
    res += '\''; \
} else { \
    res += ref.get_str(id); \
}
            res += '<';
            EXPAND_PRINT_IDX(m_idx1)
            res += ',';
            EXPAND_PRINT_IDX(m_idx2)
            res += '>';
            return res;
    }
}


std::string
ex_delta::strmap(const map_str_str &mss) const
{
    return str(internal::DEFAULT);
}


std::string
ex_delta::tex(int pflag) const
{
    if (pflag & INDEXING_ONLY) return std::string();
#ifdef EXPAND_PRINT_IDX
#undef EXPAND_PRINT_IDX
#endif
#define EXPAND_PRINT_IDX(id) \
if (id < 0) { \
    res += "\\mathrm{"; \
    res += str2tex(ref.get_str(-id), false); \
    res += "}"; \
} else { \
    res += str2tex(ref.get_str(id)); \
}
    stringhash &ref = stringhash::get_instance();
    std::string res = "\\delta^{\\langle ";
    EXPAND_PRINT_IDX(m_idx1)
    res += ',';
    EXPAND_PRINT_IDX(m_idx2)
    res += "\\rangle}";
    return res;
}


int
ex_delta::get_lag_max(bool) const
{
    return INT_MIN;
}


int
ex_delta::get_lag_min(bool) const
{
    return INT_MAX;
}


ptr_base
ex_delta::diff(const ptr_base &p) const
{
    return ex_num::zero();
}


bool
ex_delta::has(const ptr_base&, search_flag, bool) const
{
    return false;
}


bool
ex_delta::hasidx(int idx) const
{
    return ((m_idx1 == idx) || (m_idx2 == idx)) ? true : false;
}


ptr_base
ex_delta::substidx(unsigned what, int with) const
{
    int i1, i2;
    i1 = (m_idx1 == (int) what) ? with : m_idx1;
    i2 = (m_idx2 == (int) what) ? with : m_idx2;
    return create(i1, i2);
}









