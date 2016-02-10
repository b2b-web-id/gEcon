/***********************************************************
 * (c) Kancelaria Prezesa Rady Ministrów 2012-2015         *
 * Treść licencji w pliku 'LICENCE'                        *
 *                                                         *
 * (c) Chancellery of the Prime Minister 2012-2015         *
 * License terms can be found in the file 'LICENCE'        *
 *                                                         *
 * Author: Grzegorz Klima                                  *
 ***********************************************************/

/** \file ex_symbidx.cpp
 * \brief Indexed symbols.
 */

#include <ex_symbidx.h>
#include <ex_num.h>
#include <ops.h>
#include <stringhash.h>
#include <utils.h>
#include <cmp.h>
#include <climits>


using namespace symbolic;
using namespace symbolic::internal;




ex_symbidx::ex_symbidx(const std::string &n,
                       bool fixed1, const std::string &i1)
                        : ex_base(SYMBIDX | SINGLE)
{
    stringhash &ref = stringhash::get_instance();
    m_hash = ref.get_hash(n);
    m_idx1 = ref.get_hash(i1);
    if (fixed1) m_idx1 = -m_idx1;
    m_idx2 = 0;
    m_idx3 = 0;
    m_idx4 = 0;
    m_noid = 1;
    if (!fixed1) m_type |= HASIDX;
}


ex_symbidx::ex_symbidx(const std::string &n,
                       bool fixed1, const std::string &i1,
                       bool fixed2, const std::string &i2)
                        : ex_base(SYMBIDX | SINGLE)
{
    stringhash &ref = stringhash::get_instance();
    m_hash = ref.get_hash(n);
    m_idx1 = ref.get_hash(i1);
    if (fixed1) m_idx1 = -m_idx1;
    m_idx2 = ref.get_hash(i2);
    if (fixed2) m_idx2 = -m_idx2;
    m_idx3 = 0;
    m_idx4 = 0;
    m_noid = 2;
    if (!fixed1 || !fixed2) m_type |= HASIDX;
}


ex_symbidx::ex_symbidx(const std::string &n,
                       bool fixed1, const std::string &i1,
                       bool fixed2, const std::string &i2,
                       bool fixed3, const std::string &i3)
                        : ex_base(SYMBIDX | SINGLE)
{
    stringhash &ref = stringhash::get_instance();
    m_hash = ref.get_hash(n);
    m_idx1 = ref.get_hash(i1);
    if (fixed1) m_idx1 = -m_idx1;
    m_idx2 = ref.get_hash(i2);
    if (fixed2) m_idx2 = -m_idx2;
    m_idx3 = ref.get_hash(i3);
    if (fixed3) m_idx3 = -m_idx3;
    m_idx4 = 0;
    m_noid = 3;
    if (!fixed1 || !fixed2 || !fixed3) m_type |= HASIDX;
}


ex_symbidx::ex_symbidx(const std::string &n,
                       bool fixed1, const std::string &i1,
                       bool fixed2, const std::string &i2,
                       bool fixed3, const std::string &i3,
                       bool fixed4, const std::string &i4)
                        : ex_base(SYMBIDX | SINGLE)
{
    stringhash &ref = stringhash::get_instance();
    m_hash = ref.get_hash(n);
    m_idx1 = ref.get_hash(i1);
    if (fixed1) m_idx1 = -m_idx1;
    m_idx2 = ref.get_hash(i2);
    if (fixed2) m_idx2 = -m_idx2;
    m_idx3 = ref.get_hash(i3);
    if (fixed3) m_idx3 = -m_idx3;
    m_idx4 = ref.get_hash(i4);
    if (fixed4) m_idx4 = -m_idx4;
    m_noid = 4;
    if (!fixed1 || !fixed2 || !fixed3 || !fixed4) m_type |= HASIDX;
}


ptr_base
ex_symbidx::create(const std::string &n,
                   bool fixed1, const std::string &i1)
{
    return ptr_base(new ex_symbidx(n, fixed1, i1));
}


ptr_base
ex_symbidx::create(const std::string &n,
                   bool fixed1, const std::string &i1,
                   bool fixed2, const std::string &i2)
{
    return ptr_base(new ex_symbidx(n, fixed1, i1, fixed2, i2));
}


ptr_base
ex_symbidx::create(const std::string &n,
                   bool fixed1, const std::string &i1,
                   bool fixed2, const std::string &i2,
                   bool fixed3, const std::string &i3)
{
    return ptr_base(new ex_symbidx(n, fixed1, i1, fixed2, i2, fixed3, i3));
}


ptr_base
ex_symbidx::create(const std::string &n,
                   bool fixed1, const std::string &i1,
                   bool fixed2, const std::string &i2,
                   bool fixed3, const std::string &i3,
                   bool fixed4, const std::string &i4)
{
    return ptr_base(new ex_symbidx(n, fixed1, i1, fixed2, i2, fixed3, i3,
                                   fixed4, i4));
}



ex_symbidx::ex_symbidx(unsigned n, unsigned no, int i1, int i2, int i3, int i4)
                        : ex_base(SYMBIDX | SINGLE)
{
    m_hash = n;
    m_noid = no;
    m_idx1 = i1;
    m_idx2 = i2;
    m_idx3 = i3;
    m_idx4 = i4;
    if ((i1 > 0) || (i2 > 0) || (i3 > 0) || (i4 > 0)) m_type |= HASIDX;
}


ptr_base
ex_symbidx::create(unsigned n, unsigned no, int i1, int i2, int i3, int i4)
{
    return ptr_base(new ex_symbidx(n, no, i1, i2, i3, i4));
}


void
ex_symbidx::destroy(ex_base *ptr)
{
    delete ptr;
}



int
ex_symbidx::compare(const ex_symbidx &b) const
{
    int t;
    if ((t = compareT(m_hash, b.m_hash))) return t;
    if ((t = compareT(m_noid, b.m_noid))) return t;
    if ((t = compareT(m_idx1, b.m_idx1))) return ((m_idx1 < 0) && (b.m_idx1 < 0)) ? -t : t;
    if ((t = compareT(m_idx2, b.m_idx2))) return ((m_idx2 < 0) && (b.m_idx2 < 0)) ? -t : t;
    if ((t = compareT(m_idx3, b.m_idx3))) return ((m_idx3 < 0) && (b.m_idx3 < 0)) ? -t : t;
    if ((t = compareT(m_idx4, b.m_idx4))) return ((m_idx4 < 0) && (b.m_idx4 < 0)) ? -t : t;
    return 0;
}


std::string
ex_symbidx::get_name() const
{
    return stringhash::get_instance().get_str(m_hash);
}



std::string
ex_symbidx::str(int pflag) const
{
    if (pflag & INDEXING_ONLY) return std::string();
    stringhash &ref = stringhash::get_instance();
    std::string res = get_name();
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
        if (m_noid > 1) {
            res += "__";
            EXPAND_PRINT_IDX(m_idx2)
        }
        if (m_noid > 2) {
            res += "__";
            EXPAND_PRINT_IDX(m_idx3)
        }
        if (m_noid > 3) {
            res += "__";
            EXPAND_PRINT_IDX(m_idx4)
        }
        return res;
    } else if (pflag & DROP_QUOTES) {
#ifdef EXPAND_PRINT_IDX
#undef EXPAND_PRINT_IDX
#endif
#define EXPAND_PRINT_IDX(id) \
if (id < 0) { \
    res += ref.get_str(-id); \
} else { \
    res += ref.get_str(id); \
}
        res += '<';
        EXPAND_PRINT_IDX(m_idx1)
        if (m_noid > 1) {
            res += ',';
            EXPAND_PRINT_IDX(m_idx2)
        }
        if (m_noid > 2) {
            res += ',';
            EXPAND_PRINT_IDX(m_idx3)
        }
        if (m_noid > 3) {
            res += ',';
            EXPAND_PRINT_IDX(m_idx4)
        }
        res += '>';
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
        if (m_noid > 1) {
            res += ',';
            EXPAND_PRINT_IDX(m_idx2)
        }
        if (m_noid > 2) {
            res += ',';
            EXPAND_PRINT_IDX(m_idx3)
        }
        if (m_noid > 3) {
            res += ',';
            EXPAND_PRINT_IDX(m_idx4)
        }
        res += '>';
        return res;
    }
}


std::string
ex_symbidx::strmap(const map_str_str &mss) const
{
    map_str_str::const_iterator it;
    std::string name = str(CONVERT_IDX);
    it = mss.find(name);
    if (it == mss.end()) {
        return name;
    }
    return it->second;
}


std::string
ex_symbidx::tex(int pflag) const
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
    std::string res = "{" + str2tex(get_name(), true);
    res += "}^{\\langle ";
    EXPAND_PRINT_IDX(m_idx1)
    if (m_noid > 1) {
        res += ',';
        EXPAND_PRINT_IDX(m_idx2)
    }
    if (m_noid > 2) {
        res += ',';
        EXPAND_PRINT_IDX(m_idx3)
    }
    if (m_noid > 3) {
        res += ',';
        EXPAND_PRINT_IDX(m_idx4)
    }
    res += "\\rangle}";
    return res;
}


int
ex_symbidx::get_lag_max(bool) const
{
    return INT_MIN;
}


int
ex_symbidx::get_lag_min(bool) const
{
    return INT_MAX;
}


ptr_base
ex_symbidx::diff(const ptr_base &p) const
{
    if (p->type() != SYMBIDX) return ex_num::zero();
    const ex_symbidx *pp = p.get<ex_symbidx>();
    if (m_hash != pp->m_hash) return ex_num::zero();
    if (m_noid != pp->m_noid) return ex_num::zero();
    ptr_base res = ex_delta::create(m_idx1, pp->m_idx1);
    if (m_noid > 1) { res = mk_mul(res, ex_delta::create(m_idx2, pp->m_idx2)); }
    if (m_noid > 2) { res = mk_mul(res, ex_delta::create(m_idx3, pp->m_idx3)); }
    if (m_noid > 3) { res = mk_mul(res, ex_delta::create(m_idx4, pp->m_idx4)); }
    return res;
}



bool
ex_symbidx::has(const ptr_base &what, search_flag f, bool exact_idx) const
{
    if (what->type() != SYMBIDX) return false;
    const ex_symbidx *w = what.get<ex_symbidx>();
    if (compareT(m_hash, w->m_hash)) return false;
    if (compareT(m_noid, w->m_noid)) return false;
    if (exact_idx) {
        if (compareT(m_idx1, w->m_idx1)) return false;
        if (compareT(m_idx2, w->m_idx2)) return false;
        if (compareT(m_idx3, w->m_idx3)) return false;
        if (compareT(m_idx4, w->m_idx4)) return false;
    } else {
        if ((m_idx1 > 0) && (w->m_idx1 < 0)) return false;
        if ((m_idx1 < 0) && (w->m_idx1 < 0) && (m_idx1 != w->m_idx1)) return false;
        if ((m_idx2 > 0) && (w->m_idx2 < 0)) return false;
        if ((m_idx2 < 0) && (w->m_idx2 < 0) && (m_idx2 != w->m_idx2)) return false;
        if ((m_idx3 > 0) && (w->m_idx3 < 0)) return false;
        if ((m_idx3 < 0) && (w->m_idx3 < 0) && (m_idx3 != w->m_idx3)) return false;
        if ((m_idx4 > 0) && (w->m_idx4 < 0)) return false;
        if ((m_idx4 < 0) && (w->m_idx4 < 0) && (m_idx4 != w->m_idx4)) return false;
    }
    switch (f) {
        case EXACT_T:
        case ANY_T:
            return true;
        case DIFF_T:
        case LEAD_T:
        case LAG_T:
            return false;
        default:
            INTERNAL_ERROR
    }
    return false;
}


bool
ex_symbidx::hasidx(int idx) const
{
    return ((m_idx1 == idx) || (m_idx2 == idx)
            || (m_idx3 == idx) || (m_idx4 == idx)) ? true : false;
}


ptr_base
ex_symbidx::substidx(unsigned what, int with) const
{
    int i1, i2, i3, i4;
    i1 = (m_idx1 == (int) what) ? with : m_idx1;
    i2 = (m_idx2 == (int) what) ? with : m_idx2;
    i3 = (m_idx3 == (int) what) ? with : m_idx3;
    i4 = (m_idx4 == (int) what) ? with : m_idx4;
    return create(m_hash, m_noid, i1, i2, i3, i4);
}










