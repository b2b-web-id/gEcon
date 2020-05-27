/*****************************************************************************
 * This file is a part of gEcon.                                             *
 *                                                                           *
 * (c) Chancellery of the Prime Minister of the Republic of Poland 2012-2015 *
 * (c) Grzegorz Klima, Karol Podemski, Kaja Retkiewicz-Wijtiwiak 2015-2018   *
 * License terms can be found in the file 'LICENCE'                          *
 *                                                                           *
 * Author: Grzegorz Klima                                                    *
 *****************************************************************************/

/** \file ex_vartidx.cpp
 * \brief Time indexed variables with additional indices
 */

#include <error.h>
#include <ex_vartidx.h>
#include <ex_num.h>
#include <ops.h>
#include <cmp.h>
#include <stringhash.h>
#include <utils.h>
#include <climits>


using namespace symbolic;
using namespace symbolic::internal;


ex_vartidx::ex_vartidx(const std::string &n, int l,
                       bool fixed1, const std::string &i1)
    : ex_base(VARTIDX | ((l == INT_MIN) ? 0 : HAST) | SINGLE), m_lag(l)
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


ex_vartidx::ex_vartidx(const std::string &n, int l,
                       bool fixed1, const std::string &i1,
                       bool fixed2, const std::string &i2)
    : ex_base(VARTIDX | ((l == INT_MIN) ? 0 : HAST) | SINGLE), m_lag(l)
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


ex_vartidx::ex_vartidx(const std::string &n, int l,
                       bool fixed1, const std::string &i1,
                       bool fixed2, const std::string &i2,
                       bool fixed3, const std::string &i3)
    : ex_base(VARTIDX | ((l == INT_MIN) ? 0 : HAST) | SINGLE), m_lag(l)
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


ex_vartidx::ex_vartidx(const std::string &n, int l,
                       bool fixed1, const std::string &i1,
                       bool fixed2, const std::string &i2,
                       bool fixed3, const std::string &i3,
                       bool fixed4, const std::string &i4)
    : ex_base(VARTIDX | ((l == INT_MIN) ? 0 : HAST) | SINGLE), m_lag(l)
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


ex_vartidx::ex_vartidx(unsigned n, int l, unsigned no, int i1, int i2, int i3, int i4)
    : ex_base(VARTIDX | ((l == INT_MIN) ? 0 : HAST) | SINGLE), m_lag(l)
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
ex_vartidx::create(const std::string &n, int l,
                   bool fixed1, const std::string &i1)
{
    return ptr_base(new ex_vartidx(n, l, fixed1, i1));
}



ptr_base
ex_vartidx::create(const std::string &n, int l,
                   bool fixed1, const std::string &i1,
                   bool fixed2, const std::string &i2)
{
    return ptr_base(new ex_vartidx(n, l, fixed1, i1, fixed2, i2));
}



ptr_base
ex_vartidx::create(const std::string &n, int l,
                   bool fixed1, const std::string &i1,
                   bool fixed2, const std::string &i2,
                   bool fixed3, const std::string &i3)
{
    return ptr_base(new ex_vartidx(n, l, fixed1, i1, fixed2, i2, fixed3, i3));
}



ptr_base
ex_vartidx::create(const std::string &n, int l,
                   bool fixed1, const std::string &i1,
                   bool fixed2, const std::string &i2,
                   bool fixed3, const std::string &i3,
                   bool fixed4, const std::string &i4)
{
    return ptr_base(new ex_vartidx(n, l, fixed1, i1, fixed2, i2, fixed3, i3,
                                   fixed4, i4));
}


ptr_base
ex_vartidx::create(unsigned n, int l, unsigned no, int i1, int i2, int i3, int i4)
{
    return ptr_base(new ex_vartidx(n, l, no, i1, i2, i3, i4));
}

void
ex_vartidx::destroy(ex_base *ptr)
{
    delete ptr;
}


ptr_base
ex_vartidx::copy0() const
{
    return ptr_base(new ex_vartidx(m_hash, 0,
                                   m_noid, m_idx1, m_idx2, m_idx3, m_idx4));
}


int
ex_vartidx::compare(const ex_vartidx &b) const
{
    int t;
    if ((t = compareT(m_lag, b.m_lag))) return t;
    if ((t = compareT(m_hash, b.m_hash))) return t;
    if ((t = compareT(m_noid, b.m_noid))) return t;
    if ((t = compareT(m_idx1, b.m_idx1))) return ((m_idx1 < 0) && (b.m_idx1 < 0)) ? -t : t;
    if ((t = compareT(m_idx2, b.m_idx2))) return ((m_idx2 < 0) && (b.m_idx2 < 0)) ? -t : t;
    if ((t = compareT(m_idx3, b.m_idx3))) return ((m_idx3 < 0) && (b.m_idx3 < 0)) ? -t : t;
    if ((t = compareT(m_idx4, b.m_idx4))) return ((m_idx4 < 0) && (b.m_idx4 < 0)) ? -t : t;
    return 0;
}


int
ex_vartidx::compare_name(const ex_vartidx &b) const
{
    return compareT(m_hash, b.m_hash);
}



std::string
ex_vartidx::get_name() const
{
    return stringhash::get_instance().get_str(m_hash);
}



std::string
ex_vartidx::str(int pflag, bool c_style) const
{
    if (pflag & INDEXING_ONLY) return std::string();
    stringhash &ref = stringhash::get_instance();
    std::string res = get_name();
    if (pflag & DROP_IDX) {
    } else if (pflag & CONVERT_IDX) {
#ifdef EXPAND_PRINT_IDX
#undef EXPAND_PRINT_IDX
#endif
#define EXPAND_PRINT_IDX(id) \
if (id < 0) { res += ref.get_str(-id); \
} else { res += ref.get_str(id); }
            res += "__";
            EXPAND_PRINT_IDX(m_idx1)
            if (m_noid > 1) { res += "__"; EXPAND_PRINT_IDX(m_idx2) }
            if (m_noid > 2) { res += "__"; EXPAND_PRINT_IDX(m_idx3) }
            if (m_noid > 3) { res += "__"; EXPAND_PRINT_IDX(m_idx4) }
    } else if (pflag & DROP_QUOTES) {
#ifdef EXPAND_PRINT_IDX
#undef EXPAND_PRINT_IDX
#endif
#define EXPAND_PRINT_IDX(id) \
if (id < 0) { res += ref.get_str(-id); \
} else { res += ref.get_str(id); }
            res += '<';
            EXPAND_PRINT_IDX(m_idx1)
            if (m_noid > 1) { res += ','; EXPAND_PRINT_IDX(m_idx2) }
            if (m_noid > 2) { res += ','; EXPAND_PRINT_IDX(m_idx3) }
            if (m_noid > 3) { res += ','; EXPAND_PRINT_IDX(m_idx4) }
            res += '>';
    } else {
#ifdef EXPAND_PRINT_IDX
#undef EXPAND_PRINT_IDX
#endif
#define EXPAND_PRINT_IDX(id) \
if (id < 0) { res += '\''; res += ref.get_str(-id); res += '\''; \
} else { res += ref.get_str(id); }
            res += '<';
            EXPAND_PRINT_IDX(m_idx1)
            if (m_noid > 1) { res += ','; EXPAND_PRINT_IDX(m_idx2) }
            if (m_noid > 2) { res += ','; EXPAND_PRINT_IDX(m_idx3) }
            if (m_noid > 3) { res += ','; EXPAND_PRINT_IDX(m_idx4) }
            res += '>';
    }
    if (pflag & DROP_T) {
        return res;
    } else if (pflag & CONVERT_T) {
        if (m_lag == INT_MIN) return res + "__ss";
        if (m_lag < 0) return res + "__tm" + num2str(-m_lag);
        if (m_lag == 0) return res + "__t";
        return res + "__tp" + num2str(m_lag);
    } else {
        if (m_lag == INT_MIN) return res + "[ss]";
        else if (m_lag == 0) return res + "[]";
        else return res + '[' + num2str(m_lag) + ']';
    }
}


std::string
ex_vartidx::strmap(const map_str_str &mss, bool c_style) const
{
    std::string name = str(CONVERT_IDX | DROP_T, c_style);
    map_str_str::const_iterator it;
    it = mss.find(name);
    if (it == mss.end()) {
        return name;
    }
    return it->second;
}


std::string
ex_vartidx::tex(int pflag) const
{
    if (pflag & INDEXING_ONLY) return std::string();
    stringhash &ref = stringhash::get_instance();
    std::string res = "{" + str2tex(get_name());
#ifdef EXPAND_PRINT_IDX
#undef EXPAND_PRINT_IDX
#endif
#define EXPAND_PRINT_IDX(id) \
if (id < 0) { res += str2tex(ref.get_str(-id), false); \
} else { res += str2tex(ref.get_str(id)); }
    res += "}^{\\langle ";
    EXPAND_PRINT_IDX(m_idx1)
    if (m_noid > 1) { res += ','; EXPAND_PRINT_IDX(m_idx2) }
    if (m_noid > 2) { res += ','; EXPAND_PRINT_IDX(m_idx3) }
    if (m_noid > 3) { res += ','; EXPAND_PRINT_IDX(m_idx4) }
    res += "\\rangle}";
    if (pflag & DROP_T) {
        return res;
    } else {
        if (m_lag == INT_MIN) return res + "_\\mathrm{ss}";
        if (m_lag == 0) return res + "_{t}";
        if (m_lag < 0) return res + "_{t" + num2str(m_lag) + '}';
        return res + "_{t+" + num2str(m_lag) + '}';
    }
}


int
ex_vartidx::get_lag_max(bool) const
{
    return m_lag;
}


int
ex_vartidx::get_lag_min(bool) const
{
    if (m_lag == INT_MIN) return INT_MAX;
    return m_lag;
}


ptr_base
ex_vartidx::diff(const ptr_base &p) const
{
    if (p->type() != VARTIDX) return ex_num::zero();
    const ex_vartidx *pp = p.get<ex_vartidx>();
    if (m_hash != pp->m_hash) return ex_num::zero();
    if (m_noid != pp->m_noid) return ex_num::zero();
    if (m_lag != pp->m_lag) return ex_num::zero();
    ptr_base res = ex_delta::create(m_idx1, pp->m_idx1);
    if (m_noid > 1) { res = mk_mul(res, ex_delta::create(m_idx2, pp->m_idx2)); }
    if (m_noid > 2) { res = mk_mul(res, ex_delta::create(m_idx3, pp->m_idx3)); }
    if (m_noid > 3) { res = mk_mul(res, ex_delta::create(m_idx4, pp->m_idx4)); }
    return res;
}



bool
ex_vartidx::has(const ptr_base &what, search_flag f, bool exact_idx) const
{
    if (what->type() != VARTIDX) return false;
    const ex_vartidx *w = what.get<ex_vartidx>();
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
            return (m_lag == w->m_lag) ? true : false;
        case ANY_T:
            return true;
        case DIFF_T:
            return (m_lag != w->m_lag) ? true : false;
        case LEAD_T:
            return (m_lag > w->m_lag) ? true : false;
        case LAG_T:
            return (m_lag < w->m_lag) ? true : false;
        default:
            INTERNAL_ERROR
    }
    return false;
}


bool
ex_vartidx::hasidx(int idx) const
{
    return ((m_idx1 == idx) || (m_idx2 == idx)
            || (m_idx3 == idx) || (m_idx4 == idx)) ? true : false;
}



ptr_base
ex_vartidx::substidx(unsigned what, int with) const
{
    int i1, i2, i3, i4;
    i1 = (m_idx1 == (int) what) ? with : m_idx1;
    i2 = (m_idx2 == (int) what) ? with : m_idx2;
    i3 = (m_idx3 == (int) what) ? with : m_idx3;
    i4 = (m_idx4 == (int) what) ? with : m_idx4;
    return create(m_hash, m_lag, m_noid, i1, i2, i3, i4);
}


ptr_base
ex_vartidx::lag(int l) const
{
    if (l == INT_MIN) return create(m_hash, l, m_noid, m_idx1, m_idx2, m_idx3, m_idx4);
    return create(m_hash, m_lag + l, m_noid, m_idx1, m_idx2, m_idx3, m_idx4);
}








