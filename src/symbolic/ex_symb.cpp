/***********************************************************
 * (c) Kancelaria Prezesa Rady Ministrów 2012-2015         *
 * Treść licencji w pliku 'LICENCE'                        *
 *                                                         *
 * (c) Chancellery of the Prime Minister 2012-2015         *
 * License terms can be found in the file 'LICENCE'        *
 *                                                         *
 * Author: Grzegorz Klima                                  *
 ***********************************************************/

/** \file ex_symb.cpp
 * \brief Symbols.
 */

#include <ex_symb.h>
#include <ex_num.h>
#include <stringhash.h>
#include <utils.h>
#include <cmp.h>
#include <climits>


using namespace symbolic;
using namespace symbolic::internal;



ex_symb::ex_symb(const std::string &n) : ex_base(SYMB | SINGLE)
{
    m_hash = stringhash::get_instance().get_hash(n);
}

ptr_base
ex_symb::create(const std::string &n)
{
    return ptr_base(new ex_symb(n));
}


ex_symb::ex_symb(unsigned n) : ex_base(SYMB | SINGLE)
{
    m_hash = n;
}


ptr_base
ex_symb::create(unsigned n)
{
    return ptr_base(new ex_symb(n));
}


void
ex_symb::destroy(ex_base *ptr)
{
    delete ptr;
}


int
ex_symb::compare(const ex_symb &b) const
{
    return compareT(m_hash, b.m_hash);
}


std::string
ex_symb::get_name() const
{
    return stringhash::get_instance().get_str(m_hash);
}



std::string
ex_symb::str(int pflag) const
{
    if (pflag & INDEXING_ONLY) return std::string();
    return get_name();
}


std::string
ex_symb::strmap(const map_str_str &mss) const
{
    map_str_str::const_iterator it;
    std::string name = get_name();
    it = mss.find(name);
    if (it == mss.end()) {
        return name;
    }
    return it->second;
}


std::string
ex_symb::tex(int pflag) const
{
    if (pflag & INDEXING_ONLY) return std::string();
    return str2tex(get_name());
}


int
ex_symb::get_lag_max(bool) const
{
    return INT_MIN;
}


int
ex_symb::get_lag_min(bool) const
{
    return INT_MAX;
}


ptr_base
ex_symb::diff(const ptr_base &p) const
{
    if (p->type() != SYMB) return ex_num::zero();
    const ex_symb *pp = p.get<ex_symb>();
    if (!compare(*pp)) return ex_num::one();
    return ex_num::zero();
}


bool
ex_symb::has(const ptr_base &what, search_flag f, bool) const
{
    if (what->type() != SYMB) return false;
    const ex_symb *w = what.get<ex_symb>();
    if (compareT(m_hash, w->m_hash)) return false;
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
ex_symb::hasidx(int) const
{
    return false;
}




