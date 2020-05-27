/*****************************************************************************
 * This file is a part of gEcon.                                             *
 *                                                                           *
 * (c) Chancellery of the Prime Minister of the Republic of Poland 2012-2015 *
 * (c) Grzegorz Klima, Karol Podemski, Kaja Retkiewicz-Wijtiwiak 2015-2018   *
 * License terms can be found in the file 'LICENCE'                          *
 *                                                                           *
 * Author: Grzegorz Klima                                                    *
 *****************************************************************************/

/** \file idx_ex.cpp
 * \brief Indexing expression.
 */


#include <decl.h>
#include <error.h>
#include <idx_set.h>
#include <idx_ex.h>
#include <ex.h>
#include <stringhash.h>
#include <cmp.h>
#include <iostream>


using namespace symbolic;
using namespace symbolic::internal;



idx_ex::idx_ex(const string &id, idx_set is) : m_set(is)
{
    stringhash &ref = stringhash::get_instance();
    m_id = ref.get_hash(id);
    m_excl_id = 0;
}


idx_ex::idx_ex(const string &id, idx_set is, const string &exclid, bool fixed)
                : m_set(is)
{
    stringhash &ref = stringhash::get_instance();
    m_id = ref.get_hash(id);
    m_excl_id = ref.get_hash(exclid);
    if (fixed) m_excl_id = -m_excl_id;
}


idx_ex::idx_ex(const char *id, idx_set is, const char *exclid, bool fixed)
                : m_set(is)
{
    stringhash &ref = stringhash::get_instance();
    m_id = ref.get_hash(id);
    m_excl_id = ref.get_hash(exclid);
    if (fixed) m_excl_id = -m_excl_id;
}


std::string
idx_ex::str(print_flag pflag) const
{
    if (pflag != DEFAULT) USER_ERROR("invalid print flag in idx_ex::str()")
    std::string res = "<";
    stringhash &ref = stringhash::get_instance();
    res += ref.get_str(m_id);
    res += "::";
    res += m_set.name();
    if (m_excl_id) {
        res += "\\";
        if (m_excl_id < 0) {
            res += "\'";
            res += ref.get_str(-m_excl_id);
            res += "\'";
        } else {
            res += ref.get_str(m_excl_id);
        }
    }
    res += ">";
    return res;
}


std::string
idx_ex::tex(print_flag pflag) const
{
    if (pflag != DEFAULT) USER_ERROR("invalid print flag in idx_ex::str()")
    std::string res;
    stringhash &ref = stringhash::get_instance();
    res += str2tex(ref.get_str(m_id));
    res += "\\in ";
    res += str2tex(m_set.name());
    if (m_excl_id) {
        res += "\\backslash ";
        if (m_excl_id < 0) {
            res += "\\{\\mathrm{";
            res += str2tex(ref.get_str(-m_excl_id), false);
            res += "}\\}";
        } else {
            res += ref.get_str(m_excl_id);
        }
    }
    return res;
}



idx_ex
idx_ex::substidx(unsigned what, int with) const
{
    return idx_ex(m_id, m_set, (m_excl_id == (int) what) ? with : m_excl_id);
}



idx_ex
idx_ex::substidx(const string &what, const string &with, bool fixed) const
{
    stringhash &ref = stringhash::get_instance();
    return substidx(ref.get_hash(what),
                    (fixed) ? -ref.get_hash(with) : ref.get_hash(with));
}


int
idx_ex::compare(const idx_ex &ie, bool strict) const
{
    int t;
    if ((!(bool) *this) | (!(bool) ie))
        USER_ERROR("comparing void index expression(s)")
//     if ((!ie.m_set) && (!m_set)) return 0;
//     if (!m_set) return -1;
//     if (!ie.m_set) return 1;
    if ((t = compareT(m_set.m_p->m_name, ie.m_set.m_p->m_name))) return t;
    if (strict) { if ((t = compareT(m_id, ie.m_id))) return t; }
    if ((t = compareT(m_excl_id, ie.m_excl_id))) return t;
    return 0;
}



std::ostream&
symbolic::operator<<(std::ostream &os, const idx_ex &ie)
{
    os << ie.str() << '\n';
    return os;
}







