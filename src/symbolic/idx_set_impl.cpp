/*****************************************************************************
 * This file is a part of gEcon.                                             *
 *                                                                           *
 * (c) Chancellery of the Prime Minister of the Republic of Poland 2012-2015 *
 * (c) Grzegorz Klima, Karol Podemski, Kaja Retkiewicz-Wijtiwiak 2015-2018   *
 * License terms can be found in the file 'LICENCE'                          *
 *                                                                           *
 * Author: Grzegorz Klima                                                    *
 *****************************************************************************/

/** \file idx_set_impl.cpp
 * \brief Index set implementation.
 */


#include <idx_set_impl.h>
#include <stringhash.h>
#include <utils.h>


using namespace symbolic;
using namespace symbolic::internal;


idx_set_impl::idx_set_impl(const string &name) : m_name(stringhash::get_instance().get_hash(name))
{
}


bool
idx_set_impl::add(const string &s)
{
    return insert(stringhash::get_instance().get_hash(s)).second;
}


bool
idx_set_impl::add(unsigned n)
{
    return insert(n).second;
}


bool
idx_set_impl::has(const string &s)
{
    if (find(stringhash::get_instance().get_hash(s)) != end()) {
        return true;
    }
    return false;
}


bool
idx_set_impl::has(unsigned n)
{
    if (find(n) != end()) {
        return true;
    }
    return false;
}



bool
idx_set_impl::rm(const string &s)
{
    unsigned n = stringhash::get_instance().get_hash(s);
    set<unsigned>::iterator it = find(n);
    if (it != end()) {
        erase(it);
        return true;
    }
    return false;
}



bool
idx_set_impl::rm(unsigned n)
{
    set<unsigned>::iterator it = find(n);
    if (it != end()) {
        erase(it);
        return true;
    }
    return false;
}



string
idx_set_impl::name() const
{
    return stringhash::get_instance().get_str(m_name);
}


void
idx_set_impl::set_name(const string &n)
{
    m_name = stringhash::get_instance().get_hash(n);
}



string
idx_set_impl::str() const
{
    string res;
    if (m_name) {
        res = name();
        res += " = ";
    }
    res += "{ ";
    set<unsigned>::const_iterator it;
    for (it = begin(); it != end();) {
        res += '\'' + stringhash::get_instance().get_str(*it) + '\'';
        ++it;
        if (it != end()) res += ", "; else res += ' ';
    }
    res += '}';
    return res;
}


string
idx_set_impl::tex() const
{
    string res = str2tex(name());
    res += " = \\{ ";
    for (set<unsigned>::const_iterator it = begin(); it != end(); ++it) {
        if (it != begin()) res += ", ";
        res += str2tex(stringhash::get_instance().get_str(*it), false);
    }
    res += " \\}";
    return res;
}




































