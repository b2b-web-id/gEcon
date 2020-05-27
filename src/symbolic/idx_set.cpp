/*****************************************************************************
 * This file is a part of gEcon.                                             *
 *                                                                           *
 * (c) Chancellery of the Prime Minister of the Republic of Poland 2012-2015 *
 * (c) Grzegorz Klima, Karol Podemski, Kaja Retkiewicz-Wijtiwiak 2015-2018   *
 * License terms can be found in the file 'LICENCE'                          *
 *                                                                           *
 * Author: Grzegorz Klima                                                    *
 *****************************************************************************/

/** \file idx_set.cpp
 * \brief Index set.
 */


#include <idx_set.h>
#include <stringhash.h>


using namespace symbolic;
using namespace symbolic::internal;


idx_set::idx_set() : m_p(0), m_rc(0)
{
}


idx_set::idx_set(const string &name)
{
    m_rc = new int;
    *m_rc = 1;
    m_p = new idx_set_impl(name);
}


idx_set::~idx_set()
{
    dec_rc();
}


idx_set::idx_set(const idx_set &is)
{
    m_p = is.m_p;
    m_rc = is.m_rc;
    if (m_rc) ++(*m_rc);
}

idx_set::idx_set(const idx_set &is, const std::string &name)
{
    if (!is.m_p) {
        m_rc = new int;
        *m_rc = 1;
        m_p = new idx_set_impl(name);
        return;
    } else if (is.m_p->m_name) {
        m_rc = new int;
        *m_rc = 1;
        m_p = new idx_set_impl(*is.m_p);
    } else {
        m_p = is.m_p;
        m_rc = is.m_rc;
        if (m_rc) ++(*m_rc);
    }
    if (m_p) m_p->set_name(name);
}



idx_set&
idx_set::operator=(const idx_set &is)
{
    if (m_p != is.m_p) {
        dec_rc();
        m_p = is.m_p;
        m_rc = is.m_rc;
        if (m_rc) ++(*m_rc);
    }
    return *this;
}


void
idx_set::dec_rc()
{
    if (m_rc) {
        if (--(*m_rc) == 0) {
            delete m_rc;
            delete m_p;
        }
    }
}


idx_set
idx_set::prefix(const std::string &s) const
{
    idx_set res("");
    if (!*this) return res;
    if (s == "") return *this;
    for (idx_set_impl::const_iterator it = m_p->begin(); it != m_p->end(); ++it)
        res.add(s + stringhash::get_instance().get_str(*it));
    return res;
}


idx_set
idx_set::suffix(const std::string &s) const
{
    idx_set res("");
    if (!*this) return res;
    if (s == "") return *this;
    for (idx_set_impl::const_iterator it = m_p->begin(); it != m_p->end(); ++it)
        res.add(stringhash::get_instance().get_str(*it) + s);
    return res;
}



bool
idx_set::operator<=(const idx_set &b) const
{
    if (!operator bool()) return true;
    if (!b) return false;
    for (idx_set_impl::const_iterator it = m_p->begin(); it != m_p->end(); ++it)
        if (!b.has(*it)) return false;
    return true;
}


bool
idx_set::operator==(const idx_set &b) const
{
    if (size() != b.size()) return false;
    return *this <= b;
}


bool
idx_set::operator!=(const idx_set &b) const
{
    return !operator==(b);
}


idx_set
idx_set::set_intersect(const idx_set &b) const
{
    idx_set res("");
    if (!size()) return res;
    if (!b.size()) return res;
    for (idx_set_impl::const_iterator it = m_p->begin(); it != m_p->end(); ++it)
        if (b.has(*it)) res.add(*it);
    return res;
}

idx_set
idx_set::set_sum(const idx_set &b) const
{
    if (!size()) return idx_set(b, "");
    if (!b.size()) return idx_set(*this, "");
    idx_set res(b, "");
    for (idx_set_impl::const_iterator it = m_p->begin(); it != m_p->end(); ++it)
        res.add(*it);
    return res;
}


idx_set
idx_set::set_diff(const idx_set &b) const
{
    idx_set res("");
    if (!size()) return idx_set(b, "");
    if (!b.size()) return idx_set(*this, "");
    for (idx_set_impl::const_iterator it = m_p->begin(); it != m_p->end(); ++it)
        if (!b.has(*it)) res.add(*it);
    return res;
}



std::ostream&
symbolic::operator<<(std::ostream &s, const idx_set &is)
{
    s << is.str() << '\n';
    return s;
}






























