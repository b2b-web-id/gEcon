/*****************************************************************************
 * This file is a part of gEcon.                                             *
 *                                                                           *
 * (c) Chancellery of the Prime Minister of the Republic of Poland 2012-2015 *
 * (c) Grzegorz Klima, Karol Podemski, Kaja Retkiewicz-Wijtiwiak 2015-2018   *
 * License terms can be found in the file 'LICENCE'                          *
 *                                                                           *
 * Author: Grzegorz Klima                                                    *
 *****************************************************************************/

/** \file stringhash.cpp
 * \brief Hash values for strings.
 */

#include <stringhash.h>
#include <error.h>

using namespace symbolic;
using namespace symbolic::internal;

namespace {

inline
unsigned trans(char c)
{
    if (c == '\0') return 0;
    if (c == '_') return 1;
    if ((c >= '0') && (c <= '9')) return (unsigned) c - 46;
    if ((c >= 'A') && (c <= 'Z')) return (unsigned) c - 27;
    if ((c >= 'a') && (c <= 'z')) return (unsigned) c - 85;
    return 0;
}

} /* namespace */


unsigned
stringhash::get_hash(const std::string &str)
{
    if (!str.size()) return 0;

    std::map<std::string, unsigned>::const_iterator it;
    it = m_map2int.find(str);
    if (it != m_map2int.end()) {
        return it->second;
    }

    if (m_ind >= 524288) {
        throw std::bad_alloc();
    }

    unsigned ids = (trans(str[0]) << 25) + (trans(str[1]) << 19);
    ids += (++m_ind);
    m_map2int[str] = ids;
    m_map2str[ids] = str;

    return ids;
}


const std::string&
stringhash::get_str(unsigned id) const
{
    static const std::string nul(""), quest("?");
    if (id == 0) return nul;
    std::map<unsigned, std::string>::const_iterator it;
    it = m_map2str.find(id);
    if (it == m_map2str.end()) INTERNAL_ERROR
    return it->second;
}



unsigned
stringhash::append_underscore(unsigned id)
{
    unsigned idn = get_hash(get_str(id) + '_');
    m_underscore.insert(idn);
    return idn;
}


unsigned
stringhash::append(unsigned id, const std::string &s)
{
    unsigned idn = get_hash(get_str(id) + s);
    return idn;
}



bool
stringhash::has_underscore(unsigned id) const
{
    std::set<unsigned>::const_iterator it;
    it = m_underscore.find(id);
    if (it != m_underscore.end()) return true;
    return false;
}












