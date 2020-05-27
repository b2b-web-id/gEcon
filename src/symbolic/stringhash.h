/*****************************************************************************
 * This file is a part of gEcon.                                             *
 *                                                                           *
 * (c) Chancellery of the Prime Minister of the Republic of Poland 2012-2015 *
 * (c) Grzegorz Klima, Karol Podemski, Kaja Retkiewicz-Wijtiwiak 2015-2018   *
 * License terms can be found in the file 'LICENCE'                          *
 *                                                                           *
 * Author: Grzegorz Klima                                                    *
 *****************************************************************************/

/** \file stringhash.h
 * \brief Hash values for strings.
 */

#ifndef SYMBOLIC_STRINGHASH_H

#define SYMBOLIC_STRINGHASH_H

#include <string>
#include <map>
#include <set>


namespace symbolic {
namespace internal {

/// Singleton class handling hash values for strings
class stringhash
{
  public:
    /// Given string get hash values
    unsigned get_hash(const std::string &str);
    /// Given hash value retrieve reference to string
    const std::string& get_str(unsigned id) const;
    /// Given string hash value return hash value of a string with appended underscore.
    unsigned append_underscore(unsigned id);
    /// Given string hash value return hash value of a string with another string appended.
    unsigned append(unsigned id, const std::string &s);
    /// Given string hash value return hash value of a string with appended underscore.
    bool has_underscore(unsigned id) const;
    /// Get instance
    static stringhash& get_instance()
    {
        static stringhash instance;
        return instance;
    }

  private:
    // Maps
    std::map<std::string, unsigned> m_map2int;
    std::map<unsigned, std::string> m_map2str;
    // Strings with underscores at the end (used internally).
    std::set<unsigned> m_underscore;
    // Index
    unsigned m_ind;
    // Constructor is private
    stringhash() : m_ind(0) { ; }
    // Private, too. Not implemented.
    stringhash(stringhash const& copy);
    stringhash& operator=(stringhash const& copy);

}; /* class stringhash */


} /* namespace internal */
} /* namespace symbolic */

#endif /* SYMBOLIC_STRINGHASH_H */
