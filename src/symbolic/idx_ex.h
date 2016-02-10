/***********************************************************
 * (c) Kancelaria Prezesa Rady Ministrów 2012-2015         *
 * Treść licencji w pliku 'LICENCE'                        *
 *                                                         *
 * (c) Chancellery of the Prime Minister 2012-2015         *
 * License terms can be found in the file 'LICENCE'        *
 *                                                         *
 * Author: Grzegorz Klima                                  *
 ***********************************************************/

/** \file idx_ex.h
 * \brief Indexing expression.
 */


#ifndef SYMBOLIC_IDX_EX_H

#define SYMBOLIC_IDX_EX_H

#include <string>
#include <decl.h>
#include <idx_set.h>
#include <ops.h>
#include <iostream>


namespace symbolic {

using std::string;
using std::set;

/// Indexing expression
class idx_ex {
  public:
    /// Default constructor.
    idx_ex() : m_id(0), m_excl_id(0) { ; }
    /// Constructor
    idx_ex(const string &id, idx_set is);
    /// Constructor
    idx_ex(const string &id, idx_set is, const string &exclid, bool fixed = false);
    /// Constructor
    idx_ex(const char *id, idx_set is, const char *exclid, bool fixed = false);
    /// Destructor
    ~idx_ex() { ; }

    /// Is it a valid indexing expression?
    operator bool() const { return m_id; }

    /// LaTeX string representation
    std::string tex(internal::print_flag pflag = internal::DEFAULT) const;
    /// String representation
    std::string str(internal::print_flag pflag = internal::DEFAULT) const;

    /// Index substitution
    idx_ex substidx(unsigned what, int with) const;
    /// Index substitution
    idx_ex substidx(const string &what, const string &with, bool fixed = true) const;

    /// Comparison
    int compare(const idx_ex&, bool strict = true) const;

    /// Get index
    unsigned get_id() const { return m_id; }
    /// Get excluded index
    unsigned get_excl_id() const { return (m_excl_id > 0) ? (unsigned) m_excl_id : 0; }

  private:
    // Private constructor
    idx_ex(unsigned id, idx_set is, int excl_id) :
        m_id(id), m_set(is), m_excl_id(excl_id) { ; }
    // Index
    unsigned m_id;
    // Index set
    idx_set m_set;
    // Exclusion
    int m_excl_id;

    friend class internal::ex_sum;
    friend class internal::ex_prod;
    friend class internal::ex_idx;
    friend internal::ptr_base internal::add_idx(const internal::ptr_base &e, const idx_ex &ie);
    friend vec_ex expand(const symbolic::ex &e);

}; /* class idx_ex */

/// Output
std::ostream& operator<<(std::ostream&, const idx_ex&);


} /* namespace symbolic */

#endif /* SYMBOLIC_IDX_EX_H */
