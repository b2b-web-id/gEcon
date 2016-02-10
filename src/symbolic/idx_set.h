/***********************************************************
 * (c) Kancelaria Prezesa Rady Ministrów 2012-2015         *
 * Treść licencji w pliku 'LICENCE'                        *
 *                                                         *
 * (c) Chancellery of the Prime Minister 2012-2015         *
 * License terms can be found in the file 'LICENCE'        *
 *                                                         *
 * Author: Grzegorz Klima                                  *
 ***********************************************************/

/** \file idx_set.h
 * \brief Index set.
 */


#ifndef SYMBOLIC_IDX_SET_H

#define SYMBOLIC_IDX_SET_H

#include <decl.h>
#include <idx_set_impl.h>
#include <string>
#include <iostream>


namespace symbolic {

using std::string;
using std::set;

/// Set of indices
class idx_set {
  public:
    // Default constructor.
    idx_set();
    /// Constructor.
    explicit idx_set(const string &name);
    /// Copy constructor.
    idx_set(const idx_set &is);
    /// 'Copy constructor' with new name.
    idx_set(const idx_set &is, const std::string &name);
    /// Destructor.
    ~idx_set();
    /// Assignment.
    idx_set& operator=(const idx_set&);

    /// Is it a valid index set?
    operator bool() const { return m_p; }

    /// Return set with names of elements prefixed with string s.
    idx_set prefix(const std::string &s) const;
    /// Return set with names of elements suffixed with string s.
    idx_set suffix(const std::string &s) const;

    /// Is this a subset
    bool operator<=(const idx_set &b) const;
    /// Is this equal
    bool operator==(const idx_set &b) const;
    /// Is this not equal
    bool operator!=(const idx_set &b) const;

    /// Set intersection
    idx_set set_intersect(const idx_set &b) const;
    /// Set sum
    idx_set set_sum(const idx_set &b) const;
    /// Set difference
    idx_set set_diff(const idx_set &b) const;

    /// Has element
    bool has(const string &s) const { if (m_p) return m_p->has(s); else return false; }
    /// Has element
    bool has(unsigned n) const { if (m_p) return m_p->has(n); else return false; }

    /// Add element
    bool add(const string &s) { if (m_p) return m_p->add(s); else return false; }
    /// Add element
    bool add(unsigned n) { if (m_p) return m_p->add(n); else return false; }
    /// Remove element
    bool rm(const string &s) { if (m_p) return m_p->rm(s); else return false; }
    /// Remove element
    bool rm(unsigned n) { if (m_p) return m_p->rm(n); else return false; }

    /// Retrieve name
    string name() const { if (m_p) return m_p->name(); else return "??"; }

    /// Size
    unsigned size() const { if (m_p) return m_p->size(); else return 0; }

    /// String representation
    string str() const { if (m_p) return m_p->str(); else return "??"; }
    /// LaTeX string representation
    string tex() const { if (m_p) return m_p->tex(); else return "??"; }

  private:
    // Pointer to implementation
    internal::idx_set_impl *m_p;
    // Pointer to ref count.
    int *m_rc;
    // Decrease ref count.
    void dec_rc();

    /// Pointer.
    internal::idx_set_impl* operator->() { return m_p; }
    /// Const pointer.
    const internal::idx_set_impl* operator->() const { return m_p; }

    friend class idx_ex;
    friend class internal::ex_sum;
    friend class internal::ex_prod;
    friend vec_ex expand(const ex &e);

}; /* class idx_set */


/// Stream output
std::ostream& operator<<(std::ostream &s, const idx_set &is);


} /* namespace symbolic */

#endif /* SYMBOLIC_IDX_SET_H */
