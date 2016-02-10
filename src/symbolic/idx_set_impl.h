/***********************************************************
 * (c) Kancelaria Prezesa Rady Ministrów 2012-2015         *
 * Treść licencji w pliku 'LICENCE'                        *
 *                                                         *
 * (c) Chancellery of the Prime Minister 2012-2015         *
 * License terms can be found in the file 'LICENCE'        *
 *                                                         *
 * Author: Grzegorz Klima                                  *
 ***********************************************************/

/** \file idx_set_impl.h
 * \brief Index set implementation.
 */


#ifndef SYMBOLIC_IDX_SET_IMPL_H

#define SYMBOLIC_IDX_SET_IMPL_H

#include <decl.h>
#include <string>
#include <set>
#include <iostream>


namespace symbolic {
namespace internal {

using std::string;
using std::set;

/// Set of indices
class idx_set_impl : public set<unsigned> {
  public:
    /// Constructor.
    explicit idx_set_impl(const string &name);
    /// Destructor
    ~idx_set_impl() { ; }

    /// Has element
    bool has(const string &s);
    /// Has element
    bool has(unsigned n);

    /// Add element
    bool add(const string &s);
    /// Add element
    bool add(unsigned n);
    /// Remove element
    bool rm(const string &s);
    /// Remove element
    bool rm(unsigned n);

    /// Set name
    void set_name(const std::string &name);

    /// Retrieve name
    string name() const;

    /// String representation
    string str() const;
    /// LaTeX string representation
    string tex() const;

  private:
    unsigned m_name;

    // Default constructor is private.
    idx_set_impl();

    friend class symbolic::idx_ex;
    friend class symbolic::idx_set;

}; /* class idx_set_impl */



} /* namespace internal */
} /* namespace symbolic */

#endif /* SYMBOLIC_IDX_SET_IMPL_H */
