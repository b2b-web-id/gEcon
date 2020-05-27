/*****************************************************************************
 * This file is a part of gEcon.                                             *
 *                                                                           *
 * (c) Chancellery of the Prime Minister of the Republic of Poland 2012-2015 *
 * (c) Grzegorz Klima, Karol Podemski, Kaja Retkiewicz-Wijtiwiak 2015-2018   *
 * License terms can be found in the file 'LICENCE'                          *
 *                                                                           *
 * Author: Grzegorz Klima                                                    *
 *****************************************************************************/

/** \file ex_prod.h
 * \brief Product over indices.
 */

#ifndef SYMBOLIC_EX_PROD_H

#define SYMBOLIC_EX_PROD_H

#include <ex_base.h>
#include <ptr_base.h>
#include <ex.h>
#include <idx_ex.h>


namespace symbolic {
namespace internal {


/// Class representing expected value
class ex_prod : public ex_base {
  public:
    /// Constructor from indexing expression and expression
    ex_prod(const idx_ex &ie, const ptr_base &e);
    /// Destructor
    ~ex_prod() { ; }

    /// Constructor from from indexing expression and expression
    static ptr_base create(const idx_ex &ie, const ptr_base &e);
    /// Free memory (assumes that ptr is actually pointer to ex_prod)
    static void destroy(ex_base *ptr);

    /// Redindex
    ptr_base reindex(bool force = false) const;

    /// Comparison
    int compare(const ex_prod&) const;

    /// String representation
    virtual std::string str(int pflag, bool c_style) const;
    /// String representation using string 2 string map (name substitution).
    /// This ignores expected vlau operators.
    virtual std::string strmap(const map_str_str&, bool c_style) const;
    /// LaTeX string representation
    virtual std::string tex(int pflag) const;
    /// Max lag in expression
    virtual int get_lag_max(bool stop_on_E = false) const;
    /// Min lag in expression
    virtual int get_lag_min(bool stop_on_E = false) const;

    /// Derivative wrt a variable
    virtual ptr_base diff(const ptr_base&) const;

    /// Does expression have a given subexpression?
    virtual bool has(const ptr_base &what, search_flag f, bool exact_idx) const;
    /// Does expression have a given index?
    virtual bool hasidx(int idx) const;

    /// Expansion
    ptr_base expand() const;

    /// Get expression
    const ptr_base& get_e() const { return m_e; }
    /// Get indexing expression
    const idx_ex& get_ie() const { return m_ie; }

  private:
    // No default constructor
    ex_prod();
    // Indexing expression
    idx_ex m_ie;
    // Expression
    ptr_base m_e;

    friend ptr_base substidx(const ptr_base &e, unsigned what, int with);

}; /* class ex_sum */


} /* namespace internal */
} /* namespace symbolic */

#endif /* SYMBOLIC_EX_PROD_H */
