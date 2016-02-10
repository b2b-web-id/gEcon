/***********************************************************
 * (c) Kancelaria Prezesa Rady Ministrów 2012-2015         *
 * Treść licencji w pliku 'LICENCE'                        *
 *                                                         *
 * (c) Chancellery of the Prime Minister 2012-2015         *
 * License terms can be found in the file 'LICENCE'        *
 *                                                         *
 * Author: Grzegorz Klima                                  *
 ***********************************************************/

/** \file ex_pow.h
 * \brief Powers.
 */

#ifndef SYMBOLIC_EX_POW_H

#define SYMBOLIC_EX_POW_H


#include <ex_base.h>
#include <ptr_base.h>


namespace symbolic {
namespace internal {


/// Class representing power expression.
class ex_pow : public ex_base {
  public:
    /// Constructor from two arguments.
    ex_pow(const ptr_base &a, const ptr_base &b);
    /// Destructor
    virtual ~ex_pow() { ; }

    /// Constructor from two arguments.
    static ptr_base create(const ptr_base &a, const ptr_base &b);
    /// Free memory (assumes that ptr is acutally pointer to ex_pow)
    static void destroy(ex_base *ptr);

    /// Comparison
    int compare(const ex_pow&) const;

    /// String representation
    virtual std::string str(int pflag) const;
    /// String representation using string 2 string map (name substitution).
    virtual std::string strmap(const map_str_str&) const;
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

    /// Get base
    const ptr_base& get_base() const { return m_base; }
    /// Get exponent
    const ptr_base& get_exp() const { return m_exp; }

  private:
    // No default constructor
    ex_pow();
    // Ops
    ptr_base m_base, m_exp;

}; /* class ex_pow */

} /* namespace internal */
} /* namespace symbolic */

#endif /* SYMBOLIC_EX_POW_H */
