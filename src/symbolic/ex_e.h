/*****************************************************************************
 * This file is a part of gEcon.                                             *
 *                                                                           *
 * (c) Chancellery of the Prime Minister of the Republic of Poland 2012-2015 *
 * (c) Grzegorz Klima, Karol Podemski, Kaja Retkiewicz-Wijtiwiak 2015-2018   *
 * License terms can be found in the file 'LICENCE'                          *
 *                                                                           *
 * Author: Grzegorz Klima                                                    *
 *****************************************************************************/

/** \file ex_e.h
 * \brief Expected value (conditional - time indexed).
 */

#ifndef SYMBOLIC_EX_E_H

#define SYMBOLIC_EX_E_H

#include <ex_base.h>
#include <ptr_base.h>


namespace symbolic {
namespace internal {


/// Class representing expected value
class ex_e : public ex_base {
  public:
    /// Constructor from expression and lag
    ex_e(const ptr_base &e, int l);
    /// Destructor
    ~ex_e() { ; }

    /// Constructor from expression and lag
    static ptr_base create(const ptr_base &e, int l);
    /// Free memory (assumes that ptr is actually pointer to ex_e)
    static void destroy(ex_base *ptr);

    /// Comparison
    int compare(const ex_e&) const;

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

    /// Lag
    ptr_base lag(int l) const;
    /// Get lag
    int get_lag() const { return m_lag; }
    /// Get argument
    const ptr_base& get_arg() const { return m_arg; }

  private:
    // No default constructor
    ex_e();
    // Key
    ptr_base m_arg;
    // Lag
    int m_lag;

}; /* class ex_e */


} /* namespace internal */
} /* namespace symbolic */

#endif /* SYMBOLIC_EX_E_H */
