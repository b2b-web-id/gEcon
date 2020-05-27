/*****************************************************************************
 * This file is a part of gEcon.                                             *
 *                                                                           *
 * (c) Chancellery of the Prime Minister of the Republic of Poland 2012-2015 *
 * (c) Grzegorz Klima, Karol Podemski, Kaja Retkiewicz-Wijtiwiak 2015-2018   *
 * License terms can be found in the file 'LICENCE'                          *
 *                                                                           *
 * Author: Grzegorz Klima                                                    *
 *****************************************************************************/

/** \file ex_func.h
 * \brief Functions.
 */

#ifndef SYMBOLIC_EX_FUNC_H

#define SYMBOLIC_EX_FUNC_H


#include <ex_base.h>
#include <ptr_base.h>


namespace symbolic {
namespace internal {


/// Class representing number
class ex_func : public ex_base {
  public:
    /// Constructor
    ex_func(func_code c, const ptr_base &arg);
    /// Destructor
    ~ex_func() { ; }

    /// Constructor from two arguments.
    static ptr_base create(func_code c, const ptr_base &arg);
    /// Free memory (assumes that ptr is actually pointer to ex_func)
    static void destroy(ex_base *ptr);

    /// Comparison
    int compare(const ex_func&) const;

    /// String representation
    virtual std::string str(int pflag, bool c_style) const;
    /// String representation using string 2 string map (name substitution).
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

    /// Evaluate function
    static Number eval(func_code c, const Number &x);

    /// Get code
    func_code get_code() const { return m_code; };
    /// Get argument
    const ptr_base& get_arg() const { return m_arg; }

  private:
    // No default constructor
    ex_func() : ex_base(NUL), m_arg(0) { ; }
    // Function code
    func_code m_code;
    // Pointer to an argument
    ptr_base m_arg;

}; /* class ex_func */

} /* namespace internal */
} /* namespace symbolic */

#endif /* SYMBOLIC_EX_FUNC_H */
