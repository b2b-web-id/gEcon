/*****************************************************************************
 * This file is a part of gEcon.                                             *
 *                                                                           *
 * (c) Chancellery of the Prime Minister of the Republic of Poland 2012-2015 *
 * (c) Grzegorz Klima, Karol Podemski, Kaja Retkiewicz-Wijtiwiak 2015-2018   *
 * License terms can be found in the file 'LICENCE'                          *
 *                                                                           *
 * Author: Grzegorz Klima                                                    *
 *****************************************************************************/

/** \file ex_vart.h
 * \brief Time indexed variables.
 */

#ifndef SYMBOLIC_EX_VART_H

#define SYMBOLIC_EX_VART_H

#include <ex_base.h>


namespace symbolic {
namespace internal {


/// Class representing symbol
class ex_vart : public ex_base {
  public:
    /// Constructor from string and lag
    ex_vart(const std::string &n, int l);
    /// Destructor
    ~ex_vart() { ; }

    /// Return a copy with time index set to 0
    ptr_base copy0() const;

    /// Constructor from string and lag
    static ptr_base create(const std::string &n, int l);
    /// Free memory (assumes that ptr is actually pointer to ex_vart)
    static void destroy(ex_base *ptr);

    /// Comparison
    int compare(const ex_vart&) const;
    /// Name comparison
    int compare_name(const ex_vart&) const;

    /// Get name
    std::string get_name() const;

    /// String representation
    virtual std::string str(int pflag, bool c_style) const;
    /// String representation using string 2 string map (name substitution).
    /// For time indexed variables this always drops time indices.
    virtual std::string strmap(const map_str_str&, bool c_style) const;
    /// LaTeX string representation
    virtual std::string tex(int pflag) const;
    /// Max lag in expression
    virtual int get_lag_max(bool stop_on_E = false) const;
    /// Min lag in expression
    virtual int get_lag_min(bool stop_on_E = false) const;

    /// Derivative w.r.t. a variable
    virtual ptr_base diff(const ptr_base&) const;

    /// Does expression have a given subexpression?
    virtual bool has(const ptr_base &what, search_flag f, bool exact_idx) const;
    /// Does expression have a given index?
    virtual bool hasidx(int idx) const;

    /// Lag
    ptr_base lag(int l) const;
    /// Get lag
    int get_lag() const { return m_lag; }

  private:
    // No default constructor
    ex_vart() : ex_base(NUL) { ; }
    // Constructor given hash and lag
    ex_vart(unsigned i, int l);
    // Constructor given hash and lag
    static ptr_base create(unsigned i, int l);
    // Lag
    int m_lag;

    friend ptr_base append_name(const ptr_base &p, const std::string &s);
    friend ptr_base add_idx(const ptr_base &e, const idx_ex &ie);
    friend ptr_base lag0(const ptr_base &e);

}; /* class ex_vart */


} /* namespace internal */
} /* namespace symbolic */

#endif /* SYMBOLIC_EX_VART_H */
