/*****************************************************************************
 * This file is a part of gEcon.                                             *
 *                                                                           *
 * (c) Chancellery of the Prime Minister of the Republic of Poland 2012-2015 *
 * (c) Grzegorz Klima, Karol Podemski, Kaja Retkiewicz-Wijtiwiak 2015-2018   *
 * License terms can be found in the file 'LICENCE'                          *
 *                                                                           *
 * Author: Grzegorz Klima                                                    *
 *****************************************************************************/

/** \file ex_num.h
 * \brief Numbers.
 */

#ifndef SYMBOLIC_EX_NUM_H

#define SYMBOLIC_EX_NUM_H

#include <ex_base.h>
#include <ptr_base.h>


namespace symbolic {
namespace internal {


/// Expression / number
class ex_num : public ex_base {
  public:
    /// Constructor from Number
    explicit ex_num(const Number &n) : ex_base(NUM|SINGLE), m_val(n) { ; }
    /// Destructor
    ~ex_num() { ; }

    /// Constructor from Number
    static ptr_base create(const Number&);
    /// One
    static ptr_base one() { return create(1.); }
    /// Zero
    static ptr_base zero() { return create(0.); }
    /// Free memory (assumes that ptr is actually pointer to ex_num)
    static void destroy(ex_base *ptr);

    /// Value
    Number val() const { return m_val; }

    /// Comparison
    int compare(const ex_num&) const;

    /// String representation
    virtual std::string str(int pflag, bool c_style) const;
    /// String representation using string 2 string map (name substitution).
    virtual std::string strmap(const map_str_str&, bool c_style) const;
    /// LaTeX string representation
    virtual std::string tex(int) const;

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

  private:
    // No default constructor
    ex_num();
    // Value
    Number m_val;

}; /* class ex_num */

} /* namespace internal */
} /* namespace symbolic */

#endif /* SYMBOLIC_EX_NUM_H */
