/*****************************************************************************
 * This file is a part of gEcon.                                             *
 *                                                                           *
 * (c) Chancellery of the Prime Minister of the Republic of Poland 2012-2015 *
 * (c) Grzegorz Klima, Karol Podemski, Kaja Retkiewicz-Wijtiwiak 2015-2018   *
 * License terms can be found in the file 'LICENCE'                          *
 *                                                                           *
 * Author: Grzegorz Klima                                                    *
 *****************************************************************************/

/** \file ex_add.h
 * \brief Addition.
 */

#ifndef SYMBOLIC_EX_ADD_H

#define SYMBOLIC_EX_ADD_H


#include <ex_base.h>
#include <num_ex_pair_vec.h>
#include <ex.h>
#include <triplet.h>


namespace symbolic {
namespace internal {


/// Class representing addition operation
class ex_add : public ex_base {
  public:
    /// Constructor from an expression.
    explicit ex_add(const ptr_base &p);
    /// Constructor from a scalar / expression pair.
    ex_add(const Number &s, const ptr_base &p);
    /// Constructor from two arguments.
    ex_add(const ptr_base &a, const ptr_base &b);
    /// Constructor from num_ex_pair_vec
    explicit ex_add(const num_ex_pair_vec &ops, bool try_reduce = false);
    /// Destructor
    virtual ~ex_add() { ; }

    /// Comparison
    int compare(const ex_add&) const;

    /// Constructor from an expression.
    static ptr_base create(const ptr_base &p);
    /// Constructor from a scalar expression pair.
    static ptr_base create(const Number &s, const ptr_base &p);
    /// Constructor from two arguments.
    static ptr_base create(const ptr_base &a, const ptr_base &b);
    /// Constructor from num_ex_pair_vec.
    static ptr_base create(const num_ex_pair_vec &ops, bool try_reduce = false);

    /// Free memory (assumes that ptr is acutally pointer to ex_add)
    static void destroy(ex_base *ptr);

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

    /// Get operands
    const num_ex_pair_vec &get_ops() const { return m_ops; }
    /// No. of operands
    int no_ops() const { return m_ops.size(); }

  private:
    // No default constructor
    ex_add();
    // Ops
    num_ex_pair_vec m_ops;

    void update_flags();

}; /* class ex_add */

} /* namespace internal */
} /* namespace symbolic */

#endif /* SYMBOLIC_EX_ADD_H */
