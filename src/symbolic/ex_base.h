/*****************************************************************************
 * This file is a part of gEcon.                                             *
 *                                                                           *
 * (c) Chancellery of the Prime Minister of the Republic of Poland 2012-2015 *
 * (c) Grzegorz Klima, Karol Podemski, Kaja Retkiewicz-Wijtiwiak 2015-2018   *
 * License terms can be found in the file 'LICENCE'                          *
 *                                                                           *
 * Author: Grzegorz Klima                                                    *
 *****************************************************************************/

/** \file ex_base.h
 * \brief Abstract base class for expressions.
 */

#ifndef SYMBOLIC_EX_BASE_H

#define SYMBOLIC_EX_BASE_H

#include <string>
#include <stringhash.h>
#include <decl.h>
#include <number.h>


namespace symbolic {
namespace internal {



/// Abstract base class for expressions
class ex_base {

  public:
    /// Constructor
    explicit ex_base(unsigned type) : m_type(type) { ; }
    /// Destructor
    virtual ~ex_base() { ; }

    /// String representation.
    virtual std::string str(int pflag = DEFAULT, bool c_style = false) const = 0;
    /// String representation using string 2 string map (name substitution).
    virtual std::string strmap(const map_str_str&, bool c_style = false) const = 0;
    /// LaTeX string representation.
    virtual std::string tex(int pflag = DEFAULT) const = 0;

    /// Max lag in expression.
    virtual int get_lag_max(bool stop_on_E = false) const = 0;
    /// Min lag in expression.
    virtual int get_lag_min(bool stop_on_E = false) const = 0;

    /// Derivative w.r.t. a variable.
    virtual ptr_base diff(const ptr_base&) const = 0;

    /// Does expression have a given subexpression?
    virtual bool has(const ptr_base &what, search_flag f, bool exact_idx) const = 0;
    /// Does expression have a given index?
    virtual bool hasidx(int idx) const = 0;

    /// Is 0?
    bool is0() const;
    /// Is 1?
    bool is1() const;
    /// Is -1?
    bool ism1() const;
    /// Return value if number, NaN otherwise.
    Number val() const;

    /// Type
    unsigned type() const { return (m_type & 0xffffff00); }
    /// Flag
    unsigned flag() const { return (m_type & 0xff); }
    /// Hash value
    unsigned hash() const { return m_hash; }
#ifdef DEBUG_SYMBOLIC
    /// Information for debugging
    void info() const;
#endif /* DEBUG_SYMBOLIC */

  protected:
    /// No default constructor
    ex_base() { ; }
    /// Type & flags
    unsigned m_type;
    /// Hash value
    unsigned m_hash;

}; /* class ex_base */

} /* namespace internal */
} /* namespace symbolic */

#endif /* SYMBOLIC_EX_BASE_H */
