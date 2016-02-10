/***********************************************************
 * (c) Kancelaria Prezesa Rady Ministrów 2012-2015         *
 * Treść licencji w pliku 'LICENCE'                        *
 *                                                         *
 * (c) Chancellery of the Prime Minister 2012-2015         *
 * License terms can be found in the file 'LICENCE'        *
 *                                                         *
 * Author: Grzegorz Klima                                  *
 ***********************************************************/

/** \file ex_symb.h
 * \brief Symbols.
 */

#ifndef SYMBOLIC_EX_SYMB_H

#define SYMBOLIC_EX_SYMB_H

#include <ex_base.h>
#include <ptr_base.h>


namespace symbolic {
namespace internal {

/// Class representing symbol
class ex_symb : public ex_base {
  public:
    /// Constructor from string
    explicit ex_symb(const std::string &n);
    /// Destructor
    ~ex_symb() { ; }

    /// Constructor from string.
    static ptr_base create(const std::string &n);
    /// Free memory (assumes that ptr is actually pointer to ex_symb)
    static void destroy(ex_base *ptr);

    /// Comparison
    int compare(const ex_symb&) const;

    /// Get name
    std::string get_name() const;

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

  private:
    // No default constructor
    ex_symb();
    // Constructor from hash value
    ex_symb(unsigned n);
    // Constructor from hash values
    static ptr_base create(unsigned n);

    friend ptr_base append_name(const ptr_base &p, const std::string &s);
    friend ptr_base add_idx(const ptr_base &e, const idx_ex &ie);

}; /* class ex_symb */

} /* namespace internal */
} /* namespace symbolic */

#endif /* SYMBOLIC_EX_SYMB_H */
