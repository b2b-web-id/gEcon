/*****************************************************************************
 * This file is a part of gEcon.                                             *
 *                                                                           *
 * (c) Chancellery of the Prime Minister of the Republic of Poland 2012-2015 *
 * (c) Grzegorz Klima, Karol Podemski, Kaja Retkiewicz-Wijtiwiak 2015-2018   *
 * License terms can be found in the file 'LICENCE'                          *
 *                                                                           *
 * Author: Grzegorz Klima                                                    *
 *****************************************************************************/

/** \file ex_symbidx.h
 * \brief Indexed symbols.
 */

#ifndef SYMBOLIC_EX_SYMBIDX_H

#define SYMBOLIC_EX_SYMBIDX_H

#include <ex_base.h>
#include <ptr_base.h>


namespace symbolic {
namespace internal {

/// Class representing symbol
class ex_symbidx : public ex_base {
  public:
    /// Constructor from strings
    ex_symbidx(const std::string &n, bool fixed1, const std::string &i1);
    /// Constructor from strings
    ex_symbidx(const std::string &n, bool fixed1, const std::string &i1,
               bool fixed2, const std::string &i2);
    /// Constructor from strings
    ex_symbidx(const std::string &n, bool fixed1, const std::string &i1,
               bool fixed2, const std::string &i2,
               bool fixed3, const std::string &i3);
    /// Constructor from strings
    ex_symbidx(const std::string &n, bool fixed1, const std::string &i1,
               bool fixed2, const std::string &i2,
               bool fixed3, const std::string &i3,
               bool fixed4, const std::string &i4);
    /// Destructor
    ~ex_symbidx() { ; }

    /// Constructor from strings
    static ptr_base create(const std::string &n,
                           bool fixed1, const std::string &i1);
    /// Constructor from strings
    static ptr_base create(const std::string &n,
                           bool fixed1, const std::string &i1,
                           bool fixed2, const std::string &i2);
    /// Constructor from strings
    static ptr_base create(const std::string &n,
                           bool fixed1, const std::string &i1,
                           bool fixed2, const std::string &i2,
                           bool fixed3, const std::string &i3);
    /// Constructor from strings
    static ptr_base create(const std::string &n,
                           bool fixed1, const std::string &i1,
                           bool fixed2, const std::string &i2,
                           bool fixed3, const std::string &i3,
                           bool fixed4, const std::string &i4);
    /// Free memory (assumes that ptr is actually pointer to ex_symbidx)
    static void destroy(ex_base *ptr);

    /// Comparison
    int compare(const ex_symbidx&) const;

    /// Get name
    std::string get_name() const;

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

  private:
    // No default constructor
    ex_symbidx();
    // Constructor from hash values
    ex_symbidx(unsigned n, unsigned no, int i1, int i2, int i3, int i4);
    // Constructor from hash values
    static ptr_base create(unsigned n, unsigned no, int i1, int i2, int i3, int i4);
    // Index substitution
    ptr_base substidx(unsigned what, int with) const;
    // Indices
    unsigned m_noid;
    int m_idx1, m_idx2, m_idx3, m_idx4;

    friend ptr_base append_name(const ptr_base &p, const std::string &s);
    friend ptr_base add_idx(const ptr_base &e, const idx_ex &ie);
    friend ptr_base subst(const ptr_base &e, const ptr_base &what,
                          const ptr_base &with, bool all_leads_lags);
    friend ptr_base substidx(const ptr_base &e, unsigned what, int with);
    friend void collect_idx(const ptr_base &p, std::set<unsigned> &iset);

}; /* class ex_symbidx */

} /* namespace internal */
} /* namespace symbolic */

#endif /* SYMBOLIC_EX_SYMBIDX_H */
