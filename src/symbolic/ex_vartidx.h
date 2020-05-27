/*****************************************************************************
 * This file is a part of gEcon.                                             *
 *                                                                           *
 * (c) Chancellery of the Prime Minister of the Republic of Poland 2012-2015 *
 * (c) Grzegorz Klima, Karol Podemski, Kaja Retkiewicz-Wijtiwiak 2015-2018   *
 * License terms can be found in the file 'LICENCE'                          *
 *                                                                           *
 * Author: Grzegorz Klima                                                    *
 *****************************************************************************/

/** \file ex_vartidx.h
 * \brief Time indexed variables with additional indices.
 */

#ifndef SYMBOLIC_EX_VARTIDX_H

#define SYMBOLIC_EX_VARTIDX_H

#include <ex_base.h>


namespace symbolic {
namespace internal {


/// Class representing symbol
class ex_vartidx : public ex_base {
  public:
    /// Constructor from string and lag, and strings with indices
    ex_vartidx(const std::string &n, int l,
               bool fixed1, const std::string &i1);
    /// Constructor from string and lag, and strings with indices
    ex_vartidx(const std::string &n, int l,
               bool fixed1, const std::string &i1,
               bool fixed2, const std::string &i2);
    /// Constructor from string and lag, and strings with indices
    ex_vartidx(const std::string &n, int l,
               bool fixed1, const std::string &i1,
               bool fixed2, const std::string &i2,
               bool fixed3, const std::string &i3);
    /// Constructor from string and lag, and strings with indices
    ex_vartidx(const std::string &n, int l,
               bool fixed1, const std::string &i1,
               bool fixed2, const std::string &i2,
               bool fixed3, const std::string &i3,
               bool fixed4, const std::string &i4);
    /// Destructor
    ~ex_vartidx() { ; }

    /// Return a copy with time index set to 0
    ptr_base copy0() const;

    /// Constructor from string and lag, and strings with indices
    static ptr_base create(const std::string &n, int l,
                           bool fixed1, const std::string &i1);
    /// Constructor from string and lag, and strings with indices
    static ptr_base create(const std::string &n, int l,
                           bool fixed1, const std::string &i1,
                           bool fixed2, const std::string &i2);
    /// Constructor from string and lag, and strings with indices
    static ptr_base create(const std::string &n, int l,
                           bool fixed1, const std::string &i1,
                           bool fixed2, const std::string &i2,
                           bool fixed3, const std::string &i3);
    /// Constructor from string and lag, and strings with indices
    static ptr_base create(const std::string &n, int l,
                           bool fixed1, const std::string &i1,
                           bool fixed2, const std::string &i2,
                           bool fixed3, const std::string &i3,
                           bool fixed4, const std::string &i4);
    /// Free memory (assumes that ptr is actually pointer to ex_vartidx)
    static void destroy(ex_base *ptr);

    /// Comparison
    int compare(const ex_vartidx&) const;
    /// Name comparison
    int compare_name(const ex_vartidx&) const;

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
    ex_vartidx();
    // Constructor given hash and lag
    ex_vartidx(unsigned n, int l, unsigned no, int i1, int i2, int i3, int i4);
    // Constructor given hash and lag
    static ptr_base create(unsigned n, int l,
                           unsigned no, int i1, int i2, int i3, int i4);
    // Index substitution
    ptr_base substidx(unsigned what, int with) const;
    // Lag
    int m_lag;
    // Indices
    unsigned m_noid;
    int m_idx1, m_idx2, m_idx3, m_idx4;

    friend void collect_idx(const ptr_base &p, std::set<unsigned> &iset);
    friend ptr_base append_name(const ptr_base &p, const std::string &s);
    friend ptr_base add_idx(const ptr_base &e, const idx_ex &ie);
    friend ptr_base subst(const ptr_base &e, const ptr_base &what,
                          const ptr_base &with, bool all_leads_lags);
    friend ptr_base substidx(const ptr_base &e, unsigned what, int with);
    friend ptr_base lag0(const ptr_base &e);

}; /* class ex_vartidx */


} /* namespace internal */
} /* namespace symbolic */

#endif /* SYMBOLIC_EX_VARTIDX_H */
