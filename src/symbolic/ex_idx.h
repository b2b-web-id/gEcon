/***********************************************************
 * (c) Kancelaria Prezesa Rady Ministrów 2012-2015         *
 * Treść licencji w pliku 'LICENCE'                        *
 *                                                         *
 * (c) Chancellery of the Prime Minister 2012-2015         *
 * License terms can be found in the file 'LICENCE'        *
 *                                                         *
 * Author: Grzegorz Klima                                  *
 ***********************************************************/

/** \file ex_idx.h
 * \brief Indexed expression.
 */

#ifndef SYMBOLIC_EX_IDX_H

#define SYMBOLIC_EX_IDX_H

#include <ex_base.h>
#include <ptr_base.h>
#include <ex.h>
#include <idx_ex.h>


namespace symbolic {
namespace internal {


/// Class representing expected value
class ex_idx : public ex_base {
  public:
    /// Constructor from indexing expression and expression
    ex_idx(const idx_ex &ie, const ptr_base &e)
        : ex_base(IDX | HASIDX | (e->flag() & HAST) | (e->flag() & SINGLE)), m_ie(ie), m_e(e)
        { ; }
    /// Destructor
    ~ex_idx() { ; }

    /// Constructor from from indexing expression and expression
    static ptr_base create(const idx_ex &ie, const ptr_base &e);
    /// Free memory (assumes that ptr is actually pointer to ex_idx)
    static void destroy(ex_base *ptr);

    /// Comparison
    int compare(const ex_idx&) const;

    /// String representation
    virtual std::string str(int pflag) const;
    /// String representation using string 2 string map (name substitution).
    /// This ignores expected vlau operators.
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

    /// Get underlying expression (skipping all indices)
    const ptr_base& get_ptr() const;
    /// Get expression
    const ptr_base& get_e() const { return m_e; }
    /// Get indexing expression
    const idx_ex& get_ie() const { return m_ie; }

  private:
    // No default constructor
    ex_idx();
    // Indexing expression
    idx_ex m_ie;
    // Expression
    ptr_base m_e;

    // Does it have an indexing expression?
    bool has_ie(const idx_ex &ie) const;

    friend ptr_base substidx(const ptr_base &e, unsigned what, int with);
    friend void collect_idx(const ptr_base &p, std::set<unsigned> &iset,
                            bool incl_internal);
    friend void collect_fidx(const ptr_base &p, std::map<unsigned, unsigned> &im);

}; /* class ex_idx */


} /* namespace internal */
} /* namespace symbolic */

#endif /* SYMBOLIC_EX_IDX_H */
