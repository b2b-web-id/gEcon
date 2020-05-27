/*****************************************************************************
 * This file is a part of gEcon.                                             *
 *                                                                           *
 * (c) Chancellery of the Prime Minister of the Republic of Poland 2012-2015 *
 * (c) Grzegorz Klima, Karol Podemski, Kaja Retkiewicz-Wijtiwiak 2015-2018   *
 * License terms can be found in the file 'LICENCE'                          *
 *                                                                           *
 * Author: Grzegorz Klima                                                    *
 *****************************************************************************/

/** \file ex.h
 * \brief Expression.
 */

#ifndef SYMBOLIC_EX_H

#define SYMBOLIC_EX_H


#include <number.h>
#include <decl.h>
#include <ex_base.h>
#include <ptr_base.h>
#include <idx_ex.h>
#include <string>
#include <iostream>


namespace symbolic {


/// Expression class
class ex {

  public:
    /// Default constructor, sets ex to 0
    ex();
    /// Constructor from numbers
    explicit ex(const Number &v);
    /// Kronecker delta
    ex(bool fixed1, const std::string &i1, bool fixed2, const std::string &i2);
    /// Constructor from string (parameter)
    explicit ex(const std::string &n);
    /// Constructor from strings (indexed parameter)
    ex(const std::string &n, bool fixed1, const std::string &i1);
    /// Constructor from strings (indexed parameter)
    ex(const std::string &n, bool fixed1, const std::string &i1,
       bool fixed2, const std::string &i2);
    /// Constructor from strings (indexed parameter)
    ex(const std::string &n, bool fixed1, const std::string &i1,
       bool fixed2, const std::string &i2, bool fixed3, const std::string &i3);
    /// Constructor from strings (indexed parameter)
    ex(const std::string &n, bool fixed1, const std::string &i1,
       bool fixed2, const std::string &i2, bool fixed3, const std::string &i3,
       bool fixed4, const std::string &i4);
    /// Constructor from string and lag (variable)
    ex(const std::string &n, int l);
    /// Constructor from string, lag, and strings with indices
    ex(const std::string &n, int l, bool fixed1, const std::string &i1);
    /// Constructor from string, lag, and strings with indices
    ex(const std::string &n, int l, bool fixed1, const std::string &i1,
       bool fixed2, const std::string &i2);
    /// Constructor from string, lag, and strings with indices
    ex(const std::string &n, int l, bool fixed1, const std::string &i1,
       bool fixed2, const std::string &i2, bool fixed3, const std::string &i3);
    /// Constructor from string, lag, and strings with indices
    ex(const std::string &n, int l, bool fixed1, const std::string &i1,
       bool fixed2, const std::string &i2, bool fixed3, const std::string &i3,
       bool fixed4, const std::string &i4);
    /// Constructor of indexed expression
    ex(const idx_ex &ie, const ex &e);
    /// Constructor from ptr_base
    explicit ex(const internal::ptr_base &p) : m_ptr(p) { ; }
    /// Destructor
    ~ex() { ; }

    /// Assignment
    ex& operator=(const ex &rhs);

    /// Return new expression constructed by substitution.
    ex subst(const ex &what, const ex &with, bool all_leads_lags = true) const;
    /// Return new expression constructed by substitution.
    ex substidx(const std::string &what, const std::string &with, bool fix = true) const;
    /// Does expression contain a constant, a variable or and expected value?
    bool has(const ex &what, internal::search_flag f = internal::EXACT_T, bool exact_idx = true) const;
    /// Does expression contain conditional expectations?
    bool has_Es() const;

    /// String representation
    std::string str(int pflag = internal::DEFAULT, bool c_style = false) const;
    /// String representation using string 2 string map (name substitution).
    std::string strmap(const map_str_str&, bool c_style = false) const;
    /// LaTeX string representation
    std::string tex(int pflag = internal::DEFAULT) const;
    /// Max lag in expression
    int get_lag_max(bool stop_on_E = false) const;
    /// Min lag in expression
    int get_lag_min(bool stop_on_E = false) const;

    /// Comparison
    bool operator==(const ex&);
    /// Comparison
    bool operator!=(const ex&);

    /// Is not 0?
    operator bool() const { return !m_ptr->is0(); }
    /// Is this a valid number?
    bool validnum() const;
    /// Is 0?
    bool is0() const { return m_ptr->is0(); }
    /// Is 1?
    bool is1() const { return m_ptr->is1(); }
    /// Does expression contain any time indexed terms?
    bool hast() const { return m_ptr->flag() & symbolic::internal::HAST; }
    /// Is this a variable?
    bool is_var() const { return (m_ptr->type() == symbolic::internal::VART)
                            || (m_ptr->type() == symbolic::internal::VARTIDX); }

    /// Unary minus
    ex operator-() const;

    /// Retrieve const reference
    const internal::ptr_base& get_ptr_base() const { return m_ptr; }

#ifdef DEBUG_SYMBOLIC
    /// Information for debugging
    void info() const { m_ptr->info(); }
#endif /* DEBUG_SYMBOLIC */

    friend ex operator+(const ex &op1, const ex &op2);
    friend ex operator-(const ex &op1, const ex &op2);
    friend ex operator*(const ex &op1, const ex &op2);
    friend ex operator/(const ex &op1, const ex &op2);
    friend ex pow(const ex &op1, const ex &op2);
    friend ex func(symbolic::internal::func_code, const ex&);
    friend ex sum(const idx_ex &ie, const ex &e);
    friend ex sum(const ex &e);
    friend ex prod(const idx_ex &ie, const ex &e);
    friend ex prod(const ex &e);
    friend ex append_name(const ex &e, const std::string &s);
    friend ex add_idx(const ex &e, const idx_ex &ie);
    friend ex add_idx(const ex &e, const ex &ie);
    friend ex apply_idx(const ex &e, const ex &ie);
    friend ex lag(const ex &e, int l);
    friend ex lag0(const ex &e);
    friend ex ss(const ex &e);
    friend ex E(const ex &e, int l);
    friend ex diff(const ex &expression, const ex &variable);
    friend vec_ex expand(const ex &e);
    friend vec_ex expand(const vec_ex &ve);
    friend struct less_ex;
    friend triplet<bool, ex, ex> find_subst(const ex &expression,
                                            const set_ex&);
    friend triplet<bool, ex, ex> find_par_eq_num(const ex &expression);
    friend void find_Es(const ex&, set_ex&);
    friend ex drop_Es(const ex &e);
    friend void collect(const ex &e, set_ex &vars, set_ex &parms);
    friend void collect_lags(const ex &e, map_ex_int &map);
    friend void internal::collect_lags(const internal::ptr_base &p, map_ex_int &map);
    friend void collect_idx(const ex &e, std::set<unsigned> &iset);
    friend void collect_fidx(const ex &e, std::map<unsigned, unsigned> &iset);

private:
    /// Smart ptr to abstract base class ex_base
    internal::ptr_base m_ptr;

}; /* class ex */


/// Stream output
std::ostream& operator<<(std::ostream &s, const ex &e);


/// Operator +
ex operator+(const ex &op1, const ex &op2);
/// Operator +
inline ex operator+(const ex &op1, const Number &op2) { return op1 + ex(op2); }
/// Operator +
inline ex operator+(const ex &op1, const double &op2) { return op1 + ex(op2); }
/// Operator +
inline ex operator+(const ex &op1, const int &op2) { return op1 + ex(op2); }
/// Operator +
inline ex operator+(const Number &op1, const ex &op2) { return ex(op1) + op2; }
/// Operator +
inline ex operator+(const double &op1, const ex &op2) { return ex(op1) + op2; }
/// Operator +
inline ex operator+(const int &op1, const ex &op2) { return ex(op1) + op2; }

/// Operator -
ex operator-(const ex &op1, const ex &op2);
/// Operator -
inline ex operator-(const ex &op1, const Number &op2) { return op1 - ex(op2); }
/// Operator -
inline ex operator-(const ex &op1, const double &op2) { return op1 - ex(op2); }
/// Operator -
inline ex operator-(const ex &op1, const int &op2) { return op1 - ex(op2); }
/// Operator -
inline ex operator-(const Number &op1, const ex &op2) { return ex(op1) - op2; }
/// Operator -
inline ex operator-(const double &op1, const ex &op2) { return ex(op1) - op2; }
/// Operator -
inline ex operator-(const int &op1, const ex &op2) { return ex(op1) - op2; }

/// Operator *
ex operator*(const ex &op1, const ex &op2);
/// Operator *
inline ex operator*(const ex &op1, const Number &op2) { return op1 * ex(op2); }
/// Operator *
inline ex operator*(const ex &op1, const double &op2) { return op1 * ex(op2); }
/// Operator *
inline ex operator*(const ex &op1, const int &op2) { return op1 * ex(op2); }
/// Operator *
inline ex operator*(const Number &op1, const ex &op2) { return ex(op1) * op2; }
/// Operator *
inline ex operator*(const double &op1, const ex &op2) { return ex(op1) * op2; }
/// Operator *
inline ex operator*(const int &op1, const ex &op2) { return ex(op1) * op2; }

/// Operator /
ex operator/(const ex &op1, const ex &op2);
/// Operator /
inline ex operator/(const ex &op1, const Number &op2) { return op1 / ex(op2); }
/// Operator /
inline ex operator/(const ex &op1, const double &op2) { return op1 / ex(op2); }
/// Operator /
inline ex operator/(const ex &op1, const int &op2) { return op1 / ex(op2); }
/// Operator /
inline ex operator/(const Number &op1, const ex &op2) { return ex(op1) / op2; }
/// Operator /
inline ex operator/(const double &op1, const ex &op2) { return ex(op1) / op2; }
/// Operator /
inline ex operator/(const int &op1, const ex &op2) { return ex(op1) / op2; }

/// Power
ex pow(const ex &op1, const ex &op2);
/// Power
inline ex pow(const ex &op1, const Number &op2) { return pow(op1, ex(op2)); }
/// Power
inline ex pow(const ex &op1, const double &op2) { return pow(op1, ex(op2)); }
/// Power
inline ex pow(const ex &op1, const int &op2) { return pow(op1, ex(op2)); }
/// Power
inline ex pow(const Number &op1, const ex &op2) { return pow(ex(op1), op2); }
/// Power
inline ex pow(const double &op1, const ex &op2) { return pow(ex(op1), op2); }
/// Power
inline ex pow(const int &op1, const ex &op2) { return pow(ex(op1), op2); }

/// Function
ex func(symbolic::internal::func_code, const ex&);
/// Square root
inline ex sqrt(const ex &x) { return pow(x, ex(0.5)); }
/// Exponential function
inline ex exp(const ex &arg) { return func(symbolic::internal::EXP, arg); }
/// Logarithm
inline ex log(const ex &arg) { return func(symbolic::internal::LOG, arg); }
/// Sine
inline ex sin(const ex &arg) { return func(symbolic::internal::SIN, arg); }
/// Cosine
inline ex cos(const ex &arg) { return func(symbolic::internal::COS, arg); }
/// Tangent
inline ex tan(const ex &arg) { return func(symbolic::internal::TAN, arg); }
/// Arc sine
inline ex asin(const ex &arg) { return func(symbolic::internal::ASIN, arg); }
/// Arc cosine
inline ex acos(const ex &arg) { return func(symbolic::internal::ACOS, arg); }
/// Arc tangent
inline ex atan(const ex &arg) { return func(symbolic::internal::ATAN, arg); }
/// Hyperbolic sine
inline ex sinh(const ex &arg) { return func(symbolic::internal::SINH, arg); }
/// Hyperbolic cosine
inline ex cosh(const ex &arg) { return func(symbolic::internal::COSH, arg); }
/// Hyperbolic tangent
inline ex tanh(const ex &arg) { return func(symbolic::internal::TANH, arg); }
/// Error function
inline ex erf(const ex &arg) { return func(symbolic::internal::ERF, arg); }
/// Cumulative density function of the standard normal distribution
inline ex pnorm(const ex &arg) { return func(symbolic::internal::PNORM, arg); }

/// Sum over indices
ex sum(const idx_ex &ie, const ex &e);
/// Sum over indices (if argument is indexed).
ex sum(const ex &e);
/// Product over indices
ex prod(const idx_ex &ie, const ex &e);
/// Product over indices (if argument is indexed).
ex prod(const ex &e);

/// Differentiate
ex diff(const ex &expression, const ex &variable);

/// Given variable or parameter create new one by appending string to its name
ex append_name(const ex &e, const std::string &s);
/// Add indexing to a variable / parameter.
ex add_idx(const ex &e, const idx_ex &ie);
/// Add indexing to a variable / parameter.
ex add_idx(const ex &e, const ex &ie);
/// Apply indexing of one expression to another.
ex apply_idx(const ex &e, const ex &of);
/// Lag.
ex lag(const ex &e, int l);
/// Set lag of a variable to zero.
ex lag0(const ex &e);
/// Steady state
ex ss(const ex &e);
/// Conditional expectation
ex E(const ex &e, int l = 0);

/// Expand indexed expressions.
vec_ex expand(const ex &e);
/// Expand indexed expressions.
vec_ex expand(const vec_ex &ve);

/// Drop expectations
ex drop_Es(const ex &e);
/// Find expressions under expected value
void find_Es(const ex &e, set_ex&);
/// Collect variables and parameters.
void collect(const ex &e, set_ex &vars, set_ex &parms);
/// Collect variables and parameters.
void collect(const vec_ex &ve, set_ex &vars, set_ex &parms);
/// Collect variables in lag > 1
void collect_lags(const ex &e, map_ex_int &map);
/// Collect stray indices.
void collect_idx(const ex &e, std::set<unsigned> &iset);
/// Collect free indices.
void collect_fidx(const ex &e, std::map<unsigned, unsigned> &imap);
/// Collect free indices.
void collect_fidx(const idx_ex &ie, std::map<unsigned, unsigned> &imap);
/// Given an expression equal 0 try to find substitution for time indexed
/// vars with numbers / parameters
triplet<bool, ex, ex> find_subst(const ex &expression, const set_ex&);
/// Given an expression equal 0 try to find substitution for parameter with a number
triplet<bool, ex, ex> find_par_eq_num(const ex &expression);


} /* namespace symbolic */

#endif /* SYMBOLIC_EX_H */
