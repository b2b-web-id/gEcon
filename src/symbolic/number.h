/***********************************************************
 * (c) Kancelaria Prezesa Rady Ministrów 2012-2015         *
 * Treść licencji w pliku 'LICENCE'                        *
 *                                                         *
 * (c) Chancellery of the Prime Minister 2012-2015         *
 * License terms can be found in the file 'LICENCE'        *
 *                                                         *
 * Author: Grzegorz Klima                                  *
 ***********************************************************/

/** \file number.h
 * \brief Numbers.
 */

#ifndef SYMBOLIC_NUMBER_H

#define SYMBOLIC_NUMBER_H


#include <iostream>
#include <cmath>
#include <utils.h>


namespace symbolic {

/// Number. Currently just a wrapper for double.
class Number {
  public:
    /// Default constructor (0).
    Number() : m_val(0.) { ; }
    /// Constructor (double).
    Number(const double &v) { m_val = v; }
    /// Constructor (int).
    Number(int v) { m_val = (double) v; }
    /// Constructor (unsigned).
    Number(unsigned v) { m_val = (double) v; }
    /// Destructor
    ~Number() { ; }

    /// Value
    double val() const { return m_val; }

    /// Operator bool
    bool operator()(const Number &n) { return m_val != 0.; }

    /// Operator +=
    Number& operator+=(const Number &n) { m_val += n.m_val; return *this; }
    /// Operator -=
    Number& operator-=(const Number &n) { m_val -= n.m_val; return *this; }
    /// Operator *=
    Number& operator*=(const Number &n) { m_val *= n.m_val; return *this; }
    /// Operator /=
    Number& operator/=(const Number &n) { m_val /= n.m_val; return *this; }

    /// String representation
    std::string str() const;
    /// LaTeX string representation
    std::string tex() const;

  private:
    // Value
    double m_val;

    friend std::ostream& operator<<(std::ostream&, const Number&);
    friend bool operator==(const Number&, const Number&);
    friend bool operator!=(const Number&, const Number&);
    friend bool operator>(const Number&, const Number&);
    friend bool operator>=(const Number&, const Number&);
    friend bool operator<(const Number&, const Number&);
    friend bool operator<=(const Number&, const Number&);
    friend Number operator-(const Number&);
    friend Number operator+(const Number&, const Number&);
    friend Number operator-(const Number&, const Number&);
    friend Number operator*(const Number&, const Number&);
    friend Number operator/(const Number&, const Number&);
    friend Number pow(const Number&, const Number&);
    friend Number exp(const Number&);
    friend Number log(const Number&);
    friend Number sin(const Number&);
    friend Number cos(const Number&);
    friend Number tan(const Number&);
    friend Number asin(const Number&);
    friend Number acos(const Number&);
    friend Number atan(const Number&);
    friend Number sinh(const Number&);
    friend Number cosh(const Number&);
    friend Number tanh(const Number&);
    // friend Number erf(const Number&);

}; /* class Number */


/// String representation
inline
std::string
Number::str() const
{
    return symbolic::internal::num2str(m_val);
}

/// LaTeX string representation
inline
std::string
Number::tex() const
{
    return symbolic::internal::num2tex(m_val);
}



/// Stream output
inline
std::ostream&
operator<<(std::ostream &os, const Number &n)
{
    if (std::isnan(n.m_val)) return (os << "NaN");
    if (std::isinf(n.m_val)) {
        if (n.m_val == -1./0.) os << '-';
        return (os << "Inf");
    }
    return (os << n.m_val);
}


/// Operator ==
inline
bool
operator==(const Number &a, const Number &b)
{
    return (a.m_val == b.m_val);
}

/// Operator !=
inline
bool
operator!=(const Number &a, const Number &b)
{
    return !(a == b);
}

/// Operator >
inline
bool
operator>(const Number &a, const Number &b)
{
    return (a.m_val > b.m_val);
}

/// Operator <=
inline
bool
operator<=(const Number &a, const Number &b)
{
    return !(a > b);
}

/// Operator <
inline
bool
operator<(const Number &a, const Number &b)
{
    return (a.m_val < b.m_val);
}

/// Operator >=
inline
bool
operator>=(const Number &a, const Number &b)
{
    return !(a > b);
}



/// Unary minus
inline
Number
operator-(const Number &n)
{
    return Number(-n.m_val);
}

/// Operator +
inline
Number
operator+(const Number &a, const Number &b)
{
    return Number(a.m_val + b.m_val);
}

/// Operator -
inline
Number
operator-(const Number &a, const Number &b)
{
    return Number(a.m_val - b.m_val);
}

/// Operator *
inline
Number
operator*(const Number &a, const Number &b)
{
    return Number(a.m_val * b.m_val);
}

/// Operator /
inline
Number
operator/(const Number &a, const Number &b)
{
    return Number(a.m_val / b.m_val);
}

/// Powers
inline
Number
pow(const Number &a, const Number &b)
{
    return Number(std::pow(a.m_val, b.m_val));
}

/// Exp
inline Number exp(const Number &n) { return Number(std::exp(n.m_val)); }

/// Log
inline Number log(const Number &n) { return Number(std::log(n.m_val)); }

/// Sine
inline Number sin(const Number &n) { return Number(std::sin(n.m_val)); }

/// Cosine
inline Number cos(const Number &n) { return Number(std::cos(n.m_val)); }

/// Tangent
inline Number tan(const Number &n) { return Number(std::tan(n.m_val)); }

/// Arc sine
inline Number asin(const Number &n) { return Number(std::asin(n.m_val)); }

/// Arc cosine
inline Number acos(const Number &n) { return Number(std::acos(n.m_val)); }

/// Arc tangent
inline Number atan(const Number &n) { return Number(std::atan(n.m_val)); }

/// Hyperbolic sine
inline Number sinh(const Number &n) { return Number(std::sinh(n.m_val)); }

/// Hyperbolic cosine
inline Number cosh(const Number &n) { return Number(std::cosh(n.m_val)); }

/// Hyperbolic tangent
inline Number tanh(const Number &n) { return Number(std::tanh(n.m_val)); }

// /// Erf
// inline Number erf(const Number &n) { return Number(std::erf(n.m_val)); }



} /* namespace symbolic */

#endif /* SYMBOLIC_NUMBER_H */
