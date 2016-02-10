/***********************************************************
 * (c) Kancelaria Prezesa Rady Ministrów 2012-2015         *
 * Treść licencji w pliku 'LICENCE'                        *
 *                                                         *
 * (c) Chancellery of the Prime Minister 2012-2015         *
 * License terms can be found in the file 'LICENCE'        *
 *                                                         *
 * Author: Grzegorz Klima                                  *
 ***********************************************************/

/** \file triplet.h
 * \brief Triplet class template.
 */

#ifndef SYMBOLIC_TRIPLET_H

#define SYMBOLIC_TRIPLET_H


namespace symbolic
{

/// Struct holding three objects of arbitrary types.
template<class T1, class T2, class T3>
struct triplet
{
    T1 first;
    T2 second;
    T3 third;

    /// Default constructor
    triplet() : first(), second(), third() { ; }

    /// Constructor
    triplet(const T1& a, const T2& b, const T3& c)
    : first(a), second(b), third(c) { ; }

    /// Copy constructor
    template<class U1, class U2, class U3>
    triplet(const triplet<U1, U2, U3> &p)
    : first(p.first), second(p.second), third(p.third)  { ; }

    /// Comparison
    bool operator==(const triplet<T1, T2, T3> &y)
    {
        return (first == y.first) && (second == y.second) && (third == y.third);
    }
    /// Comparison
    bool operator!=(const triplet<T1, T2, T3>& y)
    {
        return !(*this == y);
    }


}; /* struct triplet */


} /* namespace symbolic */

#endif /* SYMBOLIC_TRIPLET_H */
