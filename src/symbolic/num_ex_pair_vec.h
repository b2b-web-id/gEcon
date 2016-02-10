/***********************************************************
 * (c) Kancelaria Prezesa Rady Ministrów 2012-2015         *
 * Treść licencji w pliku 'LICENCE'                        *
 *                                                         *
 * (c) Chancellery of the Prime Minister 2012-2015         *
 * License terms can be found in the file 'LICENCE'        *
 *                                                         *
 * Author: Grzegorz Klima                                  *
 ***********************************************************/

/** \file num_ex_pair_vec.h
 * \brief Vector of scalar-expression pairs used in construction
 * of classes representing addition and multiplication.
 */

#ifndef SYMBOLIC_NUM_EX_PAIR_VEC_H

#define SYMBOLIC_NUM_EX_PAIR_VEC_H

#include <number.h>
#include <ptr_base.h>
#include <vector>


namespace symbolic {
namespace internal {


typedef std::pair<Number, ptr_base> num_ex_pair;
typedef std::vector<num_ex_pair> num_ex_pair_vec_base;

/// Vector of scalar-expression pairs used in construction
/// of classes representing addition and multiplication.
class num_ex_pair_vec : public num_ex_pair_vec_base {
  public:
    /// Default constructor.
    num_ex_pair_vec() { ; }
    /// Constructor from ex_base ptr.
    num_ex_pair_vec(const ptr_base &p);
    /// Constructor from scalar and ex_base ptr.
    num_ex_pair_vec(const Number &s, const ptr_base &p);
    /// Constructor from scalar and num_ex_pair_vec (scaling).
    num_ex_pair_vec(const Number &s, const num_ex_pair_vec &p, ex_type t);
    /// Constructor from a pair of 2 scalar expression pairs
    num_ex_pair_vec(const num_ex_pair_vec &a,
                    const num_ex_pair_vec &b, ex_type t);
    /// Destructor
    ~num_ex_pair_vec() { ; }

    /// Factor out coefficient
    int compare(const num_ex_pair_vec&) const;

    /// Reduce
    void reduce(ex_type t);
    /// Max lag in expression
    int get_lag_max(bool stop_on_E = false) const;
    /// Min lag in expression
    int get_lag_min(bool stop_on_E = false) const;
    /// Does it have given subexpression?
    bool has(const ptr_base &what, search_flag f, bool exact_idx) const;
    /// Does it have given index?
    bool hasidx(int idx = 0) const;
    /// Does it have time indexed expressions?
    bool hast() const;

}; /* class num_ex_pair_vec */

} /* namespace internal */
} /* namespace symbolic */

#endif /* SYMBOLIC_NUM_EX_PAIR_VEC_H */
