/***********************************************************
 * (c) Kancelaria Prezesa Rady Ministrów 2012-2015         *
 * Treść licencji w pliku 'LICENCE'                        *
 *                                                         *
 * (c) Chancellery of the Prime Minister 2012-2015         *
 * License terms can be found in the file 'LICENCE'        *
 *                                                         *
 * Author: Grzegorz Klima                                  *
 ***********************************************************/

/** \file ops.h
 * \brief Operations.
 */

#ifndef SYMBOLIC_OPS_H

#define SYMBOLIC_OPS_H

#include <decl.h>
#include <ptr_base.h>
#include <string>


namespace symbolic {
namespace internal {

/// Expression reduction.
ptr_base reduce(const ptr_base&);
/// Indexed expression reduction.
ptr_base reduceidx(const ptr_base&, const ::symbolic::idx_ex&);
/// Unary minus.
ptr_base mk_neg(const ptr_base &lhs);
/// Addition
ptr_base mk_add(const ptr_base &lhs, const ptr_base &rhs);
/// Addition
ptr_base mk_add(const num_ex_pair_vec &ops);
/// Subtraction.
ptr_base mk_sub(const ptr_base &lhs, const ptr_base &rhs);
/// Multiplication.
ptr_base mk_mul(const ptr_base &lhs, const ptr_base &rhs);
/// Multiplication.
ptr_base mk_mul(const num_ex_pair_vec &ops);
/// Division.
ptr_base mk_div(const ptr_base &lhs, const ptr_base &rhs);
/// Power.
ptr_base mk_pow(const ptr_base &lhs, const ptr_base &rhs);
/// Conditional expected value.
ptr_base mk_E(const ptr_base &p, int l);
/// Function.
ptr_base mk_func(func_code c, const ptr_base &arg);

/// Lag expression.
ptr_base lag(const ptr_base &p, int l);
/// Steady state.
ptr_base ss(const ptr_base &p);
/// Drop expectations
ptr_base drop_Es(const ptr_base &p);
/// Substitution.
ptr_base subst(const ptr_base &e, const ptr_base &what, const ptr_base &with,
               bool all_leads_lags = true);
/// Index substitution.
ptr_base substidx(const ptr_base &e, unsigned what, int with);
/// Given variable or paramter create new one by appending string to its name
ptr_base append_name(const ptr_base &p, const std::string &s);
/// Add indexing to a variable / parameter.
ptr_base add_idx(const ptr_base &e, const idx_ex &ie);
/// Add indexing to a variable / parameter.
ptr_base add_idx(const ptr_base &e, const ptr_base &ie);
/// Apply indexing of one expression to another.
ptr_base apply_idx(const ptr_base &e, const ptr_base &of);
/// Create sum from indexed expression.
ptr_base mk_sum(const ptr_base &p);
/// Create product from indexed expression.
ptr_base mk_prod(const ptr_base &p);

/// Expand sums and products.
ptr_base expandsp(const ptr_base &p);

/// Does expression contain conditional expectations?
bool has_Es(const ptr_base &e);
/// Find expressions under expected value
void find_Es(const ptr_base&, set_ex&);

/// Collect variables and parameters.
void collect(const ptr_base&, set_ex &vars, set_ex &parms);
/// Collect variables in lag > 1
void collect_lags(const ptr_base &p, map_ex_int &map);

/// Collect indices.
void collect_idx(const ptr_base &p, std::set<unsigned> &iset);
/// Collect free indices.
void collect_fidx(const ptr_base &p, std::map<unsigned, unsigned> &im);
/// Collect free indices.
void collect_fidx(const idx_ex &ie, std::map<unsigned, unsigned> &im);

/// Differentiate
ptr_base mk_diff(const ptr_base &e, const ptr_base &v);


} /* namespace internal */
} /* namespace symbolic */

#endif /* SYMBOLIC_OPS_H */
