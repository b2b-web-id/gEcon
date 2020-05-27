/*****************************************************************************
 * This file is a part of gEcon.                                             *
 *                                                                           *
 * (c) Chancellery of the Prime Minister of the Republic of Poland 2012-2015 *
 * (c) Grzegorz Klima, Karol Podemski, Kaja Retkiewicz-Wijtiwiak 2015-2018   *
 * License terms can be found in the file 'LICENCE'                          *
 *                                                                           *
 * Author: Grzegorz Klima                                                    *
 *****************************************************************************/

/** \file model_block.h
 * \brief Class representing a block of general equilibrium model.
 */

#ifndef MODEL_MODEL_BLOCK_H

#define MODEL_MODEL_BLOCK_H

#include <string>
#include <vector>
#include <ex.h>

using symbolic::ex;
using symbolic::idx_ex;
using symbolic::idx_set;
using symbolic::set_ex;
using symbolic::vec_ex;
using symbolic::map_ex_ex;
using symbolic::map_ex_int;
using symbolic::map_ex_str;
using symbolic::map_str_str;

typedef std::pair<ex, std::string> exstr;
typedef std::pair<ex, int> exint;
typedef std::vector<exint> vec_exint;
typedef std::pair<ex, ex> expair;
typedef std::vector<expair> vec_expair;
typedef std::pair<std::string, int> strint;
typedef std::vector<std::pair<int, int> > vec_intint;
typedef symbolic::triplet<std::string, int, int> strintint;
typedef std::vector<strintint> vec_strintint;
typedef symbolic::triplet<ex, int, std::string> exintstr;
typedef std::vector<exintstr> vec_exintstr;
typedef std::vector<strint> vec_strint;
typedef std::vector<idx_ex> vec_idx_ex;

/// Class representing model block (optimising agent or set of equilibrium conditions).
class Model_block {
  public:

    /// Reference types
    enum reference {
        objective,
        constraints,
        focs,
        identities
    };

    /// Constructor.
    explicit Model_block(const std::string &n, idx_ex i1 = idx_ex(), idx_ex i2 = idx_ex())
        : m_name(n), m_i1(i1), m_i2(i2), m_static(false) { ; }

    /// Add definition.
    void add_definition(const ex &lhs, const ex &rhs, int l);

    /// Add control variables.
    void add_controls(const vec_exintstr &cl);

    /// Add objective.
    void add_objective(const ex &obj, const ex &obj_eq, const ex &lambda, int l);

    /// Add constraint.
    void add_constraint(const ex &lhs, const ex &rhs, const ex &lambda, int l);

    /// Add reference to equations in other block to constraints.
    void add_constraint_ref(std::string &blname, int ref, int l);

    /// Add identity.
    void add_identity(const ex &lhs, const ex &rhs, int l);

    /// Add shock.
    void add_shock(const ex &s, int l);
    /// Add shocks.
    void add_shocks(const vec_exint &pl);

    /// Add calibrating equation.
    void add_calibr(const ex &lhs, const ex &rhs, int l, const vec_exint &pl);

    /// Substitute definitions
    void subst_defs();

    /// Handle lags > 1.
    void lags();

    /// FOCs
    void derive_focs();
    /// FOCs in deterministic model.
    void derive_focs_deter();
    /// FOCs in static problem.
    void derive_focs_static();

    /// Collect variables and parameters, if init == TRUE, collect definitions
    /// and ignore FOCs.
    void collect_vp(set_ex &vars, set_ex &parms, bool init = false);

    /// Generate LaTeX documentation and write it to ofstream.
    void write_latex(std::ostream &os, bool static_model = false) const;

    /// Write information to logfile.
    void write_logfile(std::ostream &os, bool static_model = false) const;

    /// Get name.
    std::string get_name(bool full = false) const;

  private:
    // Block name.
    std::string m_name;
    // Indices.
    idx_ex m_i1, m_i2;
    // is optimisation problem static?
    bool m_static;
    // Definitions.
    vec_exint m_defs_lhs, m_defs_rhs;
    set_ex m_defs;
    // Controls.
    vec_exintstr m_controls;
    set_ex m_contr, m_contr_exp;
    // Objective var.
    ex m_obj_var;
    // Objective eq.
    ex m_obj_eq, m_obj_eq_in;
    int m_obj_line;
    // Lagrange multiplier for objective.
    ex m_obj_lm, m_obj_lm_in;
    // Constraints.
    vec_exint m_constraints, m_constraints_in_lhs, m_constraints_in_rhs;
    vec_strintint m_constraints_ref;
    // Calibration equations.
    vec_exint m_calibr, m_calibr_in_lhs, m_calibr_in_rhs;
    // Params list corresponding to equations.
    std::vector<vec_exint> m_calibr_pl;
    // Lagrange multipliers.
    vec_exint m_lagr_mult, m_lagr_mult_in;
    set_ex m_redlm;
    // Lagged variables.
    set_ex m_lags;
    // Shocks.
    vec_exint m_shocks;
    // E's.
    set_ex m_Es;
    // q's (= E's).
    vec_ex m_qs;
    // Lagrange multipliers on q's.
    vec_ex m_etas, m_etas_v;
    // Identities.
    vec_exint m_identities, m_identities_in_lhs, m_identities_in_rhs;
    // FOCs.
    vec_expair m_focs, m_focs_red;

    // FOCs reduction.
    void reduce();
    // Does it have referenced variable (objective, control or LM)?
    bool has_ref(const ex &v) const;

    friend class Model;

};


#endif /* MODEL_MODEL_BLOCK_H */
