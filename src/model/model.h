/*****************************************************************************
 * This file is a part of gEcon.                                             *
 *                                                                           *
 * (c) Chancellery of the Prime Minister of the Republic of Poland 2012-2015 *
 * (c) Grzegorz Klima, Karol Podemski, Kaja Retkiewicz-Wijtiwiak 2015-2018   *
 * License terms can be found in the file 'LICENCE'                          *
 *                                                                           *
 * Author: Grzegorz Klima                                                    *
 *****************************************************************************/

/** \file model.h
 * \brief Class representing general equilibrium model.
 */

#ifndef MODEL_MODEL_H

#define MODEL_MODEL_H

#include <gecon_info.h>
#include <ex.h>
#include <model_block.h>


/// Class representing general equilibrium model.
class Model {
  public:
    /// Model options
    enum option {
        silent = 0,
        verbose,
        warnings_opt,
        output_logf,
        output_latex,
        output_latex_long,
        output_latex_landscape,
        output_r,
        output_r_long,
        output_r_jacobian,
        output_r_rcpp,
        OPTIONS_LENGTH
    };

    /// Default constructor
    Model();

    /// Set model name (extracts path from name).
    void set_name(const std::string&);

    /// Add index set.
    bool add_set(const idx_set &s);

    /// Is this set name?
    bool is_set(const std::string &s);
    /// Get index set by its name.
    symbolic::idx_set get_set(const std::string &s);

    /// Add variables to be symbolically reduced.
    void add_red_vars(const vec_exint vl) {
        m_v_redvars = vl;
    }

    /// Add block.
    void add_block(const std::string &s, int ln, idx_ex i1, idx_ex i2);

    /// Get current block name.
    const std::string& get_curr_block() const { return m_blocks.back().m_name; }

    /// Check if block with a given has already been declared.
    bool block_declared(const std::string &s) const {
        return m_names.find(s) != m_names.end();
    }

    /// Add definition.
    void add_definition(const ex &lhs, const ex &rhs, int l) {
        m_blocks.back().add_definition(lhs, rhs, l);
    }

    /// Add control variable.
    void add_controls(const vec_exintstr &cl) {
        m_blocks.back().add_controls(cl);
    }

    /// Add objective.
    void add_objective(const ex &obj, const ex &obj_eq, const ex &lambda, int l) {
        m_blocks.back().add_objective(obj, obj_eq, lambda, l);
    }

    /// Add constraint.
    void add_constraint(const ex &lhs, const ex &rhs, const ex &lambda, int l) {
        m_blocks.back().add_constraint(lhs, rhs, lambda, l);
    }

    /// Add reference to equations in other block to constraints.
    void add_constraint_ref(std::string &blname, int ref, int l) {
        m_blocks.back().add_constraint_ref(blname, ref, l);
    }

    /// Add identity.
    void add_identity(const ex &lhs, const ex &rhs, int l) {
        m_blocks.back().add_identity(lhs, rhs, l);
    }

    /// Add shock.
    void add_shock(const ex &s, int l) {
        m_blocks.back().add_shock(s, l);
    }
    /// Add shocks.
    void add_shocks(const vec_exint &pl) {
        m_blocks.back().add_shocks(pl);
    }

    /// Add calibrating equation.
    void add_calibr(const ex &lhs, const ex &rhs, int l, const vec_exint &pl = vec_exint()) {
        m_blocks.back().add_calibr(lhs, rhs, l, pl);
    }

    /// Set option
    void set_option(option o, bool fl = true) { m_options[o] = fl; ++m_options_set[o]; }

    /// Get option
    bool get_option(option o) { return m_options[o]; }

    /// Do it.
    void do_it();

    /// Clear.
    void clear();

    /// Warning
    void warning(const std::string &mes, int l);
    /// Error
    void error(const std::string &mes, int l);

    /// Where there any errors?
    bool errors() const { return m_err.size(); }

    /// Where there any warnings?
    bool warnings() const { return m_warn.size(); }

    /// Retrieve error messages.
    std::string get_errs(bool log = false) const;

    /// Retrieve warning messages.
    std::string get_warns(bool log = false) const;

    /// Check number of warning messages and set writing logfile if necessary.
    void check_warns();

    /// Produce output
    void write() const;

  private:

    // Warnings & errors
    std::vector<std::string> m_warn, m_err;

    // Options
    bool m_options[OPTIONS_LENGTH];
    int m_options_set[OPTIONS_LENGTH];

    // Model path and name
    std::string m_path, m_name;

    // Index set names and map
    std::set<std::string> m_set_names;
    std::map<std::string, symbolic::idx_set> m_sets;

    // Blocks and block names
    std::vector<Model_block> m_blocks;
    std::set<std::string> m_names;
    // Vars for reduction
    vec_exint m_v_redvars;
    set_ex m_s_redvars;

    // Is deterministic? Is static?
    bool m_deter, m_static;

    // Lagrange multipliers
    set_ex m_lagr_mult;
    map_ex_str m_lagr_mult_in;
    // Objective functions
    map_ex_str m_obj;
    // Controls
    map_ex_str m_contr;
    // Max and min lag
    int m_max_lag, m_min_lag;
    // Lagged variables
    set_ex m_lags;
    // Shocks
    set_ex m_shocks;
    // Variables
    set_ex m_vars;
    // Parameters
    set_ex m_params, m_params_calibr, m_params_free;
    map_ex_ex m_params_free_set;
    // Equations
    set_ex m_eqs, m_t_eqs;
    // Steady state equations
    vec_ex m_ss, m_t_ss;
    // Calibration eq's
    set_ex m_calibr;
    // Calibration eq's
    set_ex m_calibr_init;
    // Variables, parameters / equations, calibrating equations map
    std::map<std::pair<int, int>, unsigned> m_var_eq_map;
    std::set<std::pair<int, int> > m_var_ceq_map;
    std::set<std::pair<int, int> > m_cpar_eq_map;
    std::set<std::pair<int, int> > m_cpar_ceq_map;
    std::set<std::pair<int, int> > m_fpar_eq_map;
    std::set<std::pair<int, int> > m_fpar_ceq_map;
    // Shocks / equations map
    std::set<std::pair<int, int> > m_shock_eq_map;
    // Jacobian for ss and calibration equations
    std::map<std::pair<int, int>, ex> m_jacob_ss_calibr;
    // 1st order pertubation
    std::map<std::pair<int, int>, ex> m_Atm1, m_At, m_Atp1, m_Aeps;

    // Write gEcon model info message.
    void write_model_info(const std::string &mes) const;
    // Write gEcon info message.
    void write_info(const std::string &mes) const;
    // Stop on errors.
    void terminate_on_errors();

    /// Get option name
    static std::string get_option_name(int o);
    // Check options
    void check_options();
    // Check indices
    void check_indices();
    // Check indices
    void error_sindices(const std::set<unsigned> &is, const ex &e, int lineno);
    // Check indices
    void check_findices();
    // Check indices
    void error_findices(const std::map<unsigned, unsigned> &im, const ex &e, int lineno);
    // Check list of variables for reduction (1/2)
    void check_red_vars1();
    // Check definitions
    void check_defs();
    // Check names in blocks
    void check_names();
    // Substitute definitions
    void subst_defs();
    // Check if model is deterministic
    void check_deter();
    // Check Lagrange multipliers
    void check_lagr();
    // Check references to blocks
    void check_refs();
    // Check controls
    void check_obj_contr();
    // Check / handle leads > 1
    void leads();
    // Check / handle lags > 1
    void lags();
    // Check if model is static
    void check_static();
    // Derive FOCs
    void derive_focs();
    // Collect shocks
    void collect_shocks();
    // Merge vars and parameters in blocks
    void collect_vp();
    // Collect calibration
    void collect_calibr();
    // Merge equations
    void collect_eq();
    // Check list of variables for reduction (2/2)
    void check_red_vars2();
    // Reduce equations
    void reduce();
    // Construct variables / equations map
    void var_eq_map();
    // Construct variables / calibrating equations map
    void var_ceq_map();
    // Construct calibrated parameters / equations map
    void par_eq_map();
    // Construct calibrated parameters / calibrating equations map
    void par_ceq_map();
    // Construct shocks / equations map
    void shock_eq_map();
    // Steady state
    void stst();
    // Steady state and calibration eq's Jacobian
    void ss_jacob();
    // 1st order derivatives
    void diff_eqs();

    // Write logfile
    void write_log() const;
    // Generate R code and write it to file.
    void write_r() const;
    // Generate LaTeX documentation and write it to file.
    void write_latex() const;

};


#endif /* MODEL_MODEL_H */
