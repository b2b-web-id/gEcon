/*****************************************************************************
 * This file is a part of gEcon.                                             *
 *                                                                           *
 * (c) Chancellery of the Prime Minister of the Republic of Poland 2012-2015 *
 * (c) Grzegorz Klima, Karol Podemski, Kaja Retkiewicz-Wijtiwiak 2015-2018   *
 * License terms can be found in the file 'LICENCE'                          *
 *                                                                           *
 * Author: Grzegorz Klima                                                    *
 *****************************************************************************/

/** \file model.cpp
 * \brief Class representing general equilibrium model.
 */

#include <model.h>
#include <model_parse.h>
#include <utils.h>
#include <stdexcept>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <climits>

using symbolic::internal::num2str;
using symbolic::internal::print_flag;
using symbolic::internal::DEFAULT;
using symbolic::internal::DROP_INDEXING;
using symbolic::internal::DROP_T;
using symbolic::internal::EXACT_T;
using symbolic::internal::ANY_T;
using symbolic::internal::DIFF_T;
using symbolic::internal::LEAD_T;
using symbolic::triplet;

#define INTERNAL_ERROR throw(std::runtime_error(std::string("internal error in file ") +\
                             __FILE__ + ", line " + symbolic::internal::num2str(__LINE__)));
#ifdef DEBUG
#define DEBUG_INFO(x) std::cerr << "DEBUG INFO: " << (x) << '\n';
#else
#define DEBUG_INFO(x)
#endif


namespace {

std::string
num_name_str(int n, const std::string &name)
{
    std::string mes = num2str(n);
    mes += ' ';
    mes += name;
    if (n != 1) mes += "s";
    return mes;
}

} /* namespace */



void
Model::do_it()
{
    DEBUG_INFO("preliminary check")
    terminate_on_errors();

    // options
    DEBUG_INFO("checking options")
    check_options();
    if (m_options[verbose]) {
        std::string mes = "model has " + num_name_str(m_blocks.size(), "block");
        mes += ": " + m_blocks[0].m_name;
        for (unsigned i = 1; i < m_blocks.size(); ++i)
            mes += ", " + m_blocks[i].m_name;
        write_model_info(mes);
    }

    // indices
    DEBUG_INFO("checking indices")
    check_indices();
    check_findices();
    terminate_on_errors();

    // checks
    DEBUG_INFO("checking variables selected for reduction (1/2)")
    check_red_vars1();
    DEBUG_INFO("checking definitions")
    check_defs();
    terminate_on_errors();
    DEBUG_INFO("checking names in blocks before definition substitution")
    check_names();
    terminate_on_errors();
    DEBUG_INFO("substituting definitions")
    subst_defs();
    DEBUG_INFO("checking if model is deterministic")
    check_deter();
    DEBUG_INFO("checking Lagrange multipliers")
    check_lagr();
    DEBUG_INFO("checking references")
    check_refs();
    terminate_on_errors();
    DEBUG_INFO("checking objective functions and controls")
    check_obj_contr();
    terminate_on_errors();
    DEBUG_INFO("checking / handling leads")
    leads();
    terminate_on_errors();
    DEBUG_INFO("checking / handling lags")
    lags();
    DEBUG_INFO("checking if model is static")
    check_static();
    terminate_on_errors();
    if (m_options[verbose]) {
        if (m_static) {
            write_model_info("model is static");
        } else {
            if (m_deter) {
                write_model_info("model is dynamic, deterministic");
            } else {
                write_model_info("model is dynamic, stochastic");
            }
        }
    }

    // FOC derivation
    DEBUG_INFO("deriving FOCs")
    derive_focs();

    // collecting variables, shocks & parameters
    DEBUG_INFO("collecting shocks")
    collect_shocks();
    terminate_on_errors();
    DEBUG_INFO("collecting variables and parameters")
    collect_vp();
    DEBUG_INFO("collecting model equations")
    collect_eq();
    DEBUG_INFO("collecting calibration equations")
    collect_calibr();
    terminate_on_errors();
    if (m_options[verbose]) {
        write_model_info("model has " + num_name_str(m_eqs.size(), "equation")
                         + " with " + num_name_str(m_vars.size(), "variable"));
        write_model_info("model has "
                         + num_name_str(m_calibr.size(), "calibrating equation")
                         + " and "
                         + num_name_str(m_params_calibr.size(), "non-free (calibrated) parameter"));
    }

    // variable reduction
    DEBUG_INFO("reducing model equations")
    check_red_vars2();
    terminate_on_errors();
    reduce();
    terminate_on_errors();
    if (m_options[verbose]) {
        write_model_info("after reduction the model has "
                         + num_name_str(m_eqs.size(), "equation")
                         + " with " + num_name_str(m_vars.size(), "variable"));
    }

    // maps, steady state & derivatives
    DEBUG_INFO("constructing variables / equations equations map")
    var_eq_map();
    DEBUG_INFO("constructing shocks / equations map")
    shock_eq_map();
    terminate_on_errors();
    DEBUG_INFO("determining steady state equations")
    stst();
    terminate_on_errors();
    DEBUG_INFO("constructing variables / calibrating equations map")
    var_ceq_map();
    DEBUG_INFO("constructing parameter / equations, calibrating equations maps")
    par_eq_map();
    par_ceq_map();
    if (m_options[output_r_jacobian]) {
        DEBUG_INFO("determining steady state equations Jacobian")
        ss_jacob();
    }
    DEBUG_INFO("differentiating equations for the 1st order perturbation")
    diff_eqs();
}


// options
// =============================================================================
std::string
Model::get_option_name(int o)
{
    switch (o) {
        case Model::silent: return "silent";
        case Model::verbose: return "verbose";
        case Model::warnings_opt: return "warnings";
        case Model::output_logf: return "output logfile";
        case Model::output_latex: return "output LaTeX";
        case Model::output_latex_long: return "output LaTeX long";
        case Model::output_latex_landscape: return "output LaTeX landscape";
        case Model::output_r: return "output R";
        case Model::output_r_long: return "output R long";
        case Model::output_r_jacobian: return "output R Jacobian";
        case Model::output_r_rcpp: return "output R Rcpp";
        default:
            INTERNAL_ERROR
    }
}



void
Model::check_options()
{
#ifdef R_DLL
    if (!m_options[output_r]) {
        warning("ignoring option \"output R = false;\"", 0);
    }
#endif /* R_DLL */
    if (m_options[verbose] && m_options[silent]) {
        m_options[silent] = false;
        warning("option \"verbose = true;\" overrides \"silent = true;\"", 0);
    }
    
    if (m_options_set[output_latex_long] && !m_options[output_latex]) {
        warning("ignoring option \"output LaTeX long\" when LaTeX output is turned off", 0);
    }
    
    if (m_options[output_r_rcpp] && m_options[output_r_long]) {
        m_options[output_r_rcpp] = false;
        warning("ignoring option output R rcpp when R long option is turned on", 0);
    }
    
    if (m_options_set[output_latex_landscape] && !m_options[output_latex]) {
        warning("ignoring option \"output LaTeX landscape\" when LaTeX output is turned off", 0);
    }
    
    for (int i = 0; i < OPTIONS_LENGTH; ++i) {
        if (m_options_set[i] > 1) {
            warning("option \"" + get_option_name(i) + "\" set more than once; assuming the last setting ("
                    + (m_options[i] ? "true)" : "false)"), 0);
        }
    }
}


// indices
// =============================================================================
void
Model::error_sindices(const std::set<unsigned> &is, const ex &e, int lineno)
{
    unsigned n = is.size();
    std::string mes = "stray ";
    if (n == 1) mes += "index (";
    else mes += "indices (";
    const symbolic::internal::stringhash &ref =
        symbolic::internal::stringhash::get_instance();
    std::set<unsigned>::const_iterator it = is.begin();
    mes += '\"' +  ref.get_str(*it) + '\"';
    for (++it; it != is.end(); ++it) {
        mes += ", \"" +  ref.get_str(*it) + '\"';
    }
    mes += ") in expression \"" + e.str(DROP_INDEXING) + "\"";
    error(mes, lineno);
}



void
Model::check_indices()
{
    for (std::vector<exint>::const_iterator it = m_v_redvars.begin();
         it != m_v_redvars.end(); ++it) {
        ex e = it->first;
        int line = it->second;
        std::set<unsigned> iset;
        collect_idx(e, iset);
        if (iset.size()) error_sindices(iset, e, line);
        iset.clear();
    }
    for (unsigned i = 0, n = m_blocks.size(); i < n; ++i) {
        std::set<unsigned> iset;
        idx_ex i1 = m_blocks[i].m_i1;
        idx_ex i2 = m_blocks[i].m_i2;
        ex e;
        int line;
#ifdef DO_IT
#undef DO_IT
#endif
#define DO_IT \
e = ex(i1, ex(i2, e)); \
collect_idx(e, iset); \
if (iset.size()) error_sindices(iset, e, line); \
iset.clear();
        for (unsigned j = 0, J = m_blocks[i].m_defs_lhs.size(); j < J; ++j) {
            line = m_blocks[i].m_defs_lhs[j].second;
            e = m_blocks[i].m_defs_lhs[j].first;
            DO_IT
            e = m_blocks[i].m_defs_lhs[j].first;
            collect_idx(e, iset);
            if (m_blocks[i].m_i2) {
                if (iset.find(m_blocks[i].m_i2.get_id()) == iset.end())
                    error("lhs in definition (\""+ e.str() +"\") is missing "
                          + m_blocks[i].get_name(true) + "'s block index \""
                          + symbolic::internal::stringhash::get_instance().get_str(m_blocks[i].m_i2.get_id())
                          + "\"", line);
            }
            if (m_blocks[i].m_i1) {
                if (iset.find(m_blocks[i].m_i1.get_id()) == iset.end())
                    error("lhs in definition (\""+ e.str() +"\") is missing "
                          + m_blocks[i].get_name(true) + "'s block index \""
                          + symbolic::internal::stringhash::get_instance().get_str(m_blocks[i].m_i1.get_id())
                          + "\"", line);
            }
            iset.clear();
            line = m_blocks[i].m_defs_rhs[j].second;
            e = m_blocks[i].m_defs_rhs[j].first;
            DO_IT
        }
        for (unsigned j = 0, J = m_blocks[i].m_controls.size(); j < J; ++j) {
            line = m_blocks[i].m_controls[j].second;
            e = m_blocks[i].m_controls[j].first;
            DO_IT
            e = m_blocks[i].m_controls[j].first;
            collect_idx(e, iset);
            if (m_blocks[i].m_i2) {
                if (iset.find(m_blocks[i].m_i2.get_id()) == iset.end())
                    error("control variable (\""+ e.str() +"\") is missing "
                          + m_blocks[i].get_name(true) + "'s block index \""
                          + symbolic::internal::stringhash::get_instance().get_str(m_blocks[i].m_i2.get_id())
                          + "\"", line);
            }
            if (m_blocks[i].m_i1) {
                if (iset.find(m_blocks[i].m_i1.get_id()) == iset.end())
                    error("control variable (\""+ e.str() +"\") is missing "
                          + m_blocks[i].get_name(true) + "'s block index \""
                          + symbolic::internal::stringhash::get_instance().get_str(m_blocks[i].m_i1.get_id())
                          + "\"", line);
            }
            iset.clear();
        }
        line = m_blocks[i].m_obj_line;
        e = m_blocks[i].m_obj_var;
        DO_IT
        e = m_blocks[i].m_obj_var;
        if (e) {
            collect_idx(e, iset);
            if (m_blocks[i].m_i2) {
                if (iset.find(m_blocks[i].m_i2.get_id()) == iset.end())
                    error("objective variable (\""+ e.str() +"\") is missing "
                          + m_blocks[i].get_name(true) + "'s block index \""
                          + symbolic::internal::stringhash::get_instance().get_str(m_blocks[i].m_i2.get_id())
                          + "\"", line);
            }
            if (m_blocks[i].m_i1) {
                if (iset.find(m_blocks[i].m_i1.get_id()) == iset.end())
                    error("objective variable (\""+ e.str() +"\") is missing "
                          + m_blocks[i].get_name(true) + "'s block index \""
                          + symbolic::internal::stringhash::get_instance().get_str(m_blocks[i].m_i1.get_id())
                          + "\"", line);
            }
            iset.clear();
        }
        e = m_blocks[i].m_obj_eq_in;
        DO_IT
        e = m_blocks[i].m_obj_lm_in;
        DO_IT
        e = m_blocks[i].m_obj_lm_in;
        if (e) {
            collect_idx(e, iset);
            if (m_blocks[i].m_i2) {
                if (iset.find(m_blocks[i].m_i2.get_id()) == iset.end())
                    error("Lagrange multiplier on objective equation (time aggregator) (\""+ e.str() +"\") is missing "
                          + m_blocks[i].get_name(true) + "'s block index \""
                          + symbolic::internal::stringhash::get_instance().get_str(m_blocks[i].m_i2.get_id())
                          + "\"", line);
            }
            if (m_blocks[i].m_i1) {
                if (iset.find(m_blocks[i].m_i1.get_id()) == iset.end())
                    error("Lagrange multiplier on objective equation (time aggregator) (\""+ e.str() +"\") is missing "
                          + m_blocks[i].get_name(true) + "'s block index \""
                          + symbolic::internal::stringhash::get_instance().get_str(m_blocks[i].m_i1.get_id())
                          + "\"", line);
            }
            iset.clear();
        }
        for (unsigned j = 0, J = m_blocks[i].m_constraints.size(); j < J; ++j) {
            line = m_blocks[i].m_constraints_in_lhs[j].second;
            e = m_blocks[i].m_constraints_in_lhs[j].first;
            DO_IT
            line = m_blocks[i].m_constraints_in_rhs[j].second;
            e = m_blocks[i].m_constraints_in_rhs[j].first;
            DO_IT
            line = m_blocks[i].m_lagr_mult_in[j].second;
            e = m_blocks[i].m_lagr_mult_in[j].first;
            DO_IT
            e = m_blocks[i].m_lagr_mult_in[j].first;
            if (e) {
                collect_idx(e, iset);
                if (m_blocks[i].m_i2) {
                    if (iset.find(m_blocks[i].m_i2.get_id()) == iset.end())
                        error("Lagrange multiplier (\""+ e.str() +"\") is missing "
                              + m_blocks[i].get_name(true) + "'s block index \""
                              + symbolic::internal::stringhash::get_instance().get_str(m_blocks[i].m_i2.get_id())
                              + "\"", line);
                }
                if (m_blocks[i].m_i1) {
                    if (iset.find(m_blocks[i].m_i1.get_id()) == iset.end())
                        error("Lagrange multiplier (\""+ e.str() +"\") is missing "
                              + m_blocks[i].get_name(true) + "'s block index \""
                              + symbolic::internal::stringhash::get_instance().get_str(m_blocks[i].m_i1.get_id())
                              + "\"", line);
                }
                iset.clear();
            }
        }
        for (unsigned j = 0, J = m_blocks[i].m_identities.size(); j < J; ++j) {
            line = m_blocks[i].m_identities_in_lhs[j].second;
            e = m_blocks[i].m_identities_in_lhs[j].first;
            DO_IT
            line = m_blocks[i].m_identities_in_rhs[j].second;
            e = m_blocks[i].m_identities_in_rhs[j].first;
            DO_IT
        }
        for (unsigned j = 0, J = m_blocks[i].m_shocks.size(); j < J; ++j) {
            line = m_blocks[i].m_shocks[j].second;
            e = m_blocks[i].m_shocks[j].first;
            DO_IT
        }
        for (unsigned j = 0, J = m_blocks[i].m_calibr.size(); j < J; ++j) {
            line = m_blocks[i].m_calibr_in_lhs[j].second;
            e = m_blocks[i].m_calibr_in_lhs[j].first;
            DO_IT
            line = m_blocks[i].m_calibr_in_rhs[j].second;
            e = m_blocks[i].m_calibr_in_rhs[j].first;
            DO_IT
            for (unsigned k = 0; k < m_blocks[i].m_calibr_pl[j].size(); ++k) {
                line = m_blocks[i].m_calibr_pl[j][k].second;
                e = m_blocks[i].m_calibr_pl[j][k].first;
                DO_IT
            }
        }
    }
}



void
Model::error_findices(const std::map<unsigned, unsigned> &im, const ex &e, int lineno)
{
    unsigned n = im.size();
    std::string mes = "duplicated ";
    if (n == 1) mes += "index in nested indexing expression (";
    else mes += "indices in nested indexing expression (";
    const symbolic::internal::stringhash &ref =
        symbolic::internal::stringhash::get_instance();
    std::map<unsigned, unsigned>::const_iterator it = im.begin();
    mes += '\"' +  ref.get_str(it->first) + '\"';
    for (++it; it != im.end(); ++it) {
        mes += ", \"" +  ref.get_str(it->first) + '\"';
    }
    mes += ") in expression \"" + e.str() + "\"";
    error(mes, lineno);
}



void
Model::check_findices()
{
    for (std::vector<exint>::const_iterator iit = m_v_redvars.begin();
         iit != m_v_redvars.end(); ++iit) {
        ex e = iit->first;
        int line = iit->second;
        std::map<unsigned, unsigned> imap;
        std::map<unsigned, unsigned>::iterator it;
        collect_fidx(e, imap);
        for (it = imap.begin(); it != imap.end(); ++it) {
            if (!it->second) imap.erase(it);
        }
        if (imap.size()) error_findices(imap, e, line);
        imap.clear();
    }
    for (unsigned i = 0, n = m_blocks.size(); i < n; ++i) {
        std::map<unsigned, unsigned> imap;
        std::map<unsigned, unsigned>::iterator it;
        ex e;
        int line;
#ifdef DO_IT
#undef DO_IT
#endif
#define DO_IT \
collect_fidx(m_blocks[i].m_i1, imap); \
collect_fidx(m_blocks[i].m_i2, imap); \
collect_fidx(e, imap); \
for (it = imap.begin(); it != imap.end(); ++it) { \
    if (!it->second) imap.erase(it); \
} \
if (imap.size()) error_findices(imap, ex(m_blocks[i].m_i1, ex(m_blocks[i].m_i2, e)), line); \
imap.clear();
        for (unsigned j = 0, J = m_blocks[i].m_defs_lhs.size(); j < J; ++j) {
            line = m_blocks[i].m_defs_lhs[j].second;
            e = m_blocks[i].m_defs_lhs[j].first;
            DO_IT
            line = m_blocks[i].m_defs_rhs[j].second;
            e = m_blocks[i].m_defs_rhs[j].first;
            DO_IT
        }
        for (unsigned j = 0, J = m_blocks[i].m_controls.size(); j < J; ++j) {
            line = m_blocks[i].m_controls[j].second;
            e = m_blocks[i].m_controls[j].first;
            DO_IT
        }
        line = m_blocks[i].m_obj_line;
        e = m_blocks[i].m_obj_var;
        DO_IT
        e = m_blocks[i].m_obj_eq_in;
        DO_IT
        e = m_blocks[i].m_obj_lm_in;
        DO_IT
        for (unsigned j = 0, J = m_blocks[i].m_constraints.size(); j < J; ++j) {
            line = m_blocks[i].m_constraints_in_lhs[j].second;
            e = m_blocks[i].m_constraints_in_lhs[j].first;
            DO_IT
            line = m_blocks[i].m_constraints_in_rhs[j].second;
            e = m_blocks[i].m_constraints_in_rhs[j].first;
            DO_IT
            line = m_blocks[i].m_lagr_mult_in[j].second;
            e = m_blocks[i].m_lagr_mult_in[j].first;
            DO_IT
        }
        for (unsigned j = 0, J = m_blocks[i].m_identities.size(); j < J; ++j) {
            line = m_blocks[i].m_identities_in_lhs[j].second;
            e = m_blocks[i].m_identities_in_lhs[j].first;
            DO_IT
            line = m_blocks[i].m_identities_in_rhs[j].second;
            e = m_blocks[i].m_identities_in_rhs[j].first;
            DO_IT
        }
        for (unsigned j = 0, J = m_blocks[i].m_shocks.size(); j < J; ++j) {
            line = m_blocks[i].m_shocks[j].second;
            e = m_blocks[i].m_shocks[j].first;
            DO_IT
        }
        for (unsigned j = 0, J = m_blocks[i].m_calibr.size(); j < J; ++j) {
            line = m_blocks[i].m_calibr_in_lhs[j].second;
            e = m_blocks[i].m_calibr_in_lhs[j].first;
            DO_IT
            line = m_blocks[i].m_calibr_in_rhs[j].second;
            e = m_blocks[i].m_calibr_in_rhs[j].first;
            DO_IT
            for (unsigned k = 0; k < m_blocks[i].m_calibr_pl[j].size(); ++k) {
                line = m_blocks[i].m_calibr_pl[j][k].second;
                e = m_blocks[i].m_calibr_pl[j][k].first;
                DO_IT
            }
        }
    }
}



// checks
// =============================================================================
void
Model::check_red_vars1()
{
    std::vector<exint>::const_iterator it = m_v_redvars.begin();
    for (; it != m_v_redvars.end(); ++it) {
        ex e = it->first;
        int line = it->second;
        if (e.get_lag_max()) {
            error("variable(s) \"" + lag0(e).str()
                  + "\" selected for reduction listed in lead/lag", line);
            continue;
        }
        vec_ex rv = expand(e);
        vec_ex::const_iterator iit = rv.begin();
        for (; iit != rv.end(); ++iit) {
            if (m_s_redvars.find(*iit) != m_s_redvars.end()) {
                warning("variable \"" + e.str() + "\" already selected for reduction", line);
            } else {
                m_s_redvars.insert(*iit);
            }
        }
    }
}


void
Model::check_defs()
{
    for (unsigned i = 0, n = m_blocks.size(); i < n; ++i) {
        idx_ex i1 = m_blocks[i].m_i1;
        idx_ex i2 = m_blocks[i].m_i2;
        for (unsigned d = 0, D = m_blocks[i].m_defs_lhs.size(); d < D; ++d) {
            ex lhs = m_blocks[i].m_defs_lhs[d].first;
            ex rhs = m_blocks[i].m_defs_rhs[d].first;
            if (lhs.is_var()) {
                if (lhs.get_lag_max()) {
                    error("\"" + lag0(lhs).str() + "\" defined in lead/lag",
                          m_blocks[i].m_defs_lhs[d].second);
                } else if (!rhs.hast()) {
                    error("variable \"" + lhs.str() + "\" defined as constant expression \""
                          + rhs.str() + "\"", m_blocks[i].m_defs_lhs[d].second);
                } else { // ok
                    vec_ex vdef = expand(ex(i1, ex(i2, lhs)));
                    vec_ex::const_iterator iit = vdef.begin();
                    for (; iit != vdef.end(); ++iit) {
                        if (m_s_redvars.find(*iit) != m_s_redvars.end()) {
                            warning("variable \"" + iit->str()
                                    + "\" selected for reduction appears on the LHS of definition in "
                                    + m_blocks[i].m_name + " block", m_blocks[i].m_defs_lhs[d].second);
                        }
                    }
                }
            } else {
                if (rhs.hast()) {
                    error("parameter \"" + lhs.str() + "\" defined as variable expression \""
                          + rhs.str() + "\"", m_blocks[i].m_defs_lhs[d].second);
                }
            }
            if (!m_blocks[i].m_defs.insert(lhs).second) {
                error("\"" + lhs.str() + "\" already defined",
                      m_blocks[i].m_defs_lhs[d].second);
            }
            if (rhs.has(lhs, ANY_T, false)) {
                error("\"" + lhs.str() + "\" used in its own definition",
                      m_blocks[i].m_defs_lhs[d].second);
            }
            for (unsigned dd = d + 1; dd < D; ++dd) {
                ex rhsf = m_blocks[i].m_defs_rhs[dd].first;
                if (rhsf.has(lhs, ANY_T, false)) {
                    error("\"" + lhs.str() + "\" appears in definition following "
                          "the definition of \"" + lhs.str() + "\" itself",
                          m_blocks[i].m_defs_lhs[dd].second);
                }
            }
        }
    }
}



void
Model::check_names()
{
    // NOTE: this check is within each block before substitution,
    // this is checked again after all variables and parameters are collected
    for (unsigned i = 0, n = m_blocks.size(); i < n; ++i) {
        set_ex vars, params;
        m_blocks[i].collect_vp(vars, params, true);
        set_ex::const_iterator it;
        std::set<std::string> names;
        for (it = vars.begin(); it != vars.end(); ++it) {
            std::string name = it->str(DROP_T);
            names.insert(name);
        }
        for (it = params.begin(); it != params.end(); ++it) {
            std::string name = it->str();
            if (!names.insert(name).second) {
                error("\"" + name + "\" treated both as a variable and a parameter in block "
                      + m_blocks[i].m_name + "; you might have forgotten \"[]\"", 0);
            }
        }
    }
}



void
Model::subst_defs()
{
    for (unsigned i = 0, n = m_blocks.size(); i < n; ++i) {
        m_blocks[i].subst_defs();
    }
}



void
Model::check_deter()
{
    m_deter = true;
    for (unsigned i = 0, n = m_blocks.size(); i < n; ++i) {
        if (m_blocks[i].m_shocks.size()) {
            m_deter = false;
            return;
        }
    }
    for (unsigned i = 0, n = m_blocks.size(); i < n; ++i) {
        ex eq = m_blocks[i].m_obj_eq;
        int l = m_blocks[i].m_obj_line;
        if (eq.has_Es()) {
            warning("dropping expectation in objective \""
                    + m_blocks[i].m_obj_var.str() + " = " + eq.str()
                    + "\"; model is deterministic (has no shocks)", l);
            m_blocks[i].m_obj_eq = drop_Es(m_blocks[i].m_obj_eq);
            m_blocks[i].m_obj_eq_in = drop_Es(m_blocks[i].m_obj_eq_in);
        }
        for (unsigned e = 0; e < m_blocks[i].m_constraints.size(); ++e) {
            eq = m_blocks[i].m_constraints[e].first;
            l = m_blocks[i].m_constraints[e].second;
            if (eq.has_Es()) {
                warning("dropping expectation in constraint \""
                        + eq.str() + " = 0\"; model is deterministic (has no shocks)", l);
                m_blocks[i].m_constraints[e].first =
                    drop_Es(m_blocks[i].m_constraints[e].first);
            }
        }
        for (unsigned e = 0; e < m_blocks[i].m_identities.size(); ++e) {
            eq = m_blocks[i].m_identities[e].first;
            l = m_blocks[i].m_identities[e].second;
            if (eq.has_Es()) {
                warning("dropping expectation in identity \""
                        + eq.str() + " = 0\"; model is deterministic (has no shocks)", l);
                m_blocks[i].m_identities[e].first =
                    drop_Es(m_blocks[i].m_identities[e].first);
            }
        }
    }
}



void
Model::check_lagr()
{
    for (unsigned i = 0, n = m_blocks.size(); i < n; ++i) {
        ex lmin = m_blocks[i].m_obj_lm_in;
        idx_ex i1 = m_blocks[i].m_i1;
        idx_ex i2 = m_blocks[i].m_i2;
        if (lmin) {
            if ((lmin.get_lag_max())) {
                error("Lagrange multiplier  \"" + lag0(lmin).str() + "\" declared in lead/lag",
                      m_blocks[i].m_obj_line);
            }
            vec_ex vlmin = expand(ex(i1, ex(i2, lmin)));
            vec_ex::const_iterator iit = vlmin.begin();
            for (; iit != vlmin.end(); ++iit) {
                if (!m_lagr_mult_in.insert(exstr(*iit, m_blocks[i].m_name)).second) {
                    error("Lagrange multiplier \"" + iit->str() + "\" already used",
                          m_blocks[i].m_obj_line);
                }
            }
        }
        unsigned ll = 0, LL = m_blocks[i].m_lagr_mult_in.size();
        for (; ll < LL; ++ll) {
            lmin = m_blocks[i].m_lagr_mult_in[ll].first;
            if (lmin) {
                if ((lmin.get_lag_max())) {
                    error("Lagrange multiplier  \"" + lag0(lmin).str() + "\" declared in lead/lag",
                          m_blocks[i].m_lagr_mult_in[ll].second);
                }
                vec_ex vlmin = expand(ex(i1, ex(i2, lmin)));
                vec_ex::const_iterator iit = vlmin.begin();
                for (; iit != vlmin.end(); ++iit) {
                    if (!m_lagr_mult_in.insert(exstr(*iit, m_blocks[i].m_name)).second) {
                        error("Lagrange multiplier \"" + iit->str() + "\" already used",
                              m_blocks[i].m_lagr_mult_in[ll].second);
                    }
                }
            }
        }
    }
}



void
Model::check_refs()
{
    for (unsigned i = 0, n = m_blocks.size(); i < n; ++i) {
        std::set<std::pair<std::string, int> > refs;
        for (unsigned c = 0; c < m_blocks[i].m_controls.size(); ++c) {
            if (m_blocks[i].m_controls[c].third == "") continue;
            unsigned b = 0;
            std::string bname = m_blocks[i].m_controls[c].third;
            while ((m_blocks[b].m_name != bname) && (b + 1 < n)) ++b;
            if (m_blocks[b].m_name != bname) ++b;
            if (b == n) {
                error("block \"" + bname + "\" does not exis", m_blocks[i].m_controls[c].second);
                continue;
            } else if (b > i) {
                error("reference to block \"" + bname + "\""
                      + " that is declared after \"" + m_blocks[i].m_name + "\"",
                      m_blocks[i].m_controls[c].second);
                continue;
            } else if (b == i) {
                error("reference to block \"" + bname + "\" within itself",
                      m_blocks[i].m_controls[c].second);
                continue;
            }
        }

        for (unsigned r = 0; r < m_blocks[i].m_constraints_ref.size(); ++r) {
            strintint rf = m_blocks[i].m_constraints_ref[r];
            unsigned b = 0;
            while ((m_blocks[b].m_name != rf.first) && (b + 1 < n)) ++b;
            if (m_blocks[b].m_name != rf.first) ++b;
            if (b == n) {
                error("block \"" + rf.first + "\" does not exist", rf.third);
                continue;
            } else if (b > i) {
                error("reference to block \"" + rf.first + "\""
                      + " that is declared after \"" + m_blocks[i].m_name + "\"",
                      rf.third);
                continue;
            } else if (b == i) {
                error("reference to block \"" + rf.first + "\" within itself",
                      rf.third);
                continue;
            }
            std::string what;
            if (rf.second == Model_block::objective) {
                what = "objective";
                if (!m_blocks[b].m_obj_var) {
                    error("undefined reference to " + what + " in block \""
                          + rf.first + "\"", rf.third);
                }
            } else if (rf.second == Model_block::constraints) {
                what = "constraints";
                if (!m_blocks[b].m_constraints.size()
                    && !m_blocks[b].m_constraints_ref.size())
                    error("undefined reference to " + what + " in block \""
                          + rf.first + "\"", rf.third);
            } else if (rf.second == Model_block::focs) {
                what = "first order conditions";
                if (!m_blocks[b].m_controls.size())
                    error("undefined reference to " + what + " in block \""
                          + rf.first + "\"", rf.third);
            } else if (rf.second == Model_block::identities) {
                what = "identities";
                if (!m_blocks[b].m_identities.size())
                    error("undefined reference to " + what + " in block \""
                          + rf.first + "\"", rf.third);
            } else INTERNAL_ERROR
                if (!refs.insert(std::pair<std::string, int>(rf.first, rf.second)).second) {
                    error("repeated reference to " + what + " in block \""
                          + rf.first + "\"", rf.third);
                }
        }
    }
}



void
Model::check_obj_contr()
{
    for (unsigned i = 0, n = m_blocks.size(); i < n; ++i) {
        idx_ex i1 = m_blocks[i].m_i1;
        idx_ex i2 = m_blocks[i].m_i2;
        ex obj = m_blocks[i].m_obj_var;
        if (obj) {
            if (obj.get_lag_max()) {
                error("objective variable \"" + lag0(obj).str() + "\" on LHS of "
                      + m_blocks[i].m_name + "\'s problem is in lag different than 0",
                      m_blocks[i].m_obj_line);
                obj = lag0(obj);
            }
            if (m_blocks[i].m_obj_eq.has(lag(obj, 1), DIFF_T)) {
                error("objective variable \"" + obj.str() + "\" on RHS of " + m_blocks[i].m_name
                      + "\'s problem is in lead/lag different than 1",
                      m_blocks[i].m_obj_line);
            }
            if (!m_blocks[i].m_obj_eq.has(lag(obj, 1), EXACT_T)) {
                m_blocks[i].m_static = true;
                if (m_blocks[i].m_obj_lm_in) {
                    error("Lagrange multiplier on objective (time aggregator) declared in static problem",
                          m_blocks[i].m_obj_line);
                }
            }
            for (unsigned j = 0, J = m_blocks[i].m_constraints.size(); j < J; ++j) {
                if (m_blocks[i].m_constraints[j].first.has(lag(obj, 1), DIFF_T)) {
                    error("objective variable \"" + obj.str() + "\" in " + m_blocks[i].m_name
                          + "\'s problem's constraint is in lead/lag different than 1",
                          m_blocks[i].m_constraints[j].second);
                }
            }
            vec_ex vobj = expand(ex(i1, ex(i2, obj)));
            vec_ex::const_iterator iit = vobj.begin();
            for (; iit != vobj.end(); ++iit) {
                if (!m_obj.insert(exstr(*iit, m_blocks[i].m_name)).second) {
                    error("objective variable \"" + iit->str() + "\" in " + m_blocks[i].m_name
                          + "\'s problem repeated as objective variable in " + m_obj.find(*iit)->second
                          + "\'s problem", m_blocks[i].m_obj_line);
                }
                if (m_lagr_mult_in.find(*iit) != m_lagr_mult_in.end()) {
                    error("objective variable \"" + iit->str() + "\" in " + m_blocks[i].m_name
                          + "\'s problem is Lagrange multiplier in " + m_lagr_mult_in.find(*iit)->second
                          + "\'s problem", m_blocks[i].m_obj_line);
                }
            }
        }
    }

    for (unsigned i = 0, n = m_blocks.size(); i < n; ++i) {
        for (unsigned c = 0, C = m_blocks[i].m_controls.size(); c < C; ++c) {
            ex contr = m_blocks[i].m_controls[c].first;
            bool any = false;
            if (contr.get_lag_max()) {
                contr = lag0(contr);
                error("control variable \"" + contr.str() + "\" in " + m_blocks[i].m_name
                      + "\'s problem is listed in lag different than 0",
                      m_blocks[i].m_controls[c].second);
            }
            if (m_blocks[i].m_defs.find(contr) != m_blocks[i].m_defs.end()) {
                error("control variable \"" + contr.str() + "\" in "
                      + m_blocks[i].m_name + "\'s problem was used in definitions "
                      + "(substituted)", m_blocks[i].m_controls[c].second);
            }
            if (m_blocks[i].m_obj_eq.has(contr, ANY_T, false)) any = true;
            if (m_blocks[i].m_obj_eq.has(contr, LEAD_T, false)) {
                error("control variable \"" + contr.str() + "\" in " + m_blocks[i].m_name
                      + "\'s problem is in lead", m_blocks[i].m_obj_line);
            }
            for (unsigned j = 0, J = m_blocks[i].m_constraints.size(); j < J; ++j) {
                if (m_blocks[i].m_constraints[j].first.has(contr, ANY_T, false)) any = true;
                if (m_blocks[i].m_constraints[j].first.has(contr, LEAD_T, false)) {
                    error("control variable \"" + contr.str() + "\" in " + m_blocks[i].m_name
                          + "\'s problem is in lead", m_blocks[i].m_constraints[j].second);
                }
            }
            int ln = m_blocks[i].m_controls[c].second;
            std::string ref = m_blocks[i].m_controls[c].third;
            if (!any) {
                for (unsigned r = 0; r < m_blocks[i].m_constraints_ref.size(); ++r) {
                    strintint rf = m_blocks[i].m_constraints_ref[r];
                    unsigned b = 0;
                    while (m_blocks[b].m_name != rf.first) ++b;
                    if ((rf.second == Model_block::objective)
                        || (rf.second == Model_block::focs)) {
                        if (m_blocks[b].m_obj_eq.has(contr, ANY_T, false)) any = true;
                    }
                    if ((rf.second == Model_block::constraints)
                        || (rf.second == Model_block::focs)) {
                        for (unsigned j = 0, J = m_blocks[b].m_constraints.size(); j < J; ++j) {
                            if (m_blocks[b].m_constraints[j].first.has(contr, ANY_T, false))
                                any = true;
                        }
                    }
                    if (rf.second == Model_block::identities) {
                        for (unsigned j = 0, J = m_blocks[b].m_constraints.size(); j < J; ++j) {
                            if (m_blocks[b].m_identities[j].first.has(contr, ANY_T, false))
                                any = true;
                        }
                    }
                    if (any) break;
                }
            }
            if (!any) {
                error("control variable \"" + contr.str() + "\" in " + m_blocks[i].m_name
                      + "\'s problem does not appear in objective nor in any constraint", ln);
            }
            if (!ref.empty()) {
                unsigned b = 0;
                while (m_blocks[b].m_name != ref) ++b;
                if (!m_blocks[b].has_ref(contr))
                    error("variable \"" + contr.str() + "\" referenced in " + m_blocks[i].m_name
                          + "\'s problem does not appear in " + ref
                          + "\'s problem; only objective variables, control variables "
                          + "and Langrange multipliers should be referenced", ln);
            }
            idx_ex i1 = m_blocks[i].m_i1;
            idx_ex i2 = m_blocks[i].m_i2;
            m_blocks[i].m_contr.insert(contr);
            vec_ex vcontr = expand(ex(i1, ex(i2, contr)));
            vec_ex::const_iterator iit = vcontr.begin();
            for (; iit != vcontr.end(); ++iit) {
                if (!m_blocks[i].m_contr_exp.insert(*iit).second) {
                    error("control variable \"" + iit->str() + "\" in "
                          + m_blocks[i].m_name + "\'s problem is duplicated", ln);
                }
                if (ref.empty()) {
                    if (m_obj.find(*iit) != m_obj.end()) {
                        error("control variable \"" + iit->str() + "\" in " + m_blocks[i].m_name
                              + "\'s problem is objective variable in " + m_obj.find(*iit)->second
                              + "\'s problem", ln);
                    }
                    if (m_lagr_mult_in.find(*iit) != m_lagr_mult_in.end()) {
                        error("control variable \"" + iit->str() + "\" in " + m_blocks[i].m_name
                              + "\'s problem is Lagrange multiplier in " + m_lagr_mult_in.find(*iit)->second
                              + "\'s problem", ln);
                    }
                    if (!m_contr.insert(exstr(*iit, m_blocks[i].m_name)).second) {
                        warning("control variable \"" + iit->str() + "\" in "
                                + m_blocks[i].m_name + "\'s problem is control variable in "
                                + m_contr.find(*iit)->second + "\'s problem; "
                                + "if this was intentional you might consider rewriting your model in terms "
                                + "of two different controls and adding equilibrium condition equating them, "
                                + "one of them could then be selected for reduction", ln);
                    }
                }
            }
        }
    }
}



void
Model::leads()
{
    m_max_lag = 0;
    for (unsigned i = 0, n = m_blocks.size(); i < n; ++i) {
        if (!m_deter) {
            ex eq = m_blocks[i].m_obj_eq;
            int l = m_blocks[i].m_obj_line;
            if (eq.get_lag_max(true) > 0) {
                error("forward looking variable(s) in objective \""
                      + m_blocks[i].m_obj_var.str() + " = " + eq.str()
                      + "\" outside expected value operator in a stochastic model", l);
            }
            for (unsigned e = 0; e < m_blocks[i].m_constraints.size(); ++e) {
                eq = m_blocks[i].m_constraints[e].first;
                l = m_blocks[i].m_constraints[e].second;
                if (eq.get_lag_max(true) > 0) {
                    error("forward looking variable(s) in constraint \""
                        + eq.str() + " = 0\" outside expected value operator "
                        + "in a stochastic model", l);
                }
            }
            for (unsigned e = 0; e < m_blocks[i].m_identities.size(); ++e) {
                eq = m_blocks[i].m_identities[e].first;
                l = m_blocks[i].m_identities[e].second;
                if (eq.get_lag_max(true) > 0) {
                    error("forward looking variable(s) in identity \""
                        + eq.str() + " = 0\" outside expected value operator "
                        + "in a stochastic model", l);
                }
            }
        }

        ex eq = m_blocks[i].m_obj_eq;
        int l = m_blocks[i].m_obj_line;
        int mlag;
        if ((mlag = eq.get_lag_max()) > 1) {
            error("variable(s) in objective \""
                    + m_blocks[i].m_obj_var.str() + " = " + eq.str()
                    + "\" in lead > 1", l);
        }
        if (mlag > m_max_lag) m_max_lag = mlag;
        for (unsigned e = 0; e < m_blocks[i].m_constraints.size(); ++e) {
            eq = m_blocks[i].m_constraints[e].first;
            l = m_blocks[i].m_constraints[e].second;
            if ((mlag = eq.get_lag_max()) > 1) {
                error("variable(s) in constraint \""
                      + eq.str() + " = 0\" in lead > 1", l);
            }
            if (mlag > m_max_lag) m_max_lag = mlag;
        }
        for (unsigned e = 0; e < m_blocks[i].m_identities.size(); ++e) {
            eq = m_blocks[i].m_identities[e].first;
            l = m_blocks[i].m_identities[e].second;
            if ((mlag = eq.get_lag_max()) > 1) {
                error("variable(s) in identity \""
                    + eq.str() + " = 0\" in lead > 1", l);
            }
            if (mlag > m_max_lag) m_max_lag = mlag;
        }
    }
}



void
Model::lags()
{
    m_min_lag = 0;
    for (unsigned i = 0, n = m_blocks.size(); i < n; ++i) {
        int bmlag = 0, mlag;
        ex eq = m_blocks[i].m_obj_eq;
        mlag = eq.get_lag_min();
        if (mlag < bmlag) bmlag = mlag;
        for (unsigned e = 0; e < m_blocks[i].m_constraints.size(); ++e) {
            eq = m_blocks[i].m_constraints[e].first;
            mlag = eq.get_lag_min();
            if (mlag < bmlag) bmlag = mlag;
        }
        for (unsigned e = 0; e < m_blocks[i].m_identities.size(); ++e) {
            eq = m_blocks[i].m_identities[e].first;
            mlag = eq.get_lag_min();
            if (mlag < bmlag) bmlag = mlag;
        }
        if (bmlag < m_min_lag) m_min_lag = bmlag;
        if (bmlag < -1) {
            m_blocks[i].lags();
            set_ex::const_iterator its = m_blocks[i].m_lags.begin();
            set_ex::const_iterator itse = m_blocks[i].m_lags.end();
            for (; its != itse; ++its) {
                vec_ex vlags = expand(ex(m_blocks[i].m_i1, ex(m_blocks[i].m_i2, *its)));
                vec_ex::const_iterator iit = vlags.begin();
                for (; iit != vlags.end(); ++iit) {
                    m_lags.insert(*iit);
                }
            }
        }
    }
}



void
Model::check_static()
{
    if ((m_min_lag == 0) && (m_max_lag == 0)) m_static = true;
    if (m_static && !m_deter) {
        error("model is static (has no lags and leads) but has shocks; static models must also be deterministic", 0);
    }
}



// FOC derivation
// =============================================================================
void
Model::derive_focs()
{
    for (unsigned i = 0, n = m_blocks.size(); i < n; ++i) {
        // Collect equations from referenced blocks
        if (m_blocks[i].m_constraints_ref.size()) {
            for (unsigned r = 0; r < m_blocks[i].m_constraints_ref.size(); ++r) {
                strintint rf = m_blocks[i].m_constraints_ref[r];
                unsigned b = 0;
                while (m_blocks[b].m_name != rf.first) ++b;
                Model_block &bref = m_blocks[b];
                Model_block &bcurr = m_blocks[i];
                if (rf.second == Model_block::objective) {
                    ex eq = bref.m_obj_var - bref.m_obj_eq;
                    eq = ex(bref.m_i2, eq);
                    eq = ex(bref.m_i1, eq);
                    bcurr.m_constraints.push_back(exint(eq, -1));
                    ex lam = ex("lambda__" + bcurr.m_name + "_" + num2str(1 + (unsigned) bcurr.m_lagr_mult.size()), 0);
                    lam = add_idx(lam, bref.m_i1);
                    lam = add_idx(lam, bref.m_i2);
                    lam = add_idx(lam, bcurr.m_i1);
                    lam = add_idx(lam, bcurr.m_i2);
                    lam = apply_idx(lam, eq);
                    bcurr.m_lagr_mult.push_back(exint(lam, 0));
                    bcurr.m_redlm.insert(lam);
                } else if (rf.second == Model_block::constraints) {
                    for (unsigned j = 0; j < bref.m_constraints.size(); ++j) {
                        ex eq = bref.m_constraints[j].first;
                        ex lam = ex("lambda__" + bcurr.m_name + "_" + num2str(1 + (unsigned) bcurr.m_lagr_mult.size()), 0);
                        eq = ex(bref.m_i2, eq);
                        eq = ex(bref.m_i1, eq);
                        bcurr.m_constraints.push_back(exint(eq, -1));
                        lam = add_idx(lam, bref.m_i1);
                        lam = add_idx(lam, bref.m_i2);
                        lam = add_idx(lam, bcurr.m_i1);
                        lam = add_idx(lam, bcurr.m_i2);
                        lam = add_idx(lam, eq);
                        lam = apply_idx(lam, eq);
                        bcurr.m_lagr_mult.push_back(exint(lam, 0));
                        bcurr.m_redlm.insert(lam);
                    }
                } else if (rf.second == Model_block::focs) {
                    for (unsigned j = 0; j < bref.m_focs_red.size(); ++j) {
                        ex eq = bref.m_focs_red[j].first;
                        ex lam = ex("lambda__" + bcurr.m_name + "_" + num2str(1 + (unsigned) bcurr.m_lagr_mult.size()), 0);
                        eq = ex(bref.m_i2, eq);
                        eq = ex(bref.m_i1, eq);
                        bcurr.m_constraints.push_back(exint(eq, -1));
                        lam = add_idx(lam, bref.m_i1);
                        lam = add_idx(lam, bref.m_i2);
                        lam = add_idx(lam, bcurr.m_i1);
                        lam = add_idx(lam, bcurr.m_i2);
                        lam = add_idx(lam, eq);
                        lam = apply_idx(lam, eq);
                        bcurr.m_lagr_mult.push_back(exint(lam, 0));
                        bcurr.m_redlm.insert(lam);
                    }
                } else if (rf.second == Model_block::identities) {
                    for (unsigned j = 0; j < bref.m_identities.size(); ++j) {
                        ex eq = bref.m_identities[j].first;
                        ex lam = ex("lambda__" + bcurr.m_name + "_" + num2str(1 + (unsigned) bcurr.m_lagr_mult.size()), 0);
                        eq = ex(bref.m_i2, eq);
                        eq = ex(bref.m_i1, eq);
                        bcurr.m_constraints.push_back(exint(eq, -1));
                        lam = add_idx(lam, bref.m_i1);
                        lam = add_idx(lam, bref.m_i2);
                        lam = add_idx(lam, bcurr.m_i1);
                        lam = add_idx(lam, bcurr.m_i2);
                        lam = add_idx(lam, eq);
                        lam = apply_idx(lam, eq);
                        bcurr.m_lagr_mult.push_back(exint(lam, 0));
                        bcurr.m_redlm.insert(lam);
                    }
                } else INTERNAL_ERROR
            }
        }
        // Derive FOCs
        if (m_blocks[i].m_static) {
            m_blocks[i].derive_focs_static();
        } else if (m_deter) {
            m_blocks[i].derive_focs_deter();
        } else {
            m_blocks[i].derive_focs();
        }
        for (int f = 0, F = m_blocks[i].m_focs.size(); f < F; ++f) {
            if (!m_blocks[i].m_focs[f].first)
                warning("one of your first order conditions (w.r.t. \""
                        + m_blocks[i].m_focs[f].second.str() + "\") in " + m_blocks[i].m_name +
                        "\'s problem is of the form 0 = 0", 0);
        }
    }
}


// collecting variables, shocks & parameters
// =============================================================================
void
Model::collect_shocks()
{
    for (unsigned i = 0, n = m_blocks.size(); i < n; ++i) {
        std::vector<exint>::const_iterator it, ite;
        it = m_blocks[i].m_shocks.begin();
        ite = m_blocks[i].m_shocks.end();
        for (; it != ite; ++it) {
            if (it->first.get_lag_max()) {
                error("shock \"" + lag0(it->first).str() + "\" declared in lead/lag",
                      it->second);
                continue;
            }
            if (m_blocks[i].m_defs.find(it->first) != m_blocks[i].m_defs.end()) {
                error("shock \"" + it->first.str() + "\" declared in "
                      + m_blocks[i].m_name + "\'s block was used in definitions "
                      + "(substituted)", it->second);
                continue;
            }
            idx_ex i1 = m_blocks[i].m_i1;
            idx_ex i2 = m_blocks[i].m_i2;
            vec_ex shocks = expand(ex(i1, ex(i2, it->first)));
            vec_ex::const_iterator iit = shocks.begin();
            for (; iit != shocks.end(); ++iit) {
                if (!m_shocks.insert(*iit).second) {
                    error("shock " + iit->str() + " already declared", it->second);
                }
                if (m_obj.find(*iit) != m_obj.end()) {
                    error(m_obj.find(*iit)->second + "\'s objective variable \""
                          + iit->str() + "\" declared as shock", it->second);
                }
                if (m_contr.find(*iit) != m_contr.end()) {
                    error(m_contr.find(*iit)->second + "\'s control variable \""
                          + iit->str() + "\" declared as shock", it->second);
                }
                if (m_lagr_mult_in.find(*iit) != m_lagr_mult_in.end()) {
                    error(m_lagr_mult_in.find(*iit)->second + "\'s Lagrange multiplier \""
                          + iit->str() + "\" declared as shock", it->second);
                }
            }
        }
    }
}


namespace {

// Easter eggs ;-)
const char* namescomments[] = {
    "dupa", "Oj, brzydko, brzydko ;-)",
    "Dupa", "Oj, brzydko, brzydko ;-)",
    "DUPA", "Oj, brzydko, brzydko ;-)",
    "foo", "Hey, try and be more creative ;-)",
    "Foo", "Hey, try and be more creative ;-)",
    "FOO", "Hey, try and be more creative ;-)",
    "bar", "Hey, try and be more creative ;-)",
    "Bar", "Hey, try and be more creative ;-)",
    "BAR", "Hey, try and be more creative ;-)",
    "foobar", "Hey, try and be more creative ;-)",
    "foo_bar", "Hey, try and be more creative ;-)",
    "FOOBAR", "Hey, try and be more creative ;-)",
    "FOO_BAR", "Hey, try and be more creative ;-)",
""};


void
easter_eggs(Model &model, const std::set<std::string> &names)
{
    for (unsigned i = 0; namescomments[i][0]; i += 2) {
        std::string nm = namescomments[i];
        if (names.find(nm) != names.end())
            model.warning("\"" + nm + "\"? " + namescomments[i + 1], 0);
    }
}

} /* namespace */



void
Model::collect_vp()
{
    for (unsigned i = 0, n = m_blocks.size(); i < n; ++i) {
        m_blocks[i].collect_vp(m_vars, m_params);
    }
    if (!m_vars.size()) {
        error("the model has no variables", 0);
        terminate_on_errors();
    }

    // Test names of vars and params
    set_ex::const_iterator it;
    std::set<std::string> names;
    for (it = m_vars.begin(); it != m_vars.end(); ++it) {
        std::string name = it->str(DROP_T);
        names.insert(name);
    }
    for (it = m_shocks.begin(); it != m_shocks.end(); ++it) {
        std::string name = it->str(DROP_T);
        names.insert(name);
    }
    for (it = m_params.begin(); it != m_params.end(); ++it) {
        std::string name = it->str();
        if (!names.insert(name).second) {
            error("\"" + name + "\" treated both as a variable and a parameter; you "
                  "might have forgotten \"[]\"", 0);
        }
    }
    names.insert(m_names.begin(), m_names.end());
    easter_eggs(*this, names);
    for (it = m_shocks.begin(); it != m_shocks.end(); ++it) {
        if (m_vars.find(*it) == m_vars.end()) {
            warning("shock \"" + it->str() + "\" does not appear in any model equation", 0);
        } else m_vars.erase(*it);
    }
}




namespace {

std::string
mk_diff_no_error(unsigned n1, const std::string &name1, unsigned n2, const std::string &name2)
{
    std::string mes = "There ";
    if (n1 > 1) {
        mes += "are ";
        if (n1 < n2) { mes += "only "; }
        mes += num2str(n1);
        mes += " ";
    } else if (n1 == 1) {
        mes += "is ";
        if (n1 < n2) { mes += "only "; }
        mes += "1 ";
    } else {
        mes += "are no ";
    }
    mes += name1;
    if (n1 == 1) { mes += " "; } else { mes += "s "; }
    mes += "but there ";
    if (n2 > 1) {
        mes += "are ";
        if (n1 > n2) { mes += "only "; }
        mes += num2str(n2);
        mes += " ";
    } else if (n2 == 1) {
        mes += "is ";
        if (n1 > n2) { mes += "only "; }
        mes += "1 ";
    } else {
        mes += "are no ";
    }
    mes += name2;
    if (n2 != 1) { mes += "s"; }

    return mes;
}

} /* namespace */



void
Model::collect_eq()
{
    int i, n;
    for (i = 0, n = m_blocks.size(); i < n; ++i) {
        std::vector<expair>::const_iterator itp, itep;
        idx_ex i1 = m_blocks[i].m_i1;
        idx_ex i2 = m_blocks[i].m_i2;
        itp = m_blocks[i].m_focs_red.begin();
        itep = m_blocks[i].m_focs_red.end();
        for (; itp != itep; ++itp) {
            ex e = ex(i1, ex(i2, itp->first));
            vec_ex eqs = expand(e);
            vec_ex::const_iterator iit = eqs.begin();
            for (; iit != eqs.end(); ++iit) {
                if (!m_eqs.insert(*iit).second) {
                    warning("repeating equation: \"" + itp->first.str() + " = 0\"", 0);
                }
            }
            if (eqs.size()) m_t_eqs.insert(e);
        }
        if (m_blocks[i].m_obj_eq) {
            ex lhs = m_blocks[i].m_obj_var, rhs = m_blocks[i].m_obj_eq;
            ex e = ex(i1, ex(i2, lhs - rhs));
            vec_ex eqs = expand(e);
            vec_ex::const_iterator iit = eqs.begin();
            for (; iit != eqs.end(); ++iit) {
                if (!m_eqs.insert(*iit).second) {
                    warning("equation for " + m_blocks[i].m_name + "\'s objective \""
                            + m_blocks[i].m_obj_eq.str() + "\" is duplicated", 0);
                }
            }
            if (eqs.size()) m_t_eqs.insert(e);
        }
        std::vector<exint>::const_iterator it, ite;
        it = m_blocks[i].m_constraints.begin();
        ite = m_blocks[i].m_constraints.end();
        for (; it != ite; ++it) {
            bool internal = (it->second <= 0);
            ex e = ex(i1, ex(i2, it->first));
            vec_ex eqs = expand(e);
            vec_ex::const_iterator iit = eqs.begin();
            for (; iit != eqs.end(); ++iit) {
                if ((!m_eqs.insert(*iit).second) && (!internal)) {
                    warning("repeating constraint: \"" + iit->str() + " = 0\"",
                            it->second);
                }
            }
            if (eqs.size()) m_t_eqs.insert(e);
        }
        it = m_blocks[i].m_identities.begin();
        ite = m_blocks[i].m_identities.end();
        for (; it != ite; ++it) {
            bool internal = !(it->second);
            ex e = ex(i1, ex(i2, it->first));
            vec_ex eqs = expand(e);
            vec_ex::const_iterator iit = eqs.begin();
            for (; iit != eqs.end(); ++iit) {
                if ((!m_eqs.insert(*iit).second) && (!internal)) {
                    warning("repeating identity: \"" + iit->str() + " = 0\"",
                            it->second);
                }
            }
            if (eqs.size()) m_t_eqs.insert(e);
        }
    }

    if (!m_eqs.size()) {
        error("the model has no equations", 0);
        terminate_on_errors();
    }
    unsigned nv = m_vars.size(), ne = m_eqs.size();
    if (nv != ne) {
        error(mk_diff_no_error(nv, "variable", ne, "model equation"), 0);
    }
}




void
Model::collect_calibr()
{
    int i, n;
    set_ex params_set, params_fr_add;
    for (i = 0, n = m_blocks.size(); i < n; ++i) {
        std::vector<exint>::const_iterator it, ite;
        unsigned j, m = m_blocks[i].m_calibr.size();
        for (j = 0; j < m; ++j) {
            ex ceq = m_blocks[i].m_calibr[j].first;
            int lineno = m_blocks[i].m_calibr[j].second;
            if (!ceq) {
                error("calibrating equation \"0 = 0\"", lineno);
                continue;
            }
            if (m_static) {
                set_ex v, p;
                collect(ceq, v, p);
                bool ok = true;
                for (set_ex::const_iterator it = v.begin(); it != v.end(); ++it) {
                    if (ceq.has(*it, DIFF_T)) {
                        error("variable \"" + it->str() + "\" in calibrating equation \"" + ceq.str()
                              + "\" has invalid time index (model is static)", lineno);
                        ok = false;
                    }
                }
                if (!ok) continue;
            } else {
                if (ceq.hast()) {
                    error("calibrating equation with non-steady-state values: \"" + ceq.str()
                          + " = 0\"", lineno);
                    continue;
                }
            }
            triplet<bool, ex, ex> ps = find_par_eq_num(ceq);
            if (ps.first) {
                if (m_blocks[i].m_calibr_pl[j].size()) {
                    error("equation \"" + ps.second.str() + " = " + ps.third.str()
                          + "\" sets value of a parameter and is not a proper "
                          + "calibrating equation so it cannot be followed by a parameter list",
                          lineno);
                    continue;
                }
            }
            idx_ex i1 = m_blocks[i].m_i1;
            idx_ex i2 = m_blocks[i].m_i2;
            vec_ex vceq = expand(ex(i1, ex(i2, ceq)));
            vec_ex::const_iterator iit = vceq.begin();
            for (; iit != vceq.end(); ++iit) {
                set_ex vars, pars;
                collect(*iit, vars, pars);
                set_ex::const_iterator it, ite;
                for (it = vars.begin(), ite = vars.end(); it != ite; ++it) {
                    if (m_vars.find(*it) == m_vars.end()) {
                        if (m_shocks.find(*it) != m_shocks.end()) {
                            error("shock \"" + it->str() + "\" appears in calibrating equation", lineno);
                        } else {
                            error("variable \"" + it->str() + "\" appears in calibrating equation but does not "
                                  + "appear in any model equation", lineno);
                        }
                    }
                }
                for (it = pars.begin(), ite = pars.end(); it != ite; ++it) {
                    if (m_params.find(*it) == m_params.end()) {
                        m_params.insert(*it);
                        m_params_free.insert(*it);
                        params_fr_add.insert(*it);
                    }
                }
                ps = find_par_eq_num(*iit);
                if (ps.first) {
                    if (!params_set.insert(ps.second).second) {
                        error("free parameter \"" + ps.second.str()
                              + "\" already set to \"" + m_params_free_set.find(ps.second)->second.str() + "\"",
                              lineno);
                    }
                    if (!ps.third.validnum()) {
                        error("value of parameter \"" + ps.second.str() + "\" set to \"" + ps.third.str() + "\"",
                              lineno);
                    }
                    m_params_free_set.insert(expair(ps.second, ps.third));
                    if (m_params_calibr.find(ps.second) != m_params_calibr.end()) {
                        error("parameter \"" + ps.second.str() + "\" was earlier declared as calibrated "
                              + "parameter", lineno);
                    }
                } else {
                    pars.clear();
                    for (unsigned k = 0; k < m_blocks[i].m_calibr_pl[j].size(); ++k) {
                        ex pp = m_blocks[i].m_calibr_pl[j][k].first;
                        vec_ex vpp = expand(ex(i1, ex(i2, apply_idx(pp, *iit))));
                        vec_ex::const_iterator ipp = vpp.begin();
                        for (; ipp != vpp.end(); ++ipp) {
                            if (m_params.find(*ipp) == m_params.end()) {
                                error("parameter \"" + ipp->str() + "\" declared as calibrated, but does not "
                                      + "appear in any model equation",
                                      m_blocks[i].m_calibr_pl[j][k].second);
                            }
                            if (!pars.insert(*ipp).second) {
                                error("repeating parameter \"" + ipp->str() + "\" in a parameter list",
                                      m_blocks[i].m_calibr_pl[j][k].second);
                            } else {
                                m_params_calibr.insert(*ipp);
                            }
                            map_ex_ex::const_iterator mit;
                            if ((mit = m_params_free_set.find(*ipp)) != m_params_free_set.end()) {
                                error("parameter \"" + ipp->str() + "\" declared as calibrated parameter "
                                      + "yet at the same time its value was set to \"" + mit->second.str() + "\"",
                                      m_blocks[i].m_calibr_pl[j][k].second);
                            }
                        }
                    }

                    if (!m_calibr.insert(*iit).second) {
                        warning("repeating calibration equation \"" + iit->str()
                                + " = 0\"", lineno);
                    }
                }
            }
        }
    }

    unsigned nfp = m_params_calibr.size(), ce = m_calibr.size();
    if (ce != nfp) {
        error(mk_diff_no_error(nfp, "non-free (calibrated) parameter", ce, "calibrating equation"), 0);
    }

    for (set_ex::const_iterator it = m_params.begin(), ite = m_params.end();
         it != ite; ++it) {
        if (m_params_calibr.find(*it) == m_params_calibr.end())
            m_params_free.insert(*it);
    }
    if (m_options[verbose] && params_fr_add.size()) {
        std::string plst;
        for (set_ex::const_iterator it = params_fr_add.begin(),
         ite = params_fr_add.end(); it != ite;) {
            plst += "\"" + it->str() + "\"";
            if (++it != ite) {
                plst += ", ";
            }
        }
        write_model_info("the following parameter(s) appearing in the calibrating equations only have been added to the free parameter list: " + plst);
    }
}



// variable reduction
// =============================================================================
void
Model::check_red_vars2()
{
    // Collect internally generated Lagrange multipliers for reduction
    for (unsigned i = 0, n = m_blocks.size(); i < n; ++i) {
        ex lmin = m_blocks[i].m_obj_lm_in;
        idx_ex i1 = m_blocks[i].m_i1;
        idx_ex i2 = m_blocks[i].m_i2;
        set_ex::const_iterator its = m_blocks[i].m_redlm.begin();
        set_ex::const_iterator itse = m_blocks[i].m_redlm.end();
        for (; its != itse; ++its) {
            vec_ex vlmin = expand(ex(i1, ex(i2, *its)));
            vec_ex::const_iterator iit = vlmin.begin();
            for (; iit != vlmin.end(); ++iit) {
                if (m_vars.find(*iit) != m_vars.end())
                    m_lagr_mult.insert(*iit);
            }
        }
    }
    set_ex::const_iterator it = m_s_redvars.begin();
    for (; it != m_s_redvars.end(); ++it) {
        if (m_vars.find(*it) == m_vars.end()) {
            warning("variable \"" + it->str() + "\" selected for reduction but "
                    + "does not appear in any model equation; removing from the list",
                    0);
            m_s_redvars.erase(*it);
            continue;
        }
        set_ex::const_iterator it1;
        for (it1 = m_calibr.begin(); it1 != m_calibr.end(); ++it1) {
            if (it1->has(*it, ANY_T)) {
                warning("variable \"" + it->str()
                        + "\" selected for reduction appears in calibrating equation \""
                        + it1->str() + " = 0\"", 0);
            }
        }
    }
    m_s_redvars.insert(m_lagr_mult.begin(), m_lagr_mult.end());
    m_s_redvars.insert(m_lags.begin(), m_lags.end());
}



void
Model::reduce()
{
    vec_ex eqs, eqsc;
    eqs.reserve(m_eqs.size());
    for (set_ex::iterator it = m_eqs.begin(); it != m_eqs.end(); ++it) eqs.push_back(*it);
    eqsc.reserve(m_calibr.size());
    for (set_ex::iterator it = m_calibr.begin(); it != m_calibr.end(); ++it) eqsc.push_back(*it);

    unsigned i, j, n = eqs.size(), nc = eqsc.size();
    bool try_red;
    do {
        try_red = false;
        for (i = 0; i < n; ++i) {
            triplet<bool, ex, ex> ts = find_subst(eqs[i], m_s_redvars);
            if (!ts.first) continue;
            if (ts.third.hast()) {
                ex e = ts.second, lde = lag(e, 1), lge = lag(e, -1);
                int ld = 0, lg = 0;
                for (unsigned ii = 0; ii < n; ++ii) {
                    if (eqs[ii].has(lde)) {
                        ld = 1;
                    }
                    if (eqs[ii].has(lge)) {
                        lg = -1;
                    }
                    if (ld && lg) break;
                }
                if (ts.third.get_lag_max() + ld > 1) continue;
                if (ts.third.get_lag_min() + lg < -1) continue;
                // check if we have shocks on RHS in substitution
                bool hasshock = false;
                for (set_ex::const_iterator lit = m_shocks.begin();
                     lit != m_shocks.end(); ++lit) {
                    if (ts.third.has(*lit)) {
                        hasshock = true;
                        break;
                    }
                }
                if (hasshock && (ld || lg)) continue;
            }
            ex what = ts.second, with = ts.third;
            m_lagr_mult.erase(what);
            m_lags.erase(what);
            m_vars.erase(what);
            m_s_redvars.erase(what);
            eqs[i] = ex();
            for (j = 0; j < n; ++j) {
                eqs[j] = eqs[j].subst(what, with);
            }
            for (j = 0; j < nc; ++j) {
                eqsc[j] = eqsc[j].subst(what, with);
                for (set_ex::const_iterator lit = m_shocks.begin();
                     lit != m_shocks.end(); ++lit) {
                    eqsc[j] = eqsc[j].subst(ss(*lit), ex());
                }
            }
            try_red = true;
            break;
        }
    } while (try_red && m_s_redvars.size());

    if (m_s_redvars.size()) { // Reduction of variables appearing in 1 eq only
        vec_ex rvars;
        rvars.reserve(m_s_redvars.size());
        for (set_ex::iterator it = m_s_redvars.begin(); it != m_s_redvars.end(); ++it)
            rvars.push_back(*it);
        std::vector<int> counts(rvars.size());
        for (j = 0; j < rvars.size(); ++j) {
            ex e = rvars[j];
            bool incalibr = false;
            for (i = 0; i < nc; ++i) {
                if (eqsc[i].has(e, ANY_T)) {
                    incalibr = true;
                    break;
                }
            }
            if (!incalibr) {
                for (i = 0; i < n; ++i) {
                    if (eqs[i].has(e, ANY_T)) {
                        ++counts[j];
                    }
                }
            }
        }
        for (j = 0; j < rvars.size(); ++j) {
            if (counts[j] == 1) {
                ex e = rvars[j];
                for (i = 0; i < n; ++i) {
                    if (eqs[i].has(e, ANY_T)) {
                        eqs[i] = ex();
                        m_lagr_mult.erase(e);
                        m_lags.erase(e);
                        m_vars.erase(e);
                        m_s_redvars.erase(e);
                        break;
                    }
                }
            }
        }
    }

    m_eqs.clear();
    for (i = 0; i < n; ++i) {
        ex eq = eqs[i];
        if (eq) m_eqs.insert(eq);
    }
    m_calibr.clear();
    for (i = 0; i < nc; ++i) {
        ex eq = eqsc[i];
        if (eq) m_calibr.insert(eq);
    }

    unsigned nv = m_vars.size(), ne = m_eqs.size();
    if (nv != ne) {
        error(mk_diff_no_error(nv, "variable", ne, "model equation"), 0);
    }
    unsigned nfp = m_params_calibr.size(), ce = m_calibr.size();
    if (ce != nfp) {
        error(mk_diff_no_error(nfp, "non-free (calibrated) parameter", ce, "calibrating equation"), 0);
    }

    for (set_ex::const_iterator it = m_lagr_mult.begin();
             it != m_lagr_mult.end(); ++it) {
        m_s_redvars.erase(*it);
    }
    for (set_ex::const_iterator it = m_lags.begin();
             it != m_lags.end(); ++it) {
        m_s_redvars.erase(*it);
    }
    if (m_s_redvars.size()) {
        std::string vlst;
        print_flag pflag = (m_static) ? DROP_T : DEFAULT;
        for (set_ex::const_iterator it = m_s_redvars.begin();
             it != m_s_redvars.end();) {
            vlst += "\"" + it->str(pflag) + "\"";
            if (++it != m_s_redvars.end()) {
                vlst += ", ";
            }
        }
        warning("the following variable(s) selected for reduction could not be \
symbolically reduced in the model: " + vlst, 0);
    }
    if (m_options[verbose] && (m_lagr_mult.size() || m_lags.size())) {
        std::string vlst;
        print_flag pflag = (m_static) ? DROP_T : DEFAULT;
        for (set_ex::const_iterator it = m_lagr_mult.begin();
             it != m_lagr_mult.end();) {
            vlst += "\"" + it->str(pflag) + "\"";
            if (++it != m_lagr_mult.end()) {
                vlst += ", ";
            }
        }
        if (m_lagr_mult.size() && m_lags.size()) vlst += ", ";
        for (set_ex::const_iterator it = m_lags.begin();
             it != m_lags.end();) {
            vlst += "\"" + it->str(pflag) + "\"";
            if (++it != m_lags.end()) {
                vlst += ", ";
            }
        }
        write_model_info("the following internally generated variable(s) could not be \
symbolically reduced in the model: " + vlst);
    }
}




// maps, steady state & derivatives
// =============================================================================

namespace {

enum lag_flags {
    NUL    =     0,
    LAG_M1 =   0x1,
    LAG_0  =   0x2,
    LAG_P1 =   0x4,
    LAG_SS =   0x8
};

} /* namespace */


void
Model::var_eq_map()
{
    int i, j;
    unsigned flag;
    set_ex::const_iterator it1, it2;
    ex x, e;

    if (m_static) {
        for (it1 = m_eqs.begin(), i = 1; it1 != m_eqs.end(); ++it1, ++i) {
            e = *it1;
            for (it2 = m_vars.begin(), j = 1; it2 != m_vars.end(); ++it2, ++j) {
                x = *it2;
                if (e.has(x)) m_var_eq_map.insert(std::pair<std::pair<int, int>, unsigned>(
                                            std::pair<int, int>(i, j), LAG_0));
            }
        }
    } else {
        for (it1 = m_eqs.begin(), i = 1; it1 != m_eqs.end(); ++it1, ++i) {
            e = *it1;
            for (it2 = m_vars.begin(), j = 1; it2 != m_vars.end(); ++it2, ++j) {
                flag = 0;
                x = lag(*it2, -1);
                if (e.has(x)) flag |= LAG_M1;
                x = *it2;
                if (e.has(x)) flag |= LAG_0;
                x = lag(*it2, 1);
                if (e.has(x)) flag |= LAG_P1;
                x = ss(*it2);
                if (e.has(x)) flag |= LAG_SS;
                if (flag) m_var_eq_map.insert(std::pair<std::pair<int, int>, unsigned>(
                                            std::pair<int, int>(i, j), flag));
            }
        }
    }
}


void
Model::shock_eq_map()
{
    int i, j;
    set_ex::const_iterator it1, it2;
    ex x, e;

    for (it1 = m_eqs.begin(), i = 1; it1 != m_eqs.end(); ++it1, ++i) {
        e = *it1;
        for (it2 = m_shocks.begin(), j = 1; it2 != m_shocks.end(); ++it2, ++j) {
            if (e.has(*it2, DIFF_T)) {
                error("shock \"" + it2->str() + "\" in equation \"" + it1->str()
                      + "\" has invalid time index; shocks should have time index 0", 0);
            }
            if (e.has(*it2)) m_shock_eq_map.insert(std::pair<int, int>(i, j));
        }
    }
}





void
Model::stst()
{
    set_ex ss_set;
    set_ex::const_iterator it, it2;
    unsigned vs = m_vars.size();

    if (m_static) return;

    for (it = m_eqs.begin(); it != m_eqs.end(); ++it) {
        ex ssex = ss(*it);
        for (it2 = m_shocks.begin(); it2 != m_shocks.end(); ++it2) {
            ssex = ssex.subst(ss(*it2), ex());
        }
        if (!ssex) {
            warning("steady state equation \"0 = 0\" derived from \""
                    + it->str() + "\"", 0);
        } else if (!ss_set.insert(ssex).second) {
            warning("repeating steady state equation: \"" + ssex.str() + " = 0\"", 0);
        }
        m_ss.push_back(ssex);
    }

    unsigned sse = m_ss.size();
    if (vs != sse) {
        error(mk_diff_no_error(vs, "variable", sse, "steady state equation"), 0);
    }

    for (it = m_t_eqs.begin(); it != m_t_eqs.end(); ++it) {
        m_t_ss.push_back(ss(*it));
    }
}



void
Model::var_ceq_map()
{
    int i, j;
    set_ex::const_iterator it1, it2;

    if (m_static) {
        for (it1 = m_calibr.begin(), i = 1; it1 != m_calibr.end(); ++it1, ++i) {
            for (it2 = m_vars.begin(), j = 1; it2 != m_vars.end(); ++it2, ++j) {
                if (it1->has(*it2)) {
                    m_var_ceq_map.insert(std::pair<int, int>(i, j));
                }
            }
        }
    }
    for (it1 = m_calibr.begin(), i = 1; it1 != m_calibr.end(); ++it1, ++i) {
        for (it2 = m_vars.begin(), j = 1; it2 != m_vars.end(); ++it2, ++j) {
            if (it1->has(ss(*it2))) m_var_ceq_map.insert(std::pair<int, int>(i, j));
        }
    }
}



void
Model::par_eq_map()
{
    int i, j;
    set_ex::const_iterator it1, it2;
    ex x, r, e;

    for (it1 = m_eqs.begin(), i = 1; it1 != m_eqs.end(); ++it1, ++i) {
        for (it2 = m_params_calibr.begin(), j = 1; it2 != m_params_calibr.end(); ++it2, ++j) {
            if (it1->has(*it2)) m_cpar_eq_map.insert(std::pair<int, int>(i, j));
        }
        for (it2 = m_params_free.begin(), j = 1; it2 != m_params_free.end(); ++it2, ++j) {
            if (it1->has(*it2)) m_fpar_eq_map.insert(std::pair<int, int>(i, j));
        }
    }
}



void
Model::par_ceq_map()
{
    int i, j;
    set_ex::const_iterator it1, it2;

    for (it1 = m_calibr.begin(), i = 1; it1 != m_calibr.end(); ++it1, ++i) {
        for (it2 = m_params_calibr.begin(), j = 1; it2 != m_params_calibr.end(); ++it2, ++j) {
            if (it1->has(*it2)) m_cpar_ceq_map.insert(std::pair<int, int>(i, j));
        }
        for (it2 = m_params_free.begin(), j = 1; it2 != m_params_free.end(); ++it2, ++j) {
            if (it1->has(*it2)) m_fpar_ceq_map.insert(std::pair<int, int>(i, j));
        }
    }
}






void
Model::ss_jacob()
{
    int i, j, nv = m_vars.size(), ne = m_eqs.size();
    set_ex::const_iterator it0, it1, it2;
    vec_ex::const_iterator it0v;
    ex x, e, r;

    if (m_static) {
        for (it0 = m_eqs.begin(), i = 1; it0 != m_eqs.end(); ++it0, ++i) {
            e = *it0;
            for (it2 = m_vars.begin(), j = 1; it2 != m_vars.end(); ++it2, ++j) {
                x = *it2;
                r = diff(e, x);
                if (r) {
                    m_jacob_ss_calibr.insert(std::pair<std::pair<int, int>, ex>(
                                    std::pair<int, int>(i, j), r));
                }
            }
        }
        for (it0 = m_eqs.begin(), i = 1; it0 != m_eqs.end(); ++it0, ++i) {
            e = *it0;
            for (it2 = m_params_calibr.begin(), j = 1;
                it2 != m_params_calibr.end(); ++it2, ++j) {
                x = *it2;
                r = diff(e, x);
                if (r) {
                    m_jacob_ss_calibr.insert(std::pair<std::pair<int, int>, ex>(
                                    std::pair<int, int>(i, j + nv), r));
                }
            }
        }
        for (it1 = m_calibr.begin(), i = 1; it1 != m_calibr.end(); ++it1, ++i) {
            e = *it1;
            for (it2 = m_vars.begin(), j = 1; it2 != m_vars.end(); ++it2, ++j) {
                x = *it2;
                r = diff(e, x);
                if (r) {
                    m_jacob_ss_calibr.insert(std::pair<std::pair<int, int>, ex>(
                                    std::pair<int, int>(i + ne, j), r));
                }
            }
        }
    } else {
        for (it0v = m_ss.begin(), i = 1; it0v != m_ss.end(); ++it0v, ++i) {
            e = *it0v;
            for (it2 = m_vars.begin(), j = 1; it2 != m_vars.end(); ++it2, ++j) {
                x = ss(*it2);
                r = diff(e, x);
                if (r) {
                    m_jacob_ss_calibr.insert(std::pair<std::pair<int, int>, ex>(
                                    std::pair<int, int>(i, j), r));
                }
            }
        }
        for (it0v = m_ss.begin(), i = 1; it0v != m_ss.end(); ++it0v, ++i) {
            e = *it0v;
            for (it2 = m_params_calibr.begin(), j = 1;
                it2 != m_params_calibr.end(); ++it2, ++j) {
                x = *it2;
                r = diff(e, x);
                if (r) {
                    m_jacob_ss_calibr.insert(std::pair<std::pair<int, int>, ex>(
                                    std::pair<int, int>(i, j + nv), r));
                }
            }
        }
        for (it1 = m_calibr.begin(), i = 1; it1 != m_calibr.end(); ++it1, ++i) {
            e = *it1;
            for (it2 = m_vars.begin(), j = 1; it2 != m_vars.end(); ++it2, ++j) {
                x = ss(*it2);
                r = diff(e, x);
                if (r) {
                    m_jacob_ss_calibr.insert(std::pair<std::pair<int, int>, ex>(
                                    std::pair<int, int>(i + ne, j), r));
                }
            }
        }
    }
    for (it1 = m_calibr.begin(), i = 1; it1 != m_calibr.end(); ++it1, ++i) {
        e = *it1;
        for (it2 = m_params_calibr.begin(), j = 1;
             it2 != m_params_calibr.end(); ++it2, ++j) {
            x = *it2;
            r = diff(e, x);
            if (r) {
                m_jacob_ss_calibr.insert(std::pair<std::pair<int, int>, ex>(
                                  std::pair<int, int>(i + ne, j + nv), r));
            }
        }
    }
}




void
Model::diff_eqs()
{
    int i, j;
    unsigned fl;
    set_ex::const_iterator it1, it2, it3;
    std::map<std::pair<int, int>, unsigned>::const_iterator itf;
    ex x, r, e;

    if (m_static) return;

    for (it1 = m_eqs.begin(), i = 1; it1 != m_eqs.end(); ++it1, ++i) {
        e = *it1;
        for (it2 = m_vars.begin(), j = 1; it2 != m_vars.end(); ++it2, ++j) {
            itf = m_var_eq_map.find(std::pair<int, int>(i, j));
            if (itf == m_var_eq_map.end()) continue;
            fl = itf->second;
            if (fl & LAG_M1) {
                x = lag(*it2, -1);
                r = diff(e, x);
                r = ss(r);
                if (r) {
                    for (it3 = m_shocks.begin(); (it3 != m_shocks.end()) && (r); ++it3) {
                        r = r.subst(ss(*it3), ex());
                    }
                }
                if (r) m_Atm1.insert(std::pair<std::pair<int, int>, ex>(
                            std::pair<int, int>(i, j), r));
            }
            if (fl & LAG_0) {
                x = *it2;
                r = diff(e, x);
                r = ss(r);
                if (r) {
                    for (it3 = m_shocks.begin(); (it3 != m_shocks.end()) && (r); ++it3) {
                        r = r.subst(ss(*it3), ex());
                    }
                }
                if (r) m_At.insert(std::pair<std::pair<int, int>, ex>(
                            std::pair<int, int>(i, j), r));
            }
            if (fl & LAG_P1) {
                x = lag(*it2, 1);
                r = diff(e, x);
                r = ss(r);
                if (r) {
                    for (it3 = m_shocks.begin(); (it3 != m_shocks.end()) && (r); ++it3) {
                        r = r.subst(ss(*it3), ex());
                    }
                }
                if (r) m_Atp1.insert(std::pair<std::pair<int, int>, ex>(
                            std::pair<int, int>(i, j), r));
            }
        }
        for (it2 = m_shocks.begin(), j = 1; it2 != m_shocks.end(); ++it2, ++j) {
            x = *it2;
            r = diff(e, x);
            r = ss(r);
            if (r) {
                for (it3 = m_shocks.begin(); (it3 != m_shocks.end()) && (r); ++it3) {
                    r = r.subst(ss(*it3), ex());
                }
            }
            if (r) m_Aeps.insert(std::pair<std::pair<int, int>, ex>(
                        std::pair<int, int>(i, j), r));
        }
    }
}



