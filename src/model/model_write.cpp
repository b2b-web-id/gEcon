/*****************************************************************************
 * This file is a part of gEcon.                                             *
 *                                                                           *
 * (c) Chancellery of the Prime Minister of the Republic of Poland 2012-2015 *
 * (c) Grzegorz Klima, Karol Podemski, Kaja Retkiewicz-Wijtiwiak 2015-2018   *
 * (c) Karol Podemski, Kaja Retkiewicz-Wijtiwiak 2018-2019                   *
 * License terms can be found in the file 'LICENCE'                          *
 *                                                                           *
 * Author: Grzegorz Klima, Karol Podemski                                    *
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
#include <regex>

using symbolic::internal::num2str;
using symbolic::internal::print_flag;
using symbolic::internal::DEFAULT;
using symbolic::internal::CONVERT_T;
using symbolic::internal::DROP_T;
using symbolic::internal::CONVERT_IDX;
using symbolic::internal::DROP_IDX;
using symbolic::internal::DROP_QUOTES;




namespace {

std::string
time_str()
{
    std::string ts, mo, da, ho, mi, se;
    time_t tt = time(0);
    struct tm *now = localtime(&tt);
    mo = (now->tm_mon < 9) ? '0' + num2str(now->tm_mon + 1) : num2str(now->tm_mon + 1);
    da = (now->tm_mday < 10) ? '0' + num2str(now->tm_mday) : num2str(now->tm_mday);
    ho = (now->tm_hour < 10) ? '0' + num2str(now->tm_hour) : num2str(now->tm_hour);
    mi = (now->tm_min < 10) ? '0' + num2str(now->tm_min) : num2str(now->tm_min);
    se = (now->tm_sec < 10) ? '0' + num2str(now->tm_sec) : num2str(now->tm_sec);
    ts = num2str(now->tm_year + 1900) + '-' + mo + '-' + da
        + ' ' + ho + ':' + mi + ':' + se;
    return ts;
}


std::string
info_str(std::string newline = "")
{
    return "Generated on " + time_str() + " by gEcon ver. " + gecon_ver_str()
           + '\n' + newline + gecon_web_str();
}

} /* namespace */



void
Model::write_log() const
{
    std::string full_name = m_path + m_name + ".model.log";
#ifdef DEBUG
    std::cerr << "DEBUG INFO: writing logfile to: \'" << full_name << "\'\n";
#endif /* DEBUG */
    std::ofstream logfile(full_name.c_str());
    if (!logfile.good()) {
        report_errors("(gEcon error): failed opening logfile \'" + full_name
                      + "\' for writing");
    }

    logfile << info_str() << "\n\n";
    logfile << "Model name: " << m_name << "\n\n";

    std::string errwarn;
    if (errors()) {
        errwarn = num2str((unsigned)m_err.size()) + " ERROR";
        if (m_err.size() > 1) errwarn += "S";
    }
    if (warnings()) {
        if (errors()) errwarn += ", ";
        errwarn += num2str((unsigned)m_warn.size()) + " WARNING";
        if (m_warn.size() > 1) errwarn += "S";
    }
    if (errors() || warnings()) {
        logfile << errwarn << ", see bottom of this logfile\n\n\n";
    }

    std::string tab("    ");
    if (m_sets.size()) {
        logfile << "Index sets (" << m_sets.size() << "):\n";
        std::map<std::string, symbolic::idx_set>::const_iterator it;
        for (it = m_sets.begin(); it != m_sets.end(); ++it) {
            logfile << tab << it->second.str() << '\n';
        }
        logfile << '\n';
    }

    print_flag pflag = (m_static) ? DROP_T : DEFAULT;

    if (m_v_redvars.size()) {
        logfile << "Variables selected for reduction:\n";
        vec_exint::const_iterator it = m_v_redvars.begin();
        logfile << tab << it->first.str(pflag);
        for (++it; it != m_v_redvars.end(); ++it) {
            logfile << ", " << it->first.str(pflag);
        }
        logfile << "\n\n";
    }

    for (unsigned i = 0, n = m_blocks.size(); i < n; ++i) {
        m_blocks[i].write_logfile(logfile, m_static);
    }

    set_ex::const_iterator it;
    vec_ex::const_iterator itv;

    if (m_vars.size()) {
        logfile << "Variables (" << m_vars.size() << "):\n";
        it = m_vars.begin();
        logfile << tab << it->str(pflag);
        for (++it; it != m_vars.end(); ++it) {
            logfile << ", "<< it->str(pflag);
        }
        logfile << "\n\n";
    }

    if (m_shocks.size()) {
        logfile << "Shocks (" << m_shocks.size() << "):\n";
        it = m_shocks.begin();
        logfile << tab << it->str(pflag);
        for (++it; it != m_shocks.end(); ++it) {
            logfile << ", "<< it->str(pflag);
        }
        logfile << "\n\n";
    }

    if (m_params.size()) {
        logfile << "Parameters (" << m_params.size() << "):\n";
        it = m_params.begin();
        logfile << tab << it->str(pflag);
        for (++it; it != m_params.end(); ++it) {
            logfile << ", "<< it->str(pflag);
        }
        logfile << "\n\n";
    }

    if (m_params_free.size()) {
        logfile << "Free parameters (" << m_params_free.size() << "):\n";
        it = m_params_free.begin();
        logfile << tab << it->str(pflag);
        for (++it; it != m_params_free.end(); ++it) {
            logfile << ", "<< it->str(pflag);
        }
        logfile << "\n\n";
    }

    if (m_params_calibr.size()) {
        logfile << "Calibrated parameters (" << m_params_calibr.size() << "):\n";
        it = m_params_calibr.begin();
        logfile << tab << it->str(pflag);
        for (++it; it != m_params_calibr.end(); ++it) {
            logfile << ", "<< it->str(pflag);
        }
        logfile << "\n\n";
    }

    unsigned i;
    if (m_eqs.size()) {
        logfile << "Equations (" << m_eqs.size() << "):\n";
        for (it = m_eqs.begin(), i = 1; it != m_eqs.end(); ++it, ++i) {
            logfile << " (" << i << ")  " << it->str(pflag) << " = 0\n";
        }
        logfile << "\n";
    }

    if (!m_static && m_ss.size()) {
        logfile << "Steady state equations (" << m_eqs.size() << "):\n";
        for (itv = m_ss.begin(), i = 1; itv != m_ss.end(); ++itv, ++i) {
            logfile << " (" << i << ")  " << itv->str(pflag) << " = 0\n";
        }
        logfile << "\n";
    }

    if (m_calibr.size()) {
        logfile << "Calibrating equations (" << m_calibr.size() << "):\n";
        for (it = m_calibr.begin(), i = 1; it != m_calibr.end(); ++it, ++i) {
            logfile << " (" << i << ")  " << it->str(pflag) << " = 0\n";
        }
        logfile << "\n";
    }

    if (m_params_free_set.size()) {
        logfile << "Parameter settings (" << m_params_free_set.size() << "):\n";
        map_ex_ex::const_iterator itm;
        for (itm = m_params_free_set.begin(), i = 1; itm != m_params_free_set.end(); ++itm, ++i) {
            logfile << " (" << i << ")  " << itm->first.str(pflag) << " = "
                    << itm->second.str(pflag) << "\n";
        }
        logfile << "\n";
    }

    if (!errors() && !warnings()) return;
    logfile << "\n" << errwarn << "\n";
    if (errors()) logfile << get_errs(true) << "\n";
    if (warnings()) logfile << get_warns(true) << "\n";

    if (!logfile.good()) {
        report_errors("(gEcon error): failed writing logfile \'" + full_name
                      + "\'");
    } else if (m_options[verbose]) {
        write_info("logfile written to \'" + full_name + "\'");
    }
}




namespace {

std::string
set_to_R(const std::string &s)
{
    unsigned i = 0, n = s.length();
    std::string res;
    for (; i < n; ++i) {
        if (s[i] == '{') {
            if (s[i + 2] == '}') {
                res += "c()"; break;
            } else res += "c("; ++i;
        }
        else if (s[i] == '}') { res[res.length() - 1] = ')'; break; }
        else if (s[i] == '\'') { res += '\"'; }
        else res += s[i];
    }
    return res;
}


std::string
duplslash(const std::string &s)
{
    unsigned i = 0, n = s.length();
    std::string res;
    res.reserve(2 * n);
    for (i = 0; i < n; ++i) {
        char c = s[i];
        res += c;
        if (c == '\\') res += c;
    }
    return res;
}



} /* namespace */



void
Model::write_r() const
{
    std::string full_name = m_path + m_name + ".model.R";
#ifdef DEBUG
    std::cerr << "DEBUG INFO: writing R code to file: \'" << full_name << "\'\n";
#endif /* DEBUG */
    std::ofstream R(full_name.c_str());
    if (!R.good()) {
        report_errors("(gEcon error): failed opening R file \'" + full_name
                      + "\' for writing");
    }
    
    if ((m_eqs.size() > 100) && (!m_options[output_r_rcpp])) {
        report_warns("(gEcon warning): the number of variables / equations exceeds 100 and Rcpp \
output option is set to false - this may lead to long compilation of model equations (by JIT compiler); \
if this is the case, consider using option \
\"output R Rcpp = TRUE\"");
    }
    
    set_ex::const_iterator it;
    vec_ex::const_iterator itv;
    std::string tab("    ");
    unsigned index;
    unsigned n_v = m_vars.size(), n_s = m_shocks.size(),
             n_c = m_params_calibr.size(), n_f = m_params_free.size();

    R << "# " << info_str("# ") << "\n\n# Model name: " << m_name << "\n\n";

    R << "# info\ninfo__ <- c(\"" << m_name << "\", \""
      << m_path << m_name << ".gcn\", \""
      << time_str() << "\", \""
      << (m_options[output_r_rcpp] ? "true": "false") <<"\")\n\n";

    R << "# index sets\n";
    if (m_sets.size()) {
        R << "index_sets__ <- list(";
        std::map<std::string, symbolic::idx_set>::const_iterator it;
        for (it = m_sets.begin();;) {
            R << set_to_R(it->second.str());
            if (++it != m_sets.end()) {
                R << ",\n                     ";
            } else {
                R << ")\n\n";
                break;
            }
        }
    } else {
        R << "index_sets__ <- list()\n\n";
    }

    R << "# variables\n";
    R << "variables__ <- c(";
    it = m_vars.begin();
    R << '\"' << it->str(DROP_T | CONVERT_IDX) << '\"';
    for (++it; it != m_vars.end(); ++it) {
        R << ",\n                 \"" << it->str(DROP_T | CONVERT_IDX) << '\"';
    }
    R << ")\n\n";
    R << "variables_tex__ <- c(";
    it = m_vars.begin();
    R << '\"' << duplslash(it->tex(DROP_T)) << '\"';
    for (++it; it != m_vars.end(); ++it) {
        R << ",\n                     \"" << duplslash(it->tex(DROP_T)) << '\"';
    }
    R << ")\n\n";

    R << "# shocks\n";
    if (n_s) {
        R << "shocks__ <- c(";
        it = m_shocks.begin();
        R << '\"' << it->str(DROP_T | CONVERT_IDX) << '\"';
        for (++it; it != m_shocks.end(); ++it) {
            R << ",\n              \"" << it->str(DROP_T | CONVERT_IDX) << '\"';
        }
        R << ")\n\n";
        R << "shocks_tex__ <- c(";
        it = m_shocks.begin();
        R << '\"' << duplslash(it->tex(DROP_T)) << '\"';
        for (++it; it != m_shocks.end(); ++it) {
            R << ",\n                  \"" << duplslash(it->tex(DROP_T)) << '\"';
        }
        R << ")\n\n";
    } else {
        R << "shocks__ <- character(0)\n\n";
        R << "shocks_tex__ <- character(0)\n\n";
    }

    R << "# parameters\n";
    if (m_params.size()) {
        R << "parameters__ <- c(";
        it = m_params.begin();
        R << '\"' << it->str(CONVERT_IDX) << '\"';
        for (++it; it != m_params.end(); ++it) {
            R << ",\n                  \"" << it->str(CONVERT_IDX) << '\"';
        }
        R << ")\n\n";
        R << "parameters_tex__ <- c(";
        it = m_params.begin();
        R << '\"' << duplslash(it->tex()) << '\"';
        for (++it; it != m_params.end(); ++it) {
            R << ",\n                     \"" << duplslash(it->tex()) << '\"';
        }
        R << ")\n\n";
    } else {
        R << "parameters__ <- character(0)\n\n";
        R << "parameters_tex__ <- character(0)\n\n";
    }

    R << "# free parameters\n";
    if (m_params_free.size()) {
    R << "parameters_free__ <- c(";
        it = m_params_free.begin();
        R << '\"' << it->str(CONVERT_IDX) << '\"';
        for (++it; it != m_params_free.end(); ++it) {
            R << ",\n                       \"" << it->str(CONVERT_IDX) << '\"';
        }
        R << ")\n\n";
    } else {
        R << "parameters_free__ <- character(0)\n\n";
    }


    R << "# free parameters' values\n";
    if (m_params_free.size()) {
    R << "parameters_free_val__ <- c(";
        it = m_params_free.begin();
        map_ex_ex::const_iterator mit;
        if ((mit = m_params_free_set.find(*it)) != m_params_free_set.end()) {
            R << mit->second.str();
        } else {
            R << "NA";
        }
        for (++it; it != m_params_free.end(); ++it) {
            if ((mit = m_params_free_set.find(*it)) != m_params_free_set.end()) {
                R << ",\n                           " << mit->second.str();
            } else {
                R << ",\n                           NA";
            }
        }
        R << ")\n\n";
    } else {
        R << "parameters_free_val__ <- numeric(0)\n\n";
    }

    int pflag = (m_static) ? DROP_T : DEFAULT;
    R << "# equations\n";
    if (true) {
        R << "equations__ <- c(";
        it = m_eqs.begin();
        R << '\"' << it->str(pflag) << " = 0\"";
        for (++it; it != m_eqs.end(); ++it) {
            R << ",\n                 \"" << it->str(pflag) << " = 0\"";
        }
        R << ")\n\n";
    } else {
        R << "equations__ <- character(0)\n\n";
    }

    R << "# calibrating equations\n";
    if (m_calibr.size()) {
        R << "calibr_equations__ <- c(";
        it = m_calibr.begin();
        R << '\"' << it->str(pflag) << " = 0\"";
        for (++it; it != m_calibr.end(); ++it) {
            R << ",\n                        \"" << it->str(pflag) << " = 0\"";
        }
        R << ")\n\n";
    } else {
        R << "calibr_equations__ <- character(0)\n\n";
    }

    std::map<std::pair<int, int>, unsigned>::const_iterator itim;
    int ii = 1;
    R << "# variables / equations map\n";
    R << "vareqmap__ <- sparseMatrix(i = c(";
    itim = m_var_eq_map.begin();
    R << itim->first.first;
    for (++itim; itim != m_var_eq_map.end(); ++itim, ++ii) {
        if (!(ii % 10)) R << ",\n                                 ";
        else R << ", ";
        R << itim->first.first;
    }
    R << "),\n                           j = c(";
    itim = m_var_eq_map.begin();
    ii = 1;
    R << itim->first.second;
    for (++itim; itim != m_var_eq_map.end(); ++itim, ++ii) {
        if (!(ii % 10)) R << ",\n                                 ";
        else R << ", ";
        R << itim->first.second;
    }
    R << "),\n                           x = c(";
    itim = m_var_eq_map.begin();
    ii = 1;
    R << itim->second;
    for (++itim; itim != m_var_eq_map.end(); ++itim, ++ii) {
        if (!(ii % 10)) R << ",\n                                 ";
        else R << ", ";
        R << itim->second;
    }
    R << "),\n                           dims = c(" << n_v << ", " << n_v << "))\n\n";

    std::set<std::pair<int, int> >::const_iterator itis;
    R << "# variables / calibrating equations map\n";
    if (m_var_ceq_map.size()) {
        R << "varcalibreqmap__ <- sparseMatrix(i = c(";
        itis = m_var_ceq_map.begin();
        ii = 1;
        R << itis->first;
        for (++itis; itis != m_var_ceq_map.end(); ++itis, ++ii) {
            if (!(ii % 10)) R << ",\n                                       ";
            else R << ", ";
            R << itis->first;
        }
        R << "),\n                                 j = c(";
        itis = m_var_ceq_map.begin();
        R << itis->second;
        ii = 1;
        for (++itis; itis != m_var_ceq_map.end(); ++itis, ++ii) {
            if (!(ii % 10)) R << ",\n                                       ";
            else R << ", ";
            R << itis->second;
        }
        R << "),\n                                 x = rep(1, " << m_var_ceq_map.size();
        R << "), dims = c(" << n_c << ", " << n_v << "))\n\n";
    } else {
        R << "varcalibreqmap__ <- sparseMatrix(i = NULL, j = NULL, dims = c(";
        R << n_c << ", " << n_v << "))\n\n";
    }

    R << "# calibrated parameters / equations map\n";
    if (m_cpar_eq_map.size()) {
        R << "calibrpareqmap__ <- sparseMatrix(i = c(";
        itis = m_cpar_eq_map.begin();
        R << itis->first;
        ii = 1;
        for (++itis; itis != m_cpar_eq_map.end(); ++itis, ++ii) {
            if (!(ii % 10)) R << ",\n                                       ";
            else R << ", ";
            R << itis->first;
        }
        R << "),\n                                 j = c(";
        itis = m_cpar_eq_map.begin();
        R << itis->second;
        ii = 1;
        for (++itis; itis != m_cpar_eq_map.end(); ++itis, ++ii) {
            if (!(ii % 10)) R << ",\n                                       ";
            else R << ", ";
            R << itis->second;
        }
        R << "),\n                                 x = rep(1, " << m_cpar_eq_map.size();
        R << "), dims = c(" << n_v << ", " << n_c << "))\n\n";
    } else {
        R << "calibrpareqmap__ <- sparseMatrix(i = NULL, j = NULL, dims = c(";
        R << n_v << ", " << n_c << "))\n\n";
    }

    R << "# calibrated parameters / calibrating equations map\n";
    if (m_cpar_ceq_map.size()) {
        R << "calibrparcalibreqmap__ <- sparseMatrix(i = c(";
        itis = m_cpar_ceq_map.begin();
        R << itis->first;
        ii = 1;
        for (++itis; itis != m_cpar_ceq_map.end(); ++itis, ++ii) {
            if (!(ii % 10)) R << ",\n                                             ";
            else R << ", ";
            R << itis->first;
        }
        R << "),\n                                       j = c(";
        itis = m_cpar_ceq_map.begin();
        R << itis->second;
        ii = 1;
        for (++itis; itis != m_cpar_ceq_map.end(); ++itis, ++ii) {
            if (!(ii % 10)) R << ",\n                                             ";
            else R << ", ";
            R << itis->second;
        }
        R << "),\n                                       x = rep(1, " << m_cpar_ceq_map.size();
        R << "), dims = c(" << n_c << ", " << n_c << "))\n\n";
    } else {
        R << "calibrparcalibreqmap__ <- sparseMatrix(i = NULL, j = NULL, dims = c(";
        R << n_c << ", " << n_c << "))\n\n";
    }

    R << "# free parameters / equations map\n";
    if (m_fpar_eq_map.size()) {
        R << "freepareqmap__ <- sparseMatrix(i = c(";
        itis = m_fpar_eq_map.begin();
        R << itis->first;
        ii = 1;
        for (++itis; itis != m_fpar_eq_map.end(); ++itis, ++ii) {
            if (!(ii % 10)) R << ",\n                                     ";
            else R << ", ";
            R << itis->first;
        }
        R << "),\n                               j = c(";
        itis = m_fpar_eq_map.begin();
        R << itis->second;
        ii = 1;
        for (++itis; itis != m_fpar_eq_map.end(); ++itis, ++ii) {
            if (!(ii % 10)) R << ",\n                                     ";
            else R << ", ";
            R << itis->second;
        }
        R << "),\n                               x = rep(1, " << m_fpar_eq_map.size();
        R << "), dims = c(" << n_v << ", " << n_f << "))\n\n";
    } else {
        R << "freepareqmap__ <- sparseMatrix(i = NULL, j = NULL, dims = c(";
        R << n_v << ", " << n_f << "))\n\n";
    }

    R << "# free parameters / calibrating equations map\n";
    if (m_fpar_ceq_map.size()) {
        R << "freeparcalibreqmap__ <- sparseMatrix(i = c(";
        itis = m_fpar_ceq_map.begin();
        R << itis->first;
        ii = 1;
        for (++itis; itis != m_fpar_ceq_map.end(); ++itis, ++ii) {
            if (!(ii % 10)) R << ",\n                                           ";
            else R << ", ";
            R << itis->first;
        }
        R << "),\n                                     j = c(";
        itis = m_fpar_ceq_map.begin();
        R << itis->second;
        ii = 1;
        for (++itis; itis != m_fpar_ceq_map.end(); ++itis, ++ii) {
            if (!(ii % 10)) R << ",\n                                           ";
            else R << ", ";
            R << itis->second;
        }
        R << "),\n                                     x = rep(1, " << m_fpar_ceq_map.size();
        R << "), dims = c(" << n_c << ", " << n_f << "))\n\n";
    } else {
        R << "freeparcalibreqmap__ <- sparseMatrix(i = NULL, j = NULL, dims = c(";
        R << n_c << ", " << n_f << "))\n\n";
    }

    R << "# shocks / equations map\n";
    if (m_shock_eq_map.size()) {
        R << "shockeqmap__ <- sparseMatrix(i = c(";
        itis = m_shock_eq_map.begin();
        R << itis->first;
        ii = 1;
        for (++itis; itis != m_shock_eq_map.end(); ++itis, ++ii) {
            if (!(ii % 10)) R << ",\n                                   ";
            else R << ", ";
            R << itis->first;
        }
        R << "),\n                             j = c(";
        itis = m_shock_eq_map.begin();
        R << itis->second;
        ii = 1;
        for (++itis; itis != m_shock_eq_map.end(); ++itis, ++ii) {
            if (!(ii % 10)) R << ",\n                                   ";
            else R << ", ";
            R << itis->second;
        }
        R << "),\n                             x = rep(1, " << m_shock_eq_map.size();
        R << "), dims = c(" << n_v << ", " << n_s << "))\n\n";
    } else {
        R << "shockeqmap__ <- sparseMatrix(i = NULL, j = NULL, dims = c(";
        R << n_v << ", " << n_s << "))\n\n";
    }

    // str map
    map_str_str mss;

    int index_base = 1;
    std::string line_end = "\n";
    if (m_options[output_r_rcpp]) line_end = ";" + line_end;

    if (!m_options[output_r_long]) {
        // index map
        for (it = m_vars.begin(), index = 1; it != m_vars.end(); ++it, ++index) {
            mss[it->str(CONVERT_IDX | DROP_T)] = "v[" + num2str(index) + "]";
        }
        for (it = m_params_calibr.begin(), index = 1; it != m_params_calibr.end(); ++it, ++index) {
            mss[it->str(CONVERT_IDX)] = "pc[" + num2str(index) + "]";
        }
        for (it = m_params_free.begin(), index = 1; it != m_params_free.end(); ++it, ++index) {
            mss[it->str(CONVERT_IDX)] = "pf[" + num2str(index) + "]";
        }
    }

    // replace illegal characters from model name with underscore
    std::regex non_identifier_char ("[\\W]");
    std::string fun_name_suffix = std::regex_replace(m_name, non_identifier_char, "_");

    // settings for printing R functions
    pflag = (m_static) ? DROP_T : CONVERT_T;
    pflag |= CONVERT_IDX;

    R << "# steady state equations\n";
    if (m_options[output_r_rcpp]) {
        R << "ss_eq_cpp_code <- \'\n\n";
        R << "#include <Rcpp.h>\n\n";
        R << "Rcpp::NumericVector ss_eq__" << fun_name_suffix << "(const Rcpp::NumericVector &,  const Rcpp::NumericVector &, const Rcpp::NumericVector &);\n";
        R << "extern \"C\" void ss_eq_cpp_c_" << fun_name_suffix << "(const double*, const double*, const double*, double*);\n\n";

        R << "// [[Rcpp::export]]\n";
        R << "Rcpp::NumericVector   ss_eq__" << fun_name_suffix << "(const Rcpp::NumericVector &v, const Rcpp::NumericVector &pc, const Rcpp::NumericVector &pf)\n"; 
        R << "{\n";
        R << tab << "Rcpp::NumericVector r(" << m_eqs.size() << ");\n";
        R << tab << "ss_eq_cpp_c_" << fun_name_suffix << "(&v[-1], &pc[-1], &pf[-1], &r[-1]);\n";
        R << tab << "return r;\n";
        R << "}\n\n";

        R << "extern \"C\" {\n";
        R << "#include <math.h>\n\n";
        R << "void ss_eq_cpp_c_" << fun_name_suffix << "(const double *v, const double *pc, const double *pf, double *r)\n";
        R << "{\n";
    } else {
        R << "ss_eq__ <- function(v, pc, pf)\n";
        R << "{\n";
        R << tab << "r <- numeric(" << m_eqs.size() << ")\n";
    }

    if (m_options[output_r_long]) {
        if (m_static) {
            for (it = m_vars.begin(), index = 1; it != m_vars.end(); ++it, ++index) {
                R << tab << it->str(pflag) << " = v[" << index << "]\n" ;
            }
        } else {
            for (it = m_vars.begin(), index = 1; it != m_vars.end(); ++it, ++index) {
                R << tab << ss(*it).str(pflag) << " = v[" << index << "]\n" ;
            }
        }
        if (m_vars.size()) R << '\n';
        for (it = m_params_calibr.begin(), index = 1; it != m_params_calibr.end(); ++it, ++index) {
            R << tab << it->str(pflag) << " = pc[" << index << "]\n" ;
        }
        if (m_params_calibr.size()) R << '\n';
        for (it = m_params_free.begin(), index = 1; it != m_params_free.end(); ++it, ++index) {
            R << tab << it->str(pflag) << " = pf[" << index << "]\n" ;
        }
        if (m_params_free.size()) R << '\n';
        if (m_static) {
            for (it = m_eqs.begin(), index = 1; it != m_eqs.end(); ++it, ++index) {
                R << tab << "r[" << index << "] = " << it->str(pflag) << "\n" ;
            }
        } else {
            for (itv = m_ss.begin(), index = 1; itv != m_ss.end(); ++itv, ++index) {
                R << tab << "r[" << index << "] = " << itv->str(pflag) << "\n" ;
            }
        }
    } else {
        if (m_static) {
            for (it = m_eqs.begin(), index = index_base; it != m_eqs.end(); ++it, ++index) {
                R << tab << "r[" << index << "] = " << it->strmap(mss, m_options[output_r_rcpp]) << line_end;
            }
        } else {
            for (itv = m_ss.begin(), index = index_base; itv != m_ss.end(); ++itv, ++index) {
                R << tab << "r[" << index << "] = " << itv->strmap(mss, m_options[output_r_rcpp]) << line_end;
            }
        }
    }

    if (m_options[output_r_rcpp]) {
        R << "}\n" << "}\n";
        R << "'\n\n";
        R << "sourceCpp(code = ss_eq_cpp_code)\n";
        R << "assign(x = \"ss_eq__" << fun_name_suffix << "\", value = ss_eq__" 
          << fun_name_suffix << ", envir = .gEcon_env)\n";
        R << "rm(\"ss_eq__" << fun_name_suffix << "\", envir = globalenv())\n\n";
        R << "ss_eq__ <- function(v, pc, pf)\n";
        R << "{\n";
        R << tab << "return(get(x = \"ss_eq__" << fun_name_suffix << "\", envir = .gEcon_env)(v, pc, pf))\n";
        R << "}\n\n";
    } else {
        R << "\n" << tab << "return(r)" << "\n}\n\n";
    }

    R << "# calibrating equations\n";
    if (m_options[output_r_rcpp]) {
        R << "calibr_eq_cpp_code <-  \'\n\n";
        R << "#include <Rcpp.h>\n\n";
        R << "Rcpp::NumericVector calibr_eq__" << fun_name_suffix << "(const Rcpp::NumericVector &,  const Rcpp::NumericVector &, const Rcpp::NumericVector &);\n";
        R << "extern \"C\" void calibr_eq_cpp_c_" << fun_name_suffix << "(const double*, const double*, const double*, double*);\n\n";

        R << "// [[Rcpp::export]]\n";
        R << "Rcpp::NumericVector   calibr_eq__" << fun_name_suffix << "(const Rcpp::NumericVector &v, const Rcpp::NumericVector &pc, const Rcpp::NumericVector &pf)\n";
        R << "{\n";
        R << tab << "Rcpp::NumericVector r(" << m_calibr.size() << ");\n";
        R << tab << "calibr_eq_cpp_c_" << fun_name_suffix << "(&v[-1], &pc[-1], &pf[-1], &r[-1]);\n";
        R << tab << "return r;\n";
        R << "}\n\n";

        R << "extern \"C\" {\n";
        R << "#include <math.h>\n\n";
        R << "void calibr_eq_cpp_c_" << fun_name_suffix << "(const double *v, const double *pc, const double *pf, double *r)\n";
        R << "{\n";        
    } else {
        R << "calibr_eq__ <- function(v, pc, pf)\n";
        R << "{\n";
        R << tab << "r <- numeric(" << m_calibr.size() << ")\n";
    }

    if (m_options[output_r_long]) {
        for (it = m_vars.begin(), index = 1; it != m_vars.end(); ++it, ++index) {
            R << tab << ss(*it).str(pflag) << " = v[" << index << "]\n" ;
        }
        if (m_vars.size()) R << '\n';
        for (it = m_params_calibr.begin(), index = 1; it != m_params_calibr.end(); ++it, ++index) {
            R << tab << it->str(pflag) << " = pc[" << index << "]\n" ;
        }
        if (m_params_calibr.size()) R << '\n';
        for (it = m_params_free.begin(), index = 1; it != m_params_free.end(); ++it, ++index) {
            R << tab << it->str(pflag) << " = pf[" << index << "]\n" ;
        }
        if (m_params_free.size()) R << '\n';
        for (it = m_calibr.begin(), index = 1; it != m_calibr.end(); ++it, ++index) {
            R << tab << "r[" << index << "] = " << it->str(pflag) << "\n" ;
        }
    } else {
        for (it = m_calibr.begin(), index = index_base; it != m_calibr.end(); ++it, ++index) {
            R << tab << "r[" << index << "] = " << it->strmap(mss, m_options[output_r_rcpp]) << line_end;
        }
    }

    if (m_options[output_r_rcpp]) {
        R << "}\n" << "}\n";
        R << "'\n\n";
        R << "sourceCpp(code = calibr_eq_cpp_code)\n";
        R << "assign(x = \"calibr_eq__" << fun_name_suffix << "\", value = calibr_eq__" 
          << fun_name_suffix << ", envir = .gEcon_env)\n";
        R << "rm(\"calibr_eq__" << fun_name_suffix << "\", envir = globalenv())\n\n";
        R << "calibr_eq__ <- function(v, pc, pf)\n";
        R << "{\n";
        R << tab << "return(get(x = \"calibr_eq__" << fun_name_suffix << "\", envir = .gEcon_env)(v, pc, pf))\n";
        R << "}\n\n";    
    } else {
        R << '\n' << tab << "return(r)" << "\n}\n\n";
    }

    int n = m_eqs.size() + m_calibr.size();

    std::map<std::pair<int, int>, ex>::const_iterator itm;
    if (m_options[output_r_jacobian]) {
        R << "# steady state and calibrating equations Jacobian\n";
        if (m_options[output_r_rcpp]) {
            R << "ss_calibr_eq_jacob_cpp_code <- \'\n\n";
            R << "#include <Rcpp.h>\n\n";
            R << "Rcpp::NumericVector ss_calibr_eq_jacob_cpp__" << fun_name_suffix << "(const Rcpp::NumericVector &,  const Rcpp::NumericVector &, const Rcpp::NumericVector &);\n";
            R << "extern \"C\" void ss_calibr_eq_jacob_cpp_c_" << fun_name_suffix << "(const double*, const double*, const double*, double*);\n\n";

            R << "// [[Rcpp::export]]\n";
            R << "Rcpp::NumericVector   ss_calibr_eq_jacob_cpp__" << fun_name_suffix << "(const Rcpp::NumericVector &v, const Rcpp::NumericVector &pc, const Rcpp::NumericVector &pf)\n";
            R << "{\n";
            R << tab << "Rcpp::NumericVector jac(" << m_jacob_ss_calibr.size() << ");\n";
            R << tab << "ss_calibr_eq_jacob_cpp_c_" << fun_name_suffix << "(&v[-1], &pc[-1], &pf[-1], &jac[-1]);\n";
            R << tab << "return jac;\n";
            R << "}\n\n";

            R << "extern \"C\" {\n";
            R << "#include <math.h>\n\n";
            R << "void ss_calibr_eq_jacob_cpp_c_" << fun_name_suffix << "(const double *v, const double *pc, const double *pf, double *jac)\n";
            R << "{\n";
        } else {
            R << "ss_calibr_eq_jacob__ <- function(v, pc, pf)\n";
            R << "{\n";
            R << tab << "r <- numeric(" << m_calibr.size() << ")\n";
        }


        if (m_options[output_r_long]) {
            for (it = m_vars.begin(), index = 1; it != m_vars.end(); ++it, ++index) {
                R << tab << ss(*it).str(pflag) << " = v[" << index << "]\n" ;
            }
            if (m_vars.size()) R << '\n';
            for (it = m_params_calibr.begin(), index = 1; it != m_params_calibr.end(); ++it, ++index) {
                R << tab << it->str(pflag) << " = pc[" << index << "]\n" ;
            }
            if (m_params_calibr.size()) R << '\n';
            for (it = m_params_free.begin(), index = 1; it != m_params_free.end(); ++it, ++index) {
                R << tab << it->str(pflag) << " = pf[" << index << "]\n" ;
            }
            if (m_params_free.size()) R << '\n';
            R << tab << "jacob <- Matrix(0, nrow = "
              << (m_eqs.size() + m_calibr.size()) << ", ncol = "
              << (m_vars.size() + m_params_calibr.size()) << ", sparse = TRUE)\n";
            for (itm = m_jacob_ss_calibr.begin(); itm != m_jacob_ss_calibr.end(); ++itm) {
                R << tab << "jacob[" << itm->first.first << ", " << itm ->first.second
                << "] = " << itm->second.str(pflag) << '\n';
            }
            R << "\n" << tab << "return(jacob)" << "\n}\n\n";
        } else {
            int i = index_base;
            if (!m_options[output_r_rcpp]) {
                R << tab << "jac <- numeric(" << m_jacob_ss_calibr.size() << ")\n";
            }
            for (itm = m_jacob_ss_calibr.begin(); itm != m_jacob_ss_calibr.end(); ++itm, ++i) {
                R << tab << "jac[" << i << "] = " << itm->second.strmap(mss, m_options[output_r_rcpp]) << line_end;
            }

            if (m_options[output_r_rcpp]) {
                R << "}\n" << "}\n";
                R << "\'\n\n";
                R << "sourceCpp(code = ss_calibr_eq_jacob_cpp_code)\n";
                R << "assign(x = \"ss_calibr_eq_jacob_cpp__" << fun_name_suffix << "\", value = ss_calibr_eq_jacob_cpp__" 
                  << fun_name_suffix << ", envir = .gEcon_env)\n";
                R << "rm(\"ss_calibr_eq_jacob_cpp__" << fun_name_suffix << "\", envir = globalenv())\n\n";
                R << "ss_calibr_eq_jacob_cpp__ <- function(v, pc, pf)\n";
                R << "{\n";
                R << tab << "return(get(x = \"ss_calibr_eq_jacob_cpp__" << fun_name_suffix << "\", envir = .gEcon_env)(v, pc, pf))\n";
                R << "}\n\n"; 
        
            } else {
                R << tab << "jacob <- sparseMatrix(i = c(";
                itm = m_jacob_ss_calibr.begin();
                R << itm->first.first;
                ii = 1;
                for (++itm; itm != m_jacob_ss_calibr.end(); ++itm, ++ii) {
                    if (!(ii % 10)) R << ",\n                                ";
                    else R << ", ";
                    R << itm->first.first;
                }
                R << "),\n                          j = c(";
                itm = m_jacob_ss_calibr.begin();
                R << itm->first.second;
                ii = 1;
                for (++itm; itm != m_jacob_ss_calibr.end(); ++itm, ++ii) {
                    if (!(ii % 10)) R << ",\n                                ";
                    else R << ", ";
                    R << itm->first.second;
                }
                R << "),\n                          x = jac, dims = c(" << n << ", " << n << "))\n";
                R << "\n" << tab << "return(jacob)" << "\n}\n\n";
            }
        }
    } else {
        R << "\n";
        R << "ss_calibr_eq_jacob__ <- NULL\n\n\n";
    }

    if (m_options[output_r_rcpp]) {
        if (m_options[output_r_jacobian]) {
            R << "ss_calibr_eq_jacob__ <- function(v, pc, pf)\n";
            R << "{\n";
            R << tab << "jacob <- sparseMatrix(i = c(";
            itm = m_jacob_ss_calibr.begin();
            R << itm->first.first;
            ii = 1;
            for (++itm; itm != m_jacob_ss_calibr.end(); ++itm, ++ii) {
                if (!(ii % 10)) R << ",\n                                ";
                else R << ", ";
                R << itm->first.first;
            }
            R << "),\n                          j = c(";
            itm = m_jacob_ss_calibr.begin();
            R << itm->first.second;
            ii = 1;
            for (++itm; itm != m_jacob_ss_calibr.end(); ++itm, ++ii) {
                if (!(ii % 10)) R << ",\n                                ";
                else R << ", ";
                R << itm->first.second;
            }
            R << "),\n                          x = get(x = \"ss_calibr_eq_jacob_cpp__" << fun_name_suffix 
              << "\", envir = .gEcon_env)(v, pc, pf), dims = c(" << n << ", " << n << "))\n";
            R << "\n" << tab << "return(jacob)" << "\n}\n\n";
        }
    }

    R << "# 1st order perturbation\n";
    if (m_options[output_r_rcpp]) {
        R << "pert_cpp_code <-  \'\n\n";
        R << "#include <Rcpp.h>\n\n";


        R << "Rcpp::NumericVector m_Atp1_fun_cpp_" << fun_name_suffix << "(const Rcpp::NumericVector &,  const Rcpp::NumericVector &, const Rcpp::NumericVector &);\n";
        R << "extern \"C\" void m_Atp1_fun_c_" << fun_name_suffix << "(const double*, const double*, const double*, double*);\n\n";

        R << "// [[Rcpp::export]]\n";
        R << "Rcpp::NumericVector   m_Atp1_fun_cpp_" << fun_name_suffix << "(const Rcpp::NumericVector &v, const Rcpp::NumericVector &pc, const Rcpp::NumericVector &pf)\n";
        R << "{\n";
        R << tab << "Rcpp::NumericVector Atpl(" << m_Atp1.size() << ");\n";
        R << tab << "m_Atp1_fun_c_" << fun_name_suffix << "(&v[-1], &pc[-1], &pf[-1], &Atpl[-1]);\n";
        R << tab << "return Atpl;\n";
        R << "}\n\n";

        R << "extern \"C\" {\n";
        R << "#include <math.h>\n\n";
        R << "void m_Atp1_fun_c_" << fun_name_suffix << "(const double *v, const double *pc, const double *pf, double *Atp1x)\n";
        R << "{\n";
        int i = index_base;
        for (itm = m_Atp1.begin(); itm != m_Atp1.end(); ++itm, ++i) {
            R << tab << "Atp1x[" << i << "] = " << itm->second.strmap(mss, m_options[output_r_rcpp]) << line_end;
        }
        R << "}\n" << "}\n";


        R << "Rcpp::NumericVector m_Atm1_fun_cpp_" << fun_name_suffix << "(const Rcpp::NumericVector &,  const Rcpp::NumericVector &, const Rcpp::NumericVector &);\n";
        R << "extern \"C\" void m_Atm1_fun_c_" << fun_name_suffix << "(const double*, const double*, const double*, double*);\n\n";

        R << "// [[Rcpp::export]]\n";
        R << "Rcpp::NumericVector   m_Atm1_fun_cpp_" << fun_name_suffix << "(const Rcpp::NumericVector &v, const Rcpp::NumericVector &pc, const Rcpp::NumericVector &pf)\n";
        R << "{\n";
        R << tab << "Rcpp::NumericVector Atm1(" << m_Atm1.size() << ");\n";
        R << tab << "m_Atm1_fun_c_" << fun_name_suffix << "(&v[-1], &pc[-1], &pf[-1], &Atm1[-1]);\n";
        R << tab << "return Atm1;\n";
        R << "}\n\n";

        R << "extern \"C\" {\n";
        R << "#include <math.h>\n\n";
        R << "void m_Atm1_fun_c_" << fun_name_suffix << "(const double *v, const double *pc, const double *pf, double *Atm1x)\n";
        R << "{\n";
        i = index_base;
        for (itm = m_Atm1.begin(); itm != m_Atm1.end(); ++itm, ++i) {
            R << tab << "Atm1x[" << i << "] = " << itm->second.strmap(mss, m_options[output_r_rcpp]) << line_end;
        }
        R << "}\n" << "}\n";

        R << "Rcpp::NumericVector m_At_fun_cpp_" << fun_name_suffix << "(const Rcpp::NumericVector &,  const Rcpp::NumericVector &, const Rcpp::NumericVector &);\n";
        R << "extern \"C\" void m_At_fun_c_" << fun_name_suffix << "(const double*, const double*, const double*, double*);\n\n";

        R << "// [[Rcpp::export]]\n";
        R << "Rcpp::NumericVector   m_At_fun_cpp_" << fun_name_suffix << "(const Rcpp::NumericVector &v, const Rcpp::NumericVector &pc, const Rcpp::NumericVector &pf)\n";
        R << "{\n";
        R << tab << "Rcpp::NumericVector At(" << m_At.size() << ");\n";
        R << tab << "m_At_fun_c_" << fun_name_suffix << "(&v[-1], &pc[-1], &pf[-1], &At[-1]);\n";
        R << tab << "return At;\n";
        R << "}\n\n";

        R << "extern \"C\" {\n";
        R << "#include <math.h>\n\n";
        R << "void m_At_fun_c_" << fun_name_suffix << "(const double *v, const double *pc, const double *pf, double *Atx)\n";
        R << "{\n";
        i = index_base;
        for (itm = m_At.begin(); itm != m_At.end(); ++itm, ++i) {
            R << tab << "Atx[" << i << "] = " << itm->second.strmap(mss, m_options[output_r_rcpp]) << line_end;
        }
        R << "}\n" << "}\n";
        R << "\'\n\n";
        R << "sourceCpp(code = pert_cpp_code)\n";
        R << "assign(x = " << "\"m_Atm1_fun_cpp_" << fun_name_suffix << "\", value = m_Atm1_fun_cpp_"  << fun_name_suffix 
          << ", envir = .gEcon_env)\n";
        R << "assign(x = " << "\"m_Atp1_fun_cpp_" << fun_name_suffix << "\", value = m_Atp1_fun_cpp_"  << fun_name_suffix 
          << ", envir = .gEcon_env)\n";
        R << "assign(x = " << "\"m_At_fun_cpp_" << fun_name_suffix << "\", value = m_At_fun_cpp_"  << fun_name_suffix 
          << ", envir = .gEcon_env)\n";
         R << "rm(list = c(" 
            << "\"m_Atm1_fun_cpp_" << fun_name_suffix
            << "\",\n            " << "\"m_Atp1_fun_cpp_" << fun_name_suffix
            << "\",\n            " << "\"m_At_fun_cpp_" << fun_name_suffix
            << "\"),\n   envir = globalenv())\n\n";
    }
    
    R << "pert1__ <- function(v, pc, pf)\n";
    R << "{\n";

    if (m_options[output_r_long]) {
        for (it = m_vars.begin(), index = 1; it != m_vars.end(); ++it, ++index) {
            R << tab << ss(*it).str(pflag) << " = v[" << index << "]\n" ;
        }
        if (m_vars.size()) R << '\n';
        for (it = m_params_calibr.begin(), index = 1; it != m_params_calibr.end(); ++it, ++index) {
            R << tab << it->str(pflag) << " = pc[" << index << "]\n" ;
        }
        if (m_params_calibr.size()) R << '\n';
        for (it = m_params_free.begin(), index = 1; it != m_params_free.end(); ++it, ++index) {
            R << tab << it->str(pflag) << " = pf[" << index << "]\n" ;
        }
        if (m_params_free.size()) R << '\n';
    }

    if (m_options[output_r_long] || (m_Atm1.size() == 0)) {
        R << tab << "Atm1 <- Matrix(0, nrow = " << n_v << ", ncol = " << n_v
          << ", sparse = TRUE)\n";
        for (itm = m_Atm1.begin(); itm != m_Atm1.end(); ++itm) {
            R << tab << "Atm1[" << itm->first.first << ", " << itm ->first.second
            << "] = " << itm->second.str(pflag) << '\n';
        }
    } else {
        if (!m_options[output_r_rcpp]) {
            R << tab << "Atm1x <- numeric(" << m_Atm1.size() << ")\n";
            int i = 1;
            for (itm = m_Atm1.begin(); itm != m_Atm1.end(); ++itm, ++i) {
                R << tab << "Atm1x[" << i << "] = " << itm->second.strmap(mss) << '\n';
            }
        }
        R << tab << "Atm1 <- sparseMatrix(i = c(";
        itm = m_Atm1.begin();
        R << itm->first.first;
        ii = 1;
        for (++itm; itm != m_Atm1.end(); ++itm, ++ii) {
            if (!(ii % 10)) R << ",\n                               ";
            else R << ", ";
            R << itm->first.first;
        }
        R << "),\n                         j = c(";
        itm = m_Atm1.begin();
        R << itm->first.second;
        ii = 1;
        for (++itm; itm != m_Atm1.end(); ++itm, ++ii) {
            if (!(ii % 10)) R << ",\n                               ";
            else R << ", ";
            R << itm->first.second;
        }
        if (!m_options[output_r_rcpp]) {
            R << "),\n                         x = Atm1x, dims = c(" << n_v 
              << ", " << n_v << "))\n";
        } else {
            R << "),\n                         x = get(\"m_Atm1_fun_cpp_" << fun_name_suffix 
              << "\", envir = .gEcon_env)(v, pc, pf), dims = c(" << n_v << ", " << n_v << "))\n";
        }
    }
    R << '\n';

    if (m_options[output_r_long] || (m_At.size() == 0)) {
        R << tab << "At <- Matrix(0, nrow = " << n_v << ", ncol = " << n_v
          << ", sparse = TRUE)\n";
        for (itm = m_At.begin(); itm != m_At.end(); ++itm) {
            R << tab << "At[" << itm->first.first << ", " << itm ->first.second
            << "] = " << itm->second.str(pflag) << '\n';
        }
    } else {
        if (!m_options[output_r_rcpp]) {
            R << tab << "Atx <- numeric(" << m_At.size() << ")\n";
            int i = 1;
            for (itm = m_At.begin(); itm != m_At.end(); ++itm, ++i) {
                R << tab << "Atx[" << i << "] = " << itm->second.strmap(mss) << '\n';
            }
        }
        R << tab << "At <- sparseMatrix(i = c(";
        itm = m_At.begin();
        R << itm->first.first;
        ii = 1;
        for (++itm; itm != m_At.end(); ++itm, ++ii) {
            if (!(ii % 10)) R << ",\n                             ";
            else R << ", ";
            R << itm->first.first;
        }
        R << "),\n                       j = c(";
        itm = m_At.begin();
        R << itm->first.second;
        ii = 1;
        for (++itm; itm != m_At.end(); ++itm, ++ii) {
            if (!(ii % 10)) R << ",\n                             ";
            else R << ", ";
            R << itm->first.second;
        }
        if (!m_options[output_r_rcpp]) {
            R << "),\n                         x = Atx, dims = c(" << n_v 
              << ", " << n_v << "))\n";
        } else {
            R << "),\n                         x = get(\"m_At_fun_cpp_" << fun_name_suffix 
              << "\", envir = .gEcon_env)(v, pc, pf), dims = c(" << n_v << ", " << n_v << "))\n";
        }
    }
    R << '\n';

    if (m_options[output_r_long] || (m_Atp1.size() == 0)) {
        R << tab << "Atp1 <- Matrix(0, nrow = " << n_v << ", ncol = " << n_v
          << ", sparse = TRUE)\n";
        for (itm = m_Atp1.begin(); itm != m_Atp1.end(); ++itm) {
            R << tab << "Atp1[" << itm->first.first << ", " << itm ->first.second
            << "] = " << itm->second.str(pflag) << '\n';
        }
    } else {
        if (!m_options[output_r_rcpp]) {
            R << tab << "Atp1x <- numeric(" << m_Atp1.size() << ")\n";
            int i = 1;
            for (itm = m_Atp1.begin(); itm != m_Atp1.end(); ++itm, ++i) {
                R << tab << "Atp1x[" << i << "] = " << itm->second.strmap(mss) << '\n';
            }
        }
        R << tab << "Atp1 <- sparseMatrix(i = c(";
        itm = m_Atp1.begin();
        R << itm->first.first;
        ii = 1;
        for (++itm; itm != m_Atp1.end(); ++itm, ++ii) {
            if (!(ii % 10)) R << ",\n                               ";
            else R << ", ";
            R << itm->first.first;
        }
        R << "),\n                         j = c(";
        itm = m_Atp1.begin();
        R << itm->first.second;
        ii = 1;
        for (++itm; itm != m_Atp1.end(); ++itm, ++ii) {
            if (!(ii % 10)) R << ",\n                               ";
            else R << ", ";
            R << itm->first.second;
        }
        if (!m_options[output_r_rcpp]) {
            R << "),\n                         x = Atp1x, dims = c(" << n_v 
              << ", " << n_v << "))\n";
        } else {
            R << "),\n                         x = get(\"m_Atp1_fun_cpp_" << fun_name_suffix 
              << "\", envir = .gEcon_env)(v, pc, pf), dims = c(" << n_v << ", " << n_v << "))\n";
        }
    }
    R << '\n';

    if (m_options[output_r_long] || (m_Aeps.size() == 0)) {
        R << tab << "Aeps <- Matrix(0, nrow = " << n_v << ", ncol = " << n_s
          << ", sparse = TRUE)\n";
        for (itm = m_Aeps.begin(); itm != m_Aeps.end(); ++itm) {
            R << tab << "Aeps[" << itm->first.first << ", " << itm ->first.second
            << "] = " << itm->second.str(pflag) << '\n';
        }
    } else {
        R << tab << "Aepsx <- numeric(" << m_Aeps.size() << ")\n";
        int i = 1;
        for (itm = m_Aeps.begin(); itm != m_Aeps.end(); ++itm, ++i) {
            R << tab << "Aepsx[" << i << "] = " << itm->second.strmap(mss) << '\n';
        }
        R << tab << "Aeps <- sparseMatrix(i = c(";
        itm = m_Aeps.begin();
        R << itm->first.first;
        ii = 1;
        for (++itm; itm != m_Aeps.end(); ++itm, ++ii) {
            if (!(ii % 10)) R << ",\n                               ";
            else R << ", ";
            R << itm->first.first;
        }
        R << "),\n                         j = c(";
        itm = m_Aeps.begin();
        R << itm->first.second;
        ii = 1;
        for (++itm; itm != m_Aeps.end(); ++itm, ++ii) {
            if (!(ii % 10)) R << ",\n                               ";
            else R << ", ";
            R << itm->first.second;
        }
        R << "),\n                         x = Aepsx, dims = c(" << n_v 
          << ", " << n_s << "))\n";
    }

    R << '\n' << tab << "return(list(Atm1, At, Atp1, Aeps))" << "\n}\n\n";


    if (m_options[output_r_rcpp]) {
        R << "ext__ <- list(\n";
        R << tab << "ss_eq_cpp_code = ss_eq_cpp_code,\n";
        R << tab << "calibr_eq_cpp_code = calibr_eq_cpp_code,\n";
        if (m_options[output_r_jacobian]) {
            R << tab << "ss_calibr_eq_jacob_cpp_code = ss_calibr_eq_jacob_cpp_code,\n";
        } else {
            R << tab << "ss_calibr_eq_jacob_cpp_code = NULL,\n";            
        }
        R << tab << "pert_cpp_code = pert_cpp_code\n";
        R << ")\n\n";
    } else {
        R << "ext__ <- list()\n\n";
    }

    R << "# create model object\n"
      << "gecon_model(model_info = info__,\n"
      << "            index_sets = index_sets__,\n"
      << "            variables = variables__,\n"
      << "            variables_tex = variables_tex__,\n"
      << "            shocks = shocks__,\n"
      << "            shocks_tex = shocks_tex__,\n"
      << "            parameters = parameters__,\n"
      << "            parameters_tex = parameters_tex__,\n"
      << "            parameters_free = parameters_free__,\n"
      << "            parameters_free_val = parameters_free_val__,\n"
      << "            equations = equations__,\n"
      << "            calibr_equations = calibr_equations__,\n"
      << "            var_eq_map = vareqmap__,\n"
      << "            shock_eq_map = shockeqmap__,\n"
      << "            var_ceq_map = varcalibreqmap__,\n"
      << "            cpar_eq_map = calibrpareqmap__,\n"
      << "            cpar_ceq_map = calibrparcalibreqmap__,\n"
      << "            fpar_eq_map = freepareqmap__,\n"
      << "            fpar_ceq_map = freeparcalibreqmap__,\n"
      << "            ss_function = ss_eq__,\n"
      << "            calibr_function = calibr_eq__,\n"
      << "            ss_calibr_jac_function = ss_calibr_eq_jacob__,\n"
      << "            pert = pert1__,\n"
      << "            ext = ext__)\n";


    if (!R.good()) {
        report_errors("(gEcon error): failed writing R file \'" + full_name
                      + "\'");
    } else if (m_options[verbose]) {
        write_info("R code written to \'" + full_name + "\'");
    }
}


void
Model::write_latex() const
{
    std::string full_name = m_path + m_name + ".tex";
#ifdef DEBUG
    std::cerr << "DEBUG INFO: writing LaTeX documentation to file: \'" << full_name << "\'\n";
#endif /* DEBUG */
    if ((m_eqs.size() > 100) && (m_options[output_latex_long])) {
        report_warns("(gEcon warning): the number of variables / equations exceeds 100 and LaTeX \
output option is set to long - your LaTeX document may become large; if this is the case, consider using option \
\"output LaTeX long = FALSE\"");
    }
    
    std::ofstream latex(full_name.c_str());
    if (!latex.good()) {
        report_errors("(gEcon error): failed opening LaTeX file \'" + full_name
                      + "\' for writing");
    }
    latex << "% " << info_str("% ") << "\n\n% Model name: " << m_name << "\n\n";
    latex << "\\documentclass[10pt,a4paper]{article}\n";
    latex << "\\usepackage[utf8]{inputenc}\n";
    latex << "\\usepackage{color}\n";
    latex << "\\usepackage{graphicx}\n";
    latex << "\\usepackage{epstopdf}\n";
    latex << "\\usepackage{hyperref}\n";
    latex << "\\usepackage{amsmath}\n";
    latex << "\\usepackage{amssymb}\n";
    latex << "\\DeclareMathOperator{\\erf}{erf}\n";
    latex << "\\DeclareMathOperator{\\pnorm}{pnorm}\n";
    if (m_options[output_latex_landscape]) latex << "\\usepackage{pdflscape}\n";
    latex << "\\numberwithin{equation}{section}\n";
    latex << "\\usepackage[top = 2.5cm, bottom=2.5cm, left = 2.0cm, right=2.0cm]{geometry}\n";
    latex << "\\begin{document}\n\n";
    if (m_options[output_latex_landscape]) latex << "\\begin{landscape}\n";
    latex << "\\begin{flushleft}{\\large\n";
    latex << "Generated  on " + time_str() + " by \\href{" + gecon_web_str() +
             "}{\\texttt{gEcon}} version " + gecon_ver_str() + "\\\\\n";
    latex << "Model name: \\verb+" << m_name << "+\n";
    latex << "}\\end{flushleft}\n\n";
    latex << "\\input{" << m_name << ".model.tex}\n";
    if (m_options[output_latex_landscape]) latex << "\\end{landscape}\n";
    latex << "\\input{" << m_name << ".results.tex}\n\n";
    latex << "\\end{document}\n\n";
    if (!latex.good()) {
        report_errors("(gEcon error): failed writing LaTeX file \'" + full_name
                      + "\'");
    }
    latex.close();

    full_name = m_path + m_name + ".results.tex";
    latex.open(full_name.c_str());
    if (!latex.good()) {
        report_errors("(gEcon error): failed opening LaTeX file \'" + full_name
                      + "\' for writing");
    }
    latex << "% " << info_str("% ") << "\n\n% Model name: " << m_name << "\n\n";
    if (!latex.good()) {
        report_errors("(gEcon error): failed writing LaTeX file \'" + full_name
                      + "\'");
    }
    latex.close();

    full_name = m_path + m_name + ".model.tex";
    latex.open(full_name.c_str());
    if (!latex.good()) {
        report_errors("(gEcon error): failed opening LaTeX file \'" + full_name
                      + "\' for writing");
    }
    latex << "% " << info_str("% ") << "\n\n% Model name: " << m_name << "\n\n";

    if (m_sets.size()) {
        latex << "\\section*{Index sets}\n\n";
        std::map<std::string, symbolic::idx_set>::const_iterator it;
        for (it = m_sets.begin(); it != m_sets.end(); ++it) {
            latex << "$$" << it->second.tex() << "$$\n";
        }
        latex << '\n';
    }


    for (unsigned i = 0, n = m_blocks.size(); i < n; ++i) {
        m_blocks[i].write_latex(latex, m_static);
    }

    set_ex::const_iterator it, ite;
    vec_ex::const_iterator itv, itev;
    print_flag pflag = (m_static) ? DROP_T : DEFAULT;

    if (m_sets.size()) {
        latex << "\\section{Equilibrium relationships (before expansion and reduction)}\n\n";
        it = m_t_eqs.begin();
        ite = m_t_eqs.end();
        for (; it != ite; ++it) {
            latex << "\\begin{equation}\n";
            latex << it->tex(pflag) << " = 0\n";
            latex << "\\end{equation}\n";
        }
        latex << "\n\n\n";

        if (m_options[output_latex_long]) {
            latex << "\\section{Equilibrium relationships (after expansion and reduction)}\n\n";
            it = m_eqs.begin();
            ite = m_eqs.end();
            for (; it != ite; ++it) {
                latex << "\\begin{equation}\n";
                latex << it->tex(pflag) << " = 0\n";
                latex << "\\end{equation}\n";
            }
            latex << "\n\n\n";
        }
    } else {
        latex << "\\section{Equilibrium relationships (after reduction)}\n\n";
        it = m_eqs.begin();
        ite = m_eqs.end();
        for (; it != ite; ++it) {
            latex << "\\begin{equation}\n";
            latex << it->tex(pflag) << " = 0\n";
            latex << "\\end{equation}\n";
        }
        latex << "\n\n\n";
    }

    if (!m_static) {
        if (m_sets.size()) {
            latex << "\\section{Steady state relationships (before expansion and reduction)}\n\n";
            itv = m_t_ss.begin();
            itev = m_t_ss.end();
            for (; itv != itev; ++itv) {
                latex << "\\begin{equation}\n";
                latex << itv->tex(pflag) << " = 0\n";
                latex << "\\end{equation}\n";
            }
            latex << "\n\n\n";

            if (m_options[output_latex_long]) {
                latex << "\\section{Steady state relationships (after expansion and reduction)}\n\n";
                itv = m_ss.begin();
                itev = m_ss.end();
                for (; itv != itev; ++itv) {
                    latex << "\\begin{equation}\n";
                    latex << itv->tex(pflag) << " = 0\n";
                    latex << "\\end{equation}\n";
                }
                latex << "\n\n\n";
            }
        } else {
            latex << "\\section{Steady state relationships (after reduction)}\n\n";
            itv = m_ss.begin();
            itev = m_ss.end();
            for (; itv != itev; ++itv) {
                latex << "\\begin{equation}\n";
                latex << itv->tex(pflag) << " = 0\n";
                latex << "\\end{equation}\n";
            }
            latex << "\n\n\n";
        }
    }

    if (m_calibr.size()) {
        latex << "\\section{Calibrating equations}\n\n";
        it = m_calibr.begin();
        ite = m_calibr.end();
        for (; it != ite; ++it) {
            latex << "\\begin{equation}\n";
            latex << it->tex(pflag) << " = 0\n";
            latex << "\\end{equation}\n";
        }
    }
    latex << "\n\n\n";

    if (m_params_free_set.size()) {
        latex << "\\section{Parameter settings}\n\n";
        map_ex_ex::const_iterator itm, itme;
        itm = m_params_free_set.begin();
        itme = m_params_free_set.end();
        for (; itm != itme; ++itm) {
            latex << "\\begin{equation}\n";
            latex << itm->first.tex(pflag) << " = " << itm->second.tex(pflag) << "\n";
            latex << "\\end{equation}\n";
        }
    }
    latex << "\n\n";

    if (!latex.good()) {
        report_errors("(gEcon error): failed writing LaTeX file \'" + full_name
                      + "\'");
    } else if (m_options[verbose]) {
        write_info("LaTeX documentation written to \'" + m_path + m_name
                   + ".tex\' and \'" + full_name + "\'; created file \'"
                   + m_path + m_name + ".results.tex\'");
    }
}




void
Model::write() const
{
#ifdef R_DLL
    write_r();
#else /* R_DLL */
    if (m_options[output_r]) write_r();
#endif /* R_DLL */
    if (m_options[output_logf]) write_log();
    if (m_options[output_latex]) write_latex();

}



