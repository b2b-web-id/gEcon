/***********************************************************
 * (c) Kancelaria Prezesa Rady Ministrów 2012-2015         *
 * Treść licencji w pliku 'LICENCE'                        *
 *                                                         *
 * (c) Chancellery of the Prime Minister 2012-2015         *
 * License terms can be found in the file 'LICENCE'        *
 *                                                         *
 * Author: Grzegorz Klima                                  *
 ***********************************************************/

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

using symbolic::internal::num2str;
using symbolic::internal::print_flag;
using symbolic::internal::DEFAULT;
using symbolic::internal::CONVERT_T;
using symbolic::internal::DROP_T;
using symbolic::internal::CONVERT_IDX;
using symbolic::internal::DROP_IDX;
using symbolic::internal::DROP_QUOTES;



void
Model::write_logf() const
{
    std::string full_name = m_path + m_name + ".model.log";
#ifdef DEBUG
    std::cerr << "Writing logfile to: \'" << full_name << "\'\n";
#endif /* DEBUG */
    std::ofstream logfile(full_name.c_str());
    if (!logfile.good()) {
        report_errors("(gEcon error): failed opening logfile \'" + full_name
                      + "\' for writing");
    }
    write_log(logfile);
    if (!logfile.good()) {
        report_errors("(gEcon error): failed writing logfile \'" + full_name
                      + "\'");
    } else if (m_options[verbose]) {
        write_info("logfile written to \'" + full_name + "\'");
    }

}


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
info_str()
{
    return "Generated on " + time_str() + " by gEcon (" + gecon_web_str()
           + ") ver. " + gecon_ver_str();
}

} /* namespace */



void
Model::write_log(std::ostream &logfile) const
{
    logfile << info_str() << "\n\n";
    logfile << "Model name: " << m_name << "\n\n";

    std::string tab("    ");
    if (m_sets.size()) {
        logfile << "Index sets (" << m_sets.size() << "):\n";
        std::map<std::string, symbolic::idx_set>::const_iterator it;
        for (it = m_sets.begin(); it != m_sets.end(); ++it) {
            logfile << tab << it->second.str() << '\n';
        }
        logfile << '\n';
    }

    for (unsigned i = 0, n = m_blocks.size(); i < n; ++i) {
        m_blocks[i].write_logfile(logfile, m_static);
    }

    set_ex::const_iterator it;
    print_flag pflag = (m_static) ? DROP_T : DEFAULT;

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
        for (unsigned i = 0; i < m_ss.size(); ++i) {
            logfile << " (" << (i + 1) << ")  " << m_ss[i].str(pflag) << " = 0\n";
        }
        logfile << "\n";
    }

    if (m_calibr.size()) {
        logfile << "Calibration equations (" << m_calibr.size() << "):\n";
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

} /* namespace */



void
Model::write_r() const
{
    std::string full_name = m_path + m_name + ".model.R";
#ifdef DEBUG
    std::cerr << "Writing R code to file: \'" << full_name << "\'\n";
#endif /* DEBUG */
    std::ofstream R(full_name.c_str());
    if (!R.good()) {
        report_errors("(gEcon error): failed opening R file \'" + full_name
                      + "\' for writing");
    }
    set_ex::const_iterator it;
    std::string tab("    ");
    unsigned index;
    unsigned n_v = m_vars.size(), n_s = m_shocks.size(),
             n_c = m_params_calibr.size(), n_f = m_params_free.size();

    R << "# " << info_str() << "\n# Model name: " << m_name << "\n\n";

    R << "# info\ninfo__ <- c(\"" << m_name << "\", \""
      << m_path << m_name << "\", \""
      << time_str() << "\")\n\n";

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

    R << "# shocks\n";
    if (n_s) {
        R << "shocks__ <- c(";
        it = m_shocks.begin();
        R << '\"' << it->str(DROP_T | CONVERT_IDX) << '\"';
        for (++it; it != m_shocks.end(); ++it) {
            R << ",\n              \"" << it->str(DROP_T | CONVERT_IDX) << '\"';
        }
        R << ")\n\n";
    } else {
        R << "shocks__ <- character(0)\n\n";
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
    } else {
        R << "parameters__ <- character(0)\n\n";
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
//     pflag |= DROP_QUOTES;
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
    R << "# variables / equations map\n";
    R << "vareqmap__ <- sparseMatrix(i = c(";
    itim = m_var_eq_map.begin();
    R << itim->first.first;
    for (++itim; itim != m_var_eq_map.end(); ++itim) R << ", " << itim->first.first;
    R << "),\n                           j = c(";
    itim = m_var_eq_map.begin();
    R << itim->first.second;
    for (++itim; itim != m_var_eq_map.end(); ++itim) R << ", " << itim->first.second;
    R << "),\n                           x = c(";
    itim = m_var_eq_map.begin();
    R << itim->second;
    for (++itim; itim != m_var_eq_map.end(); ++itim) R << ", " << itim->second;
    R << "),\n                           dims = c(" << n_v << ", " << n_v << "))\n\n";

    std::set<std::pair<int, int> >::const_iterator itis;

    R << "# variables / calibrating equations map\n";
    if (m_var_ceq_map.size()) {
        R << "varcalibreqmap__ <- sparseMatrix(i = c(";
        itis = m_var_ceq_map.begin();
        R << itis->first;
        for (++itis; itis != m_var_ceq_map.end(); ++itis) R << ", " << itis->first;
        R << "),\n                                 j = c(";
        itis = m_var_ceq_map.begin();
        R << itis->second;
        for (++itis; itis != m_var_ceq_map.end(); ++itis) R << ", " << itis->second;
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
        for (++itis; itis != m_cpar_eq_map.end(); ++itis) R << ", " << itis->first;
        R << "),\n                                 j = c(";
        itis = m_cpar_eq_map.begin();
        R << itis->second;
        for (++itis; itis != m_cpar_eq_map.end(); ++itis) R << ", " << itis->second;
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
        for (++itis; itis != m_cpar_ceq_map.end(); ++itis) R << ", " << itis->first;
        R << "),\n                                       j = c(";
        itis = m_cpar_ceq_map.begin();
        R << itis->second;
        for (++itis; itis != m_cpar_ceq_map.end(); ++itis) R << ", " << itis->second;
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
        for (++itis; itis != m_fpar_eq_map.end(); ++itis) R << ", " << itis->first;
        R << "),\n                               j = c(";
        itis = m_fpar_eq_map.begin();
        R << itis->second;
        for (++itis; itis != m_fpar_eq_map.end(); ++itis) R << ", " << itis->second;
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
        for (++itis; itis != m_fpar_ceq_map.end(); ++itis) R << ", " << itis->first;
        R << "),\n                                     j = c(";
        itis = m_fpar_ceq_map.begin();
        R << itis->second;
        for (++itis; itis != m_fpar_ceq_map.end(); ++itis) R << ", " << itis->second;
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
        for (++itis; itis != m_shock_eq_map.end(); ++itis) R << ", " << itis->first;
        R << "),\n                             j = c(";
        itis = m_shock_eq_map.begin();
        R << itis->second;
        for (++itis; itis != m_shock_eq_map.end(); ++itis) R << ", " << itis->second;
        R << "),\n                             x = rep(1, " << m_shock_eq_map.size();
        R << "), dims = c(" << n_v << ", " << n_s << "))\n\n";
    } else {
        R << "shockeqmap__ <- sparseMatrix(i = NULL, j = NULL, dims = c(";
        R << n_v << ", " << n_s << "))\n\n";
    }

    // str map
    map_str_str mss;
    if (!m_options[output_r_long]) {
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
    // settings for printing R functions
    pflag = (m_static) ? DROP_T : CONVERT_T;
    pflag |= CONVERT_IDX;

    R << "# steady state equations\n";
    R << "ss_eq__ <- function(v, pc, pf)\n";
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
        R << tab << "r <- numeric(" << m_eqs.size() << ")\n";
        for (index = 0; index != m_ss.size(); ++index) {
            R << tab << "r[" << (index + 1) << "] = " << m_ss[index].str(pflag) << "\n" ;
        }
    } else {
        R << tab << "r <- numeric(" << m_eqs.size() << ")\n";
        for (index = 0; index != m_ss.size(); ++index) {
            R << tab << "r[" << (index + 1) << "] = " << m_ss[index].strmap(mss) << "\n" ;
        }
    }
    R << "\n" << tab << "return(r)" << "\n}\n\n";

    R << "# calibration equations\n";
    R << "calibr_eq__ <- function(v, pc, pf)\n";
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
        R << tab << "r <- numeric(" << m_calibr.size() << ")\n";
        for (it = m_calibr.begin(), index = 1; it != m_calibr.end(); ++it, ++index) {
            R << tab << "r[" << index << "] = " << it->str(pflag) << "\n" ;
        }
    } else {
        R << tab << "r <- numeric(" << m_calibr.size() << ")\n";
        for (it = m_calibr.begin(), index = 1; it != m_calibr.end(); ++it, ++index) {
            R << tab << "r[" << index << "] = " << it->strmap(mss) << "\n" ;
        }
    }
    R << '\n' << tab << "return(r)" << "\n}\n\n";

    std::map<std::pair<int, int>, ex>::const_iterator itm;
    R << "# steady state and calibration equations Jacobian\n";
    R << "ss_calibr_eq_jacob__ <- function(v, pc, pf)\n";
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
        R << tab << "jacob <- Matrix(0, nrow = "
          << (m_eqs.size() + m_calibr.size()) << ", ncol = "
          << (m_vars.size() + m_params_calibr.size()) << ", sparse = TRUE)\n";
        for (itm = m_jacob_ss_calibr.begin(); itm != m_jacob_ss_calibr.end(); ++itm) {
            R << tab << "jacob[" << itm->first.first << ", " << itm ->first.second
            << "] = " << itm->second.str(pflag) << '\n';
        }
    } else {
        R << tab << "jacob <- Matrix(0, nrow = "
          << (m_eqs.size() + m_calibr.size()) << ", ncol = "
          << (m_vars.size() + m_params_calibr.size()) << ", sparse = TRUE)\n";
        for (itm = m_jacob_ss_calibr.begin(); itm != m_jacob_ss_calibr.end(); ++itm) {
            R << tab << "jacob[" << itm->first.first << ", " << itm ->first.second
            << "] = " << itm->second.strmap(mss) << '\n';
        }
    }
    R << "\n" << tab << "return(jacob)" << "\n}\n\n";

    R << "# 1st order perturbation\n";
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

    R << tab << "Atm1 <- Matrix(0, nrow = " << n_v << ", ncol = " << n_v
      << ", sparse = TRUE)\n";
    if (m_options[output_r_long]) {
        for (itm = m_Atm1.begin(); itm != m_Atm1.end(); ++itm) {
            R << tab << "Atm1[" << itm->first.first << ", " << itm ->first.second
            << "] = " << itm->second.str(pflag) << '\n';
        }
    } else {
        for (itm = m_Atm1.begin(); itm != m_Atm1.end(); ++itm) {
            R << tab << "Atm1[" << itm->first.first << ", " << itm ->first.second
            << "] = " << itm->second.strmap(mss) << '\n';
        }
    }
    R << '\n';

    R << tab << "At <- Matrix(0, nrow = " << n_v << ", ncol = " << n_v
      << ", sparse = TRUE)\n";
    if (m_options[output_r_long]) {
        for (itm = m_At.begin(); itm != m_At.end(); ++itm) {
            R << tab << "At[" << itm->first.first << ", " << itm ->first.second
            << "] = " << itm->second.str(pflag) << '\n';
        }
    } else {
        for (itm = m_At.begin(); itm != m_At.end(); ++itm) {
            R << tab << "At[" << itm->first.first << ", " << itm ->first.second
            << "] = " << itm->second.strmap(mss) << '\n';
        }
    }
    R << '\n';

    R << tab << "Atp1 <- Matrix(0, nrow = " << n_v << ", ncol = " << n_v
      << ", sparse = TRUE)\n";
    if (m_options[output_r_long]) {
        for (itm = m_Atp1.begin(); itm != m_Atp1.end(); ++itm) {
            R << tab << "Atp1[" << itm->first.first << ", " << itm ->first.second
            << "] = " << itm->second.str(pflag) << '\n';
        }
    } else {
        for (itm = m_Atp1.begin(); itm != m_Atp1.end(); ++itm) {
            R << tab << "Atp1[" << itm->first.first << ", " << itm ->first.second
            << "] = " << itm->second.strmap(mss) << '\n';
        }
    }
    R << '\n';

    R << tab << "Aeps <- Matrix(0, nrow = " << n_v << ", ncol = " << n_s
      << ", sparse = TRUE)\n";
    if (m_options[output_r_long]) {
        for (itm = m_Aeps.begin(); itm != m_Aeps.end(); ++itm) {
            R << tab << "Aeps[" << itm->first.first << ", " << itm ->first.second
            << "] = " << itm->second.str(pflag) << '\n';
        }
    } else {
        for (itm = m_Aeps.begin(); itm != m_Aeps.end(); ++itm) {
            R << tab << "Aeps[" << itm->first.first << ", " << itm ->first.second
            << "] = " << itm->second.strmap(mss) << '\n';
        }
    }

    R << '\n' << tab << "return(list(Atm1, At, Atp1, Aeps))" << "\n}\n\n";

    R << "# create model object\n"
      << "gecon_model(model_info = info__,\n"
      << "            index_sets = index_sets__,\n"
      << "            variables = variables__,\n"
      << "            shocks = shocks__,\n"
      << "            parameters = parameters__,\n"
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
      << "            ss_calibr_function_jac = ss_calibr_eq_jacob__,\n"
      << "            pert = pert1__)\n\n";


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
    std::cerr << "Writing LaTeX documentation to file: \'" << full_name << "\'\n";
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
    latex << "% " << info_str() << "\n% Model name: " << m_name << "\n\n";
    latex << "\\documentclass[10pt,a4paper]{article}\n";
    latex << "\\usepackage[utf8]{inputenc}\n";
    latex << "\\usepackage{color}\n";
    latex << "\\usepackage{graphicx}\n";
    latex << "\\usepackage{epstopdf}\n";
    latex << "\\usepackage{hyperref}\n";
    latex << "\\usepackage{amsmath}\n";
    latex << "\\usepackage{amssymb}\n";
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
    latex << "% " << info_str() << "\n% Model name: " << m_name << "\n\n";
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
    latex << "% " << info_str() << "\n% Model name: " << m_name << "\n\n";

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
            for (unsigned i = 0; i < m_t_ss.size(); ++i) {
                latex << "\\begin{equation}\n";
                latex << m_t_ss[i].tex(pflag) << " = 0\n";
                latex << "\\end{equation}\n";
            }
            latex << "\n\n\n";

            if (m_options[output_latex_long]) {
                latex << "\\section{Steady state relationships (after expansion and reduction)}\n\n";
                for (unsigned i = 0; i < m_ss.size(); ++i) {
                    latex << "\\begin{equation}\n";
                    latex << m_ss[i].tex(pflag) << " = 0\n";
                    latex << "\\end{equation}\n";
                }
                latex << "\n\n\n";
            }
        } else {
            latex << "\\section{Steady state relationships (after reduction)}\n\n";
            for (unsigned i = 0; i < m_ss.size(); ++i) {
                latex << "\\begin{equation}\n";
                latex << m_ss[i].tex(pflag) << " = 0\n";
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
    if (m_options[output_logf]) write_logf();
    if (m_options[output_latex]) write_latex();

}



