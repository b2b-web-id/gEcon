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

using symbolic::internal::num2str;


Model::Model() : m_deter(false), m_static(false)
{
    m_options[silent] = false;
    m_options[verbose] = false;
    m_options[warnings_opt] = true;
    m_options[output_latex] = false;
    m_options[output_latex_long] = true;
    m_options[output_latex_landscape] = false;
    m_options[output_r_long] = false;
    m_options[output_r_jacobian] = true;
    m_options[output_r_rcpp] = false;

#ifdef R_DLL
    m_options[output_r] = true;
    m_options[output_logf] = false;
#else /* R_DLL */
    m_options[output_r] = false;
    m_options[output_logf] = true;
#endif /* R_DLL */
    for (unsigned i = 0; i < OPTIONS_LENGTH; ++i)
        m_options_set[i] = 0;
}


void
Model::write_model_info(const std::string &mes) const
{
    if (!m_options[silent]) ::write_info("(gEcon model info): " + mes);
}


void
Model::write_info(const std::string &mes) const
{
    if (!m_options[silent]) ::write_info("(gEcon info): " + mes);
}


void
Model::warning(const std::string &mes, int l)
{
    if (m_options[warnings_opt]) {
        if (l) {
            m_warn.push_back(mes + "; warning near line " + num2str(l));
        } else {
            m_warn.push_back(mes);
        }
    }
}


void
Model::error(const std::string &mes, int l)
{
#ifdef DEBUG
    if (l) {
        std::cerr << "ERROR: " << mes << "; error near line " << l << '\n';
    } else {
        std::cerr << "ERROR: " << mes << '\n';
    }
#endif  /* DEBUG */
    if (l) {
        m_err.push_back(mes + "; error near line " + num2str(l));
    } else {
        m_err.push_back(mes);
    }
}


void
Model::terminate_on_errors()
{
    if (errors()) {
        write_log();
        if (warnings()) report_warns(get_warns());
        std::string errs = get_errs();
        clear();
        report_errors(errs);
    }
}


#if defined(NERRTHRESH)
#undef NERRTHRESH
#endif
#define NERRTHRESH 10


std::string
Model::get_errs(bool log) const
{
    std::string errs;
    std::vector<std::string>::const_iterator it;
    int i = 1, n = m_err.size(), nn;
    if ((n <= NERRTHRESH) || log) nn = n; else nn = NERRTHRESH - 1;
    for (it = m_err.begin(); i <= nn; ++it, ++i) {
        errs += "(gEcon model error " + num2str(i) + "): " + *it;
        if (i < nn) errs += '\n';
    }
    if (log) return errs;
    if (n > NERRTHRESH) {
        errs += "\n(gEcon model errors): " + num2str(n - NERRTHRESH + 1) + " more follow, see logfile";
    }
    time_t tt = time(0);
    struct tm *now = localtime(&tt);
    if ((now->tm_mday == 13) && (now->tm_wday == 5)) {
        errs += "\nWell, \"I Ain't Superstitious\", but it's Friday the 13th and you've got an error...";
    }
    return errs;
}


std::string
Model::get_warns(bool log) const
{
    std::string warns;
    std::vector<std::string>::const_iterator it;
    int i = 1, n = m_warn.size(), nn;
    if ((n <= NERRTHRESH) || log) nn = n; else nn = NERRTHRESH - 1;
    for (it = m_warn.begin(); i <= nn; ++it, ++i) {
        warns += "(gEcon model warning " + num2str(i) + "): " + *it;
        if (i < nn) warns += '\n';
    }
    if (log) return warns;
    if (n > NERRTHRESH) {
        warns += "\n(gEcon model warnings): " + num2str(n - NERRTHRESH + 1) + " more follow, see logfile";
    }
    return warns;
}


void
Model::check_warns()
{
    if (m_warn.size() > NERRTHRESH){
        m_options[output_logf] = true;
    }
}


void
Model::set_name(const std::string &name)
{
    size_t slash = name.find_last_of('/');
    if (slash == std::string::npos) {
        m_name = name;
        return;
    }
    m_path = name.substr(0, slash + 1);
    m_name = name.substr(slash + 1);
}


bool
Model::add_set(const idx_set &s)
{
#ifdef DEBUG
    std::cerr << "DEBUG INFO: adding index set: " << s.name() << '\n';
#endif /* DEBUG */
    if (!m_set_names.insert(s.name()).second) return false;
    m_sets.insert(std::pair<std::string, idx_set>(s.name(), s));
    return true;
}


bool
Model::is_set(const std::string &s)
{
    std::map<std::string, symbolic::idx_set>::const_iterator it;
    if ((it = m_sets.find(s)) != m_sets.end()) return true;
    return false;
}


idx_set
Model::get_set(const std::string &s)
{
    std::map<std::string, symbolic::idx_set>::const_iterator it;
    if ((it = m_sets.find(s)) != m_sets.end()) return it->second;
    return idx_set("??");
}


void
Model::add_block(const std::string &s, int ln, idx_ex i1, idx_ex i2)
{
#ifdef DEBUG
    std::cerr << "DEBUG INFO: adding block: ";
    if (i1) std::cerr << i1.str() << ' ';
    if (i2) std::cerr << i2.str() << ' ';
    std::cerr << s << '\n';
#endif /* DEBUG */
    m_names.insert(s);
    // Stray indices?
    std::set<unsigned> iset;
    unsigned idx;
    idx = i2.get_excl_id();
    if (idx) iset.insert(idx);
    iset.erase(i1.get_id());
    idx = i1.get_excl_id();
    if (idx) iset.insert(idx);
    if (iset.size()) {
        unsigned n = iset.size();
        std::string mes = "stray ";
        if (n == 1) mes += "index (";
        else mes += "indices (";
        const symbolic::internal::stringhash &ref =
            symbolic::internal::stringhash::get_instance();
        std::set<unsigned>::const_iterator it = iset.begin();
        mes += '\"' +  ref.get_str(*it) + '\"';
        for (++it; it != iset.end(); ++it) {
            mes += ", \"" +  ref.get_str(*it) + '\"';
        }
        mes += ") in block declaration \"";
        if (i1) mes += i1.str() + ' ';
        if (i2) mes += i2.str() + ' ';
        mes += s + "\"";
        error(mes, ln);
    }
    // Duplicated free indices?
    std::map<unsigned, unsigned> imap;
    collect_fidx(i1, imap);
    collect_fidx(i2, imap);
    std::map<unsigned, unsigned>::iterator it;
    for (it = imap.begin(); it != imap.end(); ++it) {
        if (!it->second) imap.erase(it);
    }
    if (imap.size()) {
        unsigned n = imap.size();
        std::string mes = "duplicated free ";
        if (n == 1) mes += "index in nested indexing expression (";
        else mes += "indices in nested indexing expression (";
        const symbolic::internal::stringhash &ref =
            symbolic::internal::stringhash::get_instance();
        it = imap.begin();
        mes += '\"' +  ref.get_str(it->first) + '\"';
        for (++it; it != imap.end(); ++it) {
            mes += ", \"" +  ref.get_str(it->first) + '\"';
        }
        mes += ") in block declaration \"";
        if (i1) mes += i1.str() + ' ';
        if (i2) mes += i2.str() + ' ';
        mes += s + "\"";
        error(mes, ln);
    }

    m_blocks.push_back(Model_block(s, i1, i2));
}


void
Model::clear()
{
    *this = Model();
}

