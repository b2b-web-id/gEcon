/*****************************************************************************
 * This file is a part of gEcon.                                             *
 *                                                                           *
 * (c) Chancellery of the Prime Minister of the Republic of Poland 2012-2015 *
 * (c) Grzegorz Klima, Karol Podemski, Kaja Retkiewicz-Wijtiwiak 2015-2018   *
 * License terms can be found in the file 'LICENCE'                          *
 *                                                                           *
 * Author: Grzegorz Klima                                                    *
 *****************************************************************************/

/** \file model_block.cpp
 * \brief Class representing a block of general equilibrium model.
 */


#include <model_block.h>
#include <utils.h>
#include <iostream>
#include <fstream>

using symbolic::internal::num2str;
using symbolic::internal::DROP_T;
using symbolic::internal::DIFF_T;
using symbolic::internal::DEFAULT;
using symbolic::internal::CONVERT_T;
using symbolic::internal::DROP_T;
using symbolic::internal::CONVERT_IDX;
using symbolic::internal::DROP_IDX;
using symbolic::internal::DROP_INDEXING;
using symbolic::internal::INDEXING_ONLY;
using symbolic::internal::str2tex2;
using symbolic::triplet;


void
Model_block::add_definition(const ex &lhs, const ex &rhs, int l)
{
    m_defs_lhs.push_back(exint(lhs, l));
    m_defs_rhs.push_back(exint(rhs, l));
}

void
Model_block::add_controls(const vec_exintstr &cl)
{
    m_controls = cl;
}


void
Model_block::add_objective(const ex &obj, const ex &obj_eq, const ex &lambda, int l)
{
    m_obj_var = obj;
    m_obj_eq = obj_eq;
    m_obj_eq_in = obj_eq;
    m_obj_lm_in = lambda;
    if (m_obj_lm_in) {
        m_obj_lm = m_obj_lm_in;
    } else {
        m_obj_lm = ex("lambda__" + m_name + "_" + obj.str(DROP_IDX | DROP_T), 0);
        m_obj_lm = add_idx(m_obj_lm, m_i1);
        m_obj_lm = add_idx(m_obj_lm, m_i2);
        m_redlm.insert(m_obj_lm);
    }
    m_obj_line = l;
}


void
Model_block::add_constraint(const ex &lhs, const ex &rhs, const ex &lambda, int l)
{
    m_constraints.push_back(exint(rhs - lhs, l));
    m_constraints_in_lhs.push_back(exint(lhs, l));
    m_constraints_in_rhs.push_back(exint(rhs, l));
    m_lagr_mult_in.push_back(exint(lambda, l));
    if (lambda) {
        m_lagr_mult.push_back(exint(lambda, l));
    } else {
        ex lam = ex("lambda__" + m_name + "_" + num2str(1 + (unsigned) m_lagr_mult.size()), 0);
        lam = add_idx(lam, m_i1);
        lam = add_idx(lam, m_i2);
        lam = add_idx(lam, lhs);
        lam = apply_idx(lam, lhs);
        m_lagr_mult.push_back(exint(lam, l));
        m_redlm.insert(lam);
    }
}


void
Model_block::add_constraint_ref(std::string &blname, int ref, int l)
{
    m_constraints_ref.push_back(strintint(blname, ref, l));
}


void
Model_block::add_identity(const ex &lhs, const ex &rhs, int l)
{
    m_identities.push_back(exint(rhs - lhs, l));
    m_identities_in_lhs.push_back(exint(lhs, l));
    m_identities_in_rhs.push_back(exint(rhs, l));
}


void
Model_block::add_shock(const ex &s, int l)
{
    m_shocks.push_back(exint(s, l));
}

void
Model_block::add_shocks(const vec_exint &pl)
{
    m_shocks = pl;
}


void
Model_block::add_calibr(const ex &lhs, const ex &rhs, int l, const vec_exint &pl)
{
    m_calibr_in_lhs.push_back(exint(lhs, l));
    m_calibr_in_rhs.push_back(exint(rhs, l));
    m_calibr.push_back(exint(lhs - rhs, l));
    m_calibr_pl.push_back(pl);
}


std::string
Model_block::get_name(bool full) const
{
    if (!full) return m_name;
    std::string res;
    if (m_i1) res += m_i1.str();
    if (m_i2) res += m_i2.str();
    return res + m_name;
}


void
Model_block::subst_defs()
{
    for (unsigned i = 0; i < m_defs_lhs.size(); ++i) {
        ex what = m_defs_lhs[i].first;
        ex with = m_defs_rhs[i].first;
        if (m_obj_var) {
            m_obj_eq_in = m_obj_eq_in.subst(what, with);
            m_obj_eq = m_obj_eq.subst(what, with);
        }
        for (unsigned j = 0; j < m_constraints.size(); ++j) {
            m_constraints_in_lhs[j].first =
                m_constraints_in_lhs[j].first.subst(what, with);
            m_constraints_in_rhs[j].first =
                m_constraints_in_rhs[j].first.subst(what, with);
            m_constraints[j].first =
                m_constraints[j].first.subst(what, with);
        }
        for (unsigned j = 0; j < m_identities.size(); ++j) {
            m_identities_in_lhs[j].first =
                m_identities_in_lhs[j].first.subst(what, with);
            m_identities_in_rhs[j].first =
                m_identities_in_rhs[j].first.subst(what, with);
            m_identities[j].first =
                m_identities[j].first.subst(what, with);
        }
        for (unsigned j = 0; j < m_calibr.size(); ++j) {
            m_calibr[j].first = m_calibr[j].first.subst(what, with);
        }
    }
}




void
Model_block::lags()
{
    map_ex_int mcontr0, mcontr, mid;
    map_ex_int::const_iterator it;

    // lags in optimization problem
    collect_lags(m_obj_eq, mcontr0);
    for (unsigned i = 0; i < m_constraints.size(); ++i) {
        collect_lags(m_constraints[i].first, mcontr0);
    }
    for (it = mcontr0.begin(); it != mcontr0.end(); ++it) {
        if (m_contr.find(it->first) != m_contr.end()) mcontr.insert(*it); else mid.insert(*it);
    }
    mcontr0.clear();
    for (it = mcontr.begin(); it != mcontr.end(); ++it) {
        int l = it->second;
        ex v = it->first;
        for (int ll = -2; ll >= l; --ll) {
            ex vl = lag(append_name(v, "__lag_"
                        + num2str(-1 - ll)), -1);
            m_lags.insert(vl);
            m_obj_eq = m_obj_eq.subst(lag(v, ll), vl, false);
            for (unsigned i = 0; i < m_constraints.size(); ++i) {
                m_constraints[i].first =
                    m_constraints[i].first.subst(lag(v, ll), vl, false);
            }
            for (unsigned i = 0; i < m_identities.size(); ++i) {
                m_identities[i].first =
                    m_identities[i].first.subst(lag(v, ll), vl, false);
            }
        }
    }
    for (it = mcontr.begin(); it != mcontr.end(); ++it) {
        int l = it->second;
        ex v = it->first;
        ex vl = append_name(v, "__lag_1"), lmvl;
        std::string vname = v.str(DROP_T | DROP_IDX | DROP_INDEXING);
        m_constraints.push_back(exint(lag(v, -1) - vl, 0));
        lmvl = ex("lambda__" + m_name + "_" + vname + "__lag_1", 0);
        lmvl = add_idx(lmvl, m_i1);
        lmvl = add_idx(lmvl, m_i2);
        lmvl = add_idx(lmvl, v);
        lmvl = apply_idx(lmvl, v);
        m_lagr_mult.push_back(exint(lmvl, 0));
        m_redlm.insert(lmvl);
        m_controls.push_back(exintstr(vl, 0, ""));
        for (int ll = -2; ll > l; --ll) {
            ex vl1 = lag(append_name(v, "__lag_"
                        + num2str(-1 - ll)), -1);
            ex vl2 = append_name(v, "__lag_" + num2str(-ll));
            vl1 = add_idx(vl1, v);
            vl1 = apply_idx(vl1, v);
            vl2 = add_idx(vl2, v);
            vl2 = apply_idx(vl2, v);
            m_constraints.push_back(exint(vl1 - vl2, 0));
            lmvl = ex("lambda__" + m_name + "_" + vname + "__lag_"
                      + num2str(-ll), 0);
            lmvl = add_idx(lmvl, m_i1);
            lmvl = add_idx(lmvl, m_i2);
            lmvl = add_idx(lmvl, v);
            lmvl = apply_idx(lmvl, v);
            m_lagr_mult.push_back(exint(lmvl, 0));
            m_redlm.insert(lmvl);
            m_controls.push_back(exintstr(vl2, 0, ""));
        }
    }

    // lags in identities
    for (unsigned i = 0; i < m_identities.size(); ++i) {
        collect_lags(m_identities[i].first, mid);
    }
    for (it = mid.begin(); it != mid.end(); ++it) {
        int l = it->second;
        ex v = it->first;
        for (int ll = -2; ll >= l; --ll) {
            ex vl = lag(append_name(v, "__lag_"
                        + num2str(-1 - ll)), -1);
            m_lags.insert(vl);
            m_obj_eq = m_obj_eq.subst(lag(v, ll), vl, false);
            for (unsigned i = 0; i < m_constraints.size(); ++i) {
                m_constraints[i].first =
                    m_constraints[i].first.subst(lag(v, ll), vl, false);
            }
            for (unsigned i = 0; i < m_identities.size(); ++i) {
                m_identities[i].first =
                    m_identities[i].first.subst(lag(v, ll), vl, false);
            }
        }
    }
    for (it = mid.begin(); it != mid.end(); ++it) {
        int l = it->second;
        ex v = it->first;
        int ls = -1;
        map_ex_int::const_iterator it2 = mcontr.find(v);
        if (it2 != mcontr.end()) {
            ls = it2->second;
        }
        if (ls == -1) {
            ex vl = append_name(v, "__lag_1");
            m_identities.push_back(exint(lag(v, -1) - vl, 0));
            --ls;
        }
        for (int ll = ls; ll > l; --ll) {
            ex vl1 = lag(append_name(v, "__lag_"
                        + num2str(-1 - ll)), -1);
            ex vl2 = append_name(v, "__lag_" + num2str(-ll));
            m_identities.push_back(exint(vl1 - vl2, 0));
        }
    }
}



void
Model_block::derive_focs()
{
    if (!m_controls.size()) return;

    // find Es
    find_Es(m_obj_eq, m_Es);
    for (unsigned i = 0; i < m_constraints.size(); ++i) {
        find_Es(m_constraints[i].first, m_Es);
    }

    // substitute qs for Es and derive FOCs w.r.t. qs
    set_ex::const_iterator it;
    unsigned j;
    for (it = m_Es.begin(), j = 0; it != m_Es.end(); ++it, ++j) {
        ex q("q__" + m_name + "_"+ num2str(j + 1), 0);
        ex eta("eta__" + num2str(j + 1), 0);
        q = add_idx(q, m_i1);
        q = add_idx(q, m_i2);
        q = add_idx(q, *it);
        q = apply_idx(q, *it);
        eta = add_idx(eta, m_i1);
        eta = add_idx(eta, m_i2);
        eta = add_idx(eta, *it);
        eta = apply_idx(eta, *it);
        m_qs.push_back(q);
        m_etas.push_back(eta);
        m_obj_eq = m_obj_eq.subst(E(*it, 0), m_qs[j]);
        for (unsigned i = 0; i < m_constraints.size(); ++i) {
            m_constraints[i].first = m_constraints[i].first.subst(E(*it, 0), m_qs[j]);
        }
        ex S;
        for (unsigned i = 0; i < m_constraints.size(); ++i) {
            S = S + sum(m_lagr_mult[i].first * m_constraints[i].first);
        }
        m_etas_v.push_back(diff(m_obj_eq + S, m_qs[j]));
    }

    // FOC w.r.t. U
    ex Z;
    bool redlmagrr = false;
    for (it = m_Es.begin(), j = 0; it != m_Es.end(); ++it, ++j) {
            Z = Z + m_etas[j] * *it;
    }
    Z = diff(Z, lag(m_obj_var, 1));
    Z = lag(Z, -1);
    for (unsigned e = 0; e < m_etas.size(); ++e) {
        Z = Z.subst(m_etas[e], m_etas_v[e]);
    }
    if ((!Z.hast()) && (!m_obj_lm_in)) redlmagrr = true;
    if (redlmagrr) {
        m_redlm.erase(m_obj_lm);
        m_obj_lm = ex();
    } else {
        m_focs.push_back(expair(Z - m_obj_lm, m_obj_var));
    }

    // FOCs w.r.t. xs
    for (unsigned n = 0; n < m_controls.size(); ++n) {
        ex L, P, M, Q, R, C, mc, eE, x = m_controls[n].first;
        for (unsigned i = 0; i < m_constraints.size(); ++i) {
            mc = sum(m_lagr_mult[i].first * m_constraints[i].first);
            L = L + mc;
            P = P + lag(mc, 1);
        }
        for (it = m_Es.begin(),  j = 0; it != m_Es.end(); ++it, ++j) {
            eE = m_etas[j] * *it;
            M = M + eE;
            Q = Q + lag(eE, 1);
        }
        R = lag(m_obj_eq, 1) + P + Q;
        if (redlmagrr) {
            C = m_obj_eq + L + M + E(Z * R, 0);
        } else {
            C = m_obj_eq + L + M + E(lag(m_obj_lm, 1) * R, 0);
        }
        C = diff(C, x);
        m_focs.push_back(expair(C, x));
    }

    // Substitute for etas from FOCs w.r.t. qs
    for (unsigned n = 0; n < m_focs.size(); ++n) {
        for (unsigned e = 0; e < m_etas.size(); ++e) {
            m_focs[n].first = m_focs[n].first.subst(m_etas[e], m_etas_v[e]);
        }
    }

    std::vector<bool> has_lags(m_qs.size(), false);
    unsigned cnstr_init_sz = m_constraints.size();
    for (it = m_Es.begin(), j = 0; j < m_qs.size(); ++j, ++it) {
        ex qq = m_qs[j];
        ex EE = E(*it, 0);
        for (unsigned i = 0; i < cnstr_init_sz; ++i) {
            if (m_constraints[i].first.has(qq, DIFF_T)) {
                has_lags[j] = true;
                break;
            }
        }
        for (unsigned i = 0; i < m_focs.size(); ++i) {
            if (m_focs[i].first.has(qq, DIFF_T)) {
                has_lags[j] = true;
                break;
            }
        }
        if (has_lags[j]) {
            m_constraints.push_back(exint(qq - EE, 0));
        }
    }
    for (it = m_Es.begin(), j = 0; j < m_qs.size(); ++j, ++it) {
        ex qq = m_qs[j];
        ex EE = E(*it, 0);

        if (!has_lags[j]) {
            for (unsigned i = 0; i < m_focs.size(); ++i) {
                m_focs[i].first = m_focs[i].first.subst(qq, EE);
            }
            for (unsigned i = 0; i < m_constraints.size(); ++i) {
                m_constraints[i].first = m_constraints[i].first.subst(qq, EE);
            }
            m_obj_eq = m_obj_eq.subst(qq, EE);
        }
    }
    // reduce
    reduce();
}



void
Model_block::derive_focs_deter()
{
    if (!m_controls.size()) return;

    // FOC w.r.t. U
    ex A = m_obj_eq;
    bool redlmagrr = false;
    for (unsigned i = 0; i < m_constraints.size(); ++i) {
        A = A + sum(m_lagr_mult[i].first * m_constraints[i].first);
    }
    A = diff(A, lag(m_obj_var, 1));
    if ((!A.hast()) && (!m_obj_lm_in)) redlmagrr = true;
    if (redlmagrr) {
        m_redlm.erase(m_obj_lm);
        m_obj_lm = ex();
    } else {
        m_focs.push_back(expair(lag(A, -1) - m_obj_lm, m_obj_var));
    }

    // FOCs w.r.t. xs
    for (unsigned n = 0; n < m_controls.size(); ++n) {
        ex L, P, C, mc, x = m_controls[n].first;
        for (unsigned i = 0; i < m_constraints.size(); ++i) {
            mc = sum(m_lagr_mult[i].first * m_constraints[i].first);
            L = L + mc;
            P = P + lag(mc, 1);
        }
        if (redlmagrr) {
            C = m_obj_eq + L + A * (lag(m_obj_eq, 1) + P);
        } else {
            C = m_obj_eq + L + lag(m_obj_lm, 1) * (lag(m_obj_eq, 1) + P);
        }
        C = diff(C, x);
        m_focs.push_back(expair(C, x));
    }
    // reduce
    reduce();
}



void
Model_block::derive_focs_static()
{
    if (!m_controls.size()) return;
    m_redlm.erase(m_obj_lm);
    m_obj_lm = ex();
    // FOCs w.r.t. xs
    for (unsigned n = 0; n < m_controls.size(); ++n) {
        ex L, P, C, mc, x = m_controls[n].first;
        for (unsigned i = 0; i < m_constraints.size(); ++i) {
            mc = sum(m_lagr_mult[i].first * m_constraints[i].first);
            L = L + mc;
        }
        C = m_obj_eq + L;
        C = diff(C, x);
        m_focs.push_back(expair(C, x));
    }
    // reduce
    reduce();
}



void
Model_block::reduce()
{
    // reduction
    vec_expair focs_red = m_focs;
    // do it
    unsigned i, n = focs_red.size();
    bool try_red;
    do {
        try_red = false;
        for (i = 0; i < n; ++i) {
            triplet<bool, ex, ex> ts = find_subst(focs_red[i].first, m_redlm);
            if (!ts.first) continue;
            if (ts.third.hast()) {
                ex e = ts.second;
                bool laglead = false;
                for (unsigned i = 0; i < n; ++i) {
                    if (focs_red[i].first.has(e, DIFF_T)) laglead = true;
                    if (laglead) break;
                }
                if (laglead) continue;
            }
            ex what = ts.second, with = ts.third;
            m_redlm.erase(what);
            focs_red[i].first = ex();
            for (unsigned j = 0; j < n; ++j) {
                focs_red[j].first = focs_red[j].first.subst(what, with);
            }
            try_red = true;
            break;
        }
    } while (try_red && m_redlm.size());
    // finalise
    m_focs_red.reserve(n);
    for (i = 0; i < n; ++i) {
        if (focs_red[i].first)
            m_focs_red.push_back(focs_red[i]);
    }
}


bool
Model_block::has_ref(const ex &v) const
{
    if (m_obj_var.has(v)) return true;
    for (unsigned i = 0; i < m_controls.size(); ++i) {
        if (m_controls[i].first.has(v)) return true;
    }
    for (unsigned i = 0; i < m_lagr_mult_in.size(); ++i) {
        if (m_lagr_mult_in[i].first.has(v)) return true;
    }
    return false;
}


void
Model_block::collect_vp(set_ex &vars, set_ex &parms, bool init)
{
    if (init) {
        collect(m_obj_var, vars, parms);
        collect(m_obj_eq, vars, parms);
        for (int d = 0, D = m_defs_lhs.size(); d < D; ++d) {
            collect(m_defs_lhs[d].first, vars, parms);
        }
        for (int d = 0, D = m_defs_rhs.size(); d < D; ++d) {
            collect(m_defs_rhs[d].first, vars, parms);
        }
        for (int c = 0, C = m_constraints.size(); c < C; ++c) {
            collect(m_constraints[c].first, vars, parms);
        }
        for (int id = 0, Id = m_identities.size(); id < Id; ++id) {
            collect(m_identities[id].first, vars, parms);
        }
    } else {
        collect(expand(ex(m_i1, ex(m_i2, m_obj_var))), vars, parms);
        collect(expand(ex(m_i1, ex(m_i2, m_obj_eq))), vars, parms);
        for (int f = 0, F = m_focs_red.size(); f < F; ++f) {
            collect(expand(ex(m_i1, ex(m_i2, m_focs_red[f].first))), vars, parms);
        }
        for (int c = 0, C = m_constraints.size(); c < C; ++c) {
            collect(expand(ex(m_i1, ex(m_i2, m_constraints[c].first))), vars, parms);
        }
        for (int id = 0, Id = m_identities.size(); id < Id; ++id) {
            collect(expand(ex(m_i1, ex(m_i2, m_identities[id].first))), vars, parms);
        }
    }
}




void
Model_block::write_latex(std::ostream &os, bool static_model) const
{
    symbolic::internal::print_flag pflag = (static_model) ?
        symbolic::internal::DROP_T : symbolic::internal::DEFAULT;
    if (m_controls.size() || m_identities.size()) {
        os << "\\section{" << str2tex2(m_name);
        if (m_i1) os << " $" << m_i1.tex() << '$';
        if (m_i2) os << " $" << m_i2.tex() << '$';
        os << "}\n\n";
    } else {
        return;
    }

    if (m_controls.size()) {
        os << "\\subsection{Optimisation problem}\n\n";
        os << "\\begin{align}\n";
        os << "&\\max_{";
        for (unsigned n = 0; n < m_controls.size() - 1; ++n) {
            os << m_controls[n].first.tex(pflag) << ", ";
        }
        os << m_controls.back().first.tex(pflag) << "\n} "
           << m_obj_var.tex(pflag) << " = "
           << m_obj_eq_in.tex(pflag) << "\\\\\n";
        if (m_constraints_in_rhs.size() || m_constraints_ref.size()) {
            os << "&\\mathrm{s.t.:}\\nonumber\\\\\n";
            if (m_constraints_in_rhs.size()) {
                for (unsigned i = 0; i < m_constraints_in_rhs.size(); ++i) {
                    os << "& "
                       << (m_constraints_in_lhs[i].first - m_constraints_in_rhs[i].first).tex(INDEXING_ONLY)
                       << m_constraints_in_lhs[i].first.tex(pflag | DROP_INDEXING) << " = "
                       << m_constraints_in_rhs[i].first.tex(pflag | DROP_INDEXING)
                       << " \\quad \\left(" << m_lagr_mult[i].first.tex(pflag | DROP_INDEXING);
                    if ((i + 1 == m_constraints_in_rhs.size()) && !m_constraints_ref.size())
                        os << "\\right)\n";
                    else os << "\\right)\\\\\n";
                }
            }
            if (m_constraints_ref.size()) {
                for (unsigned i = 0; i < m_constraints_ref.size(); ++i) {
                    os << "& \\mathrm{" << str2tex2(m_constraints_ref[i].first) << "\'s\\ ";
                    if (m_constraints_ref[i].second == objective) os << "objective";
                    else if (m_constraints_ref[i].second == constraints) os << "constraints";
                    else if (m_constraints_ref[i].second == focs) os << "first\\ order\\ conditions";
                    else if (m_constraints_ref[i].second == identities) os << "identities";
                    os << '}';
                    if (i + 1 == m_constraints_ref.size()) os << "\\nonumber\n";
                    else os << "\\nonumber\\\\\n";
                }
            }
        }
        os << "\\end{align}\n\n\n";
    }

    if (m_identities_in_lhs.size()) {
        os << "\\subsection{Identities}\n\n";
        for (unsigned id = 0; id < m_identities_in_lhs.size(); ++id) {
            os << "\\begin{equation}\n"
               << (m_identities_in_lhs[id].first - m_identities_in_rhs[id].first).tex(INDEXING_ONLY)
               << m_identities_in_lhs[id].first.tex(pflag | DROP_INDEXING) << " = "
               << m_identities_in_rhs[id].first.tex(pflag | DROP_INDEXING)
               << "\n\\end{equation}\n";
        }
        os << "\n\n";
    }

    if (m_focs.size()) {
        os << "\\subsection{First order conditions}\n\n";
        for (unsigned f = 0; f < m_focs.size(); ++f) {
            os << "\\begin{equation}\n"
               << m_focs[f].first.tex(INDEXING_ONLY)
               << m_focs[f].first.tex(pflag | DROP_INDEXING) << " = 0\n"
               << " \\quad \\left(" << m_focs[f].second.tex(pflag | DROP_INDEXING)
               << "\\right)\n\\end{equation}\n";
        }
        os << "\n\n";
        if (m_focs_red.size() < m_focs.size()) {
            os << "\\subsection{First order conditions after reduction}\n\n";
            for (unsigned f = 0; f < m_focs_red.size(); ++f) {
                os << "\\begin{equation}\n"
                   << m_focs_red[f].first.tex(INDEXING_ONLY)
                   << m_focs_red[f].first.tex(pflag | DROP_INDEXING) << " = 0\n"
                   << " \\quad \\left(" << m_focs_red[f].second.tex(pflag)
                   << "\\right)\n\\end{equation}\n";
            }
            os << "\n\n";
        }
    }

    os << "\n\n";
}



void
Model_block::write_logfile(std::ostream &os, bool static_model) const
{
    symbolic::internal::print_flag pflag = (static_model) ?
        symbolic::internal::DROP_T : symbolic::internal::DEFAULT;
    os << "Block: ";
    if (m_i1) os << m_i1.str() << ' ';
    if (m_i2) os << m_i2.str() << ' ';
    os << m_name << '\n';

    if (m_defs_lhs.size()) {
        // definitions
        os << " Definitions:\n";
        for (unsigned i = 0; i < m_defs_lhs.size(); ++i) {
            os << "    " << m_defs_lhs[i].first.str(pflag) << " = "
               << m_defs_rhs[i].first.str(pflag) << '\n';
        }
    }
    if (m_controls.size()) {
        int refctrl = 0;
        os << " Controls:\n    ";
        for (unsigned i = 0; i < m_controls.size() - 1; ++i) {
            os << m_controls[i].first.str(pflag) << ", ";
            if (m_controls[i].third != "") ++refctrl;
        }
        os << m_controls.back().first.str(pflag) << '\n';
        if (m_controls.back().third != "") ++refctrl;

        if (refctrl) {
            os << " Controls referenced from other model blocks:\n    ";
            int count = 0;
            for (unsigned i = 0; i < m_controls.size(); ++i) {
                if (m_controls[i].third != "") {
                    os << m_controls[i].first.str(pflag) << '@'
                       << m_controls[i].third;
                    ++count;
                    if (count < refctrl) os << ", "; else os << '\n';
                }
            }
        }

        os << " Objective:\n    ";
        os << m_obj_var.str(pflag) << " = " << m_obj_eq_in.str(pflag);
        if (m_obj_lm) os << "    (" << m_obj_lm.str(pflag) << ')';
        os << '\n';

        if (m_constraints.size()) {
            os << " Constraints:\n";
            for (unsigned i = 0; i < m_lagr_mult.size(); ++i) {
                os << "    " << m_constraints[i].first.str(pflag) << " = 0    ("
                   << m_lagr_mult[i].first.str(pflag) << ")\n";
            }
            for (unsigned i = m_lagr_mult.size(); i < m_constraints.size(); ++i) {
                os << "    " << m_constraints[i].first.str(pflag)
                   << " = 0    (auxiliary)\n";
            }
        }
        if (m_constraints_ref.size()) {
            os << " Constraints referencing equations in other blocks:\n";
            for (unsigned i = 0; i < m_constraints_ref.size(); ++i) {
                os << "    ";
                if (m_constraints_ref[i].second == objective) os << "objective";
                else if (m_constraints_ref[i].second == constraints) os << "constraints";
                else if (m_constraints_ref[i].second == focs) os << "first order conditions";
                else if (m_constraints_ref[i].second == identities) os << "identities";
                os << " @ " << m_constraints_ref[i].first << '\n';
            }
        }
    }

    if (m_identities.size()) {
        // identities
        os << " Identities:\n";
        for (unsigned i = 0; i < m_identities.size(); ++i) {
            os << "    " << m_identities[i].first.str(pflag) << " = 0\n";
        }
    }

    if (m_focs.size()) {
        // FOCs
        os << " First order conditions:\n";
        for (unsigned i = 0; i < m_focs.size(); ++i) {
            os << "    " << m_focs[i].first.str(pflag) << " = 0    ("
               << m_focs[i].second.str(pflag) << ")\n";
        }
        if (m_focs_red.size() < m_focs.size()) {
            // FOCs
            os << " First order conditions after reduction:\n";
            for (unsigned i = 0; i < m_focs_red.size(); ++i) {
                os << "    " << m_focs_red[i].first.str(pflag) << " = 0    ("
                << m_focs_red[i].second.str(pflag) << ")\n";
            }
        }
    }

    os << '\n';
}





