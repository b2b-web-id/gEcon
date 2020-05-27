/*****************************************************************************
 * This file is a part of gEcon.                                             *
 *                                                                           *
 * (c) Chancellery of the Prime Minister of the Republic of Poland 2012-2015 *
 * (c) Grzegorz Klima, Karol Podemski, Kaja Retkiewicz-Wijtiwiak 2015-2018   *
 * License terms can be found in the file 'LICENCE'                          *
 *                                                                           *
 * Author: Grzegorz Klima                                                    *
 *****************************************************************************/

/** \file utils.cpp
 * \brief Utilities.
 */

#include <utils.h>
#include <error.h>
#include <sstream>
#include <iostream>
#include <vector>
#include <cstring>
#include <cstdio>
#include <cmath>

using namespace symbolic::internal;


std::string
symbolic::internal::num2str(double n)
{
    if (std::isnan(n)) return "NaN";
    if (std::isinf(n)) {
        if (n == -1./0.) return "-Inf";
        else return "Inf";
    }
    char buf[30];
    std::sprintf(buf, "%g", n);
    return std::string(buf);
}


std::string
symbolic::internal::num2tex(double n)
{
    if (std::isnan(n)) return "\\frac{0}{0}";
    if (std::isinf(n)) {
        if (n == -1./0.) return "-\\infty";
        else return "\\infty";
    }
    char buf[30];
    std::sprintf(buf, "%g", n);
    std::string r0(buf), r1;
    size_t e;
    e = r0.find('e');
    if (e >= r0.size()) return r0;
    r1.assign(r0, 0, e);
    r1 += "\\cdot 10^{";
    if (r0[++e] == '-') r1 += '-';
    ++e;
    while (r0[e] == '0') ++e;
    r1.append(r0, e, r0.size());
    r1 += '}';
    return r1;
}


std::string
symbolic::internal::num2str(int n, int d)
{
    char buf[15];
    if (d) {
        sprintf(buf, "%0*d", d, n);
    } else {
        sprintf(buf, "%d", n);
    }
    return std::string(buf);
}


std::string
symbolic::internal::num2str(unsigned n, int d)
{
    char buf[15];
    if (d) {
        sprintf(buf, "%0*u", d, n);
    } else {
        sprintf(buf, "%u", n);
    }
    return std::string(buf);
}


namespace {



const char* Greeks[] = {
    "Alpha", "\\Alpha",
    "ALPHA", "\\Alpha",
    "alpha", "\\alpha",
    "BETA", "\\Beta",
    "Beta", "\\Beta",
    "beta", "\\beta",
    "GAMMA", "\\Gamma",
    "Gamma", "\\Gamma",
    "gamma", "\\gamma",
    "DELTA", "\\Delta",
    "Delta", "\\Delta",
    "delta", "\\delta",
    "EPSILON", "\\Epsilon",
    "Epsilon", "\\Epsilon",
    "epsilon", "\\epsilon",
    "varepsilon", "\\varepsilon",
    "ZETA", "\\Zeta",
    "Zeta", "\\Zeta",
    "zeta", "\\zeta",
    "ETA", "\\Eta",
    "Eta", "\\Eta",
    "eta", "\\eta",
    "THETA", "\\Theta",
    "Theta", "\\Theta",
    "theta", "\\theta",
    "vartheta", "\\vartheta",
    "IOTA", "\\Iota",
    "Iota", "\\Iota",
    "iota", "\\iota",
    "KAPPA", "\\Kappa",
    "Kappa", "\\Kappa",
    "kappa", "\\kappa",
    "LAMBDA", "\\Lambda",
    "Lambda", "\\Lambda",
    "lambda", "\\lambda",
    "MU", "\\Mu",
    "Mu", "\\Mu",
    "mu", "\\mu",
    "NU", "\\Nu",
    "Nu", "\\Nu",
    "nu", "\\nu",
    "XI", "\\Xi",
    "Xi", "\\Xi",
    "xi", "\\xi",
    "OMICRON", "\\Omicron",
    "Omicron", "\\Omicron",
    "omicron", "\\omicron",
    "PI", "\\Pi",
    "Pi", "\\Pi",
    "pi", "\\pi",
    "varpi", "\\varpi",
    "RHO", "\\Rho",
    "Rho", "\\Rho",
    "rho", "\\rho",
    "varrho", "\\varrho",
    "SIGMA", "\\Sigma",
    "Sigma", "\\Sigma",
    "sigma", "\\sigma",
    "varsigma", "\\varsigma",
    "TAU", "\\Tau",
    "Tau", "\\Tau",
    "tau", "\\tau",
    "UPSILON", "\\Upsilon",
    "Upsilon", "\\Upsilon",
    "upsilon", "\\upsilon",
    "PHI", "\\Phi",
    "Phi", "\\Phi",
    "phi", "\\phi",
    "varphi", "\\varphi",
    "CHI", "\\Chi",
    "Chi", "\\Chi",
    "chi", "\\chi",
    "PSI", "\\Psi",
    "Psi", "\\Psi",
    "psi", "\\psi",
    "OMEGA", "\\Omega",
    "Omega", "\\Omega",
    "omega", "\\omega",
""};


const char* Symbols[] = {
    "star", "\\star",
    "STAR", "\\bigstar",
    "ast", "\\ast",
    "circ", "\\circ",
    "diamond", "\\diamond",
    "clubsuit", "\\clubsuit",
    "flat", "\\flat",
    "sharp", "\\sharp",
    "natural", "\\natural",
    "spadesuit", "\\spadesuit",
    "aleph", "\\aleph",
    "dag", "\\dag",
    "ddag", "\\ddag",
    "infty", "\\infty",
    "nabla", "\\nabla",
    "backprime", "\\backprime",
    "dollar", "\\S",
    "pounds", "\\pounds",
    "in", "\\in",
    "ni", "\\ni",
""};


const char* funcs[] = {
    "sin", "\\sin",
    "sec", "\\sec",
    "csc", "\\csc",
    "cos", "\\cos",
    "tan", "\\tan",
    "cot", "\\cot",
    "arcsin", "\\arcsin",
    "arccos", "\\arccos",
    "arctan", "\\arctan",
    "arccot", "\\arccot",
    "sinh", "\\sinh",
    "cosh", "\\cosh",
    "tanh", "\\tanh",
    "coth", "\\coth",
    "exp", "\\exp",
    "log", "\\log",
    // NOTE: this requiress \DeclareMathOperator{\erf}{erf} in LaTeX source
    "erf", "\\erf",
    // NOTE: this requiress \DeclareMathOperator{\pnorm}{pnorm} in LaTeX source
    "pnorm", "\\pnorm",
""};

} /* namespace */


std::string
symbolic::internal::func2tex(const std::string &n)
{
    for (unsigned i = 0; funcs[i][0]; i += 2) {
        if (!strcmp(funcs[i], n.c_str()))
            return std::string(funcs[++i]);
    }
    INTERNAL_ERROR
    return n;
}




namespace {

std::string
try_Greeksymb(const std::string &s)
{
    for (unsigned i = 0; Greeks[i][0]; i += 2) {
        if (!strcmp(Greeks[i], s.c_str()))
            return std::string(Greeks[i + 1]);
    }
    for (unsigned i = 0; Symbols[i][0]; i += 2) {
        if (!strcmp(Symbols[i], s.c_str()))
            return std::string(Symbols[i + 1]);
    }
    return s;
}


std::string
mk_mathrm(const std::string &s)
{
    if (s.length() < 1) return s;
    return "\\mathrm{" + s + "}";
}



std::string
mk_compact(const std::string &s)
{
    if (s.length() <= 1) return s;
    unsigned i, n = s.length();
    std::string r;
    r = '{';
    r += s[0];
    for (i = 1; i < n; ++i) { r += "\\!"; r += s[i]; }
    r += '}';
    return r;
}

} /* namespace */




std::string
symbolic::internal::str2tex(const std::string &s, bool compact)
{
    unsigned i, n = s.length();
    if (n == 1) {
        if (compact) return s;
        else return mk_mathrm(s);
    }

    // decompose string
    std::vector<std::string> p;
    std::string buf;
    for (i = 0; i < n; ++i) {
        char c = s[i];
        if (c == '_') {
            p.push_back(buf);
            buf = "";
            // buf.clear();
        } else {
            buf += c;
        }
    }
    p.push_back(buf);

    // modify pieces
    n = p.size();
    buf = try_Greeksymb(p[0]);
    if (compact) {
        if (buf[0] != '\\') p[0] = mk_compact(buf); else p[0] = buf;
    } else {
        if (buf[0] != '\\') p[0] = mk_mathrm(buf); else p[0] = buf;
    }
    --n;
    i = n;
    while (p[i] == "") { p[i] = "\\prime"; --i;}
    for (i = 1; i <= n; ++i) {
        buf = try_Greeksymb(p[i]);
        if (buf[0] != '\\') p[i] = mk_mathrm(buf); else p[i] = buf;
    }

    // finalise
    ++n;
    if (n == 1) return p[0];
    buf = p[0];
    for (i = 1; i < n; ++i) {
        if ((p[i - 1] != "\\prime") && (p[i] != "")) buf += "^{";
        buf += p[i];
    }
    for (i = 1; i < n; ++i) {
        if ((p[i - 1] != "\\prime") && (p[i] != "")) buf += "}";
    }

    return buf;
}






std::string
symbolic::internal::str2tex2(const std::string &s)
{
    // return s;

    unsigned i, n = s.length();
    if (n == 1) return s;

    // decompose string
    std::vector<std::string> p;
    std::string buf;
    for (i = 0; i < n; ++i) {
        char c = s[i];
        if (c == '_') {
            p.push_back(buf);
            buf = "";
            // buf.clear();
        } else {
            buf += c;
        }
    }
    p.push_back(buf);

    // modify pieces
    n = p.size();
    for (i = 0; i < n; ++i) {
        buf = try_Greeksymb(p[i]);
        if (buf[0] != '\\') p[i] = buf; else p[i] = "$^" + buf + "$";
    }

    // finalise
    if (n == 1) return p[0];

    buf = p[0];
    for (i = 1; i < n; ++i) {
        if (p[i][0] == '$') {
            buf += p[i];
        } else {
            buf += ' ';
            buf += p[i];
        }
    }

    return buf;
}



#ifdef R_DLL

#include <R.h>
#include <Rcpp.h>

RcppExport
SEXP
str2tex(SEXP s)
{
    return Rcpp::wrap(symbolic::internal::str2tex(Rcpp::as<std::string>(s)));
}


RcppExport
SEXP
str2tex2(SEXP s)
{
    return Rcpp::wrap(symbolic::internal::str2tex2(Rcpp::as<std::string>(s)));
}

#endif /* R_DLL */


