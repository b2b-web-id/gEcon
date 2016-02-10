/***********************************************************
 * (c) Kancelaria Prezesa Rady Ministrów 2012-2015         *
 * Treść licencji w pliku 'LICENCE'                        *
 *                                                         *
 * (c) Chancellery of the Prime Minister 2012-2015         *
 * License terms can be found in the file 'LICENCE'        *
 *                                                         *
 * Author: Grzegorz Klima                                  *
 ***********************************************************/

/** \file gecon_info.cpp
 * \brief Basic information about gEcon version.
 */


#include <gecon_info.h>
#ifdef R_DLL
#include <R.h>
#include <Rcpp.h>
#endif /* R_DLL */


std::string
gecon_ver_str()
{
    return "0.9.1 (2015-05-19)";
}


std::string
gecon_bug_str()
{
    return "Grzegorz Klima <gklima@users.sourceforge.net>";
}


std::string
gecon_web_str()
{
    return "http://gecon.r-forge.r-project.org/";
}


std::string
gecon_hello_str()
{
    return "This is gEcon version " + gecon_ver_str() + '\n'
           + gecon_web_str();
}


#ifdef R_DLL
RcppExport
SEXP
welcome_R()
{
    return Rcpp::wrap(gecon_hello_str());
}
#endif /* R_DLL */
















