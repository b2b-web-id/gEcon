/***********************************************************
 * (c) Kancelaria Prezesa Rady Ministrów 2012-2015         *
 * Treść licencji w pliku 'LICENCE'                        *
 *                                                         *
 * (c) Chancellery of the Prime Minister 2012-2015         *
 * License terms can be found in the file 'LICENCE'        *
 *                                                         *
 * Author: Grzegorz Klima                                  *
 ***********************************************************/

/** \file gecon_info.h
 * \brief Basic information about gEcon version.
 */


#ifndef GECON_INFO_H

#define GECON_INFO_H

#include <string>


/// gEcon version string.
std::string gecon_ver_str();

/// Where to report gEcon bugs?
std::string gecon_bug_str();

/// gEcon web page.
std::string gecon_web_str();

/// gEcon hello string.
std::string gecon_hello_str();



#endif /* GECON_INFO_H */
