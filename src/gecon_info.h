/*****************************************************************************
 * This file is a part of gEcon.                                             *
 *                                                                           *
 * (c) Chancellery of the Prime Minister of the Republic of Poland 2012-2015 *
 * (c) Grzegorz Klima, Karol Podemski, Kaja Retkiewicz-Wijtiwiak 2015-2018   *
 * (c) Karol Podemski, Kaja Retkiewicz-Wijtiwiak 2018-2019                   *
 *                                                                           *
 * License terms can be found in the file 'LICENCE'                          *
 *                                                                           *
 * Author: Grzegorz Klima                                                    *
 *****************************************************************************/

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

/// gEcon copyrights.
std::string gecon_copyrights();

/// gEcon hello string.
std::string gecon_hello_str();



#endif /* GECON_INFO_H */
