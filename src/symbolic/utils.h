/*****************************************************************************
 * This file is a part of gEcon.                                             *
 *                                                                           *
 * (c) Chancellery of the Prime Minister of the Republic of Poland 2012-2015 *
 * (c) Grzegorz Klima, Karol Podemski, Kaja Retkiewicz-Wijtiwiak 2015-2018   *
 * License terms can be found in the file 'LICENCE'                          *
 *                                                                           *
 * Author: Grzegorz Klima                                                    *
 *****************************************************************************/

/** \file utils.h
 * \brief Utilities.
 */

#ifndef SYMBOLIC_UTILS_H

#define SYMBOLIC_UTILS_H

#include <string>


namespace symbolic {
namespace internal {

/// Convert to string
std::string num2str(double n);
/// Convert to string
std::string num2str(int n, int d = 0);
/// Convert to string
std::string num2str(unsigned n, int d = 0);
/// Convert to LaTeX string
std::string num2tex(double n);
/// Convert function name to LaTeX
std::string func2tex(const std::string&);
/// Convert string containing variable name to valid LaTeX name by removing '_'
std::string str2tex(const std::string&, bool compact = true);
/// Convert string containing variable name to valid LaTeX name by removing '_'
std::string str2tex2(const std::string&);

} /* namespace internal */
} /* namespace symbolic */


#endif /* SYMBOLIC_UTILS_H */
