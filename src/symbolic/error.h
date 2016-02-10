/***********************************************************
 * (c) Kancelaria Prezesa Rady Ministrów 2012-2015         *
 * Treść licencji w pliku 'LICENCE'                        *
 *                                                         *
 * (c) Chancellery of the Prime Minister 2012-2015         *
 * License terms can be found in the file 'LICENCE'        *
 *                                                         *
 * Author: Grzegorz Klima                                  *
 ***********************************************************/

/** \file error.h
 * \brief Handling errors.
 */

#ifndef SYMBOLIC_ERROR_H

#define SYMBOLIC_ERROR_H

#include <stdexcept>
#include <iostream>
#include <cstdlib>
#include <utils.h>

#define INTERNAL_ERROR throw(std::runtime_error(std::string("internal error in file ") +\
                             __FILE__ + ", line " + symbolic::internal::num2str(__LINE__)));
#define USER_ERROR(x) throw(std::runtime_error(std::string(x) + " in file " + __FILE__ + ", line " + symbolic::internal::num2str(__LINE__)));


#endif /* SYMBOLIC_ERROR_H */
