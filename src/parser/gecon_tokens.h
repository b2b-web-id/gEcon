/*****************************************************************************
 * This file is a part of gEcon.                                             *
 *                                                                           *
 * (c) Chancellery of the Prime Minister of the Republic of Poland 2012-2015 *
 * (c) Grzegorz Klima, Karol Podemski, Kaja Retkiewicz-Wijtiwiak 2015-2018   *
 * License terms can be found in the file 'LICENCE'                          *
 *                                                                           *
 * Author: Grzegorz Klima                                                    *
 *****************************************************************************/

/** \file gecon_tokens.h
 * \brief Lexer tokens.
 */

#ifndef PARSER_GECON_TOKENS_H

#define PARSER_GECON_TOKENS_H

/// Allocate and set table
unsigned char** mk_tnames();
/// Free table
void free_tnames(unsigned char **tnames);

#endif /* PARSER_GECON_TOKENS_H */

