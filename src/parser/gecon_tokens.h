/***********************************************************
 * (c) Kancelaria Prezesa Rady Ministrów 2012-2015         *
 * Treść licencji w pliku 'LICENCE'                        *
 *                                                         *
 * (c) Chancellery of the Prime Minister 2012-2015         *
 * License terms can be found in the file 'LICENCE'        *
 *                                                         *
 * Author: Grzegorz Klima                                  *
 ***********************************************************/

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

