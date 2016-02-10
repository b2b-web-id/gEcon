/***********************************************************
 * (c) Kancelaria Prezesa Rady Ministrów 2012-2015         *
 * Treść licencji w pliku 'LICENCE'                        *
 *                                                         *
 * (c) Chancellery of the Prime Minister 2012-2015         *
 * License terms can be found in the file 'LICENCE'        *
 *                                                         *
 * Author: Grzegorz Klima                                  *
 ***********************************************************/

/** \file model_parse.h
 * \brief Parsing model.
 */

#include <string>

#ifndef PARSER_MODEL_PARSE_H

#define PARSER_MODEL_PARSE_H

/// Parse model in file
void model_parse(const char *fname);

/// Report errors
void report_errors(const std::string&);

/// Report warnings
void report_warns(const std::string&);

/// Write information
void write_info(const std::string&);


#endif /* PARSER_MODEL_PARSE_H */

