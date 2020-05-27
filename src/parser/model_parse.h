/*****************************************************************************
 * This file is a part of gEcon.                                             *
 *                                                                           *
 * (c) Chancellery of the Prime Minister of the Republic of Poland 2012-2015 *
 * (c) Grzegorz Klima, Karol Podemski, Kaja Retkiewicz-Wijtiwiak 2015-2018   *
 * License terms can be found in the file 'LICENCE'                          *
 *                                                                           *
 * Author: Grzegorz Klima                                                    *
 *****************************************************************************/

/** \file model_parse.h
 * \brief Parsing the model.
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

