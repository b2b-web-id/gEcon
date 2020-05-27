/*****************************************************************************
 * This file is a part of gEcon.                                             *
 *                                                                           *
 * (c) Chancellery of the Prime Minister of the Republic of Poland 2012-2015 *
 * (c) Grzegorz Klima, Karol Podemski, Kaja Retkiewicz-Wijtiwiak 2015-2018   *
 * License terms can be found in the file 'LICENCE'                          *
 *                                                                           *
 * Author: Grzegorz Klima                                                    *
 *****************************************************************************/

/** \file decl.h
 * \brief Forward declarations and enums.
 */

#include <triplet.h>
#include <string>
#include <algorithm>
#include <vector>
#include <set>
#include <map>


#ifndef SYMBOLIC_DECL_H

#define SYMBOLIC_DECL_H


namespace symbolic {


// Forward decl
class ex;
class idx_set;
class idx_ex;

struct less_ex {
    bool operator()(const symbolic::ex &a, const symbolic::ex &b) const;
    typedef ex first_argument_type;
    typedef ex second_argument_type;
    typedef bool result_type;
};

typedef std::vector<ex> vec_ex;
typedef std::set<ex, less_ex> set_ex;
typedef std::map<ex, ex, less_ex> map_ex_ex;
typedef std::map<ex, int, less_ex> map_ex_int;
typedef std::map<ex, std::string, less_ex> map_ex_str;
typedef std::map<std::string, std::string> map_str_str;


namespace internal {

// Forward decl
class ex_base;
class ptr_base;
class ex_num;
class ex_delta;
class ex_symb;
class ex_symbidx;
class ex_vart;
class ex_vartidx;
class ex_pow;
class ex_func;
class ex_e;
class num_ex_pair_vec;
class ex_add;
class ex_mul;
class idx_set_impl;
class ex_sum;
class ex_prod;
class ex_idx;


/// Expression types
enum ex_type {
    // flags
    NUL     =        0,
    EMPTY   =      0x1,
    SINGLE  =      0x2,
    SCAL1   =      0x4,
    SPOW1   =      0x8,
    HAST    =     0x10,
    HASIDX  =     0x20,
    // atomic expressions
    NUM     =    0x100,
    DELTA   =    0x200,
    SYMB    =    0x400,
    SYMBIDX =    0x800,
    VART    =   0x1000,
    VARTIDX =   0x2000,
    // compound expressions
    FUN     =   0x4000,
    ADD     =   0x8000,
    MUL     =  0x10000,
    POW     =  0x20000,
    EX      =  0x40000,
    SUM     =  0x80000,
    PROD    = 0x100000,
    IDX     = 0x200000
};



/// Function codes
enum func_code {
    EXP = 0, LOG,
    SIN, COS, TAN,
    ASIN, ACOS, ATAN,
    SINH, COSH, TANH,
    ERF, PNORM
};


/// Printing flags
enum print_flag {
    DEFAULT            =    0,
    CONVERT_T          =  0x1,
    DROP_T             =  0x2,
    CONVERT_IDX        =  0x4,
    DROP_IDX           =  0x8,
    DROP_QUOTES        = 0x10,
    INDEXING_ONLY      = 0x20,
    DROP_INDEXING      = 0x40
};


/// Searching flags
enum search_flag {
    EXACT_T =   0,
    ANY_T   = 0x1,
    DIFF_T  = 0x2,
    LEAD_T  = 0x4,
    LAG_T   = 0x8
};


} /* namespace internal */
} /* namespace symbolic */

#endif /* SYMBOLIC_DECL_H */
