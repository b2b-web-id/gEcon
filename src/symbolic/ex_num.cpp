/*****************************************************************************
 * This file is a part of gEcon.                                             *
 *                                                                           *
 * (c) Chancellery of the Prime Minister of the Republic of Poland 2012-2015 *
 * (c) Grzegorz Klima, Karol Podemski, Kaja Retkiewicz-Wijtiwiak 2015-2018   *
 * License terms can be found in the file 'LICENCE'                          *
 *                                                                           *
 * Author: Grzegorz Klima                                                    *
 *****************************************************************************/

/** \file ex_num.cpp
 * \brief Numbers.
 */


#include <ex_num.h>
#include <utils.h>
#include <cmp.h>
#include <climits>

using namespace symbolic;
using namespace symbolic::internal;


ptr_base
ex_num::create(const Number &n)
{
    return ptr_base(new ex_num(n));
}


void
ex_num::destroy(ex_base *ptr)
{
    delete ptr;
}


int
ex_num::compare(const ex_num &b) const
{
    return compareT(m_val, b.m_val);
}


std::string
ex_num::str(int pflag, bool c_style) const
{
    if (pflag & INDEXING_ONLY) return std::string();
    return m_val.str();
}

std::string
ex_num::strmap(const map_str_str&, bool c_style) const
{
    return m_val.str();
}

std::string
ex_num::tex(int pflag) const
{
    if (pflag & INDEXING_ONLY) return std::string();
    return m_val.tex();
}


int
ex_num::get_lag_max(bool) const
{
    return INT_MIN;
}


int
ex_num::get_lag_min(bool) const
{
    return INT_MAX;
}



ptr_base
ex_num::diff(const ptr_base&) const
{
    return zero();
}


bool
ex_num::has(const ptr_base&, search_flag, bool) const
{
    return false;
}


bool
ex_num::hasidx(int) const
{
    return false;
}

