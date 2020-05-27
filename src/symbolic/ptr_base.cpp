/*****************************************************************************
 * This file is a part of gEcon.                                             *
 *                                                                           *
 * (c) Chancellery of the Prime Minister of the Republic of Poland 2012-2015 *
 * (c) Grzegorz Klima, Karol Podemski, Kaja Retkiewicz-Wijtiwiak 2015-2018   *
 * License terms can be found in the file 'LICENCE'                          *
 *                                                                           *
 * Author: Grzegorz Klima                                                    *
 *****************************************************************************/

/** \file ptr_base.cpp
 * \brief Reference counted pointer to abstract base class
 * for expressions.
 */


#include <ptr_base.h>
#include <ex_num.h>
#include <ex_delta.h>
#include <ex_symb.h>
#include <ex_symbidx.h>
#include <ex_vart.h>
#include <ex_vartidx.h>
#include <ex_add.h>
#include <ex_mul.h>
#include <ex_pow.h>
#include <ex_func.h>
#include <ex_e.h>
#include <ex_sum.h>
#include <ex_prod.h>
#include <ex_idx.h>
#include <error.h>


using namespace symbolic::internal;



// Constructor
ptr_base::ptr_base(ex_base *p)
{
    m_p = p;
    if (m_p) {
        m_rc = new int;
        *m_rc = 1;
    } else {
        m_rc = 0;
    }
}


// Copy constructor
ptr_base::ptr_base(const ptr_base& p)
{
    m_p = p.m_p;
    m_rc = p.m_rc;
    if (m_rc) ++(*m_rc);
}


void
ptr_base::dec_rc()
{
    if (m_rc) {
        if (--(*m_rc) == 0) {
            delete m_rc;
            switch (m_p->type()) {
#ifdef EXPAND_CASE
#undef EXPAND_CASE
#endif
#define EXPAND_CASE(TYPE, CLASS) \
                case TYPE: \
                    CLASS::destroy(m_p); \
                    break;
                EXPAND_CASE(NUM, ex_num)
                EXPAND_CASE(DELTA, ex_delta)
                EXPAND_CASE(SYMB, ex_symb)
                EXPAND_CASE(SYMBIDX, ex_symbidx)
                EXPAND_CASE(VART, ex_vart)
                EXPAND_CASE(VARTIDX, ex_vartidx)
                EXPAND_CASE(EX, ex_e)
                EXPAND_CASE(POW, ex_pow)
                EXPAND_CASE(FUN, ex_func)
                EXPAND_CASE(ADD, ex_add)
                EXPAND_CASE(MUL, ex_mul)
                EXPAND_CASE(SUM, ex_sum)
                EXPAND_CASE(PROD, ex_prod)
                EXPAND_CASE(IDX, ex_idx)
                default:
                    INTERNAL_ERROR;
            }
        }
    }
}

// Destructor
ptr_base::~ptr_base()
{
    dec_rc();
}


// Assignment
ptr_base&
ptr_base::operator=(const ptr_base& p)
{
    if (m_p != p.m_p) {
        dec_rc();
        m_p = p.m_p;
        m_rc = p.m_rc;
        if (m_rc) ++(*m_rc);
    }
    return *this;
}





