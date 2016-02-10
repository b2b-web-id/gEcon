/***********************************************************
 * (c) Kancelaria Prezesa Rady Ministrów 2012-2015         *
 * Treść licencji w pliku 'LICENCE'                        *
 *                                                         *
 * (c) Chancellery of the Prime Minister 2012-2015         *
 * License terms can be found in the file 'LICENCE'        *
 *                                                         *
 * Author: Grzegorz Klima                                  *
 ***********************************************************/

/** \file ex_base.cpp
 * \brief Abstract base class for expressions.
 */


#include <ex_base.h>
#include <ex_num.h>

using namespace symbolic;
using namespace symbolic::internal;


bool
ex_base::is0() const
{
    if (type() != NUM) return false;
    const ex_num *p = static_cast<const ex_num*>(this);
    return (p->val() == 0.);
}


bool
ex_base::is1() const
{
    if (type() != NUM) return false;
    const ex_num *p = static_cast<const ex_num*>(this);
    return (p->val() == 1.);
}


bool
ex_base::ism1() const
{
    if (type() != NUM) return false;
    const ex_num *p = static_cast<const ex_num*>(this);
    return (p->val() == -1.);
}


Number
ex_base::val() const
{
    if (type() != NUM) return 0./0.;
    return static_cast<const ex_num*>(this)->val();
}


#ifdef DEBUG_SYMBOLIC
void
ex_base::info() const
{
    std::cerr << "  pointer: " << this << '\n';
    unsigned t = type(), f = flag();
    std::cerr << "  type (" << t << "): ";
    switch (t) {
#ifdef EXPAND_CASE
#undef EXPAND_CASE
#endif
#define EXPAND_CASE(TYPE) \
        case TYPE: \
            std::cerr << #TYPE << ' '; \
            break;
        EXPAND_CASE(NUM)
        EXPAND_CASE(DELTA)
        EXPAND_CASE(SYMB)
        EXPAND_CASE(VART)
        EXPAND_CASE(EX)
        EXPAND_CASE(POW)
        EXPAND_CASE(FUN)
        EXPAND_CASE(ADD)
        EXPAND_CASE(MUL)
        EXPAND_CASE(SUM)
        EXPAND_CASE(PROD)
        EXPAND_CASE(IDX)
        default:
            std::cerr << t;
    }
    std::cerr << '\n';
    std::cerr << "  flag(s) (" << f << "): ";
#ifdef EXPAND_FLAG
#undef EXPAND_FLAG
#endif
#define EXPAND_FLAG(FLAG) \
    if (f & FLAG) std::cerr << #FLAG << ' ';
    EXPAND_FLAG(NUL)
    EXPAND_FLAG(EMPTY)
    EXPAND_FLAG(SINGLE)
    EXPAND_FLAG(SCAL1)
    EXPAND_FLAG(SPOW1)
    EXPAND_FLAG(HAST)
    EXPAND_FLAG(HASIDX)
    std::cerr << '\n';
    std::cerr << "  str(): \"" << this->str() << "\"\n";
    std::cerr << "  tex(): \"" << this->tex() << "\"\n\n";
}
#endif /* DEBUG_SYMBOLIC */













