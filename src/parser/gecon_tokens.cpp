/*****************************************************************************
 * This file is a part of gEcon.                                             *
 *                                                                           *
 * (c) Chancellery of the Prime Minister of the Republic of Poland 2012-2015 *
 * (c) Grzegorz Klima, Karol Podemski, Kaja Retkiewicz-Wijtiwiak 2015-2018   *
 * License terms can be found in the file 'LICENCE'                          *
 *                                                                           *
 * Author: Grzegorz Klima                                                    *
 *****************************************************************************/

/** \file gecon_tokens.cpp
 * \brief Lexer tokens.
 */

#include <gecon_tokens.h>
#include <cstring>
#include <cstdlib>

enum tokens {
    ACOS=4,
    AND=5,
    ASIN=6,
    AT=7,
    ATAN=8,
    BACKSLASH=9,
    BFALSE=10,
    BLOCK=11,
    BTRUE=12,
    CALIBR=13,
    CLETTER=14,
    COLON=15,
    COMMA=16,
    COMMENT=17,
    CONSTRAINTS=18,
    CONTROLS=19,
    COS=20,
    COSH=21,
    DBLCOLON=22,
    DDOT=23,
    DEFS=24,
    DELTA=25,
    DEQ=26,
    DID=27,
    DIDU=28,
    DIV=29,
    DOLLAR=30,
    DOR=31,
    DOUBLE=32,
    DQUOTE=33,
    E=34,
    EQ=35,
    EXCLAM=36,
    EXP=37,
    FOCS=38,
    ID=39,
    IDS=40,
    IDU=41,
    INF=42,
    INT=43,
    JACOBIAN=44,
    LANDSCAPE=45,
    LANGBR=46,
    LATEX=47,
    LBRACE=48,
    LBRACK=49,
    LEQ=50,
    LOG=51,
    LOGF=52,
    LONG=53,
    LPAREN=54,
    MINUS=55,
    MUL=56,
    NEQ=57,
    OBJ=58,
    OPTS=59,
    OR=60,
    OUTPUT=61,
    PLUS=62,
    PNORM=63,
    POW=64,
    PROD=65,
    QUESTION=66,
    QUOTE=67,
    R=68,
    RANGBR=69,
    RARROW=70,
    RBRACE=71,
    RBRACK=72,
    RCPP=73,
    RPAREN=74,
    SEMI=75,
    SETS=76,
    SHOCKS=77,
    SHORT=78,
    SILENT=79,
    SIN=80,
    SINH=81,
    SLETTER=82,
    SQRT=83,
    SS=84,
    SUM=85,
    TAN=86,
    TANH=87,
    TILDE=88,
    TRYREDUCE=89,
    UDID=90,
    UID=91,
    VERBOSE=92,
    WARNINGS=93,
    WS=94,
    ZERO=95,
    END_LIST
};



unsigned char**
mk_tnames()
{
    unsigned char **tnames = (unsigned char**) malloc(sizeof(unsigned char*) * END_LIST);
    for (unsigned i = 0; i < END_LIST; ++i) {
        unsigned sz;
        switch (i) {
#define EXPAND_CASE(nn, tt) \
case nn: sz = strlen(tt) + 1; tnames[i] = (unsigned char*) malloc(sz); memcpy(tnames[i], tt, sz); \
break;
            EXPAND_CASE(TILDE, "\'~\'")
            EXPAND_CASE(QUESTION, "\'?\'")
            EXPAND_CASE(EXCLAM, "\'!\'")
            EXPAND_CASE(DOLLAR, "\'$\'")
            EXPAND_CASE(AT, "\'@\'")
            EXPAND_CASE(AND, "\'&\'")
            EXPAND_CASE(OR, "\'|\'")
            EXPAND_CASE(DOR, "\'||\'")
            EXPAND_CASE(SEMI, "\';\'")
            EXPAND_CASE(COLON, "\':\'")
            EXPAND_CASE(DBLCOLON, "\'::\'")
            EXPAND_CASE(DDOT, "\'..\'")
            EXPAND_CASE(COMMA, "\',\'")
            EXPAND_CASE(RARROW, "\'->\'")
            EXPAND_CASE(PLUS, "\'+\'")
            EXPAND_CASE(MINUS, "\'-\'")
            EXPAND_CASE(MUL, "\'*\'")
            EXPAND_CASE(DIV, "\'/\'")
            EXPAND_CASE(POW, "\'^\'")
            EXPAND_CASE(EQ, "\'=\'")
            EXPAND_CASE(DEQ, "\'==\'")
            EXPAND_CASE(NEQ, "\'!=\'")
            EXPAND_CASE(LEQ, "\'<=\'")
            EXPAND_CASE(QUOTE, "\'\'\'")
            EXPAND_CASE(DQUOTE, "\'\"\'")
            EXPAND_CASE(BACKSLASH, "\'\\\'")
            EXPAND_CASE(LBRACE, "\'{\'")
            EXPAND_CASE(RBRACE, "\'}\'")
            EXPAND_CASE(LPAREN, "\'(\'")
            EXPAND_CASE(RPAREN, "\')\'")
            EXPAND_CASE(LBRACK, "\'[\'")
            EXPAND_CASE(RBRACK, "\']\'")
            EXPAND_CASE(LANGBR, "\'<\'")
            EXPAND_CASE(RANGBR, "\'>\'")
            default:
                tnames[i] = (unsigned char*) malloc(1);
                tnames[i][0] = '\0';
        }
    }
    return tnames;
}


void
free_tnames(unsigned char **tnames)
{
    for (unsigned i = 0; i < END_LIST; ++i) {
        free(tnames[i]);
    }
    free(tnames);
}


