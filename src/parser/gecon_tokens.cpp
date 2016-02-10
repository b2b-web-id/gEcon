/***********************************************************
 * (c) Kancelaria Prezesa Rady Ministrów 2012-2015         *
 * Treść licencji w pliku 'LICENCE'                        *
 *                                                         *
 * (c) Chancellery of the Prime Minister 2012-2015         *
 * License terms can be found in the file 'LICENCE'        *
 *                                                         *
 * Author: Grzegorz Klima                                  *
 ***********************************************************/

/** \file gecon_tokens.cpp
 * \brief Lexer tokens.
 */

#include <gecon_tokens.h>
#include <cstring>
#include <cstdlib>

enum tokens {
    ACOS = 4,
    AND = 5,
    ASIN = 6,
    AT = 7,
    ATAN = 8,
    BACKSLASH = 9,
    BACKWARDCOMP = 10,
    BFALSE = 11,
    BLOCK = 12,
    BTRUE = 13,
    CALIBR = 14,
    CLETTER = 15,
    COLON = 16,
    COMMA = 17,
    COMMENT = 18,
    CONSTRAINTS = 19,
    CONTROLS = 20,
    COS = 21,
    COSH = 22,
    DBLCOLON = 23,
    DDOT = 24,
    DEFS = 25,
    DELTA = 26,
    DEQ = 27,
    DID = 28,
    DIDU = 29,
    DIV = 30,
    DOLLAR = 31,
    DOR = 32,
    DOUBLE = 33,
    DQUOTE = 34,
    E = 35,
    EQ = 36,
    ERF = 37,
    EXCLAM = 38,
    EXP = 39,
    FOCS = 40,
    ID = 41,
    IDS = 42,
    IDU = 43,
    INF = 44,
    INT = 45,
    LANDSCAPE = 46,
    LANGBR = 47,
    LATEX = 48,
    LBRACE = 49,
    LBRACK = 50,
    LEQ = 51,
    LOG = 52,
    LOGF = 53,
    LONG = 54,
    LPAREN = 55,
    MINUS = 56,
    MUL = 57,
    NEQ = 58,
    OBJ = 59,
    OPTS = 60,
    OR = 61,
    OUTPUT = 62,
    PLUS = 63,
    POW = 64,
    PROD = 65,
    QUESTION = 66,
    QUOTE = 67,
    R = 68,
    RANGBR = 69,
    RARROW = 70,
    RBRACE = 71,
    RBRACK = 72,
    RPAREN = 73,
    SEMI = 74,
    SETS = 75,
    SHOCKS = 76,
    SHORT = 77,
    SIN = 78,
    SINH = 79,
    SLETTER = 80,
    SQRT = 81,
    SS = 82,
    SUM = 83,
    TAN = 84,
    TANH = 85,
    TILDE = 86,
    TRYREDUCE = 87,
    UDID = 88,
    UID = 89,
    VERBOSE = 90,
    WS = 91,
    ZERO = 92,
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


