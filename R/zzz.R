# ############################################################################
# This file is a part of gEcon.                                              #
#                                                                            #
# (c) Chancellery of the Prime Minister of the Republic of Poland 2012-2015  #
# (c) Grzegorz Klima, Karol Podemski, Kaja Retkiewicz-Wijtiwiak 2015-2018    #
# (c) Karol Podemski, Kaja Retkiewicz-Wijtiwiak 2018-2019                    #
# License terms can be found in the file 'LICENCE'                           #
#                                                                            #
# Authors: Grzegorz Klima, Karol Podemski, Kaja Retkiewicz-Wijtiwiak         #
# ############################################################################
# Intro message and working environment for storing compiled model functions
# ############################################################################

.onAttach <- function( ... )
{
    welcome_R <- .Call("welcome_R")
    packageStartupMessage(welcome_R)
}

# ############################################################################
# This creates new environment inside gEcon namespace for storing 
# compiled model functions.
# ############################################################################
.gEcon_env <- new.env()
