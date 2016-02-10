# ###################################################################
# (c) Kancelaria Prezesa Rady Ministrów 2012-2015                   #
# Treść licencji w pliku 'LICENCE'                                  #
#                                                                   #
# (c) Chancellery of the Prime Minister 2012-2015                   #
# Licence terms can be found in the file 'LICENCE'                  #
#                                                                   #
# Authors: Karol Podemski, Kaja Retkiewicz-Wijtiwiak                #
# ###################################################################
# Intro message
# ###################################################################

.onAttach <- function( ... )
{
  gEcon_lib <- dirname(system.file(package = "gEcon"))
  welcome_R <- .Call('welcome_R')
  packageStartupMessage(welcome_R)
}

