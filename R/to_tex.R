# ###################################################################
# (c) Kancelaria Prezesa Rady Ministrów 2012-2015                   #
# Treść licencji w pliku 'LICENCE'                                  #
#                                                                   #
# (c) Chancellery of the Prime Minister 2012-2015                   #
# Licence terms can be found in the file 'LICENCE'                  #
#                                                                   #
# Authors: Karol Podemski, Kaja Retkiewicz-Wijtiwiak                #
# ###################################################################
# Functions for dumping the results to LaTeX
# ###################################################################

# ###################################################################
# The to_tex  function converts gecon_model class elements
# into TeX output
# ###################################################################
# Input 
#   arg - a list, matrix or vector that has to be converted to TeX
#   filenam - the output filename
#   rownam - a list with user defined rownames, when not set
#            rownames are created based on matrix rownames
#   colnam - a list with user defined colnames, when not set
#            colnames are created based on matrix rownames
#   section - user defined section names
#   to_table - if set to TRUE, the output is printed as a table,
#              otherwise it is displayed as a matrix
# Output
#   a tex file with the output
# ###################################################################
to_tex <- function(arg, filenam, rownam = NULL, colnam = NULL,
                   section, to_table = TRUE) 
{
    math_mode <- TRUE
    if (!to_table) math_mode <- FALSE
    
    # prepares output when the input is a matrix or vector
    if ((is.matrix(arg) || is.vector(arg)) && !is.list(arg)) {
        output <- vector("list", 2)
        output[[1]] <- paste('\n\\section{', section, "}", sep = '')
        if (!is.matrix(arg)) {
            if (is.null(rownam)) {
                names_word <- names(arg)
                if (is.null(names_word)) {
                    names_word <- rep(' ', length(arg))
                }                
                rownam <- string_to_tex(names_word, math_mode = math_mode)
            }
            if (is.null(colnam)) 
                colnam <- ' '
            arg <- matrix(arg, length(arg))
        }
        
        if (is.null(rownam)) {
            rownam <- string_to_tex(rownames(arg), section = FALSE, 
                                    math_mode = math_mode)
        }
        
        if (is.null(colnam)) {
            colnam <- string_to_tex(colnames(arg), section = FALSE, 
                                    math_mode = math_mode)
        }
        
        output[[2]] <- to_tex_matrix(arg, rownam, colnam, to_table)
        output <- unlist(output)        
    }
    
    # prepares output when the input is a list
    if (is.list(arg)) {
        l_length <- length(arg)
        l_rownam <- length(rownam)
        l_colnam <- length(colnam)
        
        # when the list of names is not specified,
        # rownames and colnames from matrices are read
        if ((l_length != l_rownam) ) {
            # retriving names
            rownam <- lapply(arg,  rownames)

            
            # handling missing cases
            # finding them
            empty_rownam <- which(unlist(lapply(rownam, is.null)))


            # filling rownames
            if (length(empty_rownam)) {
                for (i in 1:length(empty_rownam)) {
                    if (is.matrix(arg[[empty_rownam[i]]])) {
                        rownam[[empty_rownam[i]]] <- 
                                rep(' ', nrow(arg[[empty_rownam[i]]]))
                    } else {
                        arg[[empty_rownam[i]]] <- 
                            matrix(arg[[empty_rownam[i]]], length(arg[[empty_rownam[i]]]))
                        rownam[[empty_rownam[i]]] <- 
                            length(arg[[empty_rownam[i]]])
                    }
                }
            }
            
            rownam <- lapply(rownam, string_to_tex, section = FALSE, math_mode)
        }
        
        # when the list of names is not specified,
        # rownames and colnames from matrices are read
        if ((l_length != l_colnam)) {
            # retriving names
            colnam <- lapply(arg, colnames)
            
            # handling missing cases
            # finding them
            empty_colnam <- which(unlist(lapply(colnam, is.null)))
           
            # filling colnames
            if (length(empty_colnam)) {
                for (i in 1:length(empty_colnam)) {
                    if (is.matrix(arg[[empty_colnam[i]]])) {
                        colnam[[empty_colnam[i]]] <- 
                                rep(' ', nrow(arg[[empty_colnam[i]]]))
                    }
                }
            }
            colnam <- lapply(colnam, string_to_tex, section = FALSE, math_mode)
        }
        ########
        output <- vector("list", (l_length * 2 + 1))
        output[[1]] <- paste('\n\\section{', 
                              section[[1]],
                              "}", sep='')  
        for (i in 1:l_length) {
            if (length(section) == (l_length + 1)) {
                output[[(2 * i) - 1 + 1]] <- paste('\n\\subsection{', 
                                             string_to_tex(section[[i + 1]], section = TRUE, math_mode),
                                             "}", sep = '')
            } else {
                output[[(2 * i) - 1 + 1]] <- paste('\n\\subsection{', 
                                             string_to_tex(names(arg)[[i]], section = TRUE, math_mode),
                                             "}", sep = '')
            }                             
            output[[(2 * i) + 1]] <- to_tex_matrix(arg[[i]], rownam[[i]],
                                                   colnam[[i]], to_table)
        }
        output <- unlist(output)
    }
    

    if (file.exists(filenam)) {
        write(output, filenam  , sep = '\n', append = TRUE)  
        cat("results written to: \n")
        cat(filenam, "\n")
    } else {
        write(output, filenam  , sep = '\n') 
        cat("results written to: \n")
        cat(filenam, "\n")        
    }
}

# ###################################################################
# The to_tex_matrix function converts matrices or vectors into arrays
# ###################################################################
# Input 
#   arg - a matrix or vector that has to be converted to TeX
#   rowname - a list with user defined rownames, when not set
#             rownames are created based on matrix rownames
#   colname - a list with user defined colnames, when not set
#             colnames are created based on matrix colnames
#   to_table - logical, if set to FALSE matrices are displayed as 
#              TeX bordermatrices
# Output
#   a vector with TeX output lines
# ###################################################################
to_tex_matrix <- function(arg, rowname, colname, to_table = TRUE) 
{
    # converts matrix into tex output
    n_rows <- dim(arg)[1]
    n_cols <- dim(arg)[2]
    if (is.numeric(arg))
        arg <- round(arg, 4)

    if (to_table) {
        # creating output
        output <- vector(length = n_rows + 5)
        # initial tags
        output[1] <- ""
        output[2] <- paste(#'\\begin{table}[h]\n',
                           #'\\centering\n',
                           '\\begin{tabular}{|',
                           paste(rep('c', n_cols + 1), collapse='|'),
                           '|}', sep="")
        output[3] <- paste("\\hline ",
                           paste(c(" ", colname), collapse=' & '),
                            "\\\\")
        # the body of the table
        for (i in 1:n_rows) {
            output[i + 3] <- paste("\\hline ", 
                                   paste(c(rowname[i], arg[i, ]), 
                                         collapse=' & '),
                                   "\\\\")
        }
        # ending tags
        output[n_rows + 4] <- '\\hline'
        output[n_rows + 5] <- paste('\\end{tabular}\n',
                                    #'\\end{table}\n',
                                    #'\\newpage', 
                                    sep = "")
    } else {
        output <- vector(length = n_rows + 4)
        output[1] <- paste(#'\\begin{table}[h]\n',
                           #'\\centering\n',
                           '$$', sep = "")
        output[2] <- "\\bordermatrix{"
        output[3] <- paste(paste(c("~", colname), collapse=' & '),
                           "\\cr")
        # the body
        for (i in 1:n_rows) {
            output[i + 3] <- paste(paste(c(rowname[i], arg[i, ]), collapse=' & '),
                                   "\\cr")
        }
        output[n_rows + 4] <-paste( '}$$',
                                    #'\\end{table}\n',
                                    #'\\newpage', 
                                    sep = "")
    }
    return(output)
}

# ###################################################################
# The string_to_tex function converts strings into valid TeX strings
# ###################################################################
# Input 
#   strng - a character or vector of characters
#   section - logical, indicates whether text should be treated as
#             title section (_ replaced with empty space) or as 
#             ordinary text
#   math_mode - a logical, if set to TRUE adds dollars before and 
#               after argument, otherwise does not add anything
# Output
#   a character or vector of characters with the correct TeX syntax
# ###################################################################
string_to_tex <- function(strng, section = FALSE, math_mode = TRUE) 
{
    if (!is.character(strng))
        stop('The str variable must be a character or vector of characters')
    
    strng <- matrix(strng)
    
    # str2tex function wrapper
    str2tex <- function(strn, math_mode = TRUE) 
    {
        if (strn != '') {
            if (math_mode)
                strn <- paste('$', .Call('str2tex', strn), '$', sep='')
            else 
                strn <- paste(.Call('str2tex', strn), sep='')
        }
        return(strn)
    }
    
    # str2tex2 function wrapper
    str2tex2 <- function(strn) 
    {
        if (strn != '')
            strn <- .Call('str2tex2', strn)
        return(strn)
    }
    
    # applying appropriate function
    if (!section) {
        strng <- apply(X = strng, MARGIN = 1, FUN = str2tex, math_mode)
    } else {
        strng <- str2tex2(strng)
    }
    
    return(as.character(strng))
}
