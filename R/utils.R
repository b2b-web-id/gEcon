# ############################################################################
# This file is a part of gEcon.                                              #
#                                                                            #
# (c) Chancellery of the Prime Minister of the Republic of Poland 2012-2015  #
# (c) Grzegorz Klima, Karol Podemski, Kaja Retkiewicz-Wijtiwiak 2015-2018    #
# License terms can be found in the file 'LICENCE'                           #
#                                                                            #
# Authors: Grzegorz Klima, Karol Podemski, Kaja Retkiewicz-Wijtiwiak         #
# ############################################################################
# Utility functions
# ############################################################################

# ############################################################################
# The add_r_ext function attaches the .model.R extension if the model file
# was specified without it
# ############################################################################
# Input
#   filename - character string, filename
# Output
#   character string with .model.R extension, if filename was specified
#   without it
# ############################################################################
add_r_ext <- function(filename)
{
    if (length(grep(pattern="(?i)\\.model\\.R$", x = filename, perl = TRUE))) {
        return (filename)
    } else {
        return (paste(filename, ".model.R", sep = ""))
    }
}

# ############################################################################
# The add_gcn_ext function attaches .gcn extension if the model file
# was specified without it
# ############################################################################
# Input
#   filename - character string, filename
# Output
#   character string with .gcn extension, if model filename was specified
#   without it
# ############################################################################
add_gcn_ext <- function(filename)
{
    if (length(grep(pattern="(?i)\\.gcn$", x = filename, perl = TRUE))) {
        return (filename)
    } else {
        return (paste(filename, ".gcn", sep = ""))
    }
}

# ############################################################################
# The rm_gcn_ext function removes the .gcn extension from string
# ############################################################################
# Input
#   filename - character string, filename
# Output
#   string without .gcn extension
# ############################################################################
rm_gcn_ext <- function(filename)
{
    filename <- gsub(pattern = "(?i)\\.gcn$", replacement = "",
                     x = filename, perl = TRUE)
    return (filename)
}

# ############################################################################
# The get_path function retrieves path without the filename from file path
# ############################################################################
# Input
#   filename - character string, filename
# Output
#   path without the filename
# ############################################################################
get_path <- function(filename)
{
    path <- gsub(pattern = "(.*)(/[^/]+)", x = filename, replacement = "\\1")
    return (path)
}


# ############################################################################
# The list2str function returns a character string with list of inidices,
# variables, parameters, etc.
# ############################################################################
# Input
#   names - character vector
#   quote - character (quotation mark)
# Output
#   character string with list of inidices, variables, parameters, etc.
# ############################################################################
list2str <- function(names, quote = "")
{
    res <- paste0(quote, names, quote)
    return (paste(res, collapse = ", "))
}


# ############################################################################
# The list2str2 function returns a character string with list of inidices,
# variables, parameters, etc.
# ############################################################################
# Input
#   names - character vector
#   quote - character (quotation mark)
#   lineln - line length for formatting
#   tab - tab length
# Output
#   character string with list of inidices, variables, parameters, etc.
# ############################################################################
list2str2 <- function(names, quote = "", lineln = 80, tab = 4)
{
    names <- paste0(quote, names, quote)
    tabs <- paste0(rep(" ", tab), collapse = "")
    res <- tabs
    while(length(names)) {
        l <- nchar(names) + nchar(quote) + 2
        cl <- cumsum(l)
        i <- which(cl <= lineln - tab)
        res <- paste0(res, paste0(names[i], collapse = ", "))
        names <- names[-i]
        if (length(names)) {
            res <- paste0(res, ",\n", tabs)
        }
    }
    return (res)
}


# ############################################################################
# The list2ind function performs check and returns a vector of indices given
# variables, parameters, or shock names or indices
# ############################################################################
# Input
#   namesind - character vector (names) or integer vector (indices)
#   modelnames - character vector, model variables, parameters, or shocks
#   what - character string ("variable", "parameter", or "shock")
# Output
#   integer vector with model variable, parameter, or shock indices
# ############################################################################
list2ind <- function(namesind, modelnames, what, what2 = paste0(what, "s"))
{
    if(is(namesind, "numeric")) {
        ind <- as.integer(namesind)
        if (any(is.na(ind))) {
            stop(paste0("invalid ", what, " index"), call. = FALSE)
        }
        N <- length(modelnames)
        inv <- (ind > N) | (ind < 1)
        if (any(inv)) {
            inv <- which(inv)
            stop(paste0("invalid ", what, " index (indices): ",
                        list2str(ind[inv])), call. = FALSE)
        }
        dupl_check <- duplicated(ind)
        if (any(dupl_check)) {
            dupl_values <- unique(ind[dupl_check])
            warning(paste0("duplicated " , what, " index (indices): ",
                           list2str(dupl_values)), call. = FALSE)
            ind <- unique(ind)
        }
    } else if(is(namesind, "character"))  {
        names <- namesind
        not_in_model <- names[which(!(names %in% modelnames))]
        if (length(not_in_model)) {
            stop(paste0("the following are not model ", what, "s: ", list2str(not_in_model, "\"")),
                 call. = FALSE)
        }
        dupl_check <- duplicated(names)
        if (any(dupl_check)) {
            dupl_values <- unique(names[dupl_check])
            warning(paste0("duplicated ", what, " name(s): ", list2str(dupl_values, "\"")),
                    call. = FALSE)
            names <- unique(names)
        }
        ind <- match(names, modelnames)
    } else {
        stop(paste0("'", what2, "' argument should be character or integer vector or value"),
             call. = FALSE)
    }

    return (ind)
}


# ############################################################################
# Function round2zero rounds results to zero with given tolerance
# ############################################################################
# Input
#   x - a vector or a matrix
#   tol -  the tolerance level (default 1e-8)
# Output
#   x with small (in absolute value) entries rounded to zero
# ############################################################################
round2zero <- function(x, tol)
{
    ind <- which(abs(x) < tol)
    if (length(ind)) {
        x[ind] <- 0
    }
    return (x)
}


# ############################################################################
# Function numnoun constructs a string of the form "1 foo" or "2 foos"
# ############################################################################
# Input
#   n - number
#   what - noun
# Output
#   character string - number plus noun in singular of plural form
# ############################################################################
numnoun <- function(n, what)
{
    if (n == 1) {
        return (paste0("1 ", what))
    } else {
        return (paste0(n, " ", what, "s"))
    }
}



# ############################################################################
# The df2textable function writes matrices or dataframes to LaTeX tables
# ############################################################################
# Input
#   arg - matrix or dataframe to be converted to TeX
#   col_sep - logical, determines whether columns are separated by vertical
#                      lines
#   row_sep - logical, determines whether rows are separated by horizontal
#                      lines
# Output
#   character string with LaTeX output
# ############################################################################
df2textable <- function(arg, col_sep = TRUE, row_sep = FALSE)
{
    cnames <- colnames(arg)
    rnames <- rownames(arg)
    n_rows <- dim(arg)[1]
    n_cols <- dim(arg)[2]
    if (col_sep) {
        output <- paste0("\n\\begin{tabular}{",
                         paste0(rep("c", n_cols + 1), collapse = "|"),
                         "|}\n")
    } else {
        output <- paste0("\n\\begin{tabular}{c|",
                         paste0(rep("c", n_cols), collapse = ""),
                         "|}\n")
    }
    cnames <- colnames(arg)
    rnames <- rownames(arg)
    n_rows <- dim(arg)[1]
    output <- paste0(output, paste(c(" ", cnames), collapse = " & "),
                     "\\\\\n\\hline\n")
    if (row_sep) {
        for (i in 1:n_rows) {
            output <- paste0(output,
                             paste(c(rnames[i], arg[i, ]), collapse = " & "),
                             " \\\\\n\\hline\n")
        }
    } else {
        for (i in 1:n_rows) {
            output <- paste0(output,
                             paste(c(rnames[i], arg[i, ]), collapse = " & "),
                             " \\\\\n")
        }
        output <- paste0(output, "\\hline\n")
    }
    output <- paste0(output, "\\end{tabular}\n\n")
    return (output)
}



# ############################################################################
# The mat2texbmat function writes matrices to LaTeX bordermatrix
# ############################################################################
# Input
#   mat - matrix
# Output
#   character string with LaTeX output
# ############################################################################
mat2texbmat <- function(mat)
{
    nr <- dim(mat)[1]
    tex <- "$$\\bordermatrix{\n"
    tex<- paste0(tex, paste0(c("~", colnames(mat)), collapse = " & "), " \\cr\n")
    rn <- rownames(mat)
    for (i in 1:nr) {
        tex <- paste0(tex, paste0(c(rn[i], mat[i, ]), collapse = " & "), " \\cr\n")
    }
    tex <- paste0(tex, "}$$\n\n")
    return (tex)
}




# ############################################################################
# The write_tex function writes LaTeX output and prints short information
# ############################################################################
# Input
#   texf - name of LaTeX file
#   tex - character string with LaTeX code
#   skip2l - logical, if TRUE 2 blank lines are inserted before information
#                     about file to LaTeX code is written
# Output
#   None
# ############################################################################
write_tex <- function(texf, tex, skip2l = FALSE)
{
    if (file.exists(texf)) {
        write(tex, texf, append = TRUE)
    } else {
        write(tex, texf)
    }
    if (skip2l) {
        cat(paste0("\n\nresults written to: \'", texf, "\'\n"))
    } else {
        cat(paste0("results written to: \'", texf, "\'\n"))
    }
}



# ############################################################################
# The eps2tex function creates LaTeX code including selected files
# ############################################################################
# Input
#   files - character vector with filenames
#   descriptions - character vector with file descriptions
#   path - character string with path to be added before filenames
#   vspace - character string determining vspace in LaTeX code
# Output
#   character string with LaTeX output
# ############################################################################
eps2tex <- function(files, descriptions, path = "plots/", vspace = "-3em")
{
    N <- length(files)
    if (N == 1) {
        i2 <- numeric()
        ir <- 1
    } else {
        i2 <- 2 * (1:(N %/% 2))
        ir <- (N %% 2) * (1 + 2 * (N %/% 2))
    }

    latex <- ""

    for (i in i2) {
        mpage <- paste0("\\begin{figure}[h]\n",
                        "\\begin{minipage}{0.5\\textwidth}\n",
                        "\\vspace*{", vspace, "}\n",
                        "\\centering\n",
                        "\\includegraphics[width=0.99\\textwidth, scale=0.55]{",
                        path, files[i - 1], "}\n",
                        "\\caption{", descriptions[i - 1] ,"}\n",
                        "\\end{minipage}\n",
                        "\\begin{minipage}{0.5\\textwidth}\n",
                        "\\vspace*{", vspace, "}\n",
                        "\\centering\n",
                        "\\includegraphics[width=0.99\\textwidth, scale=0.55]{",
                        path, files[i], "}\n",
                        "\\caption{", descriptions[i] ,"}\n",
                        "\\end{minipage}\n",
                        "\\end{figure}\n\n")
        latex <- paste0(latex, mpage)
        if (!(i %% 4)) latex <- paste0(latex, "\\pagebreak\n\n")
    }

    if (ir) {
        mpage <- paste0("\\begin{figure}[h]\n",
                        "\\centering\n",
                        "\\begin{minipage}{0.5\\textwidth}\n",
                        "\\vspace*{", vspace, "}\n",
                        "\\centering\n",
                        "\\includegraphics[width=0.99\\textwidth, scale=0.55]{",
                        path, files[ir], "}\n",
                        "\\caption{", descriptions[ir] ,"}\n",
                        "\\end{minipage}\n",
                        "\\end{figure}")
        latex <- paste0(latex, mpage)
    }

    return (latex)
}



# ############################################################################
# The unique_output_name function finds files in the specified path
# created with pattern name_[0-9]+\.ext and creates following string:
# name_ and maximum number incremented by one
# ############################################################################
# Input
#   path - path to be searched
#   name - name to be searched
#   ext - extension to be searched
# Output
#   character:
#   A type and the maximum number incremented by one
# ############################################################################
unique_output_name <- function(path, name = "plot", ext = "eps")
{
    lof <- list.files(path, pattern = paste0(name, "_[0-9]+\\.", ext))
    if (length(lof)) {
        pat <- paste0("(?<=", name, "_).*?(?=\\.", ext, ")")
        re <- regexpr(pattern = pat, text = lof, perl = TRUE)
        pl_numbers <- as.numeric(regmatches(lof, re))
        next_number <- max(pl_numbers) + 1
    } else next_number <- 1
    nname <- paste0(name, "_", next_number, ".", ext)

    return (nname)
}
