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
# Functions preprocess_model make_model and load_model
# ############################################################################

# ############################################################################
# preprocess_model function parses gEcon preprocessor extensions (describing
# model variants), selects a variant (if provided) and writes the result (final
# model formulation) to a new .gcn file.
# ############################################################################
# Input
#   filename - a name of model file (.gcn file)
#   variant - (optional) an integer value (0-99), selected variant
#   variant_name - (optional) an character value, variant name appended
#                  to the name of the resulting .gcn file
#   filename - a name of model file (.gcn file)
# Output
#   A character value, name of the resulting .gcn file.
# ############################################################################
preprocess_model <- function(filename, variant = NULL, variant_name = NULL)
{
    lines <- readLines(filename)
    varpat <- "( *[1-9]?[0-9] *| *[1-9]?[0-9] *- *[1-9]?[0-9] *)"
    varpat <- paste0("(", varpat, "(,", varpat, ")*)")
    stpat <- paste0("##", varpat, " *\\{\\{ *(//.*)?")
    endpat <- paste0("## *\\}\\} *")
    st <- which(grepl(pattern = paste0("^(", endpat, ")?", stpat, "$"), lines))
    en <- which(grepl(pattern = paste0("^", endpat, "(", stpat, ")?$"), lines))

    # parse errors
    geconerrstr <- "\n(gEcon preprocessor error): "
    preperr <- which(grepl(pattern = "^##.*((\\{\\{)|(\\}\\})).*$", lines))
    preperr <- setdiff(preperr, union(st, en))
    if (length(preperr) > 0) {
        stop(paste0(geconerrstr, "unexpected token(s); error(s) at line(s) ",
                    paste(preperr, collapse = ", ")), call. = FALSE)
    }

    # checks
    sten <- findInterval(st, c(0, en), left.open = FALSE)
    if (any(diff(sten) == 0)) {
        ii <- which.max(diff(sten) == 0)
        stop(paste0(geconerrstr, "variant opened before previous was closed; error at line ",
             st[ii + 1]), call. = FALSE)
    }
    nv <- length(en)
    if (any(sten > nv)) {
        ii <- which.max(sten > nv)
        stop(paste0(geconerrstr, "variant not closed; error at line ", st[ii]),
             call. = FALSE)
    }
    enst <- findInterval(en, c(st, length(lines)), left.open = TRUE)
    if (any(enst == 0)) {
        ii <- which.max(enst == 0)
        stop(paste0(geconerrstr, "no variant to close; error at line ", en[ii]),
             call. = FALSE)
    }
    if (any(diff(enst) == 0)) {
        ii <- which.max(diff(enst) == 0)
        stop(paste0(geconerrstr, "no variant to close; error at line ",
             en[ii + 1]), call. = FALSE)
    }

    # parse variant numbers
    vn <- sub("//.*$", "", lines[st])
    vn <- gsub("([\\{\\}#]*)", "", vn)
    vn <- gsub("\\-)", ":", vn)
    vn <- paste0("c(", vn, ")")
    vn <- lapply(vn, function(x) eval(parse(text = x)))
    vall <- unique(unlist(vn))

    # check variant value and select variant
    if (is.null(variant)) {
        if (!is.null(variant_name)) {
            stop("'variant_name' provided without 'variant'", call. = FALSE)
        }
        if (nv > 0) {
            warning("found variant marks, but no variant selected", call. = FALSE)
        }
        return (filename)
    }
    if (!is.numeric(variant) || length(variant) != 1 || variant < 0 || variant > 99) {
        stop("expected 'variant' to be a single number (0-99)", call. = FALSE)
    }
    if (!(variant %in% vall)) {
        stop(paste0("variant ", variant, " not found in the code"),
             call. = FALSE)
    }
    if (is.null(variant_name)) {
        variant_name <- as.character(variant)
    }
    if (!is.character(variant_name) || length(variant_name) != 1) {
        stop("expected 'variant_name' to be a single character value",
             call. = FALSE)
    }

    # comment unused variants, remove variant marks
    for (i in 1:nv) {
        if (en[i] - st[i] > 1) {
            ii <- (st[i] + 1L):(en[i] - 1L)
            lines[ii] <- gsub("\t", "        ", lines[ii])
            ws <- nchar(sub("^([ ]*)(([^ ])(.*))?$", "\\1", lines[ii]))
            if (sum(ws > 0) > 0) {
                indent <- min(ws[ws > 0])
            } else {
                indent <- 8
            }
            indent <- paste0(rep(" ", indent), collapse = "")
            if (!(variant %in% vn[i])) {
                lines[ii[ws > 0]] <- sub(paste0("^(", indent, ")( *[^ ].*)$"), "\\1# \\2\\",
                                 lines[ii[ws > 0]])
            }
        } else {
            indent <- paste0(rep(" ", 8), collapse = "")
        }
        if (!(variant %in% vn[i])) {
            lines[st[i]] <- sub("^([^/]*)(// ?(.*))?$", paste0(indent, "# # \\3"),
                                lines[st[i]])
        } else {
            lines[st[i]] <- sub("^([^/]*)(// ?(.*))?$", paste0(indent, "# \\3"),
                                lines[st[i]])
        }
    }
    eno <- setdiff(en, st)
    lines[eno] <- ""

    # write file and return file name
    filename <- add_gcn_ext(paste0(rm_gcn_ext(filename), "_", variant_name))
    writeLines(lines, filename)

    return (filename)
}

# ############################################################################
# The make_model function calls dynamic library, parses model file,
# generates R output and loads it into gecon_model class object
# ############################################################################
# Input
#   filename - a name of model file (.gcn file)
#   variant - (optional) an integer value (0-99), selected variant
#   variant_name - (optional) an character value, variant name appended
#                  to the name of the resulting .gcn file
# Output
#   object of gecon_model class
# ############################################################################
make_model <- function(filename, variant = NULL, variant_name = NULL)
{
    if (!is.character(filename) | length(filename) != 1 | filename == "") {
        stop("filename has to be passed as a non-empty character string")
    }
    filename <- add_gcn_ext(filename)
    res <- tryCatch(expr = { filename <- normalizePath(filename,
                                                       winslash = "/",
                                                       mustWork = TRUE)
                    TRUE
        },
        error = function(e) FALSE
    )
    if (!res) {
        stop("the specified model file does not exist")
    }
    st <- proc.time()
    filename <- preprocess_model(filename, variant, variant_name)
    ok <- .Call("parse_from_R", filename)
    en <- proc.time()
    cat(paste0("model parsed in ", round(en[[3]] - st[[3]], digits = 2), "s\n"))
    filename <- rm_gcn_ext(filename)
    return (load_model(filename))
}

# ############################################################################
# The load_model function loads the already generated .model.R file and creates
# an object of the \code{gecon_model} class
# ############################################################################
# Input
#   filename - the path to the .R file containing the model's functions and
#              variables (a .model.R file generated earlier).
# Output
#   object of gecon_model class
# ############################################################################
load_model <- function(filename)
{
    if (!is.character(filename) | length(filename) != 1 | filename == "") {
        stop("filename has to be passed as a non-empty character string")
    }
    filename <- add_r_ext(filename)
    st <- proc.time()
    res <- eval(parse(filename))
    en <- proc.time()
    cat(paste0("model loaded in ", round(en[[3]] - st[[3]], digits = 2), "s\n"))
    return (res)
}

