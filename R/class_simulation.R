# ############################################################################
# This file is a part of gEcon.                                              #
#                                                                            #
# (c) Chancellery of the Prime Minister of the Republic of Poland 2012-2015  #
# (c) Grzegorz Klima, Karol Podemski, Kaja Retkiewicz-Wijtiwiak 2015-2018    #
# License terms can be found in the file 'LICENCE'                           #
#                                                                            #
# Authors: Grzegorz Klima, Karol Podemski, Kaja Retkiewicz-Wijtiwiak         #
# ############################################################################
# Class for storing simulation results
# ############################################################################

# ############################################################################
# Class definition
# ############################################################################
setClass(
    Class = "gecon_simulation",
    representation = representation(
        sim = "array",
        shocks = "vector",
        shocks_tex = "vector",
        variables = "vector",
        variables_tex = "vector",
        sim_name = "character",
        model_info = "character",
        r_object_name = "character"
    )
)

# ############################################################################
# The function gecon_simulation is a constructor
# of an object of gecon_simulation class
# ############################################################################
# Input
#   sim - (array) simulation results or impulse response functions with 2 or 3
#       dimensions: (variables, time) or (variables, time, shocks)
#   shocks - vector of shock names
#   shocks_tex - vector of LaTeX shock names
#   variables - vector of variables' names
#   variables_tex - vector of LaTeX variables' names
#   sim_name - a character value, the simulation name
#   model_info - a character vector of length 3 containing information about
#                   the model: input file name, input file path and the date of creation
#   r_object_name - a character value, the name of an R object storing
#                   gecon_model object for which simulations were performed
# Output
#   an object of gecon_simulation class
# ############################################################################
gecon_simulation <- function(sim,
                             shocks, shocks_tex, variables, variables_tex,
                             sim_name, model_info, r_object_name)
{
    sim_obj <- new("gecon_simulation")

    if (!is.array(sim) || !(length(dim(sim)) %in% 2:3)) {
        stop("simulation results should be a 2- or 3-dimensional array")
    } else sim_obj@sim <- sim
    if (!is.vector(shocks)) {
        stop("shocks should be of vector type")
    } else sim_obj@shocks <- shocks
    if (!is.vector(shocks_tex)) {
        stop("shocks_tex should be of vector type")
    } else sim_obj@shocks_tex <- shocks_tex
    if (!is.vector(variables)) {
        stop("variables should be of vector type")
    } else sim_obj@variables <- variables
    if (!is.vector(variables_tex)) {
        stop("variables_tex should be of vector type")
    } else sim_obj@variables_tex <- variables_tex
    if (!is.vector(sim_name)) {
        stop("sim_name should be of character type")
    } else sim_obj@sim_name <- sim_name
    if (!is.character(model_info)) {
        stop("model_info should be of character type")
    } else sim_obj@model_info <- model_info
    if (!is.character(r_object_name)) {
        stop("r_object_name should be of character type")
    } else sim_obj@r_object_name <- r_object_name

    return (sim_obj)
}


# ############################################################################
# The get_simulation_results function returns the simulation results as
# a list
# ############################################################################
# Input
#   sim_obj - an object of gecon_simulation class
# Output
#   A matrix or list of matrices with the results of simulation
# ############################################################################
get_simulation_results <- function(sim_obj)
{
    if (!is(sim_obj, "gecon_simulation"))
        stop("sim_obj should be an object of gecon_simulation class")

    if (length(dim(sim_obj@sim)) == 3) {
        nv <- dim(sim_obj@sim)[1]
        nt <- dim(sim_obj@sim)[2]
        ns <- dim(sim_obj@sim)[3]
        res <- list()
        for (s in 1:ns) {
            sim <- as.matrix(sim_obj@sim[, , s], nrow = nv, ncol = nt)
            res <- c(res, list(sim))
        }
        names(res) <- sim_obj@shocks
        return (res)
    }

    return (sim_obj@sim)
}


# ############################################################################
# The show method controls how the gecon_simulation object is printed
# ############################################################################
# Input
#   object - an object of gecon_simulation class
# Output
#   Information about an object of gecon_simulation class
# ############################################################################
setMethod("show", signature(object = "gecon_simulation"),
function(object)
{
    cat(paste0(object@sim_name, "\n"))
    cat(paste0("Simulation horizon: ", dim(object@sim)[2], "\n"))
    nv <- length(object@variables)
    ns <- length(object@shocks)
    if (length(dim(object@sim)) == 3) {
        mes <- "Responses of "
        mes <- paste0(mes, numnoun(nv, "variable"))
        mes <- paste0(mes, " to ")
        mes <- paste0(mes, numnoun(ns, "shock"))
        cat(paste0(mes, "\n"))
    } else {
        mes <- "Evolution of "
        mes <- paste0(mes, numnoun(nv, "variable"))
        mes <- paste0(mes, " based on path of ")
        mes <- paste0(mes, numnoun(ns, "shock"))
        cat(paste0(mes, "\n"))
    }
})

# ############################################################################
# The print method prints information about simulations
# ############################################################################
# Input
#   x - an object of gecon_simulation class
# Output
#   Information about the simulation
# ############################################################################
setMethod("print", signature(x = "gecon_simulation"),
function(x)
{
    object <- x
    cat(paste0(object@sim_name, "\n"))
    cat(paste0("Simulation horizon: ", dim(object@sim)[2], "\n"))
    nv <- length(object@variables)
    ns <- length(object@shocks)
    cat(paste0("Variables (", nv, "):\n"))
    cat(list2str2(object@variables))
    cat("\n")
    cat(paste0("Shocks (", ns, "):\n"))
    cat(list2str2(object@shocks))
    cat("\n")
})


# ############################################################################
# The summary method displays the results of simulations
# ############################################################################
# Input
#   x - an object of gecon_simulation class
# Output
#   Prints the simulation results
# ############################################################################
setMethod("summary", signature(object = "gecon_simulation"),
function(object)
{
    print(object)
    if (length(dim(object@sim)) == 3) {
        ns <- dim(object@sim)[3]
        nt <- dim(object@sim)[2]
        nv <- dim(object@sim)[1]
        for (s in 1:ns) {
            cat(paste0("\n", paste(rep("-", 70), collapse = ""), "\n"))
            cat(paste0("Impulse responses to ", object@shocks[s], " shock\n\n"))
            sim <- as.matrix(object@sim[, , s], nrow = nv, ncol = nt)
            print(sim)
        }
    } else {
        cat(paste0("\n", paste(rep("-", 70), collapse = ""), "\n\n"))
        print(object@sim)
    }
})


# ############################################################################
# The plot_simulation function plots simulation results on the screen
# or saves them in the model's subdirectory /plots as .eps files
# ############################################################################
# Input
#   sim_obj - an object of gecon_simulation class
#   to_eps - a logical value, if TRUE, plots will be
#            saved as .eps file in the model's /plots subdirectory,
#            plots are also added to the .results.tex file
# Output
#   Plots of objects of gecon_simulation class
# ############################################################################
plot_simulation <- function(sim_obj, to_eps = FALSE)
{
    if (!is(sim_obj, "gecon_simulation"))
        stop("sim_obj should be an object of gecon_simulation class")
    if (!is.logical(to_eps))
        stop("invalid to_eps argument")

    if (to_eps) {
        path <- paste0(get_path(sim_obj@model_info[2]), "/plots")
        dir.create(path, showWarnings = FALSE)
        files <- character()
        descriptions <- character()
    }

    # number of variables on each plot
    nvaronplot <- 5

    if (length(dim(sim_obj@sim)) == 3) {
        ns <- dim(sim_obj@sim)[3]
        nv <- dim(sim_obj@sim)[1]
        for (s in 1:ns) {
            d <- as.matrix(sim_obj@sim[, , s])
            if (nv == 1) {
                d <- as.ts(d)
            } else {
                d <- as.ts(t(d))
            }
            nplots <- ((dim(d)[2] - 1) %/% nvaronplot) + 1

            for (i in 1:nplots) {
                if (i < nplots) {
                    indices <- 1:nvaronplot + (i - 1) * nvaronplot
                } else {
                    indices <- (1 + (i - 1) * nvaronplot):dim(d)[2]
                }

                plotts <- as.ts(d[, indices, drop = FALSE])
                impresp <- "Impulse response"
                if (length(indices) > 1) impresp <- paste0(impresp, "s")

                if (to_eps) {
                    file_name <- unique_output_name(path)
                    files <- c(files, file_name)
                    file_name <- paste0(path, "/", file_name)
                    descr <- paste0(sim_obj@variables_tex[indices],
                                    collapse = ", ")
                    descr <- paste0(impresp, " ($", descr, "$) to $",
                                    sim_obj@shocks_tex[s],"$ shock")
                    descriptions <- c(descriptions, descr)
                    plot_title <- NULL
                } else {
                    file_name <- NULL
                    plot_title <- paste0(impresp, " to ", sim_obj@shocks[s],
                                         " shock")
                }

                plot_gecon(plotts,
                           serieslab = sim_obj@variables[indices],
                           main = plot_title,
                           f_name = file_name)
            }
        }
    } else {
        d <- as.ts(t(sim_obj@sim))
        nplots <- ((dim(d)[2] - 1) %/% nvaronplot) + 1

        for (i in 1:nplots) {
            if (i < nplots) {
                indices <- 1:nvaronplot + (i - 1) * nvaronplot
            } else {
                indices <- (1 + (i - 1) * nvaronplot):dim(d)[2]
            }

            plotts <- as.ts(d[, indices, drop = FALSE])

            if (to_eps) {
                file_name <- unique_output_name(path)
                files <- c(files, file_name)
                file_name <- paste0(path, "/", file_name)
                descr <- paste0(sim_obj@variables_tex[indices], collapse = ", ")
                descr <- paste0(sim_obj@sim_name, " ($", descr, "$)")
                descriptions <- c(descriptions, descr)
                plot_title <- NULL
            } else {
                file_name <- NULL
                plot_title <- sim_obj@sim_name
            }

            plot_gecon(plotts,
                       serieslab = sim_obj@variables[indices],
                       main = plot_title,
                       f_name = file_name)
        }
    }

    if (to_eps) {
        tex_out <- paste0("\n\\pagebreak\n\n\\section{", sim_obj@sim_name,"}\n\n")
        tex_out <- paste0(tex_out, eps2tex(files, descriptions))
        filenam <- paste0(rm_gcn_ext(sim_obj@model_info[2]), ".results.tex")
        if (file.exists(filenam)) {
            write(tex_out, filenam, sep = "\n", append = TRUE)
        } else {
            write(tex_out, filenam  , sep = "\n")
        }
        cat(paste0("saved ", list2str(files, "\'"), " in directory \'", path, "\'\n"))
        cat(paste0("plot(s) added to \'", filenam, "\'\n"))
    }
}


