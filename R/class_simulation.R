# ###################################################################
# (c) Kancelaria Prezesa Rady Ministrów 2012-2015                   #
# Treść licencji w pliku 'LICENCE'                                  #
#                                                                   #
# (c) Chancellery of the Prime Minister 2012-2015                   #
# Licence terms can be found in the file 'LICENCE'                  #
#                                                                   #
# Authors: Karol Podemski, Kaja Retkiewicz-Wijtiwiak                #
# ###################################################################
# Class for storing simulation results
# ###################################################################

# ###################################################################
# Class definition
# ###################################################################
setClass(
    Class = "gecon_simulation",
    representation = representation(
        sim = 'array',
        shock_list = 'vector', 
        var_list = 'vector', 
        sim_type = 'character', 
        time_n = 'numeric',
        model_info = 'character',
        model_variable_name = 'character'
    ) 
)

# ###################################################################
# The function gecon_simulation is 
# a gecon_simulation class object constructor
# ###################################################################
# Input 
#   sim = (array) impulse response functions with dimensions:
#                   (variables, time, shocks)
#   shock_list = a vector of shock names
#   var_list = a vector of simulated variables names
#   sim_type = a character, type of simulation
#   time_n = a numeric with number of periods of simulation
#   model_info = (character vector, length = 3) information about
#           model: input file name, file path and the date of creation
#   model_variable_name = (character) the name of variable storing
#                   gecon_model object based on which
#                   the simulations were created
# Output
#   An object of the gecon_simulation class
# ###################################################################
gecon_simulation <- function(sim, shock_list, var_list, sim_type,
                             time_n, model_info, model_variable_name) 
{
    sim_obj <- new('gecon_simulation')
    
    if (!is.array(sim)) {
        stop('simulation results should be of array class')
    } else sim_obj@sim <- sim
    
    if (!is.vector(shock_list)) {
        stop('shock_list should be of vector class')
    } else sim_obj@shock_list <- shock_list
    
    if (!is.vector(var_list)) {
        stop('var_list should be of vector class')
    } else sim_obj@var_list <- var_list

    if (!is.vector(sim_type)) {
        stop('sim_type should be of character class')
    } else sim_obj@sim_type <- sim_type
    
    if (!is.numeric(time_n)) {
        stop('time_n should be of vector class')
    } else sim_obj@time_n <- time_n    

    if (!is.character(model_info)) {
        stop('model_info should be of character class')
    } else sim_obj@model_info <- model_info   

    if (!is.character(model_variable_name)) {
        stop('model_variable_name should be of character class')
    } else sim_obj@model_variable_name <- model_variable_name   
    
    return(sim_obj)
}

# ###################################################################
# The print method prints information about simulations
# ###################################################################
# Input 
#   x - an object of the gecon_simulation class 
# Output
#   Information about the simulation
# ###################################################################
setMethod(
    "print",
    signature(x = "gecon_simulation"), 
    function(x) 
    {
               
        cat('Simulation has been performed on model stored in variable:', 
                x@model_variable_name, '\n')
        
        cat('The model has been created based on:', 
            x@model_info[1], ' file \n')
        
        cat('Type of simulation:', 
                x@sim_type, '\n')
                
        cat("----------------------------------------------------------", "\n\n")        
        cat('Following shocks have been used to simulate the model: \n')
        cat("\n")
                shock_list <- as.data.frame(x@shock_list)
                names(shock_list) <- 'Shocks:'
                print(shock_list)
        
        cat("----------------------------------------------------------", "\n\n")       
        cat('Simulation has been performed for the following variables: \n')
        cat("\n")
                var_list <- as.data.frame(x@var_list)
                names(var_list) <- 'Variables:'
                print(var_list)
        cat("----------------------------------------------------------", "\n\n")
    }
)

# ###################################################################
# The plot_simulation function plots simulation results to the device 
# or saves them in the model's /plots subdirectory as .eps files 
# ###################################################################
# Input 
#   sim_obj - an object of class gecon_simulation
#   to_tex - a logical. if TRUE, the function generates .eps files
#            in model's /plots subdirectory and includes them
#            in tex report
#   to_eps - if TRUE, plots shall be 
#            saved as .eps file in model's /plots subdirectory
# Output
#   Plots of objects of the gecon_simulation class
# ###################################################################
plot_simulation <- function(sim_obj, 
                            to_tex = NULL,
                            to_eps = NULL) 
{
    if (!is(sim_obj, "gecon_simulation"))
        return("This function plots gecon_simulation objects only.")

    
    file_name_next <- NULL
    file_name <- NULL
    plot_title <- NULL
    
    if(is.null(to_tex))
        to_tex <- FALSE

    if (to_tex) {
        if (is.null(to_eps)) {
            to_eps <- TRUE
        } else if (!to_eps) {
            warning("to_tex option set to TRUE overrides option to_eps = FALSE.")
            to_eps <- TRUE
        }
    }
    
    if(is.null(to_eps))
        to_eps <- FALSE
                
    if (to_eps) {
        path <- gsub(pattern = paste(sim_obj@model_info[1], '$', sep=''), 
                     replacement = '', 
                     x = sim_obj@model_info[2])
        path <- paste(path, 'plots', sep='')
        dir.create(path, showWarnings = FALSE)
        file_name <- paste(path, '/',unique_output_name(path),  sep='')
        
        if (sim_obj@sim_type == "Random path simulation") 
            file_name <- paste(file_name, "_ranpath", sep="")
        else if (sim_obj@sim_type == "Simulation with user-defined shocks")
            file_name <- paste(file_name, "_sim", sep="")                
    }
    
    if (to_tex) {
        eps_tex_list <- vector(length = dim(sim_obj@sim)[3], "list") 
    }
    
    # each shock
    for (j in 1:dim(sim_obj@sim)[3]) {
        d <- as.ts(t(matrix(sim_obj@sim[, , j], nrow = length(sim_obj@var_list),
                            ncol = sim_obj@time_n)))
        
        # number of irfs on each plot
        n_pl <- 5
        # required number of plots
        numb_pl <- ((dim(d)[2] - 1) %/% n_pl) + 1
        if (to_tex) {
            shock_tex_list <- vector(length = numb_pl)
        }
        
        # loop for each of plot
        for (i in 1:numb_pl) {
            if (i < numb_pl) {
                indices <- c(1:(n_pl) + (i - 1) * n_pl)
            } else {
                indices <- c(1 + (i - 1) * n_pl):dim(d)[2]                     
            }
            
            irf_plot <- as.ts(matrix(d[, indices], ncol = length(indices),
                                            nrow = sim_obj@time_n))
            
            # computing plots            
            if (sim_obj@sim_type == "Impulse response functions") {
                if (numb_pl == 1) {
                    if (is.null(file_name)) {
                        plot_title = paste('Impulse response function for: ', 
                                        sim_obj@shock_list[j],
                                        ' shock', sep='')
                    } else {
                        file_name_next <- paste(as.character(file_name), "_irf_",
                                                sim_obj@shock_list[j], ".eps", sep = "")
                    }
                } else {
                    if (is.null(file_name)) {
                        plot_title = paste('Impulse response function for: ', 
                                        sim_obj@shock_list[j],
                                        ' shock (plot no. ', i, ')', sep='')
                    } else {
                        file_name_next <- paste(as.character(file_name), "_irf_",
                                                sim_obj@shock_list[j], i, ".eps", sep = "")
                    }                                        
                }                
            } else {
                if (is.null(file_name))
                    plot_title = sim_obj@sim_type
                    
                if (numb_pl == 1) {
                    if (!is.null(file_name)) {
                        file_name_next <- paste(as.character(file_name), 
                                                ".eps", sep = "")
                    }
                } else {
                    if (!is.null(file_name)) {
                        file_name_next <- paste(as.character(file_name), 
                                                i, ".eps", sep = "") 
                    }
                }
            }
            plot_gecon(irf_plot, 
                       serieslab = sim_obj@var_list[indices],
                       main = plot_title,
                       f_name = file_name_next)
            
            #tex list
            if (to_tex) {
                shock_tex_list[i] <- file_name_next
            }
        } # end plot loop
        if (to_tex) {
            eps_tex_list[[j]] <- shock_tex_list
        }
    } # end shock loop
    
    if (to_tex) {
        fig_numb <- (numb_pl %% 2)
        tex_out <- paste("\n\\pagebreak\n\\section{", sim_obj@sim_type,"}\n", sep = '')
     
        for (i in 1 : length(eps_tex_list)) {

            if (sim_obj@sim_type == "Random path simulation") 
                 descr = paste('Random path simulation for', sim_obj@time_n, 'periods')
            else if (sim_obj@sim_type == "Simulation with user-defined shocks")
                 descr = paste('User-defined simulation for', sim_obj@time_n, 'periods')
            else { descr = paste('Impulse response function for ', 
                               string_to_tex(sim_obj@shock_list[i]), ' shock', sep='')
                   tex_out <- paste(tex_out, paste("\n\\subsection{Shock ", 
                                    string_to_tex(sim_obj@shock_list[i]),
                                    "}\n", sep = ''), sep = '\n \n')
            }
            single_minipage <- ''
            double_minipage <- ''
            if (fig_numb) {
                single_minipage <- paste("\\begin{figure}[h]\n",
                            "\\begin{minipage}{0.5\\textwidth}\n",
                            "\\vspace*{-3em}\n",
                            "\\centering\n",
                            "\\includegraphics[width=0.99\\textwidth, scale=0.55]{", 
                            eps_tex_list[[i]][numb_pl],
                            "}\n",
                            "\\caption{", descr ,"}\n",
                            "\\end{minipage}\n",
                            "\\end{figure}", sep = '')
            }
            
            if (numb_pl > 1) {
                for (j in seq(from = 1, to = (numb_pl - fig_numb), by =2)) {
                    double_minipage_el <- paste(
                                    "\\begin{figure}[h]\n",
                                    "\\begin{minipage}{0.5\\textwidth}\n",
                                    "\\vspace*{-3em}\n",
                                    "\\centering\n",
                                    "\\includegraphics[width=0.99\\textwidth, scale=0.55]{",
                                    eps_tex_list[[i]][j],
                                    "}\n",
                                    "\\caption{", descr ,"}\n",
                                    "\\end{minipage}\n",
                                    "\\begin{minipage}{0.5\\textwidth}\n",
                                    "\\vspace*{-3em}\n",
                                    "\\centering\n",
                                    "\\includegraphics[width=0.99\\textwidth, scale=0.55]{", 
                                    eps_tex_list[[i]][j + 1],
                                    "}\n",
                                    "\\caption{", descr ,"}\n",
                                    "\\end{minipage}\n",
                                    "\\end{figure}", sep = '')
                    double_minipage <- paste(double_minipage, double_minipage_el, sep = '\n \n')
                    if ( (j + 1) %% 4 == 0  & (j != (numb_pl - fig_numb)) )
                        double_minipage <- paste(double_minipage, "\n\\pagebreak")
                }
            }
            double_minipage <- paste(double_minipage, single_minipage, sep = '\n \n')
            tex_out <- paste(tex_out, '\n', double_minipage, '\\pagebreak', '\n')
        }

        filenam <- paste(sim_obj@model_info[2],
                        '.results.tex', sep='') 
    
        if (file.exists(filenam)) {
            write(tex_out, filenam, sep = '\n', append = TRUE)  
            cat("\nPlots added to the following TeX file: \n")
            cat(filenam, "\n")
        } else {
            write(tex_out, filenam  , sep = '\n') 
            cat("\nPlots added to the following TeX file: \n")
            cat(filenam, "\n")        
        }
    }
}

# ###################################################################
# The summary method displays the results of simulations
# ###################################################################
# Input 
#   x - an object of class gecon_simulation
# Output
#   Prints the simulation results
# ###################################################################
setMethod(
    "summary",
    signature(object = "gecon_simulation"), 
    function(object, ...) 
    {
        cat("=================== SIMULATION SUMMARY ===================", "\n\n")

        cat('Simulation has been performed on model stored in variable:', 
                object@model_variable_name, '\n')
        cat("\n")                   
        cat('The model has been created based on:', 
            object@model_info[1], ' file \n')
        cat("\n")        
        cat('Time span for simulation is:', 
                object@time_n, 'periods \n')        

        cat("\n")               
        cat('Type of simulation:', 
                object@sim_type, '\n')
        cat("----------------------------------------------------------", "\n\n")        
        cat('Following shocks have been used to simulate the model: \n')
        cat("\n")
                shock_list <- as.data.frame(object@shock_list)
                names(shock_list) <- 'Shocks:'
                print(shock_list)
        
        cat("----------------------------------------------------------", "\n\n")       
        cat('Simulation has been performed for the following variables: \n')
        cat("\n")
                var_list <- as.data.frame(object@var_list)
                names(var_list) <- 'Variables:'
                print(var_list)
                
        cat("\n----------------------------------------------------------", "\n\n")        
        cat("Simulation results \n")
        cat("\n")
        if (object@sim_type == "Impulse response functions") {
            output <- vector("list", (length(object@shock_list)))
            for (i in 1:length(object@shock_list)) {
                cat("\n", object@shock_list[i], " \n")            
                d <- matrix(object@sim[, , i], length(object@var_list),
                            object@time_n)
                rownames(d) <- object@var_list
                colnames(d) <- c(1:object@time_n)
                print(d)
            }
        } else {
            # cat("\n", object@sim_type, "\n")
            d <- matrix(object@sim[, , 1], ncol = object@time_n)
            rownames(d) <- object@var_list
            colnames(d) <- c(1:object@time_n)             
            output <- vector("list", (1)) 
            print(d)             
        }
    }
)

# ###################################################################
# The get_simulation_results function returns the simulation results as 
# the list
# ###################################################################
# Input 
#   sim_obj - an object of class gecon_simulation
# Output
#   A list with the results of simulation
# ###################################################################
get_simulation_results <- function(sim_obj)
{
    if (!is(sim_obj, "gecon_simulation"))
        stop("This function accepts only arguments ",
             "of class gecon_simulation.")
        
    if (sim_obj@sim_type == "Impulse response functions") {
        output <- vector("list", (length(sim_obj@shock_list)))
        for (i in 1:length(sim_obj@shock_list)) {            
            d <- matrix(sim_obj@sim[, , i], length(sim_obj@var_list),
                        sim_obj@time_n)
            rownames(d) <- sim_obj@var_list
            colnames(d) <- c(1:sim_obj@time_n)
            output[[i]] <- d
            names(output) <- c(sim_obj@shock_list)
        }
    } else {
        d <- matrix(sim_obj@sim[, , 1], ncol = sim_obj@time_n)
        rownames(d) <- sim_obj@var_list
        colnames(d) <- c(1:sim_obj@time_n)             
        output <- vector("list", (1)) 
        output[[1]] <- d             
    }

    return(output)
}


# ###################################################################
# The show method controls how the gecon_simulation object is printed
# ###################################################################
# Input 
#   object - object of class gecon_simulation
# Output
#   Information about object of class gecon_simulation
# ################################################################### 
setMethod(
    "show", signature(object = "gecon_simulation"),
    function(object) {
        cat('\n Simulation has been performed based on model:', 
            object@model_variable_name, '\n')
        cat('\n Time span for simulation are:', 
            object@time_n, 'periods \n')

        cat('\n Simulation has been performed for', length(object@shock_list), 
            'shocks and', length(object@var_list), 'variables \n')
    }
)

# ###################################################################
# The plot_gecon function plots the matrix into console or to eps file 
# ###################################################################
# Input 
#   sim_obj - a matrix with values of gecon_simulation
#   serieslab - a vector of series names
#   main - the plot label
#   f_name - the name of eps file to be created
# Output
#   A plot in the screen or printed to eps file
# ################################################################### 
plot_gecon <- function(sim_obj,
                       serieslab = NULL,
                       main = NULL,
                       f_name = NULL)
{
    
    # labels for OX
    lab <- time(sim_obj)
    color_palette = c("red3", "black", "gray48", "darkorange2", 
                      "orangered4", "coral3", "gray73")
    
    # y limits
    yl <- ylimits(min(sim_obj), max(sim_obj), forcezero = TRUE)
    ylim = yl[1:2]
    yby = yl[3]
    xl <- xlimits(min(lab), max(lab))
    
    # for export
    if(!is.null(f_name)) {
        options(scipen = 1)
        ps.options(family = "serif", encoding = "cp1250")
        setEPS() 
        postscript(f_name, family = "serif")
    }
    
    # specifying characteristics of the plot
    layout(rbind(1, 2), heights = c(7, 1))
    
    # default margin parameters
    def_mar <- par()$mar
    
    # default axis parameters
    def_mgp <- par()$mgp
    par(mar = c(3.1, 5.1, 5.1, 3.1))
    par(mgp = c(2, 1, 0))
    
    if (dim(sim_obj)[1] > 1) {
        plot(sim_obj, plot.type = "single", yaxt = "n", xaxt = "n", 
             col = color_palette, 
             lwd = 2, 
             main = main, 
             xlab = "Periods", 
             ylab = "Deviation from the steady state",
             ylim = ylim)     
    } else  {
        plot(x = rep(1, length(sim_obj)),
             y = as.vector(sim_obj),
             col = color_palette,  
             lwd = 2, 
             main = main, 
             xlab = "Periods", 
             ylab = "Deviation from the steady state",
             ylim = ylim, xaxt = "n", yaxt = "n")     
    }
    # axis formatting
    axis(1, at = (xl), labels = xl, tcl = -0.3, cex.axis = 0.7)
    axis(2, at = seq(ylim[1], ylim[2], by = yby), 
             labels = pretty_labels(yl[1], yl[2], yl[3]), 
             cex.axis = 0.7)
    
    if (is.null(nrow(sim_obj))) {
        x1var = length(sim_obj) + 1
    } else {
        x1var = nrow(sim_obj) + 1
    }
    segments(x0 = 0, y0 = 0, x1 = x1var, col = "black", lwd = 1, lty = 2)

    # labels
    if (!is.null(serieslab) ) { 
        par(mar = c(0, 0, 0, 0))
        plot.new()
        legend("center", legend = serieslab, lty = 1, col = color_palette, 
                cex = 0.7, border = F, horiz = TRUE)
    }

    # closing export
    if(!is.null(f_name)) {
        dev.off()
        cat("plot written to: \n")
        cat(f_name, "\n")
    }
    
    # resetting default parameters
    par(mar = def_mar)
    par(mgp = def_mgp)    
}

# ###################################################################
# Function roundpow10 rounds a number to p power of 10
# ###################################################################
# Input
#   y - number to be rounded
#   p - specified power
#   up - parameter denoting whether rounding function
#         should be ceiling (TRUE) or floor (false)
# Output
#   rounded number
# ###################################################################
roundpow10 <- function(y, p, up = TRUE)
{
    if (y == 0) return(y)
    if (y < 0) {
        y = -y;
        sgn = -1;
    } else {
        sgn = 1;
    }
    if (((up) & (sgn > 0)) | ((!up) & (sgn < 0))) {
        return(sgn * ceiling(y / 10^p) * 10^p)
    } 
    return(sgn * floor(y / 10^p) * 10^p)
}


# ###################################################################
# Function ylimits finds optimal limits on plot according
# to input data
# ###################################################################
# Input
#   ymin - the minimum value of the data
#   ymax - the maximum value of the data
#   forcezero - Boolean, should 0 be always included in the range? 
#               (default FALSE)
# Output
#   Vector with three elements denoting minimum, maximum and step
#   for y axis
# ###################################################################
ylimits <- function (ymin, ymax, forcezero = FALSE)
{
    if ((ymin == 0) & (ymax == 0)) { return(c(-1, 1, 1)) }
    if (ymin > ymax) { stop("invalid range (min greater than max)") }
    if (forcezero) {
        if (ymin > 0) ymin = 0;
        if (ymax < 0) ymax = 0;
    }

    if (ymax > ymin) {
        p = floor(log10(ymax - ymin))
    } else {
        p = floor(log10(ymax));
    }

    ylo = roundpow10(ymin, p, F)
    yhi = roundpow10(ymax, p, T)

    if (max(yhi, abs(ylo)) == 10^p) p = p - 1
    yby = 10^p

    if ((ymax - ymin) / 10^p > 15) { yby = yby * 5 }
    else if ((ymax - ymin) / 10^p > 8) { yby = yby * 2 }
    else if ((ymax - ymin) / 10^p > 4) {  }
    else if ((ymax - ymin) / 10^p > 2) { yby = yby / 2 }
    else { yby = yby / 4 }

    rh = ceiling(yhi / yby) - ceiling(ymax / yby);
    rl = ceiling(abs(ymin / yby)) - ceiling(abs(ylo / yby)) ;
    if (ylo < 0) rl = -rl;
    yhi = yhi - rh * yby;
    ylo = ylo + rl * yby;
    if (yhi < ymax) yhi = yhi + yby;
    if (ylo > ymin) ylo = ylo - yby;

    if ((yhi > 0) & (ylo < 0)) {
        rh = abs(yhi / yby - round(yhi / yby))
        if (rh > 2e-16) yhi = yhi + rh * yby;
        rl = abs(abs(ylo) / yby - round(abs(ylo) / yby))
        if (rl > 2e-16) ylo = ylo - rl * yby;
    }

    return(c(ylo, yhi, yby))
}

# ###################################################################
# The xlimits function finds optimal limits on x-axis according 
# to input data
# ###################################################################
# Input 
#   ymin - the minimum value of the data
#   ymax - the maximum value of the data
# Output
#   Scale for x- axis
# ###################################################################
xlimits <- function (xmin, xmax) 
{
    if ((xmin == 1) & (xmax == 1)) 
        { return(c(1, 1)) }
    
    if ((xmin == 1) & (xmax < 10)) 
        { return(c(xmin : xmax))}
    
    if ((xmin == 1) & (xmax >= 10)) 
        {return(c(1, seq(5, 5 * floor(xmax /5), 5)))}
}

# ###################################################################
# The roundpow10 function rounds a number to the p power of 10
# ###################################################################
# Input 
#   y - a number to be rounded
#   p - a number specifying the power
#   up - parameter denoting whether rounding function
#         should be ceiling (TRUE) or floor (false) 
# Output
#   A rounded number
# ###################################################################
roundpow10 <- function(y, p, up = TRUE) 
{
    if (y == 0) return(y)
    if (up) return(ceiling(y / 10^p) * 10^p)
    else return(floor(y / 10^p) * 10^p)
}


# ###################################################################
# The unique_plot_name function finds files in the specified path 
# created with pattern type[0-9]+ and creates following string: 
# type and maximum number incremented by one
# ###################################################################
# Input 
#   path - a path to be searched
#   type - a name to be searched
# Output
#   character: 
#   A type and the maximum number incremented by one
# ###################################################################
unique_output_name <- function(path, type = "plot") 
{
    lof <- list.files(path, 
                      pattern = paste(type, "[0-9]+_.*", sep = ''))
    if (length(lof)) {
        re <- regexpr(pattern = "(?<=plot).*?(?=_)", text = lof, perl = TRUE)
        pl_numbers <- as.numeric(regmatches(lof, re))
        next_number <- max(pl_numbers) + 1
    } else next_number <- 1
    nname <- paste(type, next_number, sep ='')
    
    return(nname)
}

# ###################################################################
# Function returns prepared labels for axis from sequence parameters
#   (from result of ylimits function) 
# ###################################################################
# Input
#   yl - initial value
#   yh - limit value 
#   yb - step value 
# Output
#   character vector representing axis labels
# ###################################################################
pretty_labels <- function(yl, yh, yb)
{
    return (format(round(seq(yl, yh, yb), find_signif(yb)), scientific = FALSE))
}

# ###################################################################
# Function finds number of significant digits from specified value
# ###################################################################
# Input
#   x - numeric value
#   tol - (optional) definied tolerance 
# Output
#   not negative integer representing number of significant digits
# ###################################################################
find_signif <- function(x, tol = 1e-2)
{
    if (x == 0) { return (0) }
    p <- 0
    while (abs(1 - round(x, digits = p)/x) > tol) {
        p <- p + 1
    }
    return (p)
}
