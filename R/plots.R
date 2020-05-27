# ############################################################################
# This file is a part of gEcon.                                              #
#                                                                            #
# (c) Chancellery of the Prime Minister of the Republic of Poland 2012-2015  #
# (c) Grzegorz Klima, Karol Podemski, Kaja Retkiewicz-Wijtiwiak 2015-2018    #
# License terms can be found in the file 'LICENCE'                           #
#                                                                            #
# Authors: Grzegorz Klima, Karol Podemski, Kaja Retkiewicz-Wijtiwiak         #
# ############################################################################
# Plotting simulations
# ############################################################################

# ############################################################################
# The plot_gecon function plots simulation results
# ############################################################################
# Input
#   sim - an object of ts class
#   serieslab - a vector of series names
#   main - the plot label
#   f_name - the name of the .eps file in which plot is to be written
# Output
#   A plot on the screen or written to an .eps file
# ############################################################################
plot_gecon <- function(sim,
                       serieslab = NULL,
                       main = NULL,
                       f_name = NULL)
{
    # labels for OX
    lab <- time(sim)
    color_palette <- c("red3", "black", "gray48", "darkorange2",
                       "orangered4", "coral3", "gray73")

    # y limits
    yl <- ylimits(min(sim), max(sim), forcezero = TRUE)
    ylim <- yl[1:2]
    yby <- yl[3]
    xl <- xlimits(min(lab), max(lab))

    # for export
    if(!is.null(f_name)) {
        options(scipen = 1)
        ps.options(family = "serif", encoding = "ISOLatin1")
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

    if (dim(sim)[1] > 1) {
        plot(sim, plot.type = "single", yaxt = "n", xaxt = "n",
             col = color_palette,
             lwd = 2,
             main = main,
             xlab = "Periods",
             ylab = "Deviation from the steady state",
             ylim = ylim)
    } else  {
        plot(x = rep(1, length(sim)),
             y = as.vector(sim),
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

    if (is.null(nrow(sim))) {
        x1var = length(sim) + 1
    } else {
        x1var = nrow(sim) + 1
    }
    segments(x0 = 0, y0 = 0, x1 = x1var, col = "black", lwd = 1, lty = 2)

    if (!is.null(serieslab) ) {
        par(mar = c(0, 0, 0, 0))
        plot.new()
        legend("center", legend = serieslab, lty = 1, col = color_palette,
                cex = 0.7, border = F, horiz = TRUE)
    }

    if(!is.null(f_name)) dev.off()

    # resetting default parameters
    par(mar = def_mar)
    par(mgp = def_mgp)
}

# ############################################################################
# Function roundpow10 rounds a number to p power of 10
# ############################################################################
# Input
#   y - number to be rounded
#   p - specified power
#   up - parameter denoting whether rounding function
#         should be ceiling (TRUE) or floor (FALSE)
# Output
#   rounded number
# ############################################################################
roundpow10 <- function(y, p, up = TRUE)
{
    if (y == 0) return (y)
    if (y < 0) {
        y = -y;
        sgn = -1;
    } else {
        sgn = 1;
    }
    if (((up) & (sgn > 0)) | ((!up) & (sgn < 0))) {
        return (sgn * ceiling(y / 10^p) * 10^p)
    }
    return (sgn * floor(y / 10^p) * 10^p)
}


# ############################################################################
# Function ylimits finds optimal limits on plot according
# to input data
# ############################################################################
# Input
#   ymin - the minimum value of the data
#   ymax - the maximum value of the data
#   forcezero - Boolean, should 0 be always included in the range?
#               (default FALSE)
# Output
#   Vector with three elements denoting minimum, maximum and step
#   for y axis
# ############################################################################
ylimits <- function (ymin, ymax, forcezero = FALSE)
{
    if ((ymin == 0) & (ymax == 0)) { return (c(-1, 1, 1)) }
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

    return (c(ylo, yhi, yby))
}

# ############################################################################
# The xlimits function finds optimal limits on x-axis according
# to input data
# ############################################################################
# Input
#   ymin - the minimum value of the data
#   ymax - the maximum value of the data
# Output
#   Scale for x- axis
# ############################################################################
xlimits <- function (xmin, xmax)
{
    if ((xmin == 1) & (xmax == 1))
        { return (c(1, 1)) }

    if ((xmin == 1) & (xmax < 10))
        { return (c(xmin : xmax))}

    if ((xmin == 1) & (xmax >= 10))
        {return (c(1, seq(5, 5 * floor(xmax /5), 5)))}
}

# ############################################################################
# The roundpow10 function rounds a number to the p power of 10
# ############################################################################
# Input
#   y - a number to be rounded
#   p - a number specifying the power
#   up - parameter denoting whether rounding function
#         should be ceiling (TRUE) or floor (false)
# Output
#   A rounded number
# ############################################################################
roundpow10 <- function(y, p, up = TRUE)
{
    if (y == 0) return (y)
    if (up) return (ceiling(y / 10^p) * 10^p)
    else return (floor(y / 10^p) * 10^p)
}


# ############################################################################
# Function returns prepared labels for axis from sequence parameters
#   (from result of ylimits function)
# ############################################################################
# Input
#   yl - initial value
#   yh - limit value
#   yb - step value
# Output
#   character vector representing axis labels
# ############################################################################
pretty_labels <- function(yl, yh, yb)
{
    return (format(round(seq(yl, yh, yb), find_signif(yb)), scientific = FALSE))
}

# ############################################################################
# Function finds number of significant digits from specified value
# ############################################################################
# Input
#   x - numeric value
#   tol - (optional) definied tolerance
# Output
#   not negative integer representing number of significant digits
# ############################################################################
find_signif <- function(x, tol = 1e-2)
{
    if (x == 0) { return (0) }
    p <- 0
    while (abs(1 - round(x, digits = p) / x) > tol) {
        p <- p + 1
    }
    return (p)
}
