\name{plot_simulation}
\alias{plot_simulation}
\title{
Plot a gecon_simulation object
}
\description{
The \code{plot_simulation} function plots the simulations or saves them as \code{.eps} files in the model's subdirectory \code{/plots}.
}
\usage{
plot_simulation(sim_obj, to_tex = NULL, to_eps = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sim_obj}{
an object of the \code{gecon_simulation} class.
}
  \item{to_tex}{
logical. If TRUE, the plots are added to a .tex file.
}
 \item{to_eps}{
logical. if TRUE, plot(s) are saved as \code{.eps} file(s) in the model's subdirectory \code{/plots}.
}
}
\value{
If the number of variables of interest is greater then five, more then one plots for each impulse are created (max. 5 variables on each plot).
Separate plots are created for all the impulses, if the \code{compute_irf} function has been used for generating simulations.
}
\author{
Karol Podemski \email{karol.podemski@gmail.com}, \cr
Kaja Retkiewicz-Wijtiwiak \email{kaja.retkiewicz@gmail.com}
}
\examples{
# copy the example to the current working directory
file.copy(from = file.path(system.file("examples", package = 'gEcon'),
          'rbc.gcn'), to = getwd())

# make and load the model
rbc <- make_model('rbc.gcn')

# solve the model
rbc <- steady_state(rbc)
rbc <- solve_pert(rbc)

# compute and plot the IRFs
rbc <- set_shock_cov_mat(rbc, shock_matrix = matrix(0.01, 1, 1), shock_order = 'epsilon_Z')
irf_rbc <- compute_irf(rbc)
plot_simulation(irf_rbc)
}
