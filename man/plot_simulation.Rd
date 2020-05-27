\name{plot_simulation}
\alias{plot_simulation}
\title{
Plot simulation results stored in an object of \code{gecon_simulation} class
}
\description{
The \code{plot_simulation} function plots the simulations or saves them as \code{.eps} files in the model's subdirectory \code{/plots}.
}
\usage{
plot_simulation(sim_obj, to_eps = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sim_obj}{
an object of \code{gecon_simulation} class.
}
 \item{to_eps}{
logical. if TRUE, the plot(s) is (are) saved as \code{.eps} file(s) in the model's subdirectory \code{/plots}
and is (are) added to a .results.tex file.
}
}
\value{
If more than five variables have been selected for simulations, at least two plots are created (max. 5 variables can be depicted on one plot).
Separate plots are created for all the impulses, if the \code{compute_irf} function has been used for generating simulations.
}
\examples{
# copy the example to the current working directory
file.copy(from = file.path(system.file("examples", package = "gEcon"),
          "rbc.gcn"), to = getwd())

# make and load the model
rbc <- make_model("rbc.gcn")

# solve the model
rbc <- steady_state(rbc)
rbc <- solve_pert(rbc)

# compute and plot the IRFs
rbc <- set_shock_cov_mat(rbc, cov_matrix = matrix(0.01, 1, 1),
                         shock_order = "epsilon_Z")
irf_rbc <- compute_irf(rbc)
plot_simulation(irf_rbc)
}
