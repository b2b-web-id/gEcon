\name{gEcon-package}
\alias{gEcon-package}
\docType{package}
\title{
General Equilibrium Economic Modelling Language and Solution Framework (gEcon)
}
\description{
Package for developing and solving dynamic (stochastic) and static general equilibrium models.
}
\details{
\code{gEcon} is a framework for developing and solving large scale dynamic
(stochastic) & static general equilibrium models. It consists of model description
language and an interface with a set of solvers in R. It was developed at the Department
for Strategic Analyses at the Chancellery of the Prime Minister of the Republic
of Poland as a part of a project aiming at construction of large scale DSGE & CGE
models of the Polish economy.

Publicly available toolboxes used in RBC/DSGE modelling require
users to derive the first order conditions (FOCs) and linearisation equations by pen
& paper (e.g. Uhlig's tool-kit) or at least require manual
derivation of the FOCs (e.g. Dynare).
Derivation of FOCs is also required by GAMS and GEMPACK - probably the two most popular frameworks used in CGE modelling.
Owing to the development of an algorithm for automatic derivation
of first order conditions and implementation of a comprehensive symbolic computations library,
\code{gEcon} allows users to describe their models in terms
of optimisation problems of agents. To authors' best knowledge there is no other
publicly available framework for writing and solving DSGE & CGE models in this
natural way. Writing models in terms of optimisation problems instead of the FOCs
is far more natural to an economist, takes off the burden of tedious
differentiation, and reduces the risk of making a mistake. \code{gEcon} allows users
to focus on economic aspects of the model and makes it possible to design
large-scale (100+ variables) models. To this end, \code{gEcon} provides template mechanism
(similar to those found in CGE modelling packages), which allows to declare similar
agents (differentiated by parameters only) in a single block.
Additionally, \code{gEcon} can automatically produce a draft of LaTeX documentation
for a model.

The model description language is simple and intuitive. Given optimisation problems, constraints and identities, computer derives the FOCs, steady-state equations, and linearisation matrices automatically. Numerical solvers can be then employed to determine the steady state and approximate equilibrium laws of motion around it.
}
\author{
Karol Podemski \email{gecon.maintenance@gmail.com}  (since 10.2012), \cr
Kaja Retkiewicz-Wijtiwiak \email{kaja.retkiewicz@gmail.com} (since 10.2012), \cr
Grzegorz Klima \email{gklima@users.sourceforge.net} (lead developer 10.2012-03.2018)

Maintainer: Karol Podemski \email{gecon.maintenance@gmail.com}
}
\references{
Cf. gEcon users' guide, which can be found at \url{http://gecon.r-forge.r-project.org/}.
}
\keyword{ package }
\examples{
# copy the example to the current working directory
file.copy(from = file.path(system.file("examples", package = "gEcon"),
          "rbc.gcn"), to = getwd())

# make and load the model
rbc <- make_model("rbc.gcn")

# solve the model
rbc <- steady_state(rbc)
rbc <- solve_pert(rbc)

# compute and print correlations
rbc <- set_shock_cov_mat(rbc, cov_matrix = matrix(0.01, 1, 1),
                         shock_order = "epsilon_Z")
rbc <- compute_model_stats(rbc, ref_var = "Y")
get_model_stats(rbc)

# compute and plot the IRFs
irf_rbc <- compute_irf(rbc)
plot_simulation(irf_rbc)
summary(rbc)
}
