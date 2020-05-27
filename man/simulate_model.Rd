\name{simulate_model}
\alias{simulate_model}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulation of the model
}
\description{
The \code{simulate_model} function simulates the impact of shock paths specified by the user
on the model's variables. In particular, it allows to compute the impact of negative shocks. 
}
\usage{
simulate_model(model, variables = NULL, shocks = NULL,
               shock_path, sim_length = 40)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
an object of \code{gecon_model} class.
}
  \item{variables}{
the names or indices of variables whose paths are to be simulated. By default all variables are selected.
}
  \item{shocks}{
the names or indices of shocks corresponding to consecutive rows of the shock_path matrix. 
If missing, the \code{rownames} of the \code{shock_path} matrix are used.
}
  \item{shock_path}{
a matrix simulated paths of \code{shocks} in rows.
}
  \item{sim_length}{
the length of simulation path, the default value is 40.
}
}
\value{
An object of \code{\link{gecon_simulation}} class storing simulated paths of the model's variables.
}
\seealso{
The \code{\link{random_path}} function generates random paths of shocks for the system behaviour simulation.
This function also returns an object of \code{gecon_simulation} class. 
Generic functions such as \code{print} and \code{summary} provide information about simulations' results.
The \code{plot_simulation} function allows to visualize the impact of simulated shock paths on the model's variables.
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

# set shock path
shock_path <- matrix(0, 1, 4)
shock_path[1, 1] <- -0.05
shock_path[1, 4] <- -0.05
# simulate model and plot results
rbc_sim <- simulate_model(rbc, variables = c("K_s", "C", "Z", "I", "Y"),
                          shocks = "epsilon_Z", shock_path = shock_path)
plot_simulation(rbc_sim)
}
