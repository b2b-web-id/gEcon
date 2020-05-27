# ############################################################################
# This file is a part of gEcon.                                              #
#                                                                            #
# (c) Chancellery of the Prime Minister of the Republic of Poland 2012-2015  #
# (c) Grzegorz Klima, Karol Podemski, Kaja Retkiewicz-Wijtiwiak 2015-2018    #
# License terms can be found in the file 'LICENCE'                           #
#                                                                            #
# Authors: Karol Podemski, Kaja Retkiewicz-Wijtiwiak                         #
# ############################################################################

# ############################################################################
# The pert_solver function solves the first order perturbation
# represented by four sparse matrices.
# ############################################################################
# Input
#   matrices [class Matrix] - model matrices (corresponding to y[t-1], y[t],
#                             E[t](y[t+1]) and eps[t])
#   is_stochastic - a logical value, TRUE if model is stochastic and FALSE,
#                                   if model is deterministic
#   tol - a numeric value, the tolerance for treating column elements as zeros
#   solver - the name of RE solver, currently only "gensys" -
#            Christopher A. Sims' solver is supported
# Output
#   A list containing the following fields:
#       eig - eigenvalues of the system
#       P - matrix of impact of state variables on state variables,
#       Q - matrix of impact of exog shocks on state variables
#               (only when the model is stochastic),
#       R - matrix of impact of state variables on non-state variables,
#       S - matrix of impact of exog shocks on non-state
#               variables (only when the model is stochastic),
#       state_var_indices - indices of state variables
#       exit_code - 0 - solution exists and is unique,
#                   1 - solution exists but is not unique (sunspots),
#                   2 - solution not found,
#       exit_info - information returned by the solver
#       norm_deter - the norm of residuals of the model's deterministic part
#                        (solution check),
#       norm_stoch - the norm of residuals of the model's stochastic part
#                        (solution check)
# ############################################################################
pert_solver <- function(matrices, is_stochastic, tol = 1e-8, solver)
{
    if (!is.character(solver) || length(solver) != 1) {
        stop("invalid solver argument")
    }
    solver_list <- list("gensys")

    pos <- as.numeric(match(solver, solver_list, nomatch = 0))

    # solver not specified in solver list
    if (pos == 0) {
        stop(paste0("solver called \'", solver, "\' is not supported"), call. = FALSE)
    }
    # gensys
    if (pos == 1) {
        return (gensys_wrap(matrices, is_stochastic, tol))
    }
}




