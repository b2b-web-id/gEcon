\name{get_par_names_by_index}
\alias{get_par_names_by_index}
\title{Parameters corresponding to given indices}
\description{
The \code{get_par_names_by_index} function allows to retrieve the names of parameters with given indices.
}
\usage{
get_par_names_by_index(model, index_names)
}
\arguments{
\item{model}{an object of \code{gecon_model} class.}
\item{index_names}{a character vector of indices.}
}
\value{
The function returns a character vector of relevant parameters' names.
}
\details{
Cf. gEcon users' guide, chapter `Templates'.
}
\examples{
# copy the example to the current working directory
file.copy(from = file.path(system.file("examples", package = "gEcon"),
          "pure_exchange.gcn"), to = getwd())

# make and load the model
pure_exchange_t <- make_model("pure_exchange.gcn")

# model calibration
pure_exchange_t <- set_free_par(pure_exchange_t,
                                free_par= c("alpha__A__1" = 0.3, "alpha__A__2" = 0.4,
                                            "alpha__A__3" = 0.3, "alpha__B__1" = 0.3,
                                            "alpha__B__2" = 0.4, "alpha__B__3" = 0.3,
                                            "e_calibr__A__1" = 3, "e_calibr__B__1" = 1,
                                            "e_calibr__A__2" = 2, "e_calibr__B__2" = 1,
                                            "e_calibr__A__3" = 1, "e_calibr__B__3" = 3))

# get all parameters associated with agent A
par_names_A <- get_par_names_by_index(pure_exchange_t, index_names = "A")
par_info(pure_exchange_t, par_names_A)

# get all parameters associated with agent B
par_names_B <- get_par_names_by_index(pure_exchange_t, index_names = "B")
par_info(pure_exchange_t, par_names_B)
}
