\name{get_index_sets}
\alias{get_index_sets}
\title{List of index sets}
\description{
The \code{get_index_sets} function retrieves a list of all the index sets specified in the \code{.gcn} file.
}
\usage{
get_index_sets(model)
}
\arguments{
\item{model}{an object of \code{gecon_model} class.}
}
\value{
The function returns a list of index sets.
Each element of the list corresponds to one set and contains the set components' names as a character vector.
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

# retrieve the index sets
pure_exchange_ind_sets <- get_index_sets(pure_exchange_t)
}
