\name{get_shock_names_by_index}
\alias{get_shock_names_by_index}
\title{Shocks corresponding to given indices}
\description{
The \code{get_shock_names_by_index} function retrieves the names of shocks with given indices.
}
\usage{
get_shock_names_by_index(model, index_names)
}
\arguments{
\item{model}{an object of \code{gecon_model} class.}
\item{index_names}{a character vector of the chosen indices.}
}
\value{
The function returns a character vector of relevant shock names.
}
\details{
Cf. gEcon manual, chapter "Templates".
}
\examples{
# copy the example to the current working directory
file.copy(from = file.path(system.file("examples", package = 'gEcon'),
          'home_prod_templ.gcn'), to = getwd())

# make and load the model
home_prod_templ <- make_model('home_prod_templ.gcn')

# get shocks affecting home production technology
hp_shocks <- get_shock_names_by_index(home_prod_templ, 'H')

# print information obout the selected shocks
shock_info(home_prod_templ, hp_shocks)
}
