\name{get_model_info}
\alias{get_model_info}
\title{Accessing information about the name and the creation date of the model}
\description{
The \code{get_model_info} function returns a character vector
with information about the model.
}
\usage{
get_model_info(model)
}
\arguments{
\item{model}{an object of \code{gecon_model} class.}
}
\value{
The function returns a character vector of length 3, containing information
about the model: the input file name, the input file path, and the date of creation.
}
\seealso{
For details, see \code{\link{gecon_model-class}}.
}
\examples{
# copy the example to the current working directory
file.copy(from = file.path(system.file("examples", package = "gEcon"),
          "rbc.gcn"), to = getwd())

# make and load the model
rbc <- make_model("rbc.gcn")

# retrieve and show model information
model_info <- get_model_info(rbc)
model_info
}
