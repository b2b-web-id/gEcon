\name{load_model}
\alias{load_model}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Load model from \code{.model.R} file
}
\description{
The \code{load_model} function loads the already generated \code{.model.R} file and creates an object of the \code{gecon_model} class.
}
\usage{
load_model(filename)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{
the path to the \code{.model.R} file containing the model's functions and variables. The \code{.model.R} extension is optional.
}
}
\details{
The \code{.model.R} file with the model specification has to be created first. It can be done with the \code{make_model} function taking a \code{.gcn} file with model specification as an argument.
}
\value{
An object of the \code{gecon_model} class representing the model. Generic functions such as \code{print} and \code{summary} allow to show the model's elements.
}
\seealso{
The \code{\link{make_model}} function in order to create an .R file with the model specification.
}
\examples{
# copy the example to the current working directory
file.copy(from = file.path(system.file("examples", package = "gEcon"),
          "rbc.gcn"), to = getwd())

# make and load the model
rbc <- make_model("rbc.gcn")
print(rbc)
# load the already generated model
rbc2 <- load_model("rbc.model.R")
print(rbc2)
}
