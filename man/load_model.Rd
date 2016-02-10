\name{load_model}
\alias{load_model}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Load model from .R file
}
\description{
The \code{load_model} function loads the already generated .R file with the model and creates an object of  the \code{gecon_model} class.
}
\usage{
load_model(model_file)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model_file}{
the name of the .R file containing the model's functions and variables. It can be a name of file or a name of file ending with a .model.R extension.
}
}
\details{
The .R file with the model specification has to be created first. It can be done by using the \code{make_model} command and the \code{gcn} file model specification or manually.
}
\value{
An object of the \code{gecon_model} class representing the model. Generic functions such as \code{print} and \code{summary} allow to show the model's elements.
}
\seealso{
the \code{\link{make_model}} function in order to create an R file with the model elements based on the model specification.
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
print(rbc)
# load the already generated model
rbc2 <- load_model('rbc.model.R')
print(rbc2)
}
