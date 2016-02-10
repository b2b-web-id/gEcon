\name{make_model}
\alias{make_model}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Make model from .gcn file
}
\description{
This function calls the dynamic library, parses the .gcn model file, generates an .R file,
and loads it into a \code{gecon_model} class object.
}
\usage{
make_model(model_file)
}
\arguments{
  \item{model_file}{
the name of the .gcn file containing model formulation. It must be ended with a \code{.gcn} extension.
}
}
\details{
Cf. gEcon manual, chapters "Model description language" and "Derivation of First Order Conditions".
}
\value{
An object of the \code{gecon_model} class representing the model. Generic functions such as \code{print} and \code{summary} allow to show the model elements.
}
\seealso{
\code{\link{load_model}} function to load already created .R file with model.
}
\author{
Grzegorz Klima \email{gklima@users.sourceforge.net}\cr
Karol Podemski \email{karol.podemski@gmail.com}\cr
Kaja Retkiewicz-Wijtiwiak \email{kaja.retkiewicz@gmail.com}
}
\note{
When the function is called, an R file with the same name as the .gcn file is created in the the .gcn file directory. Additional files such as a Latex documentation files or a logfile may be created when relevant options are set in the .gcn file.
}
\examples{
# copy the example to the current working directory
file.copy(from = file.path(system.file("examples", package = 'gEcon'),
          'rbc.gcn'), to = getwd())

# make and load the model
rbc <- make_model('rbc.gcn')
print(rbc)
}
