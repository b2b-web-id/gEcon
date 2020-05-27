\name{make_model}
\alias{make_model}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Make model from \code{.gcn} file
}
\description{
This function calls the dynamic library, parses the .gcn model file, generates a .model.R file,
and loads it into a \code{gecon_model} class object.
}
\usage{
make_model(filename, variant = NULL, variant_name = NULL)
}
\arguments{
  \item{filename}{
the path to the .gcn file containing model formulation. The \code{.gcn} extension is optional.
}
  \item{variant}{
(optional) an integer value (0-99), selected variant
}
  \item{variant_name}{
(optional) an character value, variant name appended to the name of the resulting .gcn file
}
}
\details{
Cf. gEcon users' guide, chapters `Model description language' and `Derivation of First Order Conditions'.
}
\value{
An object of the \code{gecon_model} class representing the model. Generic functions such as \code{print} and \code{summary} allow to show the model elements.
}
\seealso{
\code{\link{load_model}} function to load a \code{.model.R} file created earlier.
}
\note{
When the function is called, an R file with the same name as the .gcn file is created in the .gcn file directory. Additional files such as Latex documentation files or a logfile may be created when relevant options are set in the .gcn file.
}
\examples{
# copy the example to the current working directory
file.copy(from = file.path(system.file("examples", package = "gEcon"),
          "rbc.gcn"), to = getwd())

# make and load the model
rbc <- make_model("rbc.gcn")
print(rbc)
}
