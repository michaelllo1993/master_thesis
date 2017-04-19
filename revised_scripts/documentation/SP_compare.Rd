\name{SP_compare.R}
\alias{SP_compare.R}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compare signal peptides of orthologous proteins
}
\description{
This R script can be used to compare the signal pepide sequences of organism of intrerest against one chosen orthologue
}
\usage{
SP_compare.R(prep2comp_file,organism_of_interest_exclusive_code, orthologous_exclusive_code)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{prep2comp_file}{A file name of prep2comp_<OoI>_all.csv produced by previous script \code{\link{prep2comp.pl}}}
  \item{organism_of_interest_exclusive_code}{Ensembl exculisve code of organism of inerest, e.g. mouse is ENSMUSP}
  \item{orthologous_exclusive_code}{Likewise for orthologous proteins}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
This script prints summary of the comparison on the screen
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Michal Stolarczyk
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
#cd to the dir with files mentioned above
Rscript name filename (prep2comp_<OoI>_all.csv) OoI(ensembl exclusive code) 1Ortho(ensembl exclusive code)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ compare }% use one of  RShowDoc("KEYWORDS")
\keyword{ sequences }% __ONLY ONE__ keyword per line
