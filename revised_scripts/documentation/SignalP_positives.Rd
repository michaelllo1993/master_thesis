\name{SignalP_positives.pl}
\alias{SignalP_positives.pl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Analyze the SignalP 4.1 results
}
\description{
This script analyzes the SignalP 4.1 output files and extract the positives
}
\usage{
SignalP_positives.pl(organism_name)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{organism_name}{The name of the organism of interest. Must be identical to the <organism_name>_all.out file produced by the previous script \code{\link{run_signalp.pl}}}
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{SignalP4_positives_<organism_name>.out}{This file contains information on the ENS ID of proteins with signal peptides and starting point and anding point of ones}
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
#cd to the dir with  <organism_name>_all.out
perl Signal__positives.pl homo_sapiens
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ signal peptide }% use one of  RShowDoc("KEYWORDS")
\keyword{ postitives }% __ONLY ONE__ keyword per line
