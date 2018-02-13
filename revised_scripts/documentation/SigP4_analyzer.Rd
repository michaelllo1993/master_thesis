\name{SigP4_analyzer.pl}
\alias{SigP4_analyzer.pl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extract the signal peptides
}
\description{
This script extracts the signal peptide sequences based on SignalP 4.1 results
}
\usage{
SigP4_analyzer.pl(organism_name)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{organism_name}{The name of the organism of interest. Must be indentical to the SignalP4_positives_<organism_name>.out and shortened_seq_<organism_name> files.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{sigp_ALL_<organism_name>.out}{A file with protein ENS ID, cDNA ENS ID, sequence and cleavage site of the proteins comprising signal peptides}
\item{sigpL_ALL_<organism_name>.out}{A file with protein ENS ID, cDNA ENS ID, sequence and cleavage site of signalm peptides and starts and lengths of L-SAARs the proteins comprising signal peptides with L-SAARs}
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
#cd to the dir with Signal_4_positives_<organism_name>.out & shortened_seq_<organism_name>
perl SigP4_analyzer.pl homo_sapiens
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ signal peptides }% use one of  RShowDoc("KEYWORDS")
\keyword{ L-SAARs }% __ONLY ONE__ keyword per line
