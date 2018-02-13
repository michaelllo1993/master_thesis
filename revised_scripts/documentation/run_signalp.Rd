\name{run_signalp.pl}
\alias{run_signalp.pl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Run SignalP 4.1 analysis
}
\description{
This script runs the SignalP 4.1 analysis on (trimmed) peptide sequences to detect signal pepides
}
\usage{
run_signalp.pl(organism_name, no_of_files)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{organism_name}{The name of organism of interest. Must be identical to the shortened_seq_<organism_name>
  \item{no_of_files}{The number of files produced with previous script \code{\link{ensembl_peptide_seq_cut.pl}}}
}}

\value{
{<organism_name>_all.out}. The output is a file with the results of signal peptides detection process
}

\author{
Michal Stolarczyk
}

\examples{
#cd to dir with shortened_seq_<organism_name>
perl run_signalp.pl homo_sapiens 100
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ signal peptide }% use one of  RShowDoc("KEYWORDS")
\keyword{ detection }% __ONLY ONE__ keyword per line
