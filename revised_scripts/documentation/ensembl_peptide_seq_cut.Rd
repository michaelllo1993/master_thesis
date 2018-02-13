\name{ensembl_peptide_seq_cut.pl}
\alias{ensembl_peptide_seq_cut.pl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Trimm peptide sequences
}
\description{
This script is indtended to trimm the peptide sequences to 69 amino acids long which facilitates a memory-friendly data operations in signal peptide discovery process. Signal peptides are located in the first tens of amino acid sequences.
}
\usage{
ensembl_peptide_seq_cut.pl(organism_name)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{organism_name}{Organism of interest name - identical to the input file name}
}
\details{
The script requires a .txt file with peptide sequences in FASTA format as an input. The file must be named <organism_name>.txt. 
The output of this script are multiple files (each containing a maximum of 1000 trimmed peptide sequences) for the analysis with SignalP 4.1 software.
}
\value{
\item{shortened_seq_<organism_name>}{Multiple files with trimmed peptide sequences}
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
#cd to the dir with homo_sapiens.txt file
perl ensembl_peptide_seq_cut.pl homo_sapiens
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ trimm}% use one of  RShowDoc("KEYWORDS")
\keyword{ sequence }% __ONLY ONE__ keyword per line
