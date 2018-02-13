\name{Needle_wrapper.pl}
\alias{Needle_wrapper.pl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Alignment of orthologous amino acid sequences
}
\description{
This script performs the pairwise alignment of orthologous amino acid sequences
}
\usage{
Needle_wrapper.pl(region,organism_and_orthologues_protein_ids.tsv,sigp(L)_ALL.out/organism_of_interest.out,organism_of_interest,orthologues)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{reqion}{One of following: "SP", "WHOLE" specifying the reqion of the protein to be aligned. If SP - the sigp(L)_ALL.out is needed. Otherwise - the organism_of_interest.out (the output of prep2Needle.pl)}
  \item{organism_and_orthologues_protein_ids.tsv}{A tsv file with information on the organism's of interest orthologous proteins}
  \item{sigp(L)_ALL.out}{An out of the previous script with signal peptide sequences or signal peptide with L-SAARs sequences}
  \item{organism_of_interest}{The name of the organism of interest. The one that other sequences will be aligned to}
  \item{orthologues}{Other organisms to be aligned to sequneces of organism of interest}
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{ALL_<organism_of_interest>_Ortho.needle}{The output file of all ortologous sequences aligned to sequences of organism of interest}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Michal Stolarczyk
}
\note{
<organism_and_orthologues_protein_ids.tsv> must have the same column order as <organism_of_interest>+<orthologues>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
#cd to the dir with all files described above
perl Needle_wrapper.pl <SP> <organism_and_orthologues_protein_ids.tsv> <sigp(L)_ALL.out> <organism_of_interest> <orthologues>
perl Needle_wrapper.pl <WHOLE> <organism_and_orthologues_protein_ids.tsv> <organism_of_interest.out> <organism_of_interest> <orthologues>
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ alignment }% use one of  RShowDoc("KEYWORDS")
\keyword{ signal peptides }% __ONLY ONE__ keyword per line
