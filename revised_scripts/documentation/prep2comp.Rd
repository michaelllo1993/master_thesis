\name{prep2comp.pl}
\alias{prep2comp.pl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Prepare Needle output to compare}
}
\description{
This script prepares the results of pairwise alignment by needle software
}
\usage{
prep2comp.pl(ALL_<organism_of_interest>_Ortho.needle,<organism_of_interest>_and_orthologues_protein_ids.tsv,organism_of_interest,orthologues)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ALL_<organism_of_interest>_Ortho.needle}{The output of pervious script \code{\link{Needle_wrapper.pl}}. The result of pairwise alignment of orthologous proteins}
  \item{<organism_of_interest>_and_orthologues_protein_ids.tsv}{A tsv file with information on the organism's of interest orthologous proteins}
  \item{organism_of_interest}{The name of the organism of interest. The one that other sequences will be aligned to}
  \item{orthologues}{Other organisms to be aligned to sequneces of organism of interest}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{prep2comp_<organism_of_interest>_all.csv}{A four column csv file with names and sequences to be compared in thte downstream analysis}
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


\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
#cd to the dir with files described above
perl prep2comp.pl ALL_organism_Ortho_(all/sp).needle <ortho_dataset.tsv> <organism_of_interest> <other_organisms>
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ compare }% use one of  RShowDoc("KEYWORDS")
\keyword{ alignment }% __ONLY ONE__ keyword per line
