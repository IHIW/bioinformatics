\name{hlaDataInit}
\alias{hlaDataInit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ uselies
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hlaDataInit(hlaVersion = "3.25.0")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hlaVersion}{
%%     ~~Describe \code{hlaVersion} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (hlaVersion = "3.25.0") 
{
    genfa <- fasta.index(file.path(hlaVersion, "hla_gen.fasta"), 
        seqtype = "DNA")
    genfa$name <- paste("HLA-", sapply(strsplit(genfa$desc, " "), 
        "[[", 2), sep = "")
    nucfa <- fasta.index(file.path(hlaVersion, "hla_nuc.fasta"), 
        seqtype = "DNA")
    nucfa$name <- paste("HLA-", sapply(strsplit(nucfa$desc, " "), 
        "[[", 2), sep = "")
    features <- read.table(file.path(hlaVersion, "hlaxml.tsv"), 
        sep = "\t", header = F)
    colnames(features) <- c("name", "accnum", "version", "feature", 
        "fid", "start", "end")
    reftable <- read.table(file.path(hlaVersion, "hlaRefAllele.tsv"), 
        sep = "\t", header = F)
    refAllele <- reftable$V2
    names(refAllele) = reftable$V1
    refAllele <- refAllele
    result <- list(genfa = genfa, nucfa = nucfa, features = features, 
        refAllele = refAllele)
    return(result)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
