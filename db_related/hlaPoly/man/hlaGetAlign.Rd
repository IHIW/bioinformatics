\name{hlaGetAlign}
\alias{hlaGetAlign}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{aaaa
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hlaGetAlign(hlaData, genAllele, nucAllele, conseq)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hlaData}{
%%     ~~Describe \code{hlaData} here~~
}
  \item{genAllele}{
%%     ~~Describe \code{genAllele} here~~
}
  \item{nucAllele}{
%%     ~~Describe \code{nucAllele} here~~
}
  \item{conseq}{
%%     ~~Describe \code{conseq} here~~
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
function (hlaData, genAllele, nucAllele, conseq) 
{
    attach(hlaData)
    typeIsGen = T
    genIndex = which(genfa$name == genAllele)
    nucIndex = which(genfa$name == nucAllele)
    if (length(nucIndex) == 0) {
        nucIndex = which(nucfa$name == nucAllele)
        typeIsGen = F
    }
    genFeature = features[features$name == genAllele, ]
    genFeature <- genFeature[order(genFeature$start), ]
    nucFeature = features[features$name == nucAllele, ]
    nucFeature <- nucFeature[order(nucFeature$start), ]
    genseq <- readDNAStringSet(genfa[genIndex, ])
    if (typeIsGen) {
        nucseq <- readDNAStringSet(genfa[nucIndex, ])
    }
    else {
        nucseq <- readDNAStringSet(nucfa[nucIndex, ])
    }
    pipeGenseq = DNAStringSet("+")
    for (i in 1:nrow(genFeature)) {
        pipeGenseq = xscat(pipeGenseq, subseq(genseq, start = genFeature$start[i], 
            end = genFeature$end[i]), "+")
    }
    pipeNucseq = DNAStringSet("+")
    for (i in 1:nrow(nucFeature)) {
        pipeNucseq = xscat(pipeNucseq, subseq(nucseq, start = nucFeature$start[i], 
            end = nucFeature$end[i]), "+")
    }
    conseq <- DNAStringSet(conseq)
    trio <- append(pipeGenseq, c(pipeNucseq, conseq))
    names(trio) = c("gen", "nuc", "input")
    triomsa <- msa(trio, gapExtension = 0, order = "input")
    detach(hlaData)
    hlaData$triomsa = triomsa
    hlaData$genFeature = genFeature
    hlaData$nucFeature = nucFeature
    return(hlaData)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
