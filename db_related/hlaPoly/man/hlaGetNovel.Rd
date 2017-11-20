\name{hlaGetNovel}
\alias{hlaGetNovel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{gogog
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hlaGetNovel(hlaData, locus, nucAllele, conseq)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hlaData}{
%%     ~~Describe \code{hlaData} here~~
}
  \item{locus}{
%%     ~~Describe \code{locus} here~~
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
function (hlaData, locus, nucAllele, conseq) 
{
    genAllele = as.character(hlaData$refAllele[locus])
    hlaData <- hlaGetAlign(hlaData, genAllele, nucAllele, conseq)
    hlaAlign2Novel(hlaData)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
