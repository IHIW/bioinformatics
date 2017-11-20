\name{hlaAlign2Novel}
\alias{hlaAlign2Novel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{get Novel Polymorphism from alignment
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hlaAlign2Novel(hlaData)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hlaData}{
%%     ~~Describe \code{hlaData} here~~
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
function (hlaData) 
{
    attach(hlaData)
    msamatrix <- t(as.matrix(triomsa))
    trioIndex <- c(0, 0, 0)
    names(trioIndex) = colnames(msamatrix)
    mismatches <- data.frame(reftype = factor(levels = c("gen", 
        "nuc")), refseq = character(), conseq = character(), 
        start = integer(), end = integer(), stringsAsFactors = FALSE)
    for (i in 1:nrow(msamatrix)) {
        for (j in 1:ncol(msamatrix)) {
            if (msamatrix[i, j] != "-" && msamatrix[i, j] != 
                "|") {
                trioIndex[j] = trioIndex[j] + 1
            }
        }
        if (trioIndex["input"] >= 1 && trioIndex["input"] <= 
            nchar(conseq)) {
            if (msamatrix[i, "input"] != msamatrix[i, "nuc"]) {
                if (msamatrix[i, "nuc"] == "-") {
                  if (trioIndex["nuc"] > 0 && !(trioIndex["nuc"] \%in\% 
                    nucFeature$end)) {
                    mismatches[nrow(mismatches) + 1, ] <- data.frame(reftype = "nuc", 
                      refseq = "-", conseq = msamatrix[i, "input"], 
                      start = trioIndex["nuc"], end = trioIndex["nuc"], 
                      stringsAsFactors = FALSE)
                  }
                  else {
                    if (msamatrix[i, "input"] != msamatrix[i, 
                      "gen"]) {
                      mismatches[nrow(mismatches) + 1, ] <- data.frame(reftype = "gen", 
                        refseq = msamatrix[i, "gen"], conseq = msamatrix[i, 
                          "input"], start = trioIndex["gen"], 
                        end = trioIndex["gen"], stringsAsFactors = FALSE)
                    }
                  }
                }
                else {
                  mismatches[nrow(mismatches) + 1, ] <- data.frame(reftype = "nuc", 
                    refseq = msamatrix[i, "nuc"], conseq = msamatrix[i, 
                      "input"], start = trioIndex["nuc"], end = trioIndex["nuc"], 
                    stringsAsFactors = FALSE)
                }
            }
            if (trioIndex["input"] == nchar(conseq)) {
                break
            }
        }
    }
    nucFeature$reftype = "nuc"
    genFeature$reftype = "gen"
    refFeature = rbind(nucFeature, genFeature)
    mismatches$fid = ""
    mismatches$feature_start = 0
    attach(refFeature)
    for (i in 1:nrow(mismatches)) {
        mismatches[i, c("fid", "feature_start")] = refFeature[reftype == 
            mismatches$reftype[i] & mismatches$start[i] >= start & 
            mismatches$end[i] <= end, c("fid", "start")]
    }
    detach(refFeature)
    reformMismatches <- data.frame(reftype = factor(levels = c("gen", 
        "nuc")), refseq = character(), conseq = character(), 
        start = integer(), end = integer(), fid = character(), 
        feature_start = integer(), stringsAsFactors = FALSE)
    reformMismatches[nrow(reformMismatches) + 1, ] <- mismatches[1, 
        ]
    for (i in 2:nrow(mismatches)) {
        if (mismatches$refseq[i] == "-") {
            mismatches$end[i] = mismatches$end[i] - 1
            if (reformMismatches$refseq[nrow(reformMismatches)] == 
                "-" & reformMismatches$start[nrow(reformMismatches)] == 
                mismatches$start[i] & reformMismatches$fid[nrow(reformMismatches)] == 
                mismatches$fid[i]) {
                reformMismatches$conseq[nrow(reformMismatches)] = paste(reformMismatches$conseq[nrow(reformMismatches)], 
                  mismatches$conseq[i], sep = "")
                next
            }
        }
        else if (mismatches$conseq[i] == "-") {
            if (reformMismatches$conseq[nrow(reformMismatches)] == 
                "-" & reformMismatches$end[nrow(reformMismatches)] == 
                mismatches$start[i] - 1 & reformMismatches$fid[nrow(reformMismatches)] == 
                mismatches$fid[i]) {
                reformMismatches$refseq[nrow(reformMismatches)] = paste(reformMismatches$refseq[nrow(reformMismatches)], 
                  mismatches$refseq[i], sep = "")
                reformMismatches$end[nrow(reformMismatches)] = mismatches$end[i]
                next
            }
        }
        reformMismatches[nrow(reformMismatches) + 1, ] <- mismatches[i, 
            ]
    }
    reformMismatches$start = reformMismatches$start - reformMismatches$feature_start
    reformMismatches$end = reformMismatches$end - reformMismatches$feature_start + 
        1
    detach(hlaData)
    reformMismatches
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
