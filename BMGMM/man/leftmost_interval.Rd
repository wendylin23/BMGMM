\name{leftmost_interval}
\alias{leftmost_interval}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
leftmost_interval(U, lambda, ok)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{U}{
%%     ~~Describe \code{U} here~~
}
  \item{lambda}{
%%     ~~Describe \code{lambda} here~~
}
  \item{ok}{
%%     ~~Describe \code{ok} here~~
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
function (U, lambda, ok) 
{
    H = 0.5 * log(2) + 2.5 * log(pi) - 2.5 * log(lambda) - (pi^2)/(2 * 
        lambda) + 0.5 * lambda
    lU = log(U)
    Z = 1
    X = exp((-pi^2)/(2 * lambda))
    K = lambda/pi^2
    j = 0
    j = j + 1
    Z = Z - K * X^((j^2) - 1)
    if (H + log(Z) > lU) 
        ok = 1
    else {
        j = j + 1
        Z = Z + ((j + 1)^2) * X^(((j + 1)^2) - 1)
        if (H + log(Z) < lU) 
            ok = 0
    }
    return(ok)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
