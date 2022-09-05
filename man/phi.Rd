\name{phi}
\alias{phi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ renturns a value for a measure phi


}
\description{
renturns a value for a measure phi for a given configuration and vector of parameters q. For measure definition see Pablo A. Ferrari,  Davide Gabrielli "BBS invariant measures with independent soliton components", (3.12).
}
\usage{
phi(k,q)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{k}{
a binary vector, which describes the configuration. Maximum length is 30.
}
  \item{q}{
a numeric vector of parameters. Elements of vector must have values within interval [0,1) and sum of coordinates must be finite.
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
phi(c(1,1,0,0),c(0.4,0.3,0.2,0.1))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
