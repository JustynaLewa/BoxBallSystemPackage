\name{configuration.to.walk.plot}
\alias{configuration.to.walk.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{plot corresponding walk for a given ball configuration
%%  ~~function to do ... ~~
}
\description{
function creates a ggplot of a walk, which corresponds to a configuration given as an imput. The walk starts at (0,0). The ggplot2 package is required.
}
\usage{
configuration.to.walk.plot(k)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{k}{ a binary vector, which describes the ball configuration. Maximum length is 30.
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

configuration.to.walk.plot(c(1,1,0,1,0,0,1))
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
