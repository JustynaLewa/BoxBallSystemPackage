\name{idyntify.runs}
\alias{idyntify.runs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{function to identify runs in a given configuration
%%  ~~function to do ... ~~
}
\description{
For a given configuration eta the function the runs of eta as maximal blocks of successive boxes where eta has a constant value. Function returns a numeric vector of a length of a given configuration. The value of i-th coordinate represents the length of a run that i-th box in a configuration belongs to.
}
\usage{
idyntify.runs(k)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{k}{
a binary vector, which describes the configuration. Maximum length is 30.
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

idyntify.runs(c(1,1,1,1,0,0,1,0,1,1,1,0,1,0,0,0,0,1,1,0,1,1,0,0,0,0))

#returns vector 4 4 4 4 2 2 1 1 3 3 3 1 1 4 4 4 4 2 2 1 2 2 4 4 4 4
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
