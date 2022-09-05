\name{p_x_m}
\alias{p_x_m}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 function to determine the number of maximal solitions
}
\description{
probability phi that number of maximal solitons is equel to x, given that maximal soliton size is equel to M.
}
\usage{
p_x_m(x,M,q)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
positive integer describing number of maximal solitons
}
 \item{M}{
non negative integer, M<=length(q)
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
p_x_m(2,3,c(0.4,0.3,0.2,0.1))

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
