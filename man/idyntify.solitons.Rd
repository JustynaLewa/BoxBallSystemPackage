\name{idyntify.solitons}
\alias{idyntify.solitons}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{function to idyntify solitions in a given configuration
}
\description{
idyntify.solitons is a function which return a two dimensional list which groups the cooridates of boxes that are part of solitons in a given configuration.

length(solitions) is equel to the length of a maximum soltion in a configuration.

length(solitions[[i]]) corresponds to number of solitions of a length i.

solitions[[i]] a list, which contains all coordinates of solitions of length i. At j-th place, the list contains the cooridantes of j-th solition of length i.
}
\usage{
idyntify.solitons(k)
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
idyntify.solitons(c(1,1,1,1,0,0,1,0,1,1,1,0,1,0,0,0,0,1,1,0,1,1,0,0,0,0))
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
