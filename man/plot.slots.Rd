\name{plot.slots}
\alias{plot.slots}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
plot matrix with slots
}
\description{
Funtion returns a visual of slots in a given configuration (see Pablo A. Ferrari: "BBS invariant measures with independent soliton components", Figure 5 ). Records slots are labeled with zeros. Libraries "plot.matrix" and "RColorBrewer" are required.
}
\usage{
plot.slots(k)
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
plot.slots(c(1,1,1,0,1,1,0,0,1,0,1,1,0,0,0,0))

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
