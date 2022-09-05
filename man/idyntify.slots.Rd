\name{idyntify.slots}
\alias{idyntify.slots}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
function to idyntify slots in a configuration
}
\description{
for a given configuration, function returns a matrix M x (n+2), where M is the maximum size of a soliton in configuration, and n is the length of a configuration (two columns added to represent records). M-th row in output matrix represents all one-slots in configuration (starting with record), wheras first row represents M-slots.
}
\usage{
idyntify.slots(k)
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
idyntify.slots(c(1,1,1,1,0,0,1,0,1,1,1,0,1,0,0,0,0,1,1,0,1,1,0,0,0,0))

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
