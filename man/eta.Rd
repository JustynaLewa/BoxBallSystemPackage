\name{eta}
\alias{eta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{plots a given configuartion eta after n evolutions
}
\description{
For a given configuartion eta, the functions plots configuration after n evolutions. The plot is a set of balls, black balls correspond to a full box in configuration, whereas white ball, means that box is empty.
}
\usage{
eta(k,n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{k}{a binary vector, which describes the configuration. Maximum length is 30.
}
\item {n} {a non-negative integer, which descrbes the number of evolutions. The default value is zero. For n=0, the function plots given configuration without any tranformations.}
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
eta(c(1,0,0,1,1,0,0),1)
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
