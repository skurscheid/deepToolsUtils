\name{DownloadEnsemblGTF}
\alias{DownloadEnsemblGTF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Downloads GTF file from Ensembl FTP server.
}
\description{
The function takes Ensembl version, organism and a temporary directory/path as input and attempts to download
the corresponding GTF file using curl.
}
\usage{
DownloadEnsemblGTF(version = "current",
                   url = "ftp.ensembl.org/pub/current_gtf",
                   temp_dir = "/tmp",
                   organism = "homo_sapiens")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
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
function (x)
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
