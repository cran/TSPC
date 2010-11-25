\name{tspc.train}
\alias{tspc.train}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prediction using time-course gene expression
}
\description{
Does prediction of a quantitative regression or survival outcome, using the time-course gene expression.
}
\usage{
tspc.train(data, data.test, type = c("survival", "regression"), s0.perc = 0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
Data object with components x- a list of p by n matrix of features,
one observation per column, one matrix per time point; y- n-vector of outcome measurements;
censoring.status- n-vector of censoring censoring.status (1= died or event occurred, 0=survived, or event was censored), needed
for a censored survival outcome.
}
  \item{data.test}{
Data object with components x- a list of p by n matrix of features,
one observation per column, one matrix per time point; y- n-vector of outcome measurements;
censoring.status- n-vector of censoring censoring.status (1= died or event occurred, 0=survived, or event was censored), needed
for a censored survival outcome.
}
  \item{type}{
Problem type:  "survival" for censored survival outcome, or "regression"
for simple quantitative outcome.
}
  \item{s0.perc}{
Factor for denominator of score statistic, between 0 and 1:  the percentile
of standard deviation values added to the  denominator. Default is 0.5 (the median)
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
	\item{proj.obj}{projection of training data and test data}
	\item{fit.obj}{fitted object using training data}
}
\references{
Zhang Y, Tibshirani RJ, Davis RW. Predicting patient survival from longitudinal gene expression. Stat Appl Genet Mol Biol. 2010;9(1):Article41. Epub 2010 Nov 22.
}
\author{
Yuping Zhang
}


\examples{
x = list()
for(i in 1:2){
	set.seed(i+123)
	x[[i]] = matrix(rnorm(500*100), ncol=100)
}
y = sample(c(5:100), size=100, replace=TRUE)
censoring = sample(c(0,1), size=100, replace=TRUE)

data = list(x = x, y=y, censoring.status=censoring, genenames = as.character(paste("gene", c(1:500), sep="")), geneid = 	as.character(paste("gene", c(1:500), sep="")))

x = list()
for(i in 1:2){
	set.seed(i+133)
	x[[i]] = matrix(rnorm(500*100), ncol=100)
}
y = sample(c(5:100), size=100, replace=TRUE)
censoring = sample(c(0,1), size=100, replace=TRUE)

data.test = list(x = x, y=y, censoring.status=censoring, genenames = as.character(paste("gene", c(1:500), sep="")), geneid = as.character(paste("gene", c(1:500), sep="")))

obj = tspc.train(data, data.test, type="survival")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }
\keyword{ survival }% __ONLY ONE__ keyword per line
