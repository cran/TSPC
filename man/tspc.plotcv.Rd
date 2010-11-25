\name{tspc.plotcv}
\alias{tspc.plotcv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot output from tspc.cv
}
\description{
Plots pre-validation results from plotcv, to aid in choosing best threshold
}
\usage{
tspc.plotcv(object, call.win.metafile=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
	Object returned by tspc.cv 
}
\item{call.win.metafile}{
		call.win.metafile	 Ignore, for use by Excel program
	}
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

fit = tspc.train(data, data.test, type="survival")

cv.obj = tspc.cv(fit$fit.obj, data, type="survival", n.fold=2)
tspc.plotcv(cv.obj)
}% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}% at least one, from doc/KEYWORDS
\keyword{survival}% __ONLY ONE__ keyword per line
