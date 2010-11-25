\name{tspc.cv}
\alias{tspc.cv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cross-validation
}
\description{
This function uses a form of cross-validation to estimate the optimal feature threshold in supervised principal components
}
\usage{
tspc.cv(fit, data, seed = 123, topfea = TRUE, n.topfea = 1000, n.threshold = 20, n.fold = NULL, folds = NULL, n.components = 1, min.features = 5, max.features = nrow(data$x[[1]]), type = "survival")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{
Object returned by tspc.train
}
  \item{data}{
Data object of form described in tspc.train documentation
}
  \item{seed}{A Numeric number
}
  \item{topfea}{
If it is TRUE, then the tuning paparmeter is the number of features
}
  \item{n.topfea}{
Maximum number of features used as the tuning parameter
}
  \item{n.threshold}{
Number of thresholds, when using the number of thresholds as a tuning parameter
}
  \item{n.fold}{
Number of cross-validation folds
}
  \item{folds}{
Lists of indices of cross-validation folds (optional)
}
  \item{n.components}{
Number of cross-validation components to use: 1,2 or 3.
}
  \item{min.features}{
Minimum number of features to include, in determining
range for threshold. Default 5.
}
  \item{max.features}{
Maximum number of features to include, in determining
range for threshold. Default is total number of features in the dataset. 
}
  \item{type}{
"survival" or "regression"
}
}
\details{
This function uses a form of cross-validation to estimate the optimal feature threshold.
}
\value{
 list(thresholds = thresholds, n.threshold = n.threshold, nonzero = nonzero, scor = scor, scor.lower = scor.lower, scor.upper = scor.upper, folds = folds, n.fold = n.fold, featurescores.folds = featurescores.folds, type = type)
\item{thresholds }{Vector of thresholds considered}
\item{n.threshold}{Number of thresholds}
\item{nonzero}{Number of features exceeding each value of the threshold}
\item{scor}{Full CV scores}
\item{scor.lower}{Full CV scores minus one standard error of scores}
\item{scor.upper}{Full CV scores plus one standard error of scores}
\item{folds}{Indices of CV folds used}
\item{n.fold}{Number of folds used in the cross-validation}
\item{featurescores.folds}{Feature scores for each fold}
\item{type }{problem type}
}
\author{
Yuping Zhang
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}% at least one, from doc/KEYWORDS
\keyword{survival}% __ONLY ONE__ keyword per line
