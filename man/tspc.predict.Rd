\name{tspc.predict}
\alias{tspc.predict}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Form principal components predictor from a trained tspc object
}
\description{
Computes supervised principal components, using scores from "object"
}
\usage{
tspc.predict(object, data, newdata, threshold, n.components = 3, prediction.type = c("continuous", "discrete", "nonzero"), n.class = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
Object fit.obj returned by tspc.train
}
  \item{data}{
List of projection of training data returned by tspc.train, object proj.obj$wdata.train
}
  \item{newdata}{
List of projection of test data returned by tspc.train, object proj.obj$wdata.test
}
  \item{threshold}{
Threshold for scores.
}
  \item{n.components}{
Number of principal components to compute.
Should be 1,2 or 3.
}
  \item{prediction.type}{
"continuous" for raw principal component(s);
"discrete" for principal component categorized in equal bins;
"nonzero" for indices of features that pass the threshold 
}
  \item{n.class}{
Number of classes into which predictor is binned
(for prediction.type="discrete"
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
list(v.pred = out, u = x.sml.svd$u, d = x.sml.svd$d, 
        which.features = which.features, v.pred.1df = v.pred.1df, 
        n.components = n.pc, coef = result$coef, call = this.call, 
        prediction.type = prediction.type)
\item{v.pred }{Supervised principal componients predictor}
\item{u }{U matrix from svd of weighted feature matrix}
  \item{d }{singual values  from svd of weighted feature matrix}
 \item{which.features}{Indices of features exceeding threshold}
 \item{n.components}{Number of supervised  principal components requested}
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

predict.obj<- tspc.predict(fit$fit.obj, fit$proj.obj$data.train, fit$proj.obj$data.test, threshold=1.0, n.components=1)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}% at least one, from doc/KEYWORDS
\keyword{survival}% __ONLY ONE__ keyword per line

