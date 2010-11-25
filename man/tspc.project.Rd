\name{tspc.project}
\alias{tspc.project}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Project time-course gene expression to weighted gene expression
}
\description{
Project time-course gene expression to weighted gene expression
}
\usage{
tspc.project(data, data.test, type = c("survival", "regression"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
List of training data, of form described in tspc.train documentation
}
  \item{data.test}{
List of test data, of form described in tspc.test documentation
}
  \item{type}{
Problem type:  "survival" for censored survival outcome, or "regression"
for simple quantitative outcome.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
list(data.train = wdata.train, data.test = wdata.test)
%%  If it is a LIST, use
  \item{data.train }{Projection of training data}
  \item{data.test }{Projection of test data}
%% ...
}

\author{
Yuping Zhang
}

\keyword{regression}% at least one, from doc/KEYWORDS
\keyword{survival}% __ONLY ONE__ keyword per line
