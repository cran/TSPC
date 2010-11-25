tspc.plotcv= function (object, call.win.metafile = FALSE)
{

scor = object$scor
k = nrow(scor)
if(object$type=="survival"){
  ymax = max(object$scor.upper[!is.na(object$scor.upper)])
}

if(object$type=="regression"){
  # df of denom for f is average sample size in validation fold

n.mean=0
for(i in 1:object$n.fold){
   n.mean=n.mean+length(object$folds[[i]])/object$n.fold
}
 
 denom.df=n.mean-1-nrow(scor)
 ymax = max(object$scor.upper[!is.na(object$scor.upper)], qf(0.95, nrow(scor), denom.df))

}

if(call.win.metafile){win.metafile()}

 if(object$type=="survival"){ ylab="Likelihood ratio test statistic"}
if(object$type=="regression"){ ylab="F statistic"}  ## yuping added

#ylab="Likelihood ratio test statistic"  ## yuping deleted

   # matplot(object$th, t(scor), xlab = "Threshold", ylab = ylab, ylim = c(0, ymax), lty=rep(1,k))
   # matlines(object$th, t(scor), lty=rep(1,k), ...)
     matplot(object$th, t(scor), type="b", xlab = "Threshold", ylab = ylab, ylim = c(0, ymax), lty=rep(1,k))

if(object$type=="survival"){  abline(h = qchisq(0.95, 1), lty = 2, col = 1)} ## yuping added

if(object$type=="regression"){  ## yuping added
         # df of denom for f is average sample size in validation fold
        abline(h = qf(0.95, 1, denom.df), lty = 2, col = 1)  ## yuping added
}  ## yuping added
  

for (j in 1:k) {
      delta=((-1)^j)*diff(object$th)[1]/4
      error.bars(object$th+delta*(j>1),t(object$scor.upper[j,]), t(object$scor.lower[j,]), lty=2, col=j)
   
}

if(call.win.metafile){dev.off()}
return(TRUE)

}

error.bars <-function(x, upper, lower, width = 0.005, ...) {
  xlim <- range(x)
  barw <- diff(xlim) * width
  segments(x, upper, x, lower, ...)
  segments(x - barw, upper, x + barw, upper, ...)
  segments(x - barw, lower, x + barw, lower, ...)
  range(upper, lower)
}
jitter<-
function(x)
{
        x + 0.03 * abs(x) * sign(rnorm(length(x)))
}

