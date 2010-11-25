
tspc.train = function(data, data.test, type=c("survival", "regression"), s0.perc=0.5){
	library(superpc)
	proj.obj = tspc.project(data, data.test, type=type)
	data.train = proj.obj$data.train;
	fit.obj = superpc.train(data.train, type=type, s0.perc=s0.perc);
	return(train.obj = list(proj.obj = proj.obj, fit.obj = fit.obj))
}
