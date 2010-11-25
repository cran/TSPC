tspc.project = function(data, data.test, type=c("survival", "regression")){
	X = as.matrix(data$x[[1]]);
	for(i in 2:length(data$x)){
		X = cbind(X, data$x[[i]])
	}
	if(type=="survival"){ par = t(apply(X, 1, cox.coef.func, data$y, data$censoring.status))}
	else{par = t(apply(X, 1, lm.coef.func, data$y))}
	
	wx.train = multiply.func(data$x, par)
	wx.test = multiply.func(data.test$x, par)
	
	if(type=="survival"){
		wdata.train = list(x = wx.train, y = data$y, censoring.status=data$censoring.status, featurenames=data$featurenames);
		wdata.test = list(x = wx.test, y= data.test$y, censoring.status=data.test$censoring.status, featurenames=data.test$featurenames);
	}else{
		wdata.train = list(x = wx.train, y = data$y, featurenames = data$featurenames);
		wdata.test = list(x = wx.test, y = data.test$y, featurenames = data.test$featurenames);
	}
	return(proj.obj = list(data.train = wdata.train, data.test = wdata.test))
}
