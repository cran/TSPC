tspc.cv <- function (fit, data, seed=123, topfea=TRUE, n.topfea=1000, n.threshold = 20, n.fold= NULL, folds= NULL, n.components = 1, min.features = 5, max.features = nrow(data$x[[1]]), type="survival") 
{
	set.seed(seed)
    if (n.components > 5) {
		cat("Max # of components is 5", fill = TRUE);
	}
	n.components <- min(5, n.components)
	#n.scor.dim <- n.components + n.components -1  
	n.scor.dim = n.components
	if(is.null(max.features)){
		max.features = nrow(data$x[[1]])
	}
	if (is.null(n.fold) & is.null(folds)) {
		n.fold = 2;
	}
	n = length(data$y);
	cur.tt <- as.numeric(fit$feature.scores);
	if(topfea==TRUE){
		max.features = min(n.topfea, max.features)
		n.threshold = max.features-min.features+1
		cur.tt.sort = sort.int(abs(cur.tt), decreasing=TRUE, index.return=TRUE)
		thresholds = cur.tt.sort$x[min.features:max.features]
		lower = cur.tt.sort$x[max.features]
		upper = cur.tt.sort$x[min.features]
	}else{
		lower <- quantile(abs(cur.tt), 1 - (max.features/nrow(data$x[[1]])));
		upper <- quantile(abs(cur.tt), 1 - (min.features/nrow(data$x[[1]])));
		thresholds <- seq(from = lower, to = upper, length = n.threshold)
	}
	
	if (!is.null(folds)) {
		n.fold = length(folds)
	}
	if (is.null(folds)) {
		if (sum(data$censoring.status == 0) > 0) {
			folds = balanced.folds(data$censoring.status,nfolds = n.fold)
			n.fold = length(folds)
		}else {
			folds <- vector("list", n.fold)
			breaks <- round(seq(from = 1, to = (n + 1), length = (n.fold + 1)))
			cv.order <- sample(1:n)
			for (j in 1:n.fold) {
				folds[[j]] <- cv.order[(breaks[j]):(breaks[j + 1] - 1)]
			}
		}
	}
	featurescores.folds <- matrix(NA, nrow = nrow(data$x[[1]]), ncol = n.fold)
	
	nonzero <- rep(0, n.threshold)
	scor = array(NA, c(n.scor.dim, n.threshold, n.fold)) 
    
	scor.preval <- matrix(NA, nrow = n.components, ncol = n.threshold)
	scor.lower = NULL
	scor.upper = NULL
	v.preval <- array(NA, c(n, n.components, n.threshold))
	
	first = 1
	last = n.fold
	for (j in first:last) {
	  if(type == "survival"){
		par.train = matrix(nrow=nrow(data$x[[1]]), ncol=length(data$x));
		y.train.tmp = data$y[-folds[[j]]];
		y.test.tmp = data$y[folds[[j]]];
		censoring.status.tmp = data$censoring.status[-folds[[j]]];
		censoring.status.test.tmp = data$censoring.status[folds[[j]]];
		
		x.train.tmp =  matrix(unlist(data$x[[1]][, -folds[[j]]]), nrow=nrow(data$x[[1]]));
		
		xlist.train = list();
		xlist.train[[1]] = data$x[[1]][, -folds[[j]]];
		xlist.test = list();
		xlist.test[[1]] = data$x[[1]][, folds[[j]]];
		for(kkk in 2:length(data$x)){
			xtmp.kkk = matrix(unlist(data$x[[kkk]][, -folds[[j]]]), nrow=nrow(data$x[[1]]));
			x.train.tmp = cbind(x.train.tmp, xtmp.kkk);
			xlist.train[[kkk]] = data$x[[kkk]][, -folds[[j]]];
			xlist.test[[kkk]] = data$x[[kkk]][, folds[[j]]];
		}
		
		
		par.train = t(apply(x.train.tmp, 1, cox.coef.func, y.train.tmp, censoring.status.tmp)); 
		
		wx.train.tmp = multiply.func(xlist.train,  par.train);
		
		wx.test.tmp = multiply.func(xlist.test,  par.train);
		
		data.temp = list(x = wx.train.tmp, y = y.train.tmp, censoring.status = censoring.status.tmp);
		data.test.temp = list(x = wx.test.tmp, y = y.test.tmp, censoring.status = censoring.status.test.tmp);
		
		train.obj.temp = superpc.train(data.temp, type = "survival", s0.perc = fit$s0.perc);
		cur.tt <- train.obj.temp$feature.scores;
		featurescores.folds[, j] <- cur.tt;
	  }else{
		  par.train = matrix(nrow=nrow(data$x[[1]]), ncol=length(data$x));
		  y.train.tmp = data$y[-folds[[j]]];
		  y.test.tmp = data$y[folds[[j]]];
		  x.train.tmp = matrix(unlist(data$x[[1]][, -folds[[j]]]), nrow=nrow(data$x[[1]]));
		  xlist.train = list();
		  xlist.train[[1]] = data$x[[1]][, -folds[[j]]];
		  xlist.test = list();
		  xlist.test[[1]] = data$x[[1]][, folds[[j]]];
		  for(kkk in 2:length(data$x)){
			  xtmp.kkk = matrix(unlist(data$x[[kkk]][, -folds[[j]]]), nrow=nrow(data$x[[1]]));
			  x.train.tmp = cbind(x.train.tmp, xtmp.kkk);
			  xlist.train[[kkk]] = data$x[[kkk]][, -folds[[j]]];
			  xlist.test[[kkk]] = data$x[[kkk]][, folds[[j]]];
		  }
		  par.train = t(apply(x.train.tmp, 1, lm.coef.func, y.train.tmp));
		  wx.train.tmp = multiply.func(xlist.train, par.train);
		  wx.test.tmp = multiply.func(xlist.test, par.train);
		  data.temp = list(x = wx.train.tmp, y = y.train.tmp);
		  data.test.temp = list(x = wx.test.tmp, y=y.test.tmp);
		  train.obj.temp = superpc.train(data.temp, type="regression", s0.perc = fit$s0.perc);
		  cur.tt = train.obj.temp$feature.scores;
		  featurescores.folds[, j] = cur.tt;
		}
		for (i in 1:n.threshold) {
			cur.features <- (abs(cur.tt) > thresholds[i]);
			if (sum(cur.features) > 1) {
				nonzero[i] <- nonzero[i] + sum(cur.features)/n.fold
				cur.svd = mysvd(data.temp$x[cur.features, ], n.components = n.components); 
				xtemp = data.test.temp$x[cur.features, , drop = FALSE];
				xtemp <- t(scale(t(xtemp), center = cur.svd$feature.means, 
								 scale = FALSE))
				cur.v.all <- scale(t(xtemp) %*% cur.svd$u, center = FALSE, 
								   scale = cur.svd$d)
				n.components.eff <- min(sum(cur.features), n.components)
				cur.v <- matrix(cur.v.all[, 1:n.components.eff],  ncol=n.components.eff)
				v.preval[folds[[j]], 1:n.components.eff, i] <- cur.v
				for (k in 1:ncol(cur.v)) {
						if (type == "survival") {
							require(survival)
				#			junk0 = coxph(Surv(data.temp$y, data.temp$censoring.status) ~ cur.svd$v[, k], control = coxph.control(iter.max = 10))$coefficients ;
				#			junk01 = cur.v[, k] * junk0
				#			junk02 =  summary(coxph(Surv(data.test.temp$y, data.test.temp$censoring.status) ~ junk01))$waldtest[1]         
				#			scor[k, i, j] <-junk02;
				#			if(k > 1){
								junk1 = coxph(Surv(data.temp$y, data.temp$censoring.status) ~ cur.svd$v[, 1:k], control = coxph.control(iter.max = 10))$coefficients;
								junk2 = cur.v[, 1:k] %*% matrix(junk1, ncol=1)
								junk3 = summary(coxph(Surv(data.test.temp$y, data.test.temp$censoring.status) ~ junk2))$waldtest[1]
				#				scor[(k+n.components-1), i, j] <-  junk3;
								scor[k, i, j] <-  junk3;
				#			}              
						}
						else {
					#		junk0 = lm(data.temp$y ~ cur.svd$v[,k])$coefficients;
					#		junk01 = cur.v[, k] * junk0;
					#		junk02 = summary(lm(data.temp$y ~ junk01))$fstat[1];
					#		scor[k, i, j] = junk02;
					#		if(k > 1){
								junk1 = lm(data.temp$y ~ cur.svd$v[, 1:k])$coefficients;
								junk2 = cur.v[, 1:k] %*% matrix(junk1, ncol=1);
								junk3 = summary(lm(data.test.temp$y ~ junk2));
								scor[k, i, j] <- junk3$fstat[1]
					#		}
						}
					}
			}
		}
    }
    lscor = apply(log(scor), c(1, 2), mean.na)
    se.lscor = apply(log(scor), c(1, 2), se.na)
    scor.lower = exp(lscor - se.lscor)
    scor.upper = exp(lscor + se.lscor)
    scor <- exp(lscor)
   	
    junk <- list(thresholds = thresholds, n.threshold = n.threshold, 
				 nonzero = nonzero, scor = scor, scor.lower = scor.lower, scor.upper = scor.upper, folds = folds, 
				 n.fold = n.fold, featurescores.folds = featurescores.folds, type = type);
    class(junk) <- "tspc.cv";
    return(junk);
}

