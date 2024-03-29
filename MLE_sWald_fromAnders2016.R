###########################################################
## Section 1: Directions / Installation 
###########################################################

# 1. Install any missing libraries by using the following commands:

install.packages("SuppDists"); install.packages("plyr"); install.packages("Rcpp")

# 2. Source all the code of Section 2 
#    (copy and paste all of Section 2 into your R session)

# 3. Navigate to Section 3 for example applications / walk-throughs

###########################################################
## 2. Section 2: Code to Sourcea
###########################################################

## To simplify the code, we load these libraries in the R session:
library(SuppDists); library(plyr)  

## swfit() fits one Shifted Wald to a vector of observations 
swfit <- function(x,type=1,outlier=FALSE,dtype=1,search=FALSE,ql=.01,qu=.99,outl=10,outu=10,criterion=1){
  # For thorough performance, searches fit performance across a number of ranges
  if(search==TRUE){
    quantl <- list(
      # Insert here to add additional options (quantiles, distance types) to search
      # quants is the range of quantiles, dtype = 1 for L1-norm, 2 for L2-norm distance  
      list(quants = seq(.01,.99,len=100), dtype = 1),
      list(quants = seq(.01,.98,len=100), dtype = 1),
      list(quants = seq(.01,.97,len=100), dtype = 1),
      list(quants = seq(.02,.99,len=100), dtype = 1),
      list(quants = seq(.02,.98,len=100), dtype = 1),
      list(quants = seq(.03,.97,len=100), dtype = 1),
      list(quants = seq(.001,.999,len=100), dtype = 1),
      list(quants = seq(.001,.99,len=100), dtype = 1),
      list(quants = seq(.01,.999,len=100), dtype = 1)
    )}; 
  # For speedier performance, calculates a fit on only one given range
  if(search==FALSE){quantl <- list(list(quants = seq(ql,qu,len=100), dtype = 1))}
  
  calcparam <- function(beta,x){
    #calcparam(), given a beta value, calculates the other SW parameters
    # Compute a_0 and an initial estimate of theta
    a_0 <- ((mean(x) - min(x))^3 / mean((x-mean(x))^2))^.5
    thetaest <- function(x,a,beta) { n <- length(x);
    tmpf <- function(i,beta,n) (1-pinvGauss(i,nu=beta,1))^n
    min(x) - (a^2)*integrate(tmpf, lower=0, upper=Inf, beta,n)$value   }
    t_init <- thetaest(x,a_0,beta)
    # Then compute an initial estimate of alpha
    alphaest <- function(x,t) { (mean((x-t)^-1) - (mean(x)-t)^-1)^-.5 }
    a_init <- alphaest(x,t_init)
    # Compute a final estimate of theta
    t_final <- thetaest(x,a_init,beta)
    # Compute a final estimate of alpha
    a_final <- alphaest(x,t_final)
    # Compute a final estimate of gamma
    g_final <- 1/(beta*a_final)  
    return(c(gamma=g_final, alpha=a_final, theta=t_final,mu=beta*(a_final^2), lambda=a_final^2,tau=t_final,beta=beta))
  }
  optfunc <- function(beta,x,qnts,dtype=1){
    # optfunc() is the function that is minimized by the deviance criterion
    # Given a beta, calcparam is called, then the deviance criterion is calculated
    # returns the quantile deviation statistic. 
    vars <- as.list(calcparam(beta=beta,x=x))
    if(dtype==1){return(sum(abs((qinvGauss(qnts,vars$mu,vars$lambda)+vars$tau)-x)))}
    if(dtype==2){return(sum(((qinvGauss(qnts,vars$mu,vars$lambda)+vars$tau)-x)^2))}
  }
  tmp <- sol <- NULL;
  for(i in 1:length(quantl)){
    qnts <- quantl[[i]]$quants; dtype <- quantl[[i]]$dtype; 
    xq <- quantile(x, qnts,type=7) 
    sol <- optimize(optfunc, x=xq, qnts=qnts, dtype=dtype, interval=c(.0001,10000), maximum=F)  
    tmp <- rbind(tmp,c(calcparam(beta=sol$minimum,x=xq),res=sol$objective))
  }; vars <- data.frame(tmp)
  
  # Criterions that will decide which fit (from quantl) to return.
  quantl2 <- list(
    # Insert here to add additional criterions, access by e.g. criterion=2 in swapply() 
    seq(.05,.95,len=91)
  )
  xqnts2 <- quantile(x, quantl2[[criterion]],type=7);
  pxqnts2 <- data.frame(t(sapply(1:dim(vars)[1],function(x) qinvGauss(quantl2[[criterion]],vars$mu[x],vars$lambda[x])+vars$tau[x]))); 
  return(unlist(vars[order(rank(
    rank(apply(abs(sweep(pxqnts2,2,xqnts2)[,1:5]),1,sum))+  # which fits best the first percentiles (5th-10th)
      rank(apply(abs(sweep(pxqnts2,2,xqnts2)[,1:91]),1,sum)) # which fits most consistently over all percentiles (5th-95th)
    ,ties.method="random"))[1],])) # return these parameters
}

## swapply() fits a Shifted Wald to a full data set (e.g. data frame), by each factor specified
## Two factors minimum, e.g. "condition" and "subject" for one factor, use swfit() with tapply()
swapply <- function(dat,facs,obsvar = "RT",exclude=NA,type="1",doplot=FALSE,ql=.01,qu=.99,search=TRUE,criterion=1,outlier=NA,outl=3,outu=6){
  if(!is.na(outlier)){
    dat2 <- dat; dat <- NULL
    for(i in unique(dat2[[outlier]])){dat <- rbind(dat,dat2[which(dat2[[outlier]]==i & dat2[[obsvar]] >= median(dat2[[obsvar]][which(dat2[[outlier]]==i)])-(outl*mad(dat2[[obsvar]][which(dat2[[outlier]]==i)])) &  dat2[[obsvar]] <= median(dat2[[obsvar]][which(dat2[[outlier]]==i)])+(outu*mad(dat2[[obsvar]][which(dat2[[outlier]]==i)]))),])}
  }
  qnts <- seq(.1,.9,len=9); set.seed(1)
  vars <- with(dat, aggregate(eval(parse(text=obsvar)),
                              by = lapply((lapply(facs,function(x) parse(text=x))),function(x) as.factor(eval(x))),
                              FUN = function(x) return(swfit(x,type=type,search=search,ql=ql,qu=qu,criterion=1)))); 
  
  dimnames(vars)[[2]][1:length(facs)] <- facs
  vars <- data.frame(vars[,length(facs)+1],vars[,1:length(facs)]);
  xqnts <- with(dat, aggregate(eval(parse(text=obsvar)),
                               by = lapply((lapply(facs,function(x) parse(text=x))),function(x) as.factor(eval(x))), 
                               FUN = function(x) quantile(x, qnts,type=7))); dimnames(xqnts)[[2]][1:length(facs)] <- facs; xqnts <- data.frame(xqnts[,length(facs)+1],xqnts[,1:length(facs)]); dimnames(xqnts)[[2]][1:length(qnts)] <- round(qnts,3); xqnts <- xqnts[,1:length(qnts)]
  pxqnts <- data.frame(t(sapply(1:dim(vars)[1],function(x) qinvGauss(qnts,vars$mu[x],vars$lambda[x])+vars$tau[x]
  )));  dimnames(pxqnts)[[2]][1:length(qnts)] <- round(qnts,3)
  
  if(doplot==TRUE){
    swplot(vars=vars,dat=dat,obsvar=obsvar,exclude=exclude) 
  }
  return(list(vars=vars,xqnts=xqnts,pxqnts=pxqnts,dat=dat))
}

## swplot() plots the parameters and model fit checks, works for within or between subject designs
# example usage: swplot(vars,dat,exclude="subject");  
# (exclude="subject" excludes plotting individual subject parameters)
swplot <- function(vars,dat,obsvar,exclude=NA,xqnts=NA,pxqnts=NA){
  facs <- dimnames(vars)[[2]][9:dim(vars)[2]]
  fplot <- facs[!facs%in%exclude]
  qnts <- seq(.1,.9,.1)
  if(is.na(xqnts) || is.na(pxqnts)){
    xqnts <- with(dat, aggregate(eval(parse(text=obsvar)),
                                 by = lapply((lapply(facs,function(x) parse(text=x))),function(x) as.factor(eval(x))),
                                 FUN = function(x) quantile(x, qnts,type=7))); dimnames(xqnts)[[2]][1:length(facs)] <- facs; xqnts <- data.frame(xqnts[,length(facs)+1],xqnts[,1:length(facs)]); dimnames(xqnts)[[2]][1:length(qnts)] <- round(qnts,3); xqnts <- xqnts[,1:length(qnts)]
                                 pxqnts <- data.frame(t(sapply(1:dim(vars)[1],function(x) qinvGauss(qnts,vars$mu[x],vars$lambda[x])+vars$tau[x]
                                 )));  dimnames(pxqnts)[[2]][1:length(qnts)] <- round(qnts,3)
  }
  nf <- layout(matrix(c(1,1,1,4,
                        2,2,2,5,
                        3,3,3,6), 3, 4, byrow=TRUE), respect=F)
  par(oma = c(0, 0, 3, 0),mar=c(4,4,0,1),mgp=c(2.75,.75,0)); param <- c("gamma","alpha","theta")
  ylabs <- c(expression(paste("Acc. Rate   ",gamma)),expression(paste("Threshold   ",alpha)),expression(paste("TEA   ",theta)))
  # Parameter Values Plot
  for(k in 1:3){ me <- NULL; se <- NULL; ln <- c(0); for(i in 1:length(fplot)){
    me <- c(me,with(vars, tapply(eval(parse(text=param[k])), lapply((lapply(fplot[i],function(x) parse(text=x))),function(x) as.factor(eval(x))), mean)))
    se <- c(se,with(vars, tapply(eval(parse(text=param[k])), lapply((lapply(fplot[i],function(x) parse(text=x))),function(x) as.factor(eval(x))), stderr)))    
    ln <- c(ln,length(me))     }; xlocs <- 1:length(names(me)); Z <- 1
    plot(xlocs,me,type="p",xlab="",xaxt="n",pch=21,ylim=c(min(me-(se*Z)),max(me+(se*Z))),bg="black",las=1,ylab=ylabs[k])
    for(i in 1:length(fplot)){lines(x=xlocs[(ln[i]+1):ln[i+1]],me[(ln[i]+1):ln[i+1]])}
    segments(as.numeric(xlocs[1:length(me)]),me-(se*Z), as.numeric(xlocs[1:length(me)]), me+(se*Z))
    arrows(as.numeric(xlocs[1:length(me)]),me-(se*Z), as.numeric(xlocs[1:length(me)]), me+(se*Z),code=3,angle=90,length=.025)
    axis(1,at=as.numeric(xlocs[1:length(me)]),labels=names(me))
    title(xlab="Factor Levels") }
  
  swchecks(vars=vars,dat=dat,method=1,xqnts=xqnts,pxqnts=pxqnts)
}

## swplotws() plots the parameters and model fit checks for within-subject designs
## the standard error bars are corrected for within subject error 
# example usage: swplotws(vars,exclude="subject")
# (exclude="subject" provides corrected error bars based on within-subject variation)
swplotws <- function(vars,exclude=NA,xticklab=NULL,frame=TRUE){
  capwords <- function(s, strict = FALSE) {cap <- function(s) paste(toupper(substring(s, 1, 1)),{s <- substring(s, 2); if(strict) tolower(s) else s},sep = "", collapse = " " );sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))}
  if(frame==TRUE){nf <- layout(matrix(c(1,1,1,4,
                                        2,2,2,5,
                                        3,3,3,6), 3, 4, byrow=TRUE), respect=F)
  par(oma = c(0, 0, 3, 0),mar=c(4,4,0,1),mgp=c(2.75,.75,0))}
  param <- c("gamma","alpha","theta")
  facs <- dimnames(vars)[[2]][9:dim(vars)[2]]
  fplot <- facs[!facs%in%exclude]
  if(is.null(xticklab)){xticklab <- unlist(sapply(which(names(vars)%in%fplot), function(ind) 1:length(levels(vars[,ind])) ))}
  
  ylabs <- c(expression(paste("Acc. Rate   ",gamma)),expression(paste("Threshold   ",alpha)),expression(paste("TEA   ",theta)))
  len <- length(fplot); lfacs <- length(facs)
  for(k in 1:3){
    me <- NULL; se <- se2 <- lv <- lvmid <- meraw <- NULL; ln <- c(0)
    for(i in which(!facs%in%exclude)){              
      sm <- ddply(vars, facs[c(i,lfacs)], .fun = function(x) mean(x[[param[k]]])); names(sm)[3] <- "value"
      me <- c(me,ddply(sm, facs[i], .fun = function(x) mean(x[["value"]]))[,2])
      meraw <- c(meraw,with(vars, tapply(eval(parse(text=param[k])), lapply((lapply(facs[i],function(x) parse(text=x))),function(x) as.factor(eval(x))), mean)))
      se <- c(se,ddply(sm,facs[i], .fun = function(x) sqrt(var(x[["value"]],na.rm=TRUE)/length(na.omit(x[["value"]]))))[,2])
      if(length(unique(vars[[facs[i]]]))<3){
        tmp <- apply(t(do.call(cbind, lapply(split(sm, sm[[facs[lfacs]]]), function(x) x$value))),1,diff)
        se2 <- c(se2,sqrt(var(tmp,na.rm=TRUE)/length(na.omit(tmp)))); rm(tmp)
      }else{
        se2 <- c(se2,apply(t(apply(t(do.call(cbind, lapply(split(sm, sm[[facs[lfacs]]]), function(x) x$value))),1,diff)),2,function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))))
      }
      Z <- 1.96; Z <- 1;ln <- c(ln,length(me));lv <- rbind(lv,data.frame(lvl=1:length(levels(vars[[facs[i]]])),fac=rep(i,times=length(levels(vars[[facs[i]]]))) )); #lv <- c(lv,1:length(levels(vars[[facs[i]]])))
    }
    x <- barplot(me,space=as.numeric(1:dim(lv)[1] %in% which(lv[["lvl"]]==1)[-1]),add=FALSE,ylab=ylabs[k],ylim=c(min(me-(se*Z)),max(me+(se*Z))),xpd=F,axes=F,xaxt='n')
    axis(1,at=x,labels=xticklab,cex.axis=.925); axis(2,las=1); box(); inds <- which(lv[,"lvl"] %% 2 == 1) #inds <- which(lv %% 2 == 1)
    lv <- data.frame(lv,ind =1:dim(lv)[1],x=x); 
    lv <- data.frame(lv,ind2=unlist(with(lv, tapply(ind,fac,function(y) y!=max(y)))))
    inds <- lv$x[lv$ind2]+.5
    mid <- 0.5; me <- apply(rbind(me[which(lv$ind2)],me[which(lv$ind2)+1]),2,min) # to plot se bars between
    arrows(lwd=2,x0=x[lv$ind2]+mid,y0=me-se2, x1=x[lv$ind2]+mid, y1=me+se2,code=3,angle=90,length=.025)
    #if(btwsub==1){
    #  arrows(x0=x[1:2],y0=meraw[1:2]-se[1:2], x1=x[1:2], y1=meraw[1:2]+se[1:2],code=3,angle=90,length=.025)
    #}
    mtext(capwords(fplot),cex=.6,line = 2, adj = .5,side=1,at=ddply(lv,.(fac),summarize,inds= median(x))[,"inds"])
  }
}

## swchecks() provides model fit checks (without plotting parameters)
# example usage: plot the parameters with swplotws(), then swchecks(vars,dat,obsvar="rt")
swchecks <- function(vars,dat,obsvar="rt",xqnts=NA,pxqnts=NA,method=1,qqplot=3,xliml=7,frceaxis=FALSE,doqq=TRUE,dodist=TRUE,ylimlqq=0,ylimuqq=0,sim,xlimudec=0,ylimudec=0,ylimlcel=0,ylimucel=0,frame=FALSE){
  if(frame==TRUE){nf <- layout(matrix(c(1,4,4,4,
                                        2,5,5,5,
                                        3,6,6,6), 3, 4, byrow=TRUE), respect=F)
  par(oma = c(0, 0, 3, 0),mar=c(4,4,0,1),mgp=c(2.75,.75,0))}
  qnts <- seq(.1,.9,len=9);  ndat <- dim(vars)[1];
  facs <- dimnames(vars)[[2]][9:dim(vars)[2]]
  if(is.na(xqnts) || is.na(pxqnts)){
    xqnts <- with(dat, aggregate(eval(parse(text=obsvar)),
                                 by = lapply((lapply(facs,function(x) parse(text=x))),function(x) as.factor(eval(x))), 
                                 FUN = function(x) quantile(x, qnts,type=7))); dimnames(xqnts)[[2]][1:length(facs)] <- facs; xqnts <- data.frame(xqnts[,length(facs)+1],xqnts[,1:length(facs)]); dimnames(xqnts)[[2]][1:length(qnts)] <- round(qnts,3); xqnts <- xqnts[,1:length(qnts)]
                                 pxqnts <- data.frame(t(sapply(1:dim(vars)[1],function(x) qinvGauss(qnts,vars$mu[x],vars$lambda[x])+vars$tau[x]
                                 )));  dimnames(pxqnts)[[2]][1:length(qnts)] <- round(qnts,3)
  }
  ## Standardized Residual
  if(method==1){
    xse <- pxse <- sapply(1:ndat, function(ind)  sqrt(vars$alpha[ind]/(vars$gamma[ind]^3)))
    resid <- abs((xqnts/xse)-(pxqnts/pxse))
  }
  ## Non-standardized Residual
  if(method==2){
    xse <- pxse <- sapply(1:ndat, function(ind)  sqrt(vars$alpha[ind]/(vars$gamma[ind]^3)))
    resid <- abs(xqnts-pxqnts)
  }  
  ##### BEGIN QQ-Plot
  if(frceaxis==TRUE){ plot(NA,xlim=c(min(ylimlqq),max(ylimuqq)),ylim=c(min(ylimlqq),max(ylimuqq)),las=1,ylab="Reaction Time (Fit)",xlab="",main="")
  }else{plot(NA,xlim=c(max(min(xqnts,pxqnts),ylimlqq),max(xqnts,pxqnts,ylimuqq)),ylim=c(max(min(xqnts,pxqnts),ylimlqq),max(xqnts,pxqnts,ylimuqq)),las=1,ylab="Reaction Time (Fit)",xlab="",main="")
  }
  if(qqplot==1 || qqplot==3){points(as.matrix(xqnts),as.matrix(pxqnts),col="light grey")}
  if(qqplot==2 || qqplot==3){mxqnts <- apply(xqnts,2,mean)
  mpxqnts <- apply(pxqnts,2,mean)
  residraw <- abs(xqnts-pxqnts); stdev <- apply(residraw,2,sd)
  arrows(mxqnts,mpxqnts-stdev,mxqnts,mpxqnts+stdev,angle=90,length=.025,code=3)
  points(mxqnts,mpxqnts,pch=19,bg="dark grey",col="dark grey")
  }
  box()
  par(new=TRUE);plot(function(x) x=x,las=1,ylab="",xlab="",axes="F")
  title(sub="Reaction Time (data)",line=2)
  ##### END QQ-Plot  
  ##### BEGIN Decile Plot 
  dens <- lapply(1:9,function(i) density(resid[,i]))
  xmax <- max(sapply(1:9,function(i) quantile(resid[,i],p=.5)  ))
  ymax <- max(sapply(1:9,function(i) max(density(resid[,i])$y)  ))
  plot(NA,ylim=c(0,max(ymax,ylimudec)),xlim=c(0,max(xmax,xlimudec)),las=1,main="",xlab="",ylab="Density")
  invisible(sapply(1:9,function(i){
    points(dens[[i]],type="l")
    text(dens[[i]]$x[which(dens[[i]]$y==max(dens[[i]]$y))], .9*max(dens[[i]]$y),labels=paste("",i,sep="")) 
  })); title(sub="Std. Residual (per Decile)",line=2)
  ##### End Decile Plot 
  ##### BEGIN Cell Plot 
  cresid <- apply(abs(resid),1,sum)
  plot(NA, xlab="",ylab=expression(paste("Std. Residual")),las=1,ylim=c(min(0,cresid,ylimlcel),max(cresid*1.25,ylimucel)),xlim=c(0,ndat),col=c("dark grey"),xaxt='n')
  segments(x0=1:length(cresid),x1=1:length(cresid),y0=rep(0,times=length(cresid)),y1=cresid,col="dark grey")
  axis(1); box()
  abline(b=0,a=quantile(cresid,c(.05,.95))[2],col="black",lty=2)
  abline(b=0,a=quantile(cresid,c(.05,.95))[1],col="black",lty=2)
  text(ndat/2,quantile(cresid,c(.01,.978))[2],labels="95%")
  text(ndat/2,quantile(cresid,c(.001,.978))[1],labels="5%")
  points(x=ndat/2,y=mean(cresid),col="black",bg="black",pch=19)
  arrows(ndat/2,mean(cresid)-stderr(cresid),ndat/2,mean(cresid)+stderr(cresid),angle=90,length=.1,code=3)
  legend("topright",legend=sapply(c(bquote(bar(Delta) == .(round(mean(cresid),2))),
                                    bquote(bar(sigma)[X] == .(round(mean(xse),0))), 
                                    bquote(rho[Delta][sigma] == .(round(cor(cresid,xse),2)))
  ),as.expression) ,bty="n" ,xjust=1)
  title(sub="Cell",line=2) 
  ##### END Cell Plot  
} 

# stderr() is a supplementary function that calculates the standard error
stderr <- function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))

###########################################################
## End of Section 2 
###########################################################

setwd("add repository path here")

d <- read.csv(file="RP4.csv") 
# d$loc <- NULL
d$suj <- as.numeric(factor(d$suj))
d = d[complete.cases(d),]
unique(d$suj)
d = d[d$soa == "-100" | d$soa == "200",]
# d$soa[d$soa==-550] = 1
d$soa[d$soa==-100] = 1
d$soa[d$soa==200] = 2
# d$soa[d$soa==500] = 4


size <- NA
## compute mean number of trials and sd
for (i in 1:20){

  size[1 + 4*(i-1)] <- length(d$rt[d$suj == i & d$soa == 1 & d$cg == 1]) 
  size[2 + 4*(i-1)] <- length(d$rt[d$suj == i & d$soa == 2 & d$cg == 1]) 
  size[3 + 4*(i-1)] <- length(d$rt[d$suj == i & d$soa == 1 & d$cg == 0]) 
  size[4 + 4*(i-1)] <- length(d$rt[d$suj == i & d$soa == 2 & d$cg == 0]) 

}
meansize <- mean(size)
sdsize <- sd(size)



## 2. Fit the model to the data, by searching a variety of fitting ranges 

f <- swapply(dat=d,search=TRUE, obsvar="rt", #name of RT variable
             facs=c("soa","cg","suj"))  #design cell combos to fit

## 4. Plot the variables with within-subject error corrected bars 

swplotws(f$vars,c("suj")) 

## 5. Plot the model fit diagnostics

swchecks(f$vars,d,xqnts=f$xqnts,pxqnts=f$pxqnts)

## 6. Use str(f) to look deeper into the results

f$vars   # gives the fitted parameters
f$data   # gives the observed data
f$xqnts  # gives the 9 deciles of the data
f$pxqnts # gives the 9 deciles by the model fit
library("ez")
datANoVA <- data.frame(1:80)
datANoVA$gamma <- f$vars$gamma
datANoVA$alpha <- f$vars$alpha
datANoVA$theta <- f$vars$theta
datANoVA$suj <- f$vars$suj
datANoVA$soa <- f$vars$soa
datANoVA$cg <- f$vars$cg

# write.csv(datANoVA, file = "sWald_output.csv")

ezANOVA(data = datANoVA,dv = gamma, wid = suj,within =.(soa,cg))#, within_covariates = loc)
ezANOVA(data = datANoVA,dv = alpha, wid = suj,within =.(soa,cg))#, within_covariates = loc)
ezANOVA(data = datANoVA,dv = theta, wid = suj,within =.(soa,cg))#, within_covariates = loc)



