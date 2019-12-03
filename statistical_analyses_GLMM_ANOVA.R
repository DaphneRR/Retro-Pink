# conduct glmer analyses on retro-pink data

# d prime cannot be the object of a formal glmer analysis, however we can try looking into
# RTs using a Normal or inverse.gaussian link function. The IES can also be explored this way.

# load data

setwd("D:/RetroPink/")
library("lme4")
d = read.csv("RP4.csv",na.strings=NaN)
d = d[complete.cases(d),]
unique(d$suj)
# in case we need to rescale
# d$soa[d$soa == -600] = 1
# d$soa[d$soa == -150] = 2
# d$soa[d$soa == 150] = 3
# d$soa[d$soa == 450] = 4
d$soa = as.factor(d$soa)
d$suj = as.factor(d$suj)
d$congr = as.factor(d$congr)

# perform glmm using inverse.gaussian link. Notice the interaction is significant
fit = glmer(rt ~ soa*congr + (1|suj), data = d, family =inverse.gaussian(link = "identity"))
r.squaredGLMM(fit)

drop1(fit,test="Chisq") # here we have a small convergence issue.
fit2 = glmer(rt ~ soa + congr + (1|suj), data = d, family =inverse.gaussian(link = "identity"),glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
drop1(fit2,test="Chisq")
# here, we have a small convergence issue. 


## diagnostic plot
require(ggplot2)
ggplot()+geom_point(aes(fitted(fit),y=residuals(fit, type = "pearson")),pch=21)+stat_smooth(aes(fitted(fit),y=residuals(fit, type = "pearson")),method = "loess",se = T, n=5,col="red")+theme_bw()+labs(x="fitted", y = "residuals")


# the reference for this methodology is Lo and Andrews 2015
# we have an interaction, we'll perform pairwise comparisons

# here we use the vignette for Tukey's tests on interactions in a glm framework as a tutorial
# to set up our contrasts and perform our planned comparisons
# link:
# https://cran.r-project.org/web/packages/multcomp/vignettes/multcomp-examples.pdf

library (multcomp)
summary(glht(fit, linfct = mcp(soa <= "Tukey")))

tmp <- expand.grid(congr = unique(d$congr),soa = unique(d$soa))
X = model.matrix(~soa*congr,data = tmp)
glht(fit,linfct = X)
Tukey = contrMat(table(d$congr), "Tukey")
K1 = cbind(Tukey,matrix(0,nrow = nrow(Tukey),ncol = ncol(Tukey)),matrix(0,nrow = nrow(Tukey),ncol = ncol(Tukey)),matrix(0,nrow = nrow(Tukey),ncol = ncol(Tukey)))
rownames(K1) = paste(levels(d$soa)[1],rownames(K1), sep = ":")
K2 = cbind(matrix(0,nrow = nrow(Tukey),ncol = ncol(Tukey)),Tukey,matrix(0,nrow = nrow(Tukey),ncol = ncol(Tukey)),matrix(0,nrow = nrow(Tukey),ncol = ncol(Tukey)))
rownames(K2) = paste(levels(d$soa)[2],rownames(K2), sep = ":")
K3 = cbind(matrix(0,nrow = nrow(Tukey),ncol = ncol(Tukey)),matrix(0,nrow = nrow(Tukey),ncol = ncol(Tukey)),Tukey,matrix(0,nrow = nrow(Tukey),ncol = ncol(Tukey)))
rownames(K3) = paste(levels(d$soa)[3],rownames(K3), sep = ":")
K4 = cbind(matrix(0,nrow = nrow(Tukey),ncol = ncol(Tukey)),matrix(0,nrow = nrow(Tukey),ncol = ncol(Tukey)),matrix(0,nrow = nrow(Tukey),ncol = ncol(Tukey)),Tukey)
rownames(K4) = paste(levels(d$soa)[4],rownames(K4), sep = ":")
K = rbind(K1,K2,K3,K4)
colnames(K) = c(colnames(Tukey),colnames(Tukey),colnames(Tukey),colnames(Tukey))
summary(glht(fit,linfct = K %*% X, test = adjusted("holm")))
confint(glht(fit,linfct = K %*% X, test = adjusted("holm")))
# this says there is an effect of congruency for the -150 and +150ms SOAs



# We can repeat the procedure for the IES. The distribution does not change shape from a simple division so we'll use the inverse Gaussian parent.
# reload data and calculate average accuracy by condition
d = read.csv("RP4.csv",na.strings=NaN)
Sujs = unique(d$suj)
for(i in Sujs){
  d$meanAcc[d$soa == -600 & d$congr == 0 & d$suj == i] <- mean(d$acc[!is.na(d$soa) & d$soa == -600 & !is.na(d$congr) & d$congr == 0 & d$suj == i]) 
  d$meanAcc[d$soa == -600 & d$congr == 1 & d$suj == i] <- mean(d$acc[!is.na(d$soa) & d$soa == -600 & !is.na(d$congr) & d$congr == 1 & d$suj == i]) 
  d$meanAcc[d$soa == -150 & d$congr == 0 & d$suj == i] <- mean(d$acc[!is.na(d$soa) & d$soa == -150 & !is.na(d$congr) & d$congr == 0 & d$suj == i]) 
  d$meanAcc[d$soa == -150 & d$congr == 1 & d$suj == i] <- mean(d$acc[!is.na(d$soa) & d$soa == -150 & !is.na(d$congr) & d$congr == 1 & d$suj == i]) 
  d$meanAcc[d$soa == +150 & d$congr == 0 & d$suj == i] <- mean(d$acc[!is.na(d$soa) & d$soa == +150 & !is.na(d$congr) & d$congr == 0 & d$suj == i]) 
  d$meanAcc[d$soa == +150 & d$congr == 1 & d$suj == i] <- mean(d$acc[!is.na(d$soa) & d$soa == +150 & !is.na(d$congr) & d$congr == 1 & d$suj == i]) 
  d$meanAcc[d$soa == +450 & d$congr == 0 & d$suj == i] <- mean(d$acc[!is.na(d$soa) & d$soa == +450 & !is.na(d$congr) & d$congr == 0 & d$suj == i]) 
  d$meanAcc[d$soa == +450 & d$congr == 1 & d$suj == i] <- mean(d$acc[!is.na(d$soa) & d$soa == +450 & !is.na(d$congr) & d$congr == 1 & d$suj == i]) 
}
# d$meanAcc = d$meanAcc*100
d$ies <- d$rt/d$meanAcc
d = d[complete.cases(d),]
d$soa[d$soa == -600] = 1
d$soa[d$soa == -150] = 2
d$soa[d$soa == 150] = 3
d$soa[d$soa == 450] = 4
d$soa = as.factor(d$soa)
d$suj = as.factor(d$suj)
d$congr = as.factor(d$congr)

# perform the glmm on IES. We add more iterations to avoid a failure to converge
fit = glmer(ies ~ soa*congr + (1|suj), data = d, family =inverse.gaussian(link = "identity"),glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))

## these lines help determine which factors are significant
drop1(fit,test="Chisq")
fit2 = glmer(ies ~ soa + congr + (1|suj), data = d, family =inverse.gaussian(link = "identity"),glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
drop1(fit2,test="Chisq")
# we find an effect of all fixed factors.
r.squaredGLMM(fit)



## diagnostic plot
require(ggplot2)
ggplot()+geom_point(aes(fitted(fit),y=residuals(fit, type = "pearson")),pch=21)+stat_smooth(aes(fitted(fit),y=residuals(fit, type = "pearson")),method = "loess",se = T, n=5,col="red")+theme_bw()+labs(x="fitted", y = "residuals")


# the reference for this methodology is Lo and Andrews 2015
# now that we do have an interaction, we'll perform pairwise comparisons

# here we use the vignette for Tukey's tests on interactions in a glm framework as a tutorial
# to set up our contrasts and perform our planned comparisons
# link:
# https://cran.r-project.org/web/packages/multcomp/vignettes/multcomp-examples.pdf

library (multcomp)
summary(glht(fit, linfct = mcp(soa <= "Tukey")))

tmp <- expand.grid(congr = unique(d$congr),soa = unique(d$soa))
X = model.matrix(~soa*congr,data = tmp)
glht(fit,linfct = X)
Tukey = contrMat(table(d$congr), "Tukey")
K1 = cbind(Tukey,matrix(0,nrow = nrow(Tukey),ncol = ncol(Tukey)),matrix(0,nrow = nrow(Tukey),ncol = ncol(Tukey)),matrix(0,nrow = nrow(Tukey),ncol = ncol(Tukey)))
rownames(K1) = paste(levels(d$soa)[1],rownames(K1), sep = ":")
K2 = cbind(matrix(0,nrow = nrow(Tukey),ncol = ncol(Tukey)),Tukey,matrix(0,nrow = nrow(Tukey),ncol = ncol(Tukey)),matrix(0,nrow = nrow(Tukey),ncol = ncol(Tukey)))
rownames(K2) = paste(levels(d$soa)[2],rownames(K2), sep = ":")
K3 = cbind(matrix(0,nrow = nrow(Tukey),ncol = ncol(Tukey)),matrix(0,nrow = nrow(Tukey),ncol = ncol(Tukey)),Tukey,matrix(0,nrow = nrow(Tukey),ncol = ncol(Tukey)))
rownames(K3) = paste(levels(d$soa)[3],rownames(K3), sep = ":")
K4 = cbind(matrix(0,nrow = nrow(Tukey),ncol = ncol(Tukey)),matrix(0,nrow = nrow(Tukey),ncol = ncol(Tukey)),matrix(0,nrow = nrow(Tukey),ncol = ncol(Tukey)),Tukey)
rownames(K4) = paste(levels(d$soa)[4],rownames(K4), sep = ":")
K = rbind(K1,K2,K3,K4)
colnames(K) = c(colnames(Tukey),colnames(Tukey),colnames(Tukey),colnames(Tukey))
summary(glht(fit,linfct = K %*% X, test = adjusted("holm")))
confint(glht(fit,linfct = K %*% X, test = adjusted("holm")))
# this says there is an effect of congruency for the -150 and +150ms SOAs


# load d prime data
dprime = read.csv("dprime_data.csv",sep=";")
dprime$suj = as.factor(dprime$suj)
dprime$cg = as.factor(dprime$cg)
dprime$soa = as.factor(dprime$soa)


ezANOVA(data = dprime,dv = val, wid = suj,within =.(soa,cg))#, within_covariates = loc)

