## ---- echo = FALSE------------------------------------------------------------
knitr::opts_chunk$set(
 fig.width  = 5 ,
 fig.height = 3.5,
 fig.align  = 'center'
)
oldpar <- list(mar = par()$mar, mfrow = par()$mfrow)

## -----------------------------------------------------------------------------
library("e1071")
library("classmap")

## -----------------------------------------------------------------------------
set.seed(1)
X = matrix(rnorm(200*2),ncol=2)
X[1:100,] = X[1:100,]+2
X[101:150,] = X[101:150,]-2
y = as.factor(c(rep("blue",150),rep("red",50))) # two classes
table(y)
dim(X) # 200 2
length(y) # 200
cols = c("deepskyblue3","red")
col = cols[as.numeric(y)]
plot(X,col=col,pch=19) 

## -----------------------------------------------------------------------------
dat = data.frame(X=X,y=y)
set.seed(1) 
svmfit = svm(y~.,data=dat,scale=F,kernel="radial",cost=10,
             gamma=1, probability=T)

## -----------------------------------------------------------------------------
names(svmfit)

plot(svmfit$decision.values,col=cols[as.numeric(y)]); abline(h=0)

plot(svmfit,data=dat,X.2~X.1, col = cols)

## -----------------------------------------------------------------------------
vcr.train = vcr.svm.train(X, y, svfit=svmfit)

## -----------------------------------------------------------------------------
names(vcr.train)
vcr.train$predint 
vcr.train$pred    
vcr.train$altint  
vcr.train$altlab  

## -----------------------------------------------------------------------------
vcr.train$PAC[1:3] 
summary(vcr.train$PAC)

## -----------------------------------------------------------------------------
vcr.train$fig[1:5,] 

## -----------------------------------------------------------------------------
vcr.train$farness[1:5]
summary(vcr.train$farness)

## -----------------------------------------------------------------------------
summary(vcr.train$ofarness)

## -----------------------------------------------------------------------------
confmat.vcr(vcr.train, cutoff = 0.975)

## -----------------------------------------------------------------------------
confmat.vcr(vcr.train)

## -----------------------------------------------------------------------------
confmat.vcr(vcr.train,showClassNumbers=T)

## -----------------------------------------------------------------------------
stackedplot(vcr.train, classCols = cols, separSize = 1.5,
            minSize=1, showLegend=T, cutoff = 0.975)

stplot = stackedplot(vcr.train, classCols = cols, 
                     separSize = 1.5, minSize = 1)
stackedplot(vcr.train,classCols=cols)

## -----------------------------------------------------------------------------
par(mfrow = c(1,1))
classmap(vcr.train,1,classCols=cols)
classmap(vcr.train,2,classCols=cols)
par(oldpar)

## -----------------------------------------------------------------------------
Xnew = X[c(1:25,151:175),]; dim(Xnew) # 50 200
ynew = y
ynew[c(21:25,171:175)] = NA
(ynew = ynew[c(1:25,151:175)]) 
length(ynew) # 50

## -----------------------------------------------------------------------------
vcr.test = vcr.svm.newdata(Xnew,ynew,vcr.train)

## -----------------------------------------------------------------------------
vcr.test$predint 
length(vcr.test$predint) 
vcr.test$altint 
length(vcr.test$altint) 
vcr.test$PAC 
length(vcr.test$PAC) 
vcr.test$farness 
length(vcr.test$farness) 
vcr.test$ofarness 
length(vcr.test$ofarness) 

## -----------------------------------------------------------------------------
plot(vcr.test$yintnew[c(1:20,26:45)],
     vcr.train$yint[c(1:20,151:170)]); abline(0,1)
plot(vcr.test$predint[c(1:20,26:45)],
     vcr.train$predint[c(1:20,151:170)]); abline(0,1)
plot(vcr.test$altint[c(1:20,26:45)],
     vcr.train$altint[c(1:20,151:170)]); abline(0,1)
plot(vcr.test$PAC[c(1:20,26:45)],
     vcr.train$PAC[c(1:20,151:170)]); abline(0,1)
plot(as.vector(vcr.test$fig[c(1:20,26:45),]),
     as.vector(vcr.train$fig[c(1:20,151:170),])); abline(0,1)
plot(vcr.test$farness[c(1:20,26:45)],
     vcr.train$farness[c(1:20,151:170)]); abline(0,1)
plot(vcr.test$ofarness[c(1:20,26:45)],
     vcr.train$ofarness[c(1:20,151:170)]); abline(0,1)

## -----------------------------------------------------------------------------
Kxx = makeKernel(X,svfit=svmfit)
dim(Kxx) 

## -----------------------------------------------------------------------------
outFV = makeFV(Kxx)
Xf = outFV$Xf 

## -----------------------------------------------------------------------------
max(abs(as.vector(Kxx - Xf%*%t(Xf)))) # tiny, this is good.

## -----------------------------------------------------------------------------
explvar = cumsum(apply(Xf,2,var))
plot(explvar/explvar[ncol(Xf)]); abline(h=0.95)
min(which(explvar > 0.95*explvar[ncol(Xf)])) 

## -----------------------------------------------------------------------------
range(rowSums(Xf^2)) 

## -----------------------------------------------------------------------------
pairs(Xf[,1:5], col=col)
plot(Xf[,1],Xf[,5],col=col,pch=19)

## -----------------------------------------------------------------------------
Xfdat = data.frame(X=Xf,y=y)
set.seed(1)
svmfit1 = svm(y~.,data=Xfdat,scale=F,kernel="linear",cost=10,
              probability=T)
plot(svmfit1$decision.values,svmfit$decision.values); abline(0,1)
vcr.trainf = vcr.svm.train(X=Xf, y, svfit=svmfit1)

## -----------------------------------------------------------------------------
plot(vcr.train$yint,vcr.trainf$yint); abline(0,1)
plot(vcr.train$predint,vcr.trainf$predint); abline(0,1)
plot(vcr.train$altint,vcr.trainf$altint); abline(0,1)
plot(vcr.train$PAC,vcr.trainf$PAC); abline(0,1)
plot(vcr.train$fig,vcr.trainf$fig); abline(0,1)
vcr.train$figparams$farscales
plot(vcr.train$farness,vcr.trainf$farness); abline(0,1)
plot(vcr.train$ofarness,vcr.trainf$ofarness); abline(0,1)

## -----------------------------------------------------------------------------
Xnew = Xf[c(1:25,151:175),]; dim(Xnew) 

## -----------------------------------------------------------------------------
vcr.testf = vcr.svm.newdata(Xnew,ynew,vcr.trainf)

## -----------------------------------------------------------------------------
plot(vcr.testf$yintnew,vcr.test$yintnew); abline(0,1)
plot(vcr.testf$predint,vcr.test$predint); abline(0,1)
plot(vcr.testf$altint,vcr.test$altint); abline(0,1)
plot(vcr.testf$PAC,vcr.test$PAC); abline(0,1)
plot(vcr.testf$fig,vcr.test$fig); abline(0,1)
plot(vcr.testf$farness,vcr.test$farness); abline(0,1)

## -----------------------------------------------------------------------------
data(data_bookReviews)
X = data_bookReviews[1:500, 1]; y = data_bookReviews[1:500, 2] 
Xnew  = data_bookReviews[501:1000, 1]; ynew = data_bookReviews[501:1000, 2] 
length(X); length(Xnew); length(y); length(ynew) 

## -----------------------------------------------------------------------------
cols = c("deepskyblue3","red")
col = cols[as.numeric(y)]
X[5]

## -----------------------------------------------------------------------------
strdot = kernlab::stringdot(type = "spectrum", length = 7)
strdot
ptm = proc.time()
Kr_train = kernlab::kernelMatrix(kernel=strdot, x=X)
(proc.time()-ptm)[3] 
dim(Kr_train) 

## -----------------------------------------------------------------------------
outFV = makeFV(Kr_train)
Xf = outFV$Xf 
dim(Xf) 

## -----------------------------------------------------------------------------
inprod = Xf%*%t(Xf) 
max(abs(as.vector(Kr_train-inprod))) 

## -----------------------------------------------------------------------------
range(rowSums(Xf^2)) # all 1

## -----------------------------------------------------------------------------
explvar = cumsum(apply(Xf,2,var))
plot(explvar/explvar[ncol(Xf)]); abline(h=0.95)
min(which(explvar > 0.95*explvar[ncol(Xf)]))

## -----------------------------------------------------------------------------
Xfdat = data.frame(X = Xf, y = y)
set.seed(1) # we have to fix the seed!
svmfit = svm(y~.,data=Xfdat,scale=F,kernel="linear",cost=2,
             probability=T)

## -----------------------------------------------------------------------------
plot(svmfit$decision.values) 

## -----------------------------------------------------------------------------
vcr.train = vcr.svm.train(Xf, y, svfit=svmfit)

## -----------------------------------------------------------------------------
confmat.vcr(vcr.train, showOutliers = F)

## -----------------------------------------------------------------------------
ptm = proc.time()
Kr_test = kernlab::kernelMatrix(kernel=strdot, x=Xnew, y=X)
(proc.time()-ptm)[3] 

## -----------------------------------------------------------------------------
FVtest = makeFV(kmat=Kr_test, transfmat=outFV$transfmat)
Zf = FVtest$Xf 
dim(Zf) 

## -----------------------------------------------------------------------------
inprod = Zf%*%t(Xf) 
max(abs(as.vector(Kr_test-inprod))) # tiny, good.

## -----------------------------------------------------------------------------
vcr.test = vcr.svm.newdata(Zf,ynew,vcr.train)

## -----------------------------------------------------------------------------
confmat.vcr(vcr.test)

## -----------------------------------------------------------------------------
stackedplot(vcr.test,classCols=cols,separSize=1.5,minSize=1.5)

## -----------------------------------------------------------------------------
# To identify the atypical points, uncomment the next line:
# classmap(vcr.test,1,classCols=cols,identify=T)
# Identified point(s): 
# [1]  67 127 118 145  94  45
#
# pdf("Bookreviews_500_classmap_test_negative.pdf")
par(mar = c(3.5, 3.5, 2.0, 0.2))
coords = classmap(vcr.test, 1, classCols = cols, cex = 1.5,
                  main = "predictions of the negative book reviews",
                  cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5,
                  maxfactor = 1.05, identify = FALSE)
indstomark = c(67, 127, 118, 145, 94, 45) # from identify=T above
labs  = c(  "a",  "b", "x",  "c",   "d",  "e")
xvals = coords[indstomark, 1] +
   c(-0.08, 0, 0, 0.09, 0.09, 0.09) # visual fine tuning
yvals = coords[indstomark, 2] +
   c(0, 0.033, 0.030, 0, 0, 0)
text(x = xvals, y = yvals, labels = labs, cex = 1.5)
legend("topright", fill = cols,
       legend = c("negative", "positive"), cex = 1.4,
       ncol = 1, bg = "white")
# dev.off()
par(oldpar)

## -----------------------------------------------------------------------------
(Xnew[which(ynew == 1)])[67] # (a) in the paper.

(Xnew[which(ynew == 1)])[127] # (b) in the paper.

(Xnew[which(ynew == 1)])[118] # (x) 
# Negative review, but highly recommends another book!

(Xnew[which(ynew == 1)])[145] # (c) in the paper.
# Contains some grammatical errors
# Is negative, but expressed very indirectly.

(Xnew[which(ynew == 1)])[94] # (d) in the paper. 
# Definitely negative, but quite short.

(Xnew[which(ynew == 1)])[45] # (e) in the paper.

## -----------------------------------------------------------------------------
#
# To identify the atypical points:
# classmap(vcr.test,2,classCols=cols,identify=T)
# Identified point(s): 
# [1] 48  32 18 168 113
#
# pdf("Bookreviews_500_classmap_test_positive.pdf")
par(mar = c(3.5, 3.5, 2.0, 0.2))
coords = classmap(vcr.test, 2, classCols = cols, cex = 1.5,
                  main = "predictions of the positive book reviews",
                  cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5,
                  maxfactor = 1.05, identify = FALSE)
indstomark = c(48, 32, 18, 168, 113) # from identify=T above
coords[indstomark,]
labs  = c( "f",   "g",  "h",  "i",  "j")
xvals = coords[indstomark, 1] +
   c(0, 0, 0.09, 0.09, 0.09) # visual fine tuning
yvals = coords[indstomark, 2] +
   c(0.032, 0.032, 0, 0, 0)
text(x = xvals, y = yvals, labels = labs, cex = 1.5)
legend("left", fill = cols,
       legend = c("negative", "positive"), cex = 1.4,
       ncol = 1, bg = "white")
# dev.off()
par(oldpar)

## -----------------------------------------------------------------------------
(Xnew[which(ynew == 2)])[48] # (f) in the paper.

(Xnew[which(ynew == 2)])[32] # (g) in the paper.

(Xnew[which(ynew == 2)])[18] # (h) in the paper.

(Xnew[which(ynew == 2)])[168] # (i) in the paper.

(Xnew[which(ynew == 2)])[113] # (j) in the paper.
# very short (and spelling error).

