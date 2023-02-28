library(glmnet)

######Model builing and Testing ########
modellingLasso <- read.csv("Train_Data.csv", header = T, stringsAsFactors = T)
modellingLasso$Response

lassoModel <- glmnet(
  x=data.matrix(modellingLasso[,-1]),
  y=modellingLasso$Response,
  standardize=TRUE,
  alpha=0.5,
  family="binomial")
plot(lassoModel, xvar="lambda")

cv.lassoModel <- cv.glmnet(
  x=data.matrix(modellingLasso[,-1]),
  y=modellingLasso$Response,
  standardize=TRUE,
  alpha=0.5,
  nfolds=10,
  family="binomial",
  parallel=TRUE)


plot(cv.lassoModel)

idealLambda <- cv.lassoModel$lambda.min
idealLambda1se <- cv.lassoModel$lambda.1se
print(idealLambda); print(idealLambda1se)

# derive coefficients for each gene
co <- coef(cv.lassoModel, s=idealLambda, exact=TRUE)
co

co.se <- coef(cv.lassoModel, s=idealLambda1se, exact=TRUE)
co.se

#identify predictors for each sub-type
x <- rownames(co.se)[which(co.se != 0)]
rownamesco.se$RD[whichco.se$RD != 0]


which(co.se$pCR)
rownames(co.se)

finalLasso <- glm(modellingLasso$Response ~ 
                   
                   	  ,
                  data=modellingLasso,
                  family=binomial(link="logit"), maxit = 100)
summary(finalLasso)



require(pROC)

roc <- roc(modellingLasso$Response, fitted(finalLasso), smooth=FALSE)
ci(roc)
plot.roc(
  roc,
  grid=FALSE,
  grid.lwd=2,
  col="red",
  main="")
text(0.3, 0.45,
     col="red",
     paste("AUC (AUC):",
           paste(round(ci.auc(roc)[1:3], 3), sep=" ", collapse=","), sep=" "),
     cex=0.7)

#####New Prediction ######
new <- read.csv('new_data.csv', header = T, row.names = 1)
new <-data.matrix(new)
nrow(new)
results<-predict(lassoModel, newx=new, s = idealLambda)
results
