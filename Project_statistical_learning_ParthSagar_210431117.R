##########################################################################
#                      PARTH SAGAR                                       #        
#                       210431117                                        #
#               Project in statistical learning in r                     #
##########################################################################

#install.packages("caTools")
#install.packages("tidyverse")


library(tidyverse)
library(mlbench)
library(caTools)
library(leaps)
#loading the data
data("BreastCancer")
#checking the size
dim(BreastCancer)

#printing the first 6 rows of the data
head(BreastCancer)
#removing the rows with null values
BreastCancer=BreastCancer[complete.cases(BreastCancer),]
#checking population distribution in each category
table(BreastCancer$Class)
#numerical exploratory data analysis
summary(BreastCancer)
#graphical plot of exploratory data analysis
pairs(BreastCancer[,2:10])

#picking predictor variables and converting them as numeric
(X1=BreastCancer[,2:10])
X1$Cl.thickness=as.numeric(X1$Cl.thickness)
X1$Cell.size=as.numeric(X1$Cell.size)
X1$Cell.shape=as.numeric(X1$Cell.shape)
X1$Marg.adhesion=as.numeric(X1$Marg.adhesion)
X1$Epith.c.size=as.numeric(X1$Epith.c.size)
X1$Bare.nuclei=as.numeric(X1$Bare.nuclei)
X1$Bl.cromatin=as.numeric(X1$Bl.cromatin)
X1$Normal.nucleoli=as.numeric(X1$Normal.nucleoli)
X1$Mitoses=as.numeric(X1$Mitoses)
#scaling the variable
(X1= scale(X1))
#Pick out response variable
y=BreastCancer[,11]
y = as.integer(y)
y=y-1
#combining to make a new data frame
breast_cancer=data.frame(X1,y)

#storing the value of n and p
n=nrow(breast_cancer)
p=ncol(breast_cancer)-1


#doing logistic regression
#logreg_fit=glm(y~., data = breast_cancer, family = "binomial")
#summary(logreg_fit)

## Set the seed to make the analysis reproducible
set.seed(1)
#installing the library for best subset selection
library(bestglm)
#Apply best subset selection
bss_fit_AIC= bestglm(breast_cancer, family=binomial, IC="AIC")
bss_fit_BIC= bestglm(breast_cancer, family=binomial, IC="BIC")
#examine the results
bss_fit_AIC$Subsets
bss_fit_BIC$Subsets

#identify best-fitting models
(best_AIC=bss_fit_AIC$ModelReport$Bestk)
(best_BIC=bss_fit_BIC$ModelReport$Bestk)
bss_fit_AIC$Subsets
bss_fit_BIC$Subsets

coef(bss_fit_AIC,5)

#create multi panel plotting device
par(mfrow=c(1,2))
#Produce plots,highlighting optimal value of K
plot(0:p, bss_fit_AIC$Subsets$AIC, xlab = "Number of predictors", ylab = "AIC", type="b")
points(best_AIC, bss_fit_AIC$Subsets$AIC[best_AIC+1], col="red", pch=16)
plot(0:p, bss_fit_BIC$Subsets$BIC, xlab = "Number of predictors", ylab = "BIC", type="b")
points(best_BIC, bss_fit_BIC$Subsets$BIC[best_BIC+1], col="red", pch=16)




########################################## BIC #########################################
set.seed(10000)
pstar=5
bss_fit_BIC$Subsets[pstar+1,]
(indices = as.logical(bss_fit_BIC$Subsets[pstar+1, 2:(p+1)]))
BreastCancer_dataBIC = data.frame(X1[,indices], y)


## Obtain regression coefficients for this model
logreg1_fit = glm(y ~ ., data=BreastCancer_dataBIC, family="binomial")
summary(logreg1_fit)

#predicting values through logistic regression
predict_logisticreg= predict(logreg1_fit, BreastCancer_dataBIC, type = "response")
#formatting the output variables into categories
yhat_log=ifelse(predict_logisticreg > 0.5, 1, 0)
#obtaining confusion matrix
(confusion_matrix_log=table(Observed=BreastCancer_dataBIC$y,Predicted=yhat_log))
#computing the error
(log_tst_err = (1-mean(BreastCancer_dataBIC$y==yhat_log)))

################################
## 10-fold cross validation
nfolds = 10
## Sample fold-assignment index
fold_index = sample(nfolds, n, replace=TRUE)
## Print first few fold-assignments
head(fold_index)

#cross validation function
log_cv = function(X1, y, fold_ind) {
  Xy = data.frame(X1, y=y)
  nfolds = max(fold_ind)
  if(!all.equal(sort(unique(fold_ind)), 1:nfolds)) stop("Invalid fold partition.")
  cv_errors = numeric(nfolds)
  for(fold in 1:nfolds) {
    #applying logistic regression
    tmp_fit = glm(y ~ ., data=Xy[fold_ind!=fold,],family = "binomial")
    #predicting the response variable
    yhat_log = predict(tmp_fit, Xy[fold_ind==fold,], type="response")
    #categorizing the data 
    yhat = ifelse(yhat_log > 0.5, 1, 0)
    yobs = y[fold_ind==fold]
    #calculating the test error for BIC
    cv_errors[fold] = 1-mean((yobs == yhat))
  }
  fold_sizes = numeric(nfolds)
  for(fold in 1:nfolds) fold_sizes[fold] = length(which(fold_ind==fold))
  test_error = weighted.mean(cv_errors, w=fold_sizes)
  return(test_error)
}
(log_cv_error = log_cv(BreastCancer_dataBIC[,1:5], BreastCancer_dataBIC$y, fold_index))

## cross validation error=0.033


################################## AIC ####################################
set.seed(1234)

pstar=7
bss_fit_AIC$Subsets[pstar+1,]
(indices = as.logical(bss_fit_AIC$Subsets[pstar+1, 2:(p+1)]))
BreastCancer_dataAIC = data.frame(X1[,indices], y)

## Obtain regression coefficients for this model
logreg2_fit = glm(y ~ ., data=BreastCancer_dataAIC, family="binomial")
summary(logreg2_fit)


predict_logisticreg= predict(logreg1_fit, BreastCancer_dataAIC, type = "response")
yhat_log=ifelse(predict_logisticreg > 0.5, 1, 0)
(confusion_matrix_log=table(Observed=BreastCancer_dataAIC$y,Predicted=yhat_log))
(log_tst_err = (1-mean(BreastCancer_dataAIC$y==yhat_log)))

################################
## 10-fold cross validation
nfolds = 10
## Sample fold-assignment index
fold_index = sample(nfolds, n, replace=TRUE)
## Print first few fold-assignments
head(fold_index)

#cross validation function
log_cv = function(X1, y, fold_ind) {
  Xy = data.frame(X1, y=y)
  nfolds = max(fold_ind)
  if(!all.equal(sort(unique(fold_ind)), 1:nfolds)) stop("Invalid fold partition.")
  cv_errors = numeric(nfolds)
  for(fold in 1:nfolds) {
    #applying logistic regression
    tmp_fit = glm(y ~ ., data=Xy[fold_ind!=fold,],family = "binomial")
    #predicting the response variable
    yhat_log = predict(tmp_fit, Xy[fold_ind==fold,], type="response")
    #categorizing the data 
    yhat = ifelse(yhat_log > 0.5, 1, 0)
    yobs = y[fold_ind==fold]
    #calculating the test error for BIC
    cv_errors[fold] = 1-mean((yobs == yhat))
  }
  fold_sizes = numeric(nfolds)
  for(fold in 1:nfolds) fold_sizes[fold] = length(which(fold_ind==fold))
  test_error = weighted.mean(cv_errors, w=fold_sizes)
  return(test_error)
}
(log_cv_error = log_cv(BreastCancer_dataAIC[,1:7], BreastCancer_dataBIC$y, fold_index))

##cross validation error==0.03074671

















#################################### ridge #############################################


#load the glmnet package
library(glmnet)
set.seed(0111)
#choose grid of values for the tuning parameter
grid=10^seq(-4,2, length.out=100)

#fit a model with ridge penalty for each value of the tuning parameter
ridge_fit = glmnet(X1, y, family="binomial", alpha=0, standardize=FALSE, lambda=grid)

#Examine the effect of the tuning parameter on the parameter estimates
plot(ridge_fit, xvar="lambda",col=rainbow(p), label=TRUE)

##summary of ridge
(summ_ridge=summary(ridge_fit))

#k fold cross validation for optimal value of lambda
ridge_cv_fit= cv.glmnet(X1, y, family="binomial", alpha=0, standardize=FALSE, lambda=grid, foldid = fold_index)

#value of minimum lambda
(lamb_min = ridge_cv_fit$lambda.min)

#coefficient value at optimum value of tuning parameter
coef(ridge_fit, s=lamb_min)

#displaying all coefficients
coef(ridge_cv_fit)

##test error
ridge_cv_fit$cvm[which(ridge_cv_fit$lambda == lamb_min)]









########################### LASSO ################################



#load the glmnet package
library(glmnet)
set.seed(0111)
#choose grid of values for the tuning parameter
grid=10^seq(-4,2, length.out=100)

#fit a model with ridge penalty for each value of the tuning parameter
LASSO_fit = glmnet(X1, y, family="binomial", alpha=1, standardize=FALSE, lambda=grid)

#Examine the effect of the tuning parameter on the parameter estimates
plot(LASSO_fit, xvar="lambda",col=rainbow(p), label=TRUE)

##summary of ridge
summary(LASSO_fit)

#k fold cross validation for optimal value of lambda
LASSO_cv_fit= cv.glmnet(X1, y, family="binomial", alpha=1, standardize=FALSE, lambda=grid, foldid = fold_index)

#value of minimum lambda
(lamb_min = LASSO_cv_fit$lambda.min)

#coefficient value at optimum value of tuning parameter
coef(LASSO_fit, s=lamb_min)

#displaying all coefficients
coef(LASSO_cv_fit)

##test error
LASSO_cv_fit$cvm[which(LASSO_cv_fit$lambda == lamb_min)]





################################## LDA ############################################


# Load MASS and nclSLR package
library(MASS)
library(nclSLR)
linDA(variables= as.matrix(breast_cancer[-10]),  group=breast_cancer$y)

## 10-fold cross validation
nfolds = 10
## Sample fold-assignment index
fold_index = sample(nfolds, n, replace=TRUE)
## Print first few fold-assignments
head(fold_index)

# define the function for cross validation.
reg_cv_lda = function(X1, y, fold_ind) {
  Xy = data.frame(X1, y=y)
  nfolds = max(fold_ind)
  if(!all.equal(sort(unique(fold_ind)), 1:nfolds)) stop("Invalid fold partition.")
  cv_errors = numeric(nfolds)
  for(fold in 1:nfolds) {
    
    # Fitting LDA
    tmp_fit_lda = lda(Xy[fold_ind!=fold,]$y ~ ., data=Xy[fold_ind!=fold,][,-10])
    # Predicting the response variable
    phat_lda = predict(tmp_fit_lda, Xy[fold_ind==fold,][-10])
    # convert prediction into numeric from probability
    yhat_lda = phat_lda$class
    yobs = y[fold_ind==fold]
    # Calculating the test error and storing it
    cv_errors[fold] = 1 - mean((yobs == yhat_lda))
  }
  
  fold_sizes = numeric(nfolds)
  for(fold in 1:nfolds) fold_sizes[fold] = length(which(fold_ind==fold))
  test_error = weighted.mean(cv_errors, w=fold_sizes)
  return(test_error)
}

# 0.03953148
(lda_error = reg_cv_lda(breast_cancer[,-10],breast_cancer$y, fold_index))












###################################### QDA #######################################
 set.seed(78889)


## 10-fold cross validation
nfolds = 10
## Sample fold-assignment index
fold_index = sample(nfolds, n, replace=TRUE)
## Print first few fold-assignments
head(fold_index)

# define the function for cross validation.
reg_cv_qda = function(X1, y, fold_ind) {
  Xy = data.frame(X1, y=y)
  nfolds = max(fold_ind)
  if(!all.equal(sort(unique(fold_ind)), 1:nfolds)) stop("Invalid fold partition.")
  cv_errors = numeric(nfolds)
  for(fold in 1:nfolds) {
    
    # Fitting QDA
    tmp_fit_qda = qda(Xy[fold_ind!=fold,]$y ~ ., data=Xy[fold_ind!=fold,][,-10])
    # Predicting the response variable
    phat_qda = predict(tmp_fit_qda, Xy[fold_ind==fold,][-10])
    # convert prediction into numeric from probability
    yhat_qda = phat_qda$class
    yobs = y[fold_ind==fold]
    # Calculating the test error and storing it
    cv_errors[fold] = 1 - mean((yobs == yhat_qda))
  }
  
  fold_sizes = numeric(nfolds)
  for(fold in 1:nfolds) fold_sizes[fold] = length(which(fold_ind==fold))
  test_error = weighted.mean(cv_errors, w=fold_sizes)
  return(test_error)
}

# 0.04978038
(qda_error = reg_cv_qda(breast_cancer[,-10],breast_cancer$y, fold_index))





