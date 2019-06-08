library(randomForest)
library(readr)

dir.create(file.path("output"), showWarnings = FALSE)


data_MS <- read_delim("MS_Data.csv", ";", escape_double = FALSE, 
                      locale = locale(decimal_mark = ",", grouping_mark = "."), 
                      trim_ws = TRUE, skip = 1)

data <- data_MS[,c(6,9:102)]


#############################################
#       Regression Tree Random Forest       #
#############################################


#Create vector with as many numbers as variables, used for amount of predictors considered at tree split
#Also create empty vector for storing the R^2 adjusted and set counter to 1
mtrys <- seq(1, dim(data)[2] - 1)
performance <- vector(mode ="numeric", length = length(mtrys))
i = 1

#Create randomForest models for different m-values and determine the R^2_adj using OOB error for all values of mtry
for (value in mtrys){
  set.seed(1)
  tree.mtry <- randomForest(N_Content ~., data = data, mtry = value, importance =TRUE)
  r2 <- 1 - sum((tree.mtry$y - predict(tree.mtry)) ^ 2) / sum((tree.mtry$y - mean(tree.mtry$y))^2)
  performance[i] <- r2
  i <- i + 1
}

#Plot the R^2  values and associated m-values
jpeg("output/m-value.jpg")
plot(mtrys, performance, main = "R-squared for different m-values in random forest", ylab = "R-squared (based on OOB error)",
     xlab = "m-value")
lines(mtrys, performance, type = "b", col = "red", lwd = 2)
dev.off()

#Pick m-value with highest R^2  
optimal.mtry <- mtrys[which.max(performance)]

#Random forest trees are created for optimal mtry
set.seed(1)
tree.mtry <- randomForest(N_Content ~., data = data, mtry = optimal.mtry, importance =TRUE)

#The trees are plotted in order to find optimal number for n-tree
jpeg("output/n-tree.jpg")
plot(tree.mtry, main = "OOB error vs amount of trees")
dev.off()

#Tree is created for ntree = 250, data is predicted and R^2 adjusted is calculated
set.seed(1)
tree.n <- randomForest(Fresh_Biomass ~ ., data = data, mtry = optimal.mtry, ntree = 200, importance =TRUE)
OOB.n <- mean((predict(tree.n) - data$Fresh_Biomass)^2)
r2_n <- 1 - sum((tree.n$y - predict(tree.n))^2) / sum((tree.n$y - mean(tree.n$y))^2)

#Check importance of variables
importance(tree.n, sort = T)
varImpPlot(tree.n)

#Remove least important (< IncNodePurity) feature, and observe increase in R^2 adj.
set.seed(1)

tree.optimal <- randomForest(N_Content ~ .-b565-b570-b590-b540-b520-b560-b715-b720-b555-b480-b705-b530-b700-b475-b695-b725-b575-b635-b595-b545-b610-b525-b510-b600-b535-b450-b495-b690-b585-b485-b615-b580-b710-b550-b515, data = data, mtry = optimal.mtry, ntree = 200, importance =TRUE)
OOB.optimal <- mean((data$N_Content - predict(tree.optimal))^2)
RMSE.optimal <- sqrt(OOB.optimal)
r2_optimal <- 1 - sum((tree.optimal$y - predict(tree.optimal))^2)/sum((tree.optimal$y - mean(tree.optimal$y))^2)
which.min(importance(tree.optimal)[,1])


#Observe a decrease in performance. The best performing Random Forest model is tree.optimal
4.078, 4.257, 4.013, 4.103, 4.094, 4.286, 4.341, 4.411, 4.143, 4.160, 4.155, 4.267, 4.202, 4.241, 4.199, 4.270, 4.317, -,-,4.345,-,4.32,4.344,-,-,-,4.304,-,4.320,-,-,-,4.363,4.446,-,-,4.400,-,-,-,-,-,-,-,-,-,-,-,-,-,-,



############################# PLS #############################
library(pls)
refl.pls <- plsr(N_Content ~ ., ncomp=20, validation = "LOO", data = data)
validationplot(refl.pls) 
aa <- (RMSEP(refl.pls, estimate = "CV"))
LV = which.min(aa$val)-1

LV = 3

R2.1 <- R2(refl.pls, ncomp=3)
R2.1 <- unlist(R2.1)
R2.1 <- round(R2.1$val2, digits=3)
R2.1

## calculate the RPD 
pred.1 <-as.data.frame(fitted(refl.pls))
pred.1 <- as.numeric(as.matrix(pred.1[1]))

RMSE_pls.1 <- sqrt(mean((data$N_Content - pred.1)^2))
RMSE_pls.1 <- round(RMSE_pls.1, digits=3)
RPD.1 <- round((sd(data$N_Content)/RMSE_pls.1), digits=3)
RPD.1
coefplot(refl.pls, ncom=3, labels = data[,-1])




########## Ridge ##########
library(glmnet)
rr.lm.cv <- lm(N_Content ~ ., data = data)
#First save predictors as a matrix, make sequence for lambda values
preds.cv <- as.matrix(data[, 2:95])
lambdas.cv <- 10 ^ seq(3, -10, by = -.01)

#Cross validate the ridge model to find optimal lambda
rr.fit.cv <- cv.glmnet(preds.cv, as.matrix(data[,1]), alpha = 0, lambda = lambdas.cv)

#Create plot with cros-validation error vs lambda
jpeg("output/optimal_lambda_rr.jpg")
plot(log(rr.fit.cv$lambda), rr.fit.cv$cvm, main = "CV error versus log(lambda)", xlab = "log(lambda)", ylab = "CV error", col = "red")
dev.off()

#Find optimal value for lambda
lambda.cv <- rr.fit.cv$lambda.min
lambda.cv

#Create a model fit for ridge regression, using all features as predictors
#Predict using the ridge regression model, with the optimal lambda
rr.fit.cv <- glmnet(preds.cv, as.matrix(data[,1]), alpha = 0, lambda = lambda.cv)
rr.pred.cv <- predict(rr.fit.cv, s = lambda.cv, newx = preds.cv)

#Calculate the errors and adjusted R squared to check performace
rr.sst.cv <- sum((as.matrix(data[,1]) - mean(as.matrix(data[,1])))^2)
rr.sse.cv <- sum((rr.pred.cv - data[,1])^2)
r2.rr.cv <- 1 - rr.sse.cv/rr.sst.cv
rr.RMSE.cv <- sqrt(rr.sse.cv/56)








### LOOCV ###
library(devtools)
#devtools::install_github("m-Py/minDiff")
library(minDiff)
library(glmnet)

set.seed(123)
data2 <- create_groups(as.data.frame(data), criteria_scale = "N_Content", sets_n = 8, repetitions = 100000)
table(data2$newSet)
tapply(data2$N_Content, data2$newSet, mean)
order.folds <- order(data2$newSet)
data2 <- data2[order.folds, ]
folds <- data.frame(1:7, 8:14, 15:21, 22:28, 29:35, 36:42, 43:49, 50:56)

SD_folds <- sd(c(8.760429, 8.708857, 9.006143, 9.141429, 9.156000, 8.960429, 9.069143, 9.116286))
Mean_folds<- mean(c(8.760429, 8.708857, 9.006143, 9.141429, 9.156000, 8.960429, 9.069143, 9.116286))

summary(data2$N_Content[data2$newSet == 5])

SSE.pls <- vector(length = dim(folds)[1], "numeric")
SSE.rf <- vector(length = dim(folds)[1], "numeric") 
SSE.rr <- vector(length = dim(folds)[1], "numeric") 

R2.pls <- vector(length = dim(folds)[1], "numeric") 
R2.rf <- vector(length = dim(folds)[1], "numeric") 
R2.rr <- vector(length = dim(folds)[1], "numeric") 


set.seed(123)
i=1
for (k in folds) {
  #First split the data, with the one observation as test set, rest as training set
  data.train <- data2[-k, ]
  data.test <- data2[k, ]
  
  # PLS
  pls.fit <- plsr(N_Content ~ ., ncomp=20, data = data.train, validation = "LOO")
  LV.fit <- (RMSEP(pls.fit, estimate = "CV"))
  LV.cv = which.min(LV.fit$val)-1
  pls.fit <- plsr(N_Content ~ ., ncomp=LV.cv, data = data.train)
  pls.pred <- predict(pls.fit, newdata = data.test[,-1])
  
  
  #RIDGE REGRESSION
  #First save predictors as a matrix, make sequence for lambda values
  preds.cv <- as.matrix(data.train[, 2:95])
  lambdas.cv <- 10 ^ seq(3, -10, by = -.01)
  
  #Cross validate the ridge model to find optimal lambda
  rr.fit.cv <- cv.glmnet(preds.cv, as.matrix(data.train[,1]), alpha = 0, lambda = lambdas.cv)
  lambda.cv <- rr.fit.cv$lambda.min
  
  rr.pred <- predict(rr.fit.cv, s = lambda.cv, newx = as.matrix(data.test[,2:95]))
  
  
  ##RANDOM FOREST
  mtrys <- seq(1, dim(data)[2] - 1)
  performance <- vector(mode ="numeric", length = length(mtrys))
  f = 1
  
  #Create randomForest models for different m-values and determine the R^2 using OOB error for all values of mtry
  for (value in mtrys){
    tree.mtry <- randomForest(N_Content ~., data = data, mtry = value, importance =TRUE)
    r2 <- 1 - sum((tree.mtry$y - predict(tree.mtry)) ^ 2) / sum((tree.mtry$y - mean(tree.mtry$y))^2)
    performance[f] <- r2
    f <- f + 1
  }
  
  optimal.mtry <- mtrys[which.max(performance)]
  
  rf.fit <- randomForest(N_Content ~ ., data = data, mtry = optimal.mtry, ntree = 250, importance =TRUE)
  rf.pred <- predict(rf.fit, newdata = data.test[,-1])
  

  #Determine the MSE for the two predictions and store in MSE vectors
  SSE.pls[i] <- sum((data.test[,1] - pls.pred)^2)
  SSE.rf[i] <- sum((data.test[,1] - rf.pred)^2)
  SSE.rr[i] <- sum((data.test[,1] - rr.pred)^2)
  
  R2.pls[i] <- 1 - (SSE.pls[i] / sum((data.test$N_Content - mean(data.train$N_Content))^2))
  R2.rf[i] <- 1 - (SSE.rf[i] / sum((data.test$N_Content - mean(data.train$N_Content))^2))
  R2.rr[i] <- 1 - (SSE.rr[i] / sum((data.test$N_Content - mean(data.train$N_Content))^2))
  
  
  i = i+1
}
hist(SSE.pls, main = "Histogram of SSE of all folds for PLS", xlab = "Sum of squared errors")
hist(SSE.rf, breaks = 5, main = "Histogram of SSE of all folds for random forest", xlab = "Sum of squared errors")
hist(SSE.rr, main = "Histogram of SSE of all folds for ridge regression", xlab = "Sum of squared errors")

sd(SSE.pls[-5])

boxplot(SSE.pls, SSE.rf, SSE.rr, names = c("PLS", "RF", "RR"), boxwex = 0.3, main = "Boxplot sum of squared errors 8-fold CV", ylab = "SSE")

CV.pls <- mean(SSE.pls)
CV.rf<- mean(SSE.rf)
CV.rr <- mean(SSE.rr)

RMSEP.pls <- sqrt(CV.pls)
RMSEP.rf <- sqrt(CV.rf)
RMSEP.rr <- sqrt(CV.rr)

SDE.pls <- sd(SSE.pls)
SDE.rf <- sd(SSE.rf)
SDE.rr <- sd(SSE.rr)

SSE.pls

plot(SSE.pls)

mean(R2.pls)
mean(R2.rf)
R2.pls
R2.rf
sd(R2.pls)
sd(R2.rf)
median(R2.pls)
median(R2.rf)

hist()

hist(data$N_Content)

SSE.pls
SSE.rf
