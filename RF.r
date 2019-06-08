library(randomForest)
library(readr)

dir.create(file.path("output"), showWarnings = FALSE)


data_MS <- read_delim("MS_Data.csv", ";", escape_double = FALSE, 
                      locale = locale(decimal_mark = ",", grouping_mark = "."), 
                      trim_ws = TRUE, skip = 1)

data <- data_MS[,c(6,9:102)]

data <- data[data$Cal_Val == "Calibration",]
data <- data[,-1]

data_t <-data_MS[,c(3,6,9:102)]
data_t <-data_t[data_t$Cal_Val == "Validation",]
data_t <-data_t[,-1]

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
tree.n <- randomForest(N_Content ~ ., data = data, mtry = optimal.mtry, ntree = 250, importance =TRUE)
OOB.n <- mean((predict(tree.n) - data$N_Content)^2)
r2_n <- 1 - sum((tree.n$y - predict(tree.n))^2) / sum((tree.n$y - mean(tree.n$y))^2)

#Check importance of variables
importance(tree.n, sort = T)
varImpPlot(tree.n)

#Remove least important (< IncNodePurity) feature, and observe increase in R^2 adj.
set.seed(1)
tree.optimal <- randomForest(N_Content ~ .-b565-b570-b590-b540-b520-b560-b715-b720-b555-b480-b705-b530-b700-b475-b695-b725-b575-b635-b595-b545-b610-b525-b510-b600-b535-b450-b495-b690-b585-b485-b615-b580-b710-b550-b515, data = data, mtry = optimal.mtry, ntree = 100, importance =TRUE)
OOB.optimal <- mean((data$N_Content - predict(tree.optimal))^2)
r2_optimal <- 1 - sum((tree.optimal$y - predict(tree.optimal))^2)/sum((tree.optimal$y - mean(tree.optimal$y))^2)
which.min(importance(tree.optimal)[,1])
importance(tree.optimal)
varImpPlot(tree.optimal)


#Observe a decrease in performance. The best performing Random Forest model is tree.optimal
4.078, 4.257, 4.013, 4.103, 4.094, 4.286, 4.341, 4.411, 4.143, 4.160, 4.155, 4.267, 4.202, 4.241, 4.199, 4.270, 4.317, -,-,4.345,-,4.32,4.344,-,-,-,4.304,-,4.320,-,-,-,4.363,4.446,-,-,4.400,-,-,-,-,-,-,-,-,-,-,-,-,-,-,



############################# PLS #############################
library(pls)
refl.pls <- plsr(N_Content ~ ., ncomp=20, validation = "LOO", data = data)
validationplot(refl.pls) 

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

### LOOCV ###
hist()

folds <- data.frame(1:8, 9:16, 17:24, 25:32, 33:40, 41:48, 49:56)
folds <- data.frame(1:7, 8:14, 15:21, 22:28, 29:35, 36:42, 43:49, 50:56)
data <- data[sample(nrow(data)),]
folds <- seq(1,56,1)

SSE.pls <- vector(length = dim(folds)[1], "numeric")
SSE.rf <- vector(length = dim(folds)[1], "numeric") 
R2.pls <- vector(length = dim(folds)[1], "numeric") 
R2.rf <- vector(length = dim(folds)[1], "numeric") 

set.seed(1)
i=1
for (k in folds) {
  #First split the data, with the one observation as test set, rest as training set
  data.train <- data[-k, ]
  data.test <- data[k, ]
  
  #Split the data for the backwards stepwise selection model in the same way
  pls.fit <- plsr(N_Content ~ ., ncomp=3, data = data.train)
  pls.pred <- predict(pls.fit, newdata = data.test[,-1])

  #pls.pred <- predict(pls.fit, newdata = data.test)
  #pls.pred <-as.data.frame(fitted(pls.fit))
  #pls.pred <- as.numeric(as.matrix(pls.pred[1]))
  
  #Fit linear model using all predictors of the training set, test on the test data
  rf.fit <- randomForest(N_Content ~ .-b565-b570-b590-b540-b520-b560-b715-b720-b555-b480-b705-b530-b700-b475-b695-b725-b575-b635-b595-b545-b610-b525-b510-b600-b535-b450-b495-b690-b585-b485-b615-b580-b710-b550-b515, data = data, mtry = optimal.mtry, ntree = 100, importance =TRUE)
  rf.pred <- predict(rf.fit, newdata = data.test[,-1])
  
  #Fit linear model using the backwards selection predictors of training set and predict test data

  #Determine the MSE for the two predictions and store in MSE vectors
  SSE.pls[i] <- sum((data.test[,1] - pls.pred)^2)
  SSE.rf[i] <- sum((data.test[,1] - rf.pred)^2)
  
  R2.pls[i] <- 1 - (SSE.pls[i] / sum((data.test$N_Content - mean(data.train$N_Content))^2))
  R2.rf[i] <- 1 - (SSE.rf[i] / sum((data.test$N_Content - mean(data.train$N_Content))^2))
  
  i = i+1
  
  #Determine MSE for ridge regression by applying the shrinkage term: the lambda's found multiplied by the coefficients of the models
  #MSE.lm.rr[k] <- MSE.lm[k] + lambda * sum(k.fit$coefficients^2)
  #MSE.bs.rr[k] <- MSE.bs[k] + lambda.bs * sum(k.fit.bs$coefficients^2) 
}

CV.pls <- mean(SSE.pls)
CV.rf<- mean(SSE.rf)

RMSEP.pls <- sqrt(CV.pls)
RMSEP.rf <- sqrt(CV.rf)


mean(R2.pls)
mean(R2.rf)
R2.pls
R2.rf
sd(R2.pls)
sd(R2.rf)
median(R2.pls)
median(R2.rf)
hist(R2.rf)
hist(data$N_Content)

MSE.pls
MSE.rf
