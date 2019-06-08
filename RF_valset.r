library(randomForest)
library(readr)

rm(list=ls())

dir.create(file.path("output"), showWarnings = FALSE)


data_MS <- read_delim("MS_Data.csv", ";", escape_double = FALSE, 
                      locale = locale(decimal_mark = ",", grouping_mark = "."), 
                      trim_ws = TRUE, skip = 1)

###############################
# splitting the data          #
###############################
data <- data_MS[,c(3,6,9:102)]
# data$VS <- ifelse(data$Treatment == "Vs", 1, 0)
# data$LPTR  <- ifelse(data$Treatment == "Lp+Tr", 1, 0)
# data$RSVS <- ifelse(data$Treatment == "Rs+Vs", 1, 0)
# data$TR <- ifelse(data$Treatment == "Tr", 1, 0)
# data$LP <- ifelse(data$Treatment == "Lp", 1, 0)
# data$FA <- ifelse(data$Treatment == "Fa", 1, 0)

data_t <-data[data$Cal_Val == "Validation",]
data_t <-data_t[,-1]

data <- data[data$Cal_Val == "Calibration",]
data <- data[,-1]

MN <- as.vector(data[,1])
MN <- sum(MN)/28

################################
# random forest validation set #
################################

mtrys <- seq(1, dim(data)[2] - 1)
performance <- vector(mode ="numeric", length = length(mtrys))
performance_2 <- vector(mode ="numeric", length = length(mtrys))
OOBerror <- vector(mode ="numeric", length = length(mtrys))


i = 1
for (value in mtrys){
  set.seed(1)
  tree.mtry <- randomForest(N_Content ~., data = data, mtry = value, importance =TRUE)
  OOBerror[i] <- mean((tree.mtry$y - predict(tree.mtry)) ^ 2)
  r2 <- 1 - sum((tree.mtry$y - predict(tree.mtry)) ^ 2) / sum((tree.mtry$y - mean(tree.mtry$y))^2)
  r2.adj <- 1 - (1-r2)  * (( (dim(data)[1] - 1))/ (dim(data)[1] - dim(data)[2] - 2))
  performance[i] <- r2
  performance_2[i] <- r2.adj
  i <- i + 1
}


plot(mtrys, OOBerror, main = "R-squared for different m-values in random forest", ylab = "R-squared (based on OOB error)",
     xlab = "m-value")
lines(mtrys, performance, type = "b", col = "red", lwd = 2)

optimal.mtry <- mtrys[which.max(performance)]


#Random forest trees are created for optimal mtry
set.seed(1)
tree.optimal <- randomForest(N_Content ~., data = data, mtry = optimal.mtry, ntree = 200, importance =TRUE)

OOB.rf.val <- sum((data[,1] - predict(tree.optimal))^2)
r2.rf.val <- 1 - OOB.rf.val/pls.sst
RMSE.rf.val <- sqrt(OOB.rf.val/28) 

pls.sst <- sum((as.matrix(data[,1]) - mean(as.matrix(data[,1])))^2)
pls.sse <- sum((pred.cal - data[,1])^2)
r2.pls.cal <- 1 - pls.sse/pls.sst
pls.RMSE.cal <- sqrt(pls.sse/28)


#The trees are plotted in order to find optimal number for n-tree
jpeg("output/n-tree.jpg")
plot(tree.mtry, main = "OOB error vs amount of trees")
dev.off()

#tree.optimal <- randomForest(N_Content ~ .-b715-b530-b585-b560-b720-b515-b575-b560-b650-b645-b535, data = data, mtry = optimal.mtry, ntree = 200, importance =TRUE)

rf.pred.val <- predict(tree.optimal, newdata = data_t)
rf.residuals <- data_t[,1] - rf.pred.val
rf.MSEP.val <- mean((data_t[,1] - rf.pred.val)^2)
RMSEP.rf.val <- sqrt(rf.MSEP.val)

plot(rf.pred.val, rf.residuals[,1], main = "Prediction nitrogen content vs. residuals \nfor validation set random forest", xlab = expression(paste("nitrogen predicted [% g"^"-1", "]")), ylab = "residual")
plot (rf.pred.val ~ data_t$N_Content, xlab=expression(paste("nitrogen observed [% g"^"-1", "]")), ylab=expression(paste("nitrogen predicted [% g"^"-1", "]")), main="Nitrogen observed vs predicted")

plot(tree.optimal, log="y")
varImpPlot(tree.optimal, n.var = 10)
varImp(tree.optimal)
hist(rf.residuals$N_Content)

varImp(tree.optimal)

#############################################
#      PLSR validation set                  #
#############################################

library(pls)
refl.pls <- plsr(N_Content ~ as.matrix(data[,-1]), ncomp=20, validation = "LOO", data = data)
validationplot(refl.pls) 

LV = 11

# r2.cal.pls <- R2(refl.pls, ncomp=11)
# r2.cal.pls <- unlist(r2.cal.pls)
# r2.cal.pls <- round(r2.cal.pls$val2, digits=3)
# r2.cal.pls
pred.cal <- predict(refl.pls, ncomp = 11, newdata = data[,-1])

pls.sst <- sum((as.matrix(data[,1]) - mean(as.matrix(data[,1])))^2)
pls.sse <- sum((pred.cal - data[,1])^2)
r2.pls.cal <- 1 - pls.sse/pls.sst
pls.RMSE.cal <- sqrt(pls.sse/28)

r2.pls.val <- 1 - sum((data$N_Content - pred.cal)^2)/sum((data$N_Content - mean(data$N_Content))^2)
pred.1 <-as.data.frame(fitted(refl.pls))
pred.1 <- as.numeric(as.matrix(pred.1[1]))

RMSE.pls.cal <- sqrt(mean((data$N_Content - pred.1)^2))
RMSE.pls.cal <- round(RMSE.pls.cal, digits=3)

coefplot(refl.pls, ncomp=11, labels = data[,-1])

### Test set
pls.fit.val <- plsr(N_Content ~ ., ncomp=11, data = data)


pls.pred.val <- predict(pls.fit.val, ncomp = 11, newdata = data_t[,-1])
pls.residuals <- data_t$N_Content - pls.pred.val
pls.residuals <- as.matrix(pls.residuals)
pls.RMSEP.val <- RMSEP(pls.fit.val, newdata = data_t, ncomp = 11, estimate = "test")
pls.RMSEP.val

varImp(pls.fit.val, scale = T)

plot(pls.pred.val, pls.residuals[,1], main = "Prediction nitrogen content vs. residuals \nfor validation set PLS", xlab = expression(paste("nitrogen predicted [% g"^"-1", "]")), ylab = "residual")
plot(pls.pred.val ~ data_t$N_Content, xlab=expression(paste("nitrogen observed [% g"^"-1", "]")), ylab=expression(paste("nitrogen predicted [% g"^"-1", "]")), main="Nitrogen observed vs predicted")
hist(pls.residuals)



####################
# Ridge regression #
####################
library(glmnet)
library(broom)
library(caret)

rr.lm <- lm(N_Content ~ ., data = data)
#First save predictors as a matrix, make sequence for lambda values
preds <- as.matrix(data[, 2:95])
lambdas <- 10 ^ seq(3, -10, by = -.01)

#Cross validate the ridge model to find optimal lambda
set.seed(1)
rr.fit.cal <- cv.glmnet(preds, as.matrix(data[,1]), alpha = 0, lambda = lambdas)

#Create plot with cros-validation error vs lambda
jpeg("output/optimal_lambda_rr.jpg")
plot(log(rr.fit.cal$lambda), rr.fit.cal$cvm, main = "CV error versus log(lambda)", xlab = "log(lambda)", ylab = "CV error", col = "red")
dev.off()

#Find optimal value for lambda
lambda <- rr.fit.cal$lambda.min
lambda

#Create a model fit for ridge regression, using all features as predictors
#Predict using the ridge regression model, with the optimal lambda
rr.pred.cal <- predict(rr.fit.cal, s = lambda, newx = preds)

#Calculate the errors and adjusted R squared to check performace
rr.sst <- sum((as.matrix(data[,1]) - mean(as.matrix(data[,1])))^2)
rr.sse <- sum((rr.pred.cal - data[,1])^2)
r2.rr.cal <- 1 - rr.sse/rr.sst
rr.RMSE.cal <- sqrt(rr.sse/28)

rr.pred.val <- predict(rr.fit.cal, s = lambda, newx = as.matrix(data_t[,-1]))
rr.sse.val <- sum((rr.pred.val - data_t[,1])^2)
rr.RMSEP.val <- sqrt(rr.sse.val/28)

rr.fit.cal1 <- glmnet(preds, as.matrix(data[,1]), alpha = 0, lambda = lambda)
varImp(rr.fit.cal1, lambda = lambda, scale =T, nonpara = T)

rr.residuals <- data_t[,1] - rr.pred.val
plot(rr.pred.val, rr.residuals[,1], main = "Prediction nitrogen content vs. residuals \nfor validation set ridge regression", xlab = expression(paste("nitrogen predicted [% g"^"-1", "]")), ylab = "residual")
hist(rr.residuals$N_Content)
plot(rr.pred.val ~ data_t$N_Content, xlab=expression(paste("nitrogen observed [% g"^"-1", "]")), ylab=expression(paste("nitrogen predicted [% g"^"-1", "]")), main="Nitrogen observed vs predicted")

rr.fit.cal$glmnet.fit

rr.fit.cal$df

plot(rr.fit.cal$beta, preds)

get( getOption( "device" ) )()
par(mfrow=c(3,2), cex = 0.53)
plot(pls.pred.val, pls.residuals[,1],cex.lab = 1.1, xlab = expression(paste("N predicted [% g"^"-1", "]")), ylab = "residual")
abline(h = 0, lty = 2)
plot(pls.pred.val ~ data_t$N_Content,cex.lab = 1.1, xlab=expression(paste("N observed [% g"^"-1", "]")), ylab=expression(paste("N predicted [% g"^"-1", "]")))
abline(0,1, lty=2)
plot(rr.pred.val, rr.residuals[,1],cex.lab = 1.1, xlab = expression(paste("N predicted [% g"^"-1", "]")), ylab = "residual")
abline(h = 0, lty = 2)
plot(rr.pred.val ~ data_t$N_Content,cex.lab = 1.1, xlab=expression(paste("N observed [% g"^"-1", "]")), ylab=expression(paste("N predicted [% g"^"-1", "]")))
abline(0,1, lty=2)
plot(rf.pred.val, rf.residuals[,1], cex.lab = 1.1, xlab = expression(paste("N predicted [% g"^"-1", "]")), ylab = "residual")
abline(h = 0, lty = 2)
plot (rf.pred.val ~ data_t$N_Content,cex.lab = 1.1, xlab=expression(paste("N observed [% g"^"-1", "]")), ylab=expression(paste("N predicted [% g"^"-1", "]")))
abline(0,1, lty=2)


m <- matrix(c(1,2,3,4,5,6,7,7,7),nrow = 3,ncol = 3,byrow = TRUE)

layout(mat = m,heights = c(0.4,0.4,0.2))

for (i in 1:6){
  par(mar = c(2,2,1,1))
  plot(runif(5),runif(5),xlab = "",ylab = "")
}


