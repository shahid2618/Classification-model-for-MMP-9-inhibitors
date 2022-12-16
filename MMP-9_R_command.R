
library(e1071)
library(caret)
library(dplyr)

options(repos = c(CRAN = "http://cran.rstudio.com"))
dat <- read.csv("Pubchem_MMP-9_5_Zinc.csv", header = T)

may <- read.csv("pubchem_finger_maybridge.csv", header = T)
ev_large <- read.csv("EV_mmp9_Pubchem_final.csv", header = T)
ev_decoy <- read.csv("Active_683_Decoy_EV_MMP9_Pubchem.csv", header = T)

########################### NORMALIZE########
normalize <- function(x) {
  num <- x - min(x)
  denom <- max(x) - min(x)
  return(num/denom)
}

dat_norm <- as.data.frame(lapply(dat[2:882], normalize))


ev_norm <- as.data.frame(lapply(ev[2:882], normalize))
names(dat)[168]<-"class"

### Zero Variable ####
nsv <- nearZeroVar(dat_norm)
dat_norm <- dat_norm[,-nsv]

nsv <- nearZeroVar(ev_norm)
ev_norm<- ev_norm[,-nsv]

dat1 <- cbind(dat_norm, dat$class)
datset <- cbind(dat1, dat$Name)

datset <- cbind(dat_norm, dat$CLASS)
levels(datset$`dat$CLASS`)[levels(datset$`dat$CLASS`)=="ACTIVE"] <- "1"
levels(datset$`dat$CLASS`)[levels(datset$`dat$CLASS`)=="INACTIVE"] <- "0"


names(datset)[340]<-"class"
datset <- cbind(datset, dat$Name)

### Data partition ###
library(caTools)
set.seed(1234)
split <- sample.split(datset, SplitRatio =0.75)
training_set <- subset(datset, split == TRUE)                                                                                                                                                                                                                                          
test_set <- subset(datset, split == FALSE)


## New RF code using for modelling 


###### Random Forest ###
library(randomForest)
customRF <- list(type = "Classification",
                 library = "randomForest",
                 loop = NULL)

customRF$parameters <- data.frame(parameter = c("mtry", "ntree"),
                                  class = rep("numeric", 2),
                                  label = c("mtry", "ntree"))

customRF$grid <- function(x, y, len = NULL, search = "grid") {}

customRF$fit <- function(x, y, wts, param, lev, last, weights, classProbs) {
  randomForest(x, y,
               mtry = param$mtry,
               ntree=param$ntree)
}

#Predict label
customRF$predict <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
  predict(modelFit, newdata)

#Predict prob
customRF$prob <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
  predict(modelFit, newdata, type = "prob")

customRF$sort <- function(x) x[order(x[,1]),]
customRF$levels <- function(x) x$classes
library(doParallel)
cores <- makeCluster(detectCores()-1)
registerDoParallel(cores = cores)

control <- trainControl(method="repeatedcv", 
                        number=10, 
                        repeats=3,
                        allowParallel = TRUE)
tunegrid <- expand.grid(.mtry=c(1:15),.ntree=c(100,200,300,500))
set.seed(1234)
metric <- "Accuracy"
class.weights= c("0" = 1, "1" = 3)
custom <- train(class~., data=training_set, 
                method=customRF, 
                metric=metric, 
                tuneGrid=tunegrid, 
                trControl=control,
                class.weights=class.weights)
custom
rf<- predict(custom, test_set)
confusionMatrix(rf, test_set$class)
train.pred <- predict(custom, training_set)
table(training_set$class, train.pred)

pred_rf <- prediction(as.numeric(rf), as.numeric(test_set$class))
perf_rf <- performance(pred_rf, 'tpr','fpr')
plot(perf_rf,main='ROC Curve for random forest_Pubchem',xlab ="False positive rate (1-specificity)", 
     ylab = "True positive rate (sensitivity)", col= "blue")
abline(a=0, b= 1,lty=2)


pred_may_1 <- predict(custom, may[2:882])
summary(pred_may_1)
write.table(pred_may_1, file = "VS_pubchem_MMP9_Zinc.csv")
VS_pubchem_MMP9_Zinc <- read.csv("VS_pubchem_MMP9_Zinc.csv", header = T)

write.table(pred_may_1, file = "MACCS_zinc_rf_new.csv")
MACCS_zinc_rf_new <- read.csv("MACCS_zinc_rf_new.csv", header = T)

##### EV ###
ev_large<- ev_large[2:883]
ev_decoy_2<- ev_decoy_2[2:883]

ev_test_large <- predict(custom, ev_large)
table(ev$CLASS, ev_test_large)

EV.Pubchem_3_MMP9 <- predict(custom, ev)
table(ev$CLASS, EV.Pubchem_3_MMP9)

EV_test_decoy <- predict(custom, ev_decoy_2[2:883])
table(ev_decoy_2$CLASS, EV_test_decoy)


#################### SVM ################

#### SVM ###
# Tuning Parameters #
yTrain <- training_set$class
xTrain <- training_set[1:332]
tune_svm <- tune(svm, train.x=xTrain, train.y=yTrain, 
                 kernel="radial", ranges=list(cost=10^(-3:3), gamma= 10^(-6:-3)),class.weights= c("0" = 1, "1" = 3))

tune_svm
trctrl <- trainControl(method = "cv", number = 10)
set.seed(1234)
grid <- data.frame(C=seq(0.01,5,0.5))
final_svm <- svm(class ~ ., data= training_set, kernel="radial",trControl=trctrl,tuneGrid=grid, cost=10, gamma= 0.001, class.weights= c("0" = 1, "1" = 3))
final_svm
pred <- predict(final_svm, test_set)
confusionMatrix(pred, test_set$class)

train.pred <- predict(final_svm, training_set)
table(training_set$class, train.pred)


pred_may_1 <- predict(final_svm, may[2:882])
summary(pred_may_1)
write.table(training_set, file = "trainingset_5_MMP9.csv")
trainingset_5_MMP9 <- read.csv("trainingset_5_MMP9.csv", header = T)

write.table(test_set, file = "testset_5_MMP9.csv")
testset_5_MMP9 <- read.csv("testset_5_MMP9.csv", header = T)



pred_svm_radial <- prediction(as.numeric(pred), as.numeric(test_set$class))
perf_svm_radial <- performance(pred_svm_radial, 'tpr','fpr')
plot(perf_svm_radial,main='ROC Curve for SVM_radial_Pubchem',xlab ="False positive rate (1-specificity)", 
     ylab = "True positive rate (sensitivity)", col= "red")
abline(a=0, b= 1,lty=2)



## EV####

ev<- ev[2:883]
ev.Pubchem_3_MMP9 <- predict(final_svm, ev)
table(ev$CLASS, ev.Pubchem_3_MMP9)

ev_test_decoy <- predict(final_svm, ev_decoy_2[2:883])
table(ev_decoy_2$CLASS, ev_test_decoy)



ev_test_large <- predict(custom, ev_large)
table(ev$CLASS, ev_test_large)

### NB #####

class.weights= c("0" = 1, "1" = 3)
model.nb_wt <- naiveBayes(class ~., data = training_set, trControl=control,class.weights= class.weights) 
pred_nb <- predict(model.nb_wt, test_set)
confusionMatrix(pred_nb, test_set$class)

train.pred <- predict(model.nb_wt, training_set)
table(training_set$class, train.pred)

prediction_nb <- prediction(as.numeric(pred_nb), as.numeric(test_set$class))
perf_nb <- performance(prediction_nb, 'tpr','fpr')
plot(perf_nb,main='ROC Curve for NB_Pubchem',xlab ="False positive rate (1-specificity)", 
     ylab = "True positive rate (sensitivity)", col= "green")
abline(a=0, b= 1,lty=2)

nbauc<-performance(pred_nb,"auc")@y.values[[1]]
svmauc<-performance(pred_svm_radial,"auc")@y.values[[1]]
## EV ##
EV_test_naive <- predict(model.nb_wt, ev)
table(ev$CLASS,EV_test_naive)

ev.Pubchem_3_MMP9 <- predict(model.nb_wt, ev)
table(ev$CLASS, ev.Pubchem_3_MMP9)

ev_test_decoy <- predict(model.nb_wt, ev_decoy_2[2:883])
table(ev_decoy_2$CLASS, ev_test_decoy)



##### Radial SVM ###
rbf.tune <- tune.svm(class~., data=training_set, kernel="radial", gamma=c(0.0001,0.001,0.01,0.1,0.5,1))


rbf.tune <- tune.svm(class~., data=training_set, kernel="radial",gamma= 0.001,trControl = control,
                     class.weights= c("0" = 1, "1" = 10))
rbf.tune

best.rbf <- rbf.tune$best.model
rbf.test <- predict(best.rbf, newdata=may[2:1026])
confusionMatrix(rbf.test, test_set$class, positive="1")
confusionMatrix(rbf.test, may$CLASS, positive="1")


train.pred <- predict(best.rbf, training_set)
table(training_set$class, train.pred)
test.pred <- predict(best.rbf, test_set)
table(test_set$class,test.pred)

pred_svm_radial <- prediction(as.numeric(rbf.test), as.numeric(test_set$class))

perf_svm_radial <- performance(pred_svm_radial, 'tpr','fpr')
plot(perf_svm_radial)
plot(perf_svm_radial,col='green',lty=1, main='ROC Curve for SVM_radial Extended_P fingerprints')
L<-list(bquote("SVM"== .(perf_svm_radial)))
abline(a=0, b= 1,lty=2)
plot(perf,colorize = TRUE, curve=TRUE)

svmauc<-performance(pred_svm_radial,"auc")@y.values[[1]]
nbauc<-performance(pred_nb,"auc")@y.values[[1]]

## EV ##
may.test <- predict(best.linear, may)
table(may$CLASS,may.test)


##### Combine All SVM algorithm ####

plot(perf_nb,col='green',lty=1,main='The Comparison of ROC Curve for dataset-B')
plot(perf_svm_radial,col='red',add=TRUE,lty=1)
plot(perf_rf,col='blue',add=TRUE,lty=1)
abline(a=0, b= 1,lty=2)

L<-list(bquote("NB"== .(perf_nb)), bquote("SVM"== .(perf_svm_radial)), bquote("RF"== .(perf_rf)))
legend(x=.73, y=.30, legend=c("SVM (1B)", "RF (2B)", "NB (3B)"), col=c("red", "blue", "green"), lty=1)  


#####
plot(rf.auc, colorize = TRUE,main='ROC Curve for RF')
L<-list(bquote("AUC"== .(rfauc)))
legend("bottomright",legend=L,col=c('green'),lwd=2,bg="gray",pch=18,text.font=2.5,cex=0.9)
abline(a=0, b= 1,lty=2)

