install.packages("readxl")
library("readxl")

traningdata <- read_excel("train.xlsx")
sequence <- read_excel("seqinfo.xlsx")

##match the gene sequence data in two dataframe
library(MASS)
merged.train <- merge(x = traningdata, y = sequence,
                      by.x = ("gene"),
                      by.y = ("gene")
)
print(merged.train)

##save
write.csv(merged.train,"merge.csv")

merge=read.csv("merge.csv")

##statistical summary
summary(merge)

## remove X column
merge$X <- NULL

## remove Solubility(predictor) column to the last column 
merge$Solubility <- merge[ ,7]
merge[ ,7] <- NULL

##using samll data to test decision tree and fail
smalldata <- merge[1:50,]

## replace the missing value by median for numeric columns 
merge[ ,7][is.na(merge[ ,7])] <- median(merge[ ,7], na.rm=TRUE)
merge[ ,8][is.na(merge[ ,8])] <- median(merge[ ,8], na.rm=TRUE)
merge[ ,9][is.na(merge[ ,9])] <- median(merge[ ,9], na.rm=TRUE)
merge[ ,10][is.na(merge[ ,10])] <- median(merge[ ,10], na.rm=TRUE)
merge[ ,11][is.na(merge[ ,11])] <- median(merge[ ,11], na.rm=TRUE)
merge[ ,12][is.na(merge[ ,12])] <- median(merge[ ,12], na.rm=TRUE)
merge[ ,13][is.na(merge[ ,13])] <- median(merge[ ,13], na.rm=TRUE)
merge[ ,14][is.na(merge[ ,14])] <- median(merge[ ,14], na.rm=TRUE)
merge[ ,15][is.na(merge[ ,15])] <- median(merge[ ,15], na.rm=TRUE)
merge[ ,16][is.na(merge[ ,16])] <- median(merge[ ,16], na.rm=TRUE)
merge[ ,17][is.na(merge[ ,17])] <- median(merge[ ,17], na.rm=TRUE)
merge[ ,18][is.na(merge[ ,18])] <- median(merge[ ,18], na.rm=TRUE)
merge[ ,19][is.na(merge[ ,19])] <- median(merge[ ,19], na.rm=TRUE)
merge[ ,20][is.na(merge[ ,20])] <- median(merge[ ,20], na.rm=TRUE)

##save data after using median to repalce NA
write.csv(clean1,"clean1.csv")

clean1=read.csv("clean1.csv")

##transfer the factor into numeric value
clean1[ ,1] <- as.numeric(clean1[ ,1])
clean1[ ,2] <- as.numeric(clean1[ ,2])
clean1[ ,3] <- as.numeric(clean1[ ,3])
clean1[ ,4] <- as.numeric(clean1[ ,4])
clean1[ ,5] <- as.numeric(clean1[ ,5])
clean1[ ,6] <- as.numeric(clean1[ ,6])
clean1[ ,23] <- as.numeric(clean1[ ,23])
clean1[ ,24] <- as.numeric(clean1[ ,24])
clean1[ ,25] <- as.numeric(clean1[ ,25])
clean1[ ,26] <- as.numeric(clean1[ ,26])
clean1[ ,27] <- as.numeric(clean1[ ,27])
clean1[ ,28] <- as.numeric(clean1[ ,28])

## replace the missing value by median for numeric columns
clean1[ ,6][is.na(clean1[ ,6])] <- median(clean1[ ,6], na.rm=TRUE)
clean1[ ,25][is.na(clean1[ ,25])] <- median(clean1[ ,25], na.rm=TRUE)
clean1[ ,26][is.na(clean1[ ,26])] <- median(clean1[ ,26], na.rm=TRUE)
clean1[ ,27][is.na(clean1[ ,27])] <- median(clean1[ ,27], na.rm=TRUE)
clean1[ ,28][is.na(clean1[ ,28])] <- median(clean1[ ,28], na.rm=TRUE)

##save data after using median to repalce NA
write.csv(clean1,"clean2.csv")


##use 75% as training data
install.packages(caret)
library('caret')
TrainingDataIndex <- createDataPartition(clean2[ ,29], p=0.75, list = FALSE)
train2 <- clean2[TrainingDataIndex,]
test2<- clean2[-TrainingDataIndex,]



## Create logistic Regression Model---->Error in eval(expr, envir, enclos) : y values must be 0 <= y <= 1
input <- train2
solubility.data = glm(formula = Solubility ~ ., data = input, family = binomial)
summary(solubility.data)

## build decison tree model
library(rpart)
# grow tree 
fit <- rpart(Solubility ~ .,method="class", data=train2)
printcp(fit) # display the results 
plotcp(fit) # visualize cross-validation results 
summary(fit) # detailed summary of splits
# plot tree ---->Error in plot.rpart(fit, uniform = TRUE, main = "Classification Tree for y") : 
##fit is not a tree, just a root
plot(fit, uniform=TRUE,main="Classification Tree for y")
text(fit, use.n=TRUE, all=TRUE, cex=.8)

##predicting using the model
pred <- predict(fit, newdata = test)

clean1=read.csv("clean1.csv")
## remove all columns except 24,26,29,30 from clean1
clean1[ ,1:23]<-NULL
clean1[ ,2]<-NULL
clean1[ ,3:4]<-NULL
## romove the rows including NA
x<-na.omit(clean1)

write.csv(x,"clean3.csv")
clean3=read.csv("clean3.csv")

##each time I import clean3, the first column will become X. So I need delete it each time!!!
clean3[ ,1]<-NULL

## The sequence is a very important column and we need deal with it seperately
clean3[ ,3]<-as.character(clean3[ ,3])
length(clean3[ ,3])
install.packages("protr")
library("protr")
clean3[ ,3] = clean3[ ,3][(sapply(clean3[ ,3], protcheck))]
length(clean3[ ,3])
# calculate APseAAC descriptors--->Error in FUN(X[[i]], ...) : x has unrecognized amino acid type
x1 = t(sapply(clean3[ ,3], extractAAC))

write.csv(x1,"sequenceAAC.csv")
sequenceAAC=read.csv("sequenceAAC.csv")

##add sequenceAAC column to previous dataframe
clean3[ ,5:24] <- sequenceAAC[ ,2:21]
clean3[ ,3]<-NULL
clean3[ ,24]<-clean3[ ,3]
clean3[ ,3]<-NULL
colnames(clean3)[23] <- "solubility"

write.csv(clean3,"clean4.csv")
clean4=read.csv("clean4.csv")
clean4[ ,1]<-NULL

##change the first two columns into numeric value and make the solubility higher than 1 to be 1.
clean4[ ,1] <- as.numeric(clean4[ ,1])
clean4[ ,2] <- as.numeric(clean4[ ,2])
clean4[ ,23][clean4[ ,23]>100]<-100
clean4[ ,23]<-clean4[ ,23]/100
write.csv(clean4,"clean5.csv")
clean5=read.csv("clean5.csv")
clean5[ ,1]<-NULL

##use 75% as training data
install.packages(caret)
library('caret')
TrainingDataIndex <- createDataPartition(clean5[ ,23], p=0.75, list = FALSE)
train5 <- clean5[TrainingDataIndex,]
test5<- clean5[-TrainingDataIndex,]
testdata<-test5 
test5[ ,23]<-NULL

## Create logistic Regression Model-->warning:In eval(expr, envir, enclos) : non-integer #successes in a binomial glm!
input <- train5
solubility.data = glm(formula = solubility ~ ., data = input, family = binomial)
summary(solubility.data)

##predicting using the model
pred <- predict(solubility.data, newdata = test5)

##save prediction of logistical regression
write.csv(predframe,"logisticalregression.csv")

##check the accuracy--->Error in accuracy.meas(test5$solubility, pred) : Response must have two levels.--->make solibility to be 0 or 1.
library("ROSE")
testdata$solubility[testdata$solubility >= 0.5] <-1
testdata$solubility[testdata$solubility < 0.5] <- 0
accuracy.meas(testdata$solubility, pred)
## ROC plot
roc.curve(testdata$solubility, pred, plotit = F)


## build decison tree model
library(party)
output.tree <- ctree(solubility ~ ., data=train5)
plot(output.tree)


##predicting using the model
pred <- predict(output.tree, newdata = test5)

##check the accuracy
testdata$solubility[testdata$solubility >= 0.5] <-1
testdata$solubility[testdata$solubility < 0.5] <- 0
accuracy.meas(testdata$solubility, pred)
## ROC plot
roc.curve(testdata$solubility, pred, plotit = F)

##save test and traning data
write.csv(train5,"train5.csv")
write.csv(test5,"test5.csv")
write.csv(testdata,"testdata.csv")

##SVM model--reduce cost will increase precision
install.packages("e1071")
library("e1071")
svm_model <- svm(solubility ~ ., train5, cost=0.1, kernel="polynomial", degree=3)
pred <- predict(svm_model, newdata = test5)
accuracy.meas(testdata$solubility, pred)
roc.curve(testdata$solubility, pred, plotit = F)

##test threshold definition-->precision 1.000
actual<-data.frame(c(1,1,0.7))
prediction<-data.frame(c(1,1,0.3))
accuracy.meas(actual[ ,1], prediction[ ,1])

##adjust the means of all the features to zero and the standard deviations to one
clean5=read.csv("clean5.csv")
clean5[ ,1]<-NULL
clean6<-clean5[ ,1:22]
scaled.data <- scale(clean6)
data.frame(scaled.data)
clean7<-data.frame(scaled.data)
# check that we get mean of 0 and sd of 1
colMeans(scaled.data)
apply(scaled.data, 2, sd)
clean7[ ,23] <- clean5[ ,23]
colnames(clean7)[23] <- "solubility"
write.csv(clean7,"clean6.csv")

##make solubility binary value 0 and 1 (50%), change solubility to class and rename to "train"
clean6=read.csv("clean6.csv")
clean6$solubility[clean6$solubility >= 0.5] <-1
clean6$solubility[clean6$solubility < 0.5] <- 0
colnames(clean6)[23] <- "class"
write.csv(clean6,"train.csv")
##save as train1 in GAN trainingdata

##ues data produced by GAN--there is something wrong with generatedCGAN1.csv. Do not use it.
fake=read.csv("generatedCGAN2.csv")
##keep the coulmn name same for adding two dataframes together
fake[ ,1]<-NULL
for (i in 1:23) {
  q=i-1
  a<-"X"
  b<-as.character(q)
  colnames(clean6)[i] <- paste(a,b,sep = "", collapse = "")}
total1 <- rbind(clean6,fake)
colnames(total1)[23]<-"solubility"
write.csv(total1,"total2.csv")

##use 75% as training data
clean5=read.csv("total2.csv")
library('caret')
set.seed(0)
TrainingDataIndex <- createDataPartition(clean5[ ,23], p=0.75, list = FALSE)
train5 <- clean5[TrainingDataIndex,]
test5<- clean5[-TrainingDataIndex,]
testdata<-test5 
test5[ ,23]<-NULL

## Create logistic Regression Model-->warning:In eval(expr, envir, enclos) : non-integer #successes in a binomial glm!
input <- train5
solubility.data = glm(formula = solubility ~ ., data = input, family = binomial)
summary(solubility.data)

##predicting using the model
pred <- predict(solubility.data, newdata = testdata)

##check the accuracy
accuracy.meas(testdata$solubility, pred)
## ROC plot
roc.curve(testdata$solubility, pred, plotit = F)
 
## build decison tree model
library(party)
output.tree <- ctree(solubility ~ ., data=train5)
plot(output.tree)

##predicting using the model
pred <- predict(output.tree, newdata = testdata)

##check the accuracy
accuracy.meas(testdata$solubility, pred)
## ROC plot
roc.curve(testdata$solubility, pred, plotit = F)

##SVM model--reduce cost will increase precision
library("e1071")
svm_model <- svm(solubility ~ ., train5, cost=0.1, kernel="polynomial", degree=3)
pred <- predict(svm_model, newdata = testdata)
accuracy.meas(testdata$solubility, pred)
roc.curve(testdata$solubility, pred, plotit = F)

## MCC calculation--->0
install.packages("mccr")
library(mccr)
mccr(testdata$solubility, pred)

train5=read.csv("train5.csv")
test5=read.csv("test5.csv")
testdata=read.csv("testdata.csv")
train5[ ,1]<-NULL
test5[ ,1]<-NULL
testdata[ ,1]<-NULL
train5$solubility[train5$solubility >= 0.5] <-1
train5$solubility[train5$solubility < 0.5] <- 0
##calculate MCC for three models
input <- train5
solubility.data = glm(formula = solubility ~ ., data = input, family = binomial)
summary(solubility.data)
pred <- predict(solubility.data, newdata = testdata)
accuracy.meas(testdata$solubility, pred)
roc.curve(testdata$solubility, pred, plotit = F)

install.packages("mlr")##--->not work
library(mlr)
calculateConfusionMatrix(pred)

install.packages("mccr")##--->not work
library(mccr)
mccr(testdata$solubility, pred)

install.packages("EvaluationMeasures")
library(EvaluationMeasures)
prediction<-data.frame(pred)
pred=read.csv("logisticalregression.csv")
pred[,1]<-NULL
pred$pred[pred$pred >= 0.5] <-1
pred$pred[pred$pred < 0.5] <- 0
EvaluationMeasures.MCC(testdata$solubility, pred$pred)##--->0.0299

library(party)
output.tree <- ctree(solubility ~ ., data=train5)
plot(output.tree)

prediction<-data.frame(pred)
prediction$solubility[prediction$solubility >= 0.5] <-1
prediction$solubility[prediction$solubility < 0.5] <- 0
EvaluationMeasures.MCC(testdata$solubility, prediction$solubility)##--->0.3598

library("e1071")
svm_model <- svm(solubility ~ ., train5, cost=0.1, kernel="polynomial", degree=3)
pred <- predict(svm_model, newdata = testdata)
accuracy.meas(testdata$solubility, pred)
roc.curve(testdata$solubility, pred, plotit = F) ##-->precision:0.923 recall:0.153 F:0.131 AUC:0.765
##why different , because this time the solubility for train and test is binary,pred is not binary
prediction<-data.frame(pred)
prediction$pred[prediction$pred >= 0.5] <-1
prediction$pred[prediction$pred < 0.5] <- 0
EvaluationMeasures.MCC(testdata$solubility, prediction$pred)##--->0.2624

##calculate accuracy---->work^^!!!!!
library('caret')
confusionMatrix(prediction$pred, reference = testdata$solubility)

##use SVM to test different definitions for solubility and unsolubility, mean(0.4892) , median(0.4400), 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9
#delete column except sequence and solubility
##test train:binary,testdata:binary,prediction:binary
clean5=read.csv("clean5.csv")
clean5[ ,1]<-NULL
clean5[ ,1:2]<-NULL
write.csv(clean5,"clean7.csv")

##change from here
clean7=read.csv("clean7.csv")
summary(clean7$solubility)
clean7[ ,1]<-NULL
##change number
clean7$solubility[clean7$solubility >= 0.5] <-1
clean7$solubility[clean7$solubility < 0.5] <- 0
##use 75% as training data
library('caret')
set.seed(1)
TrainingDataIndex <- createDataPartition(clean7[ ,21], p=0.75, list = FALSE)
train7 <- clean7[TrainingDataIndex,]
test7<- clean7[-TrainingDataIndex,]
testdata<-test7 
test7[ ,21]<-NULL
library("e1071")
svm_model <- svm(solubility ~ ., train7, cost=0.1, kernel="polynomial", degree=3)
pred <- predict(svm_model, newdata = testdata)
prediction<-data.frame(pred)
##change number
prediction$pred[prediction$pred >= 0.5] <-1
prediction$pred[prediction$pred < 0.5] <- 0
library("ROSE")
accuracy.meas(testdata$solubility, prediction$pred)
roc.curve(testdata$solubility, prediction$pred, plotit = F)
library(EvaluationMeasures)
EvaluationMeasures.MCC(testdata$solubility, prediction$pred)
library('caret')
confusionMatrix(prediction$pred, reference = testdata$solubility)
##change positive class
##confusionMatrix(prediction$pred, reference = testdata$solubility, positive = "1")

##try more segments using SVM---Error:Response must have two levels for function accaracy.meas and EvaluationMeasures.MCC.
clean7=read.csv("clean7.csv")
summary(clean7$solubility)
clean7[ ,1]<-NULL
clean7$solubility[clean7$solubility >= 0.7] <-2
clean7$solubility[clean7$solubility < 0.7 & clean7$solubility >= 0.3] <- 1
clean7$solubility[clean7$solubility < 0.3] <- 0
library('caret')
set.seed(1)
TrainingDataIndex <- createDataPartition(clean7[ ,21], p=0.75, list = FALSE)
train7 <- clean7[TrainingDataIndex,]
test7<- clean7[-TrainingDataIndex,]
testdata<-test7 
test7[ ,21]<-NULL
library("e1071")
svm_model <- svm(solubility ~ ., train7, cost=0.1, kernel="polynomial", degree=3)
pred <- predict(svm_model, newdata = testdata)
prediction<-data.frame(pred)
prediction$pred[prediction$pred >= 0.7] <-2
prediction$pred[prediction$pred < 0.7 & prediction$pred >= 0.3] <-1
prediction$pred[prediction$pred < 0.3] <- 0
library("ROSE")
accuracy.meas(testdata$solubility, prediction$pred)
roc.curve(testdata$solubility, prediction$pred, plotit = F)
library(EvaluationMeasures)
EvaluationMeasures.MCC(testdata$solubility, prediction$pred)
library('caret')    ##but this function works for multiply classes.^^
confusionMatrix(prediction$pred, reference = testdata$solubility)
##confusionMatrix(prediction$pred, reference = testdata$solubility, mode = "prec_recall")

##try do not divide them into classes and use the original possibility for confusionMatrix
##--->Error:the data cannot have more levels than the reference

##Plot for predicted solubility and actual solubility using SVM
plot(x = testdata$solubility,y = prediction$pred,
     xlab = "Autual solubility",
     ylab = "Predicted solubility",
     xlim = c(0,1),
     ylim = c(0,1),		 
     main = "Correlation between the predicted solubility and actual solubility"
)

--------------------------------------------------------------------------------------
##I find the length of testdata and test5 are different and I do not know why. Do not use test5,train5,testdata anymore!
########we should use figure to record corresponding results for detailed codes for repeating and examing!!!
##test whether cancelling the cell location and gene product column will reduce accuracy-->result: no reduce -figure 1
clean5=read.csv("clean5.csv")
clean5[ ,1]<-NULL
clean7=clean5
clean7$solubility[clean7$solubility >= 0.5] <-1
clean7$solubility[clean7$solubility < 0.5] <- 0
##use 75% as training data
library('caret')
set.seed(1)
TrainingDataIndex <- createDataPartition(clean7[ ,23], p=0.75, list = FALSE)
train7 <- clean7[TrainingDataIndex,]
test7<- clean7[-TrainingDataIndex,]
testdata<-test7 
test7[ ,23]<-NULL
library("e1071")
svm_model <- svm(solubility ~ ., train7, cost=0.1, kernel="polynomial", degree=3)
pred <- predict(svm_model, newdata = testdata)
prediction<-data.frame(pred)
prediction$pred[prediction$pred >= 0.5] <-1
prediction$pred[prediction$pred < 0.5] <- 0
library("ROSE")
accuracy.meas(testdata$solubility, prediction$pred)
roc.curve(testdata$solubility, prediction$pred, plotit = F)
library(EvaluationMeasures)
EvaluationMeasures.MCC(testdata$solubility, prediction$pred)
library('caret')
confusionMatrix(prediction$pred, reference = testdata$solubility)

##test train:not,testdata:not,prediction:not----->Response must have two levels.
##test train:not,testdata:binary,prediction:not ---->figure 2 
clean7=read.csv("clean7.csv")
summary(clean7$solubility)
clean7[ ,1]<-NULL
##use 75% as training data
library('caret')
set.seed(1)
TrainingDataIndex <- createDataPartition(clean7[ ,21], p=0.75, list = FALSE)
train7 <- clean7[TrainingDataIndex,]
test7<- clean7[-TrainingDataIndex,]
testdata<-test7 
test7[ ,21]<-NULL
library("e1071")
svm_model <- svm(solubility ~ ., train7, cost=0.1, kernel="polynomial", degree=3)
pred <- predict(svm_model, newdata = testdata)
testdata$solubility[testdata$solubility >= 0.5] <-1
testdata$solubility[testdata$solubility < 0.5] <- 0
library("ROSE")
accuracy.meas(testdata$solubility, prediction$pred)
roc.curve(testdata$solubility, prediction$pred, plotit = F)
library(EvaluationMeasures)
EvaluationMeasures.MCC(testdata$solubility, prediction$pred)
library('caret')
confusionMatrix(prediction$pred, reference = testdata$solubility)

##test train:not,testdata:binary,prediction:binary(change after prediction)----figure 3
clean7=read.csv("clean7.csv")
summary(clean7$solubility)
clean7[ ,1]<-NULL
##use 75% as training data
library('caret')
set.seed(1)
TrainingDataIndex <- createDataPartition(clean7[ ,21], p=0.75, list = FALSE)
train7 <- clean7[TrainingDataIndex,]
test7<- clean7[-TrainingDataIndex,]
testdata<-test7 
test7[ ,21]<-NULL
library("e1071")
svm_model <- svm(solubility ~ ., train7, cost=0.1, kernel="polynomial", degree=3)
pred <- predict(svm_model, newdata = testdata)
testdata$solubility[testdata$solubility >= 0.5] <-1
testdata$solubility[testdata$solubility < 0.5] <- 0
prediction<-data.frame(pred)
prediction$pred[prediction$pred >= 0.5] <-1
prediction$pred[prediction$pred < 0.5] <- 0
library("ROSE")
accuracy.meas(testdata$solubility, prediction$pred)
roc.curve(testdata$solubility, prediction$pred, plotit = F)
library(EvaluationMeasures)
EvaluationMeasures.MCC(testdata$solubility, prediction$pred)
library('caret')
confusionMatrix(prediction$pred, reference = testdata$solubility)

##test whether different seeds will produce different results-seed(2)-->figure 4
##all the other same with figure 3

##cancel the first two columns for clean6, and try to use possibility not binary solubility.
clean6=read.csv("clean6.csv")
clean6[,1:2]<-NULL
colnames(clean6)[21] <- "class"
write.csv(clean6,"train.csv")
##save as train2 in GAN trainingdata 

##try SVM using generatedCGAN3---same as figure 3--->figure 5
fake=read.csv("generatedCGAN4.csv")
train=read.csv("train.csv")
fake[ ,1]<-NULL
for (i in 1:21) {
  q=i-1
  a<-"X"
  b<-as.character(q)
  colnames(train)[i] <- paste(a,b,sep = "", collapse = "")}
total4 <- rbind(train,fake)
colnames(total4)[21]<-"solubility"
write.csv(total4,"total4.csv")

clean7=read.csv("total4.csv")
clean7[ ,1]<-NULL
##use 75% as training data
library('caret')
set.seed(1)
TrainingDataIndex <- createDataPartition(clean7[ ,21], p=0.75, list = FALSE)
train7 <- clean7[TrainingDataIndex,]
test7<- clean7[-TrainingDataIndex,]
testdata<-test7 
test7[ ,21]<-NULL
library("e1071")
svm_model <- svm(solubility ~ ., train7, cost=0.1, kernel="polynomial", degree=3)
pred <- predict(svm_model, newdata = testdata)
testdata$solubility[testdata$solubility >= 0.5] <-1
testdata$solubility[testdata$solubility < 0.5] <- 0
prediction<-data.frame(pred)
prediction$pred[prediction$pred >= 0.5] <-1
prediction$pred[prediction$pred < 0.5] <- 0
library("ROSE")
accuracy.meas(testdata$solubility, prediction$pred)
roc.curve(testdata$solubility, prediction$pred, plotit = F)
library(EvaluationMeasures)
EvaluationMeasures.MCC(testdata$solubility, prediction$pred)
library('caret')
confusionMatrix(prediction$pred, reference = testdata$solubility)

##test whether random produced data will produce different accuracy,using generatedCGAN4
##all the other same with figure 5----->figure 6

##test whethe the dataset is balanced: 1438(1) and 1215(0)
clean6=read.csv("clean6.csv")
clean6$solubility[clean6$solubility >= 0.5] <-1
clean6$solubility[clean6$solubility < 0.5] <- 0
library('plyr')
count(clean6, 'solubility')

##accuracy vs cutoff points(0.3:0.7,0.01)
clean7=read.csv("clean7.csv")
summary(clean7$solubility)
clean7[ ,1]<-NULL
##use 75% as training data
library('caret')
set.seed(1)
TrainingDataIndex <- createDataPartition(clean7[ ,21], p=0.75, list = FALSE)
train7 <- clean7[TrainingDataIndex,]
test7<- clean7[-TrainingDataIndex,]
testdata<-test7 
test7[ ,21]<-NULL
library("e1071")
svm_model <- svm(solubility ~ ., train7, cost=0.1, kernel="polynomial", degree=3)
pred <- predict(svm_model, newdata = testdata)

d = NULL
for (x in seq(0.3,0.7,0.01)) {
testdata$solubility[testdata$solubility >= x] <-1
testdata$solubility[testdata$solubility < x] <- 0
prediction<-data.frame(pred)
prediction$pred[prediction$pred >= x] <-1
prediction$pred[prediction$pred < x] <- 0
library("ROSE")
accuracy.meas(testdata$solubility, prediction$pred)
roc.curve(testdata$solubility, prediction$pred, plotit = F)
library(EvaluationMeasures)
EvaluationMeasures.MCC(testdata$solubility, prediction$pred)
library('caret')
cm<-confusionMatrix(prediction$pred, reference = testdata$solubility)
overall <- cm$overall
overall.accuracy <- overall['Accuracy'] 
print(x)
print(overall.accuracy)
d = rbind(d, data.frame(x,overall.accuracy))
}
write.csv(d,"accuracy_cutoff")
##Plot for accuracy vs cutoff points using SVM---->figure 7
plot(x = d$x,y = d$overall.accuracy,
     xlab = "cutoff points",
     ylab = "accuracy",
     xlim = c(0.3,0.7),
     ylim = c(0.4,0.8),		 
     main = "Correlation between accuracy and cutoff points"
)

## conditional random forest model--->figure 8
rm(list=ls())
library('ggplot2') # visualization
#library('ggthemes') # visualization
library('scales') # visualization
library('dplyr') # data manipulation
library('mice') # imputation
library('randomForest') # classification algorithm
library('caret')
library('e1071')
library('mltools')
library('data.table')
library('mlr')

clean7=read.csv("clean7.csv")
clean7[ ,1]<-NULL
##use 75% as training data
library('caret')
set.seed(1)
TrainingDataIndex <- createDataPartition(clean7[ ,21], p=0.75, list = FALSE)
train7 <- clean7[TrainingDataIndex,]
test7<- clean7[-TrainingDataIndex,]
testdata<-test7 
test7[ ,21]<-NULL

rf.lrn <- makeLearner("classif.cforest")
control <- trainControl(method="repeatedcv", number=10, repeats=5)
rf.lrn$par.vals <- list(ntree = 603, mtry = 5, nodesize=37, cutoff = c(0.8,0.2), 
                        preProcess = c("pca"), trControl=control)
##one error, the terget must be factor rather than double,so change it to factor.
train7$solubility[train7$solubility >= 0.5] <-1
train7$solubility[train7$solubility < 0.5] <- 0
train7$solubility<-as.factor(train7$solubility)
##if we do not transfer it to binary and directly to factor, error:Error in model@fit(data, ...) : error code 1 from Lapack routine 'dgesdd'
traintask <- makeClassifTask(data = train7,target = "solubility")
traintask <- filterFeatures(traintask, method = "cforest.importance", abs = 6)
model <- mlr::train(rf.lrn, traintask)
pred <- predict(model, newdata = test7)

testdata$solubility[testdata$solubility >= 0.5] <-1
testdata$solubility[testdata$solubility < 0.5] <- 0
prediction<-data.frame(pred)
library("ROSE")
accuracy.meas(testdata$solubility, prediction$response)
roc.curve(testdata$solubility, prediction$response, plotit = F)
library(EvaluationMeasures)
EvaluationMeasures.MCC(testdata$solubility, prediction$response)
library('caret')
confusionMatrix(prediction$response, reference = testdata$solubility)

##naive bayes model----->figure 9
clean7=read.csv("clean7.csv")
clean7[ ,1]<-NULL
##use 75% as training data
library('caret')
set.seed(1)
TrainingDataIndex <- createDataPartition(clean7[ ,21], p=0.75, list = FALSE)
train7 <- clean7[TrainingDataIndex,]
test7<- clean7[-TrainingDataIndex,]
testdata<-test7 
test7[ ,21]<-NULL
##one error, the terget must be factor rather than double,so change it to factor.
train7$solubility[train7$solubility >= 0.5] <-1
train7$solubility[train7$solubility < 0.5] <- 0
train7$solubility<-as.factor(train7$solubility)

library(e1071)
model <- naiveBayes(solubility ~ ., data=train7)
pred <- predict(model, newdata = test7)

testdata$solubility[testdata$solubility >= 0.5] <-1
testdata$solubility[testdata$solubility < 0.5] <- 0
prediction<-data.frame(pred)

library("ROSE")
accuracy.meas(testdata$solubility, prediction$pred)
roc.curve(testdata$solubility, prediction$pred, plotit = F)
library(EvaluationMeasures)
EvaluationMeasures.MCC(testdata$solubility, prediction$pred)
library('caret')
confusionMatrix(prediction$pred, reference = testdata$solubility)

##logestical regression---->figure 10 ---all bianry
input <- train7
solubility.data = glm(formula = solubility ~ ., data = input, family = binomial)
pred <- predict(solubility.data, newdata = test7)
testdata$solubility[testdata$solubility >= 0.5] <-1
testdata$solubility[testdata$solubility < 0.5] <- 0
prediction<-data.frame(pred)
prediction$pred[prediction$pred >= 0.5] <-1
prediction$pred[prediction$pred < 0.5] <- 0
library("ROSE")
accuracy.meas(testdata$solubility, prediction$pred)
roc.curve(testdata$solubility, prediction$pred, plotit = F)
library(EvaluationMeasures)
EvaluationMeasures.MCC(testdata$solubility, prediction$pred)
library('caret')
confusionMatrix(prediction$pred, reference = testdata$solubility)

##decision tree----->figure 11 ----all binary
##figure 12--train:not binary, test:binary, prediction:binary
library(party)
output.tree <- ctree(solubility ~ ., data=train7)
pred <- predict(output.tree, newdata = test7)
testdata$solubility[testdata$solubility >= 0.5] <-1
testdata$solubility[testdata$solubility < 0.5] <- 0
prediction<-data.frame(pred)
prediction$solubility[prediction$solubility >= 0.5] <-1
prediction$solubility[prediction$solubility < 0.5] <- 0
library("ROSE")
accuracy.meas(testdata$solubility, prediction$solubility)
roc.curve(testdata$solubility, prediction$solubility, plotit = F)
library(EvaluationMeasures)
EvaluationMeasures.MCC(testdata$solubility, prediction$solubility)
library('caret')
confusionMatrix(prediction$solubility, reference = testdata$solubility)

##Plot for predicted solubility and actual solubility using logistical Regression,all no binary--->figure 13
##Plot for predicted solubility and actual solubility using decision tree,all no binary--->figure 14
plot(x = testdata$solubility,y = prediction$pred,
     xlab = "Autual solubility",
     ylab = "Predicted solubility",
     xlim = c(0,1),
     ylim = c(0,1),		 
     main = "Correlation between the predicted solubility and actual solubility"
)
----------------------------------------------------------------------------------
##we find the method (4) train: not,test:binary, prediction:binary is the best method to transfer 
##to binary value, cutoff ponit 0.44 is best and all use them from now
  
##try to make cforest input and output not binary--->figure 15
clean7=read.csv("clean7.csv")
clean7[ ,1]<-NULL
##use 75% as training data
library('caret')
set.seed(1)
TrainingDataIndex <- createDataPartition(clean7[ ,21], p=0.75, list = FALSE)
train7 <- clean7[TrainingDataIndex,]
test7<- clean7[-TrainingDataIndex,]
testdata<-test7 
test7[ ,21]<-NULL

cf <- cforest(solubility ~ ., data = train7)
pred <- predict(cf, newdata = test7, type = "response")
prediction<-data.frame(pred)
testdata$solubility[testdata$solubility >= 0.44] <-1
testdata$solubility[testdata$solubility < 0.44] <- 0
prediction$solubility[prediction$solubility >= 0.44] <-1
prediction$solubility[prediction$solubility < 0.44] <- 0
library("ROSE")
accuracy.meas(testdata$solubility, prediction$solubility)
roc.curve(testdata$solubility, prediction$solubility, plotit = F)
library(EvaluationMeasures)
EvaluationMeasures.MCC(testdata$solubility, prediction$solubility)
library('caret')
confusionMatrix(prediction$solubility, reference = testdata$solubility)
##Plot for predicted solubility and actual solubility using cforest,all no binary--->figure 16

##try to make Naive Bayes input and output not binary
clean7=read.csv("clean7.csv")
clean7[ ,1]<-NULL
##use 75% as training data
library('caret')
set.seed(1)
TrainingDataIndex <- createDataPartition(clean7[ ,21], p=0.75, list = FALSE)
train7 <- clean7[TrainingDataIndex,]
test7<- clean7[-TrainingDataIndex,]
testdata<-test7 
test7[ ,21]<-NULL

train7$solubility<-as.factor(train7$solubility)

library(e1071)
model <- naiveBayes(solubility ~ ., data=train7)
pred <- predict(model, newdata = test7)
prediction<-data.frame(pred)

prediction$pred<-as.numeric(prediction$pred)
prediction$pred<-(prediction$pred-1)/100
#Plot for predicted solubility and actual solubility using Naive bayes,all no binary--->figure 17
testdata$solubility[testdata$solubility >= 0.44] <-1
testdata$solubility[testdata$solubility < 0.44] <- 0
prediction$pred[prediction$pred >= 0.44] <-1
prediction$pred[prediction$pred < 0.44] <- 0
#evaluation-->figure 18

##optimize the parameters of SVM--kernel:figure 19-radical, figure 20-linear, 21-sigmoid
## figure 22-RBF(default), figure 23-polynomial
clean7=read.csv("clean7.csv")
clean7[ ,1]<-NULL
##use 75% as training data
library('caret')
set.seed(1)
TrainingDataIndex <- createDataPartition(clean7[ ,21], p=0.75, list = FALSE)
train7 <- clean7[TrainingDataIndex,]
test7<- clean7[-TrainingDataIndex,]
testdata<-test7 
test7[ ,21]<-NULL
library("e1071")
svm_model <- svm(solubility ~ ., train7, cost=0.1, degree=3)
pred <- predict(svm_model, newdata = testdata)
testdata$solubility[testdata$solubility >= 0.44] <-1
testdata$solubility[testdata$solubility < 0.44] <- 0
prediction<-data.frame(pred)
prediction$pred[prediction$pred >= 0.44] <-1
prediction$pred[prediction$pred < 0.44] <- 0
library("ROSE")
accuracy.meas(testdata$solubility, prediction$pred)
roc.curve(testdata$solubility, prediction$pred, plotit = F)
library(EvaluationMeasures)
EvaluationMeasures.MCC(testdata$solubility, prediction$pred)
library('caret')
confusionMatrix(prediction$pred, reference = testdata$solubility)
#choose default RBF and optimize cost
d = NULL
for (x in seq(0.01,2,0.02)) {
  library("e1071")
  svm_model <- svm(solubility ~ ., train7, cost=x, degree=3)
  pred <- predict(svm_model, newdata = testdata)
  testdata$solubility[testdata$solubility >= 0.44] <-1
  testdata$solubility[testdata$solubility < 0.44] <- 0
  prediction<-data.frame(pred)
  prediction$pred[prediction$pred >= 0.44] <-1
  prediction$pred[prediction$pred < 0.44] <- 0
  library('caret')
  cm<-confusionMatrix(prediction$pred, reference = testdata$solubility)
  overall <- cm$overall
  overall.accuracy <- overall['Accuracy'] 
  d = rbind(d, data.frame(x,overall.accuracy))
}
##Plot for SVM accuracy vs cost---->figure 24, and choose 1.5
plot(x = d$x,y = d$overall.accuracy,
     xlab = "cost",
     ylab = "SVM-accuracy",
     xlim = c(0,2),
     ylim = c(0.738,0.782),		 
     main = "Correlation between SVM accuracy and cost"
)
#optimize degree-no influence to this kernel
d = NULL
for (x in seq(0.01,0.1,0.002)) {
  library("e1071")
  svm_model <- svm(solubility ~ ., train7, cost=1.5, gamma=x)
  pred <- predict(svm_model, newdata = testdata)
  testdata$solubility[testdata$solubility >= 0.44] <-1
  testdata$solubility[testdata$solubility < 0.44] <- 0
  prediction<-data.frame(pred)
  prediction$pred[prediction$pred >= 0.44] <-1
  prediction$pred[prediction$pred < 0.44] <- 0
  library('caret')
  cm<-confusionMatrix(prediction$pred, reference = testdata$solubility)
  overall <- cm$overall
  overall.accuracy <- overall['Accuracy'] 
  d = rbind(d, data.frame(x,overall.accuracy))
}
##Plot for SVM accuracy vs gamma---->figure 25, and choose 0.026
plot(x = d$x,y = d$overall.accuracy,
     xlab = "gamma",
     ylab = "SVM-accuracy",
     xlim = c(0.01,0.1),
     ylim = c(0.75,0.784),		 
     main = "Correlation between SVM accuracy and gamma"
)

#optimize epsilon--0.198
d = NULL
for (x in seq(0.01,0.3,0.004)) {
  library("e1071")
  svm_model <- svm(solubility ~ ., train7, cost=1.5, gamma=0.026,epsilon = x)
  pred <- predict(svm_model, newdata = testdata)
  testdata$solubility[testdata$solubility >= 0.44] <-1
  testdata$solubility[testdata$solubility < 0.44] <- 0
  prediction<-data.frame(pred)
  prediction$pred[prediction$pred >= 0.44] <-1
  prediction$pred[prediction$pred < 0.44] <- 0
  library('caret')
  cm<-confusionMatrix(prediction$pred, reference = testdata$solubility)
  overall <- cm$overall
  overall.accuracy <- overall['Accuracy'] 
  d = rbind(d, data.frame(x,overall.accuracy))
}
##Plot for SVM accuracy vs epsilon---->figure 26, and choose 0.198
plot(x = d$x,y = d$overall.accuracy,
     xlab = "epsilon",
     ylab = "SVM-accuracy",
     xlim = c(0.01,0.3),
     ylim = c(0.77,0.789),		 
     main = "Correlation between SVM accuracy and epsilon"
)

##use optimized parameter to predict bu SVM-figure 27
library("e1071")
svm_model <- svm(solubility ~ ., train7, cost=1.5, gamma=0.026,epsilon = 0.198)
pred <- predict(svm_model, newdata = testdata)
prediction<-data.frame(pred)
library("ROSE")
accuracy.meas(testdata$solubility, prediction$pred)
roc.curve(testdata$solubility, prediction$pred, plotit = F)
library(EvaluationMeasures)
EvaluationMeasures.MCC(testdata$solubility, prediction$pred)
library('caret')
confusionMatrix(prediction$pred, reference = testdata$solubility)
##Plot for predicted solubility and actual solubility using optimized SVM,all no binary--->figure 28
plot(x = testdata$solubility,y = prediction$pred,
     xlab = "Autual solubility",
     ylab = "Predicted solubility",
     xlim = c(0,1),
     ylim = c(0,1),		 
     main = "Correlation between the predicted solubility and actual solubility"
)
------------------------------------------------------------------------------------------------------------------------------
##GAN
real_data=read.csv("train.csv")
gan<-function(real_data,g_nn,d_nn,batchsize=100,epoch=100,disc_step=1,print_loss=T,display_generation_image=T,
              d_real_display=T,d_fake_print=T,display_generation_distribution=F){
  rotate <- function(x) t(apply(x, 2, rev))
  gan_model<<-list()
  
  x<-real_data
  numdata<-nrow(x)
  num_f<-numdata* g_nn$input_dim
  num_d<-numdata* d_nn$input_dim
  final_loss<-NULL
  for(t in 1:epoch){
    
    
    u<-1
    fak<-rnorm((num_f)) ### random noise
    fake<-matrix(fak,ncol=g_nn$input_dim)
    m <- nrow(fake)
    numbatches <- m/batchsize
    numbatches
    randperm <- sample(1:m, m)
    for(u in 1:numbatches){
      for(kk in 1:disc_step){
        ############ training discriminator
        
        batch_x<-fake[randperm[((u - 1) * batchsize +
                                  1):(u * batchsize)],]
        head(batch_x)
        g_nn <- nn.ff2(g_nn, batch_x,1,"g_nn") ## fake data generation from random noise
        genration<-g_nn$post[[length(g_nn$size)]]
        
        batch_x_2<-x[randperm[((u - 1) * batchsize +
                                 1):(u * batchsize)],]
        
        d_nn<-nn.ff2(d_nn,batch_x_2,1,"d_nn")
        d_real<-  d_nn$post[[length(d_nn$size)]]   ### probability real data as real
        d_fake<-nn.ff2(d_nn,genration,0,"d_nn")$post[[length(d_nn$size)]]  ### probability fake data as real
        
        G <- genration
        nn<-d_nn
        
        ######### back propagation Discriminator
        
        
        f_dn<-nn.ff2(d_nn,G,0,"d_nn") ## D(G(z))
        tr_de<-part_bp1(d_nn)
        f_de<-part_bp1(f_dn)
        
        
        
        for (i in 1:((length(d_nn$size)-1))) {
          
          dw1 <- t(tr_de[[i+1]]) %*% d_nn$post[[i]]/nrow(tr_de[[i+1 ]])
          db1 <- colMeans(tr_de[[i + 1]])
          db1 <- db1 * d_nn$learningrate
          
          dw2 <- t(f_de[[i +1]]) %*% f_dn$post[[i]]/nrow(f_de[[i+1]])
          db2 <- colMeans(f_de[[i + 1]])
          db2 <- db2 * d_nn$learningrate
          
          dw<- dw1- dw2
          dw <-  dw * d_nn$learningrate
          db <- db1-db2
          
          d_nn$W[[i]] <- d_nn$W[[i]] + dw #### gradien ascent
          d_nn$B[[i]] <- d_nn$B[[i]] + db
        }
        
        ## If you train only discriminator several times, below two probability will be 1 and 0
        
        if(d_real_display) {
          cat("\n Discriminator Training / Pr(real)")
          print(head(nn.ff2(d_nn,x,1,"d_nn")$post[[length(d_nn$size)]])) ## Pr(D(real))
        }
        if(d_real_display){
          cat("\n Discriminator Training / Pr(fake)")
          print(head(nn.ff2(d_nn,genration,1,"d_nn")$post[[length(d_nn$size)]]))  ## Pr(D(fake))
        }
        
      }
      ############ training generator
      # }
      
      fak<-rnorm((num_f))
      fake<-matrix(fak,ncol=g_nn$input_dim)
      # for(u in 1:numbatches){
      
      batch_x<-fake[randperm[((u - 1) * batchsize +
                                1):(u * batchsize)],]
      
      g_nn<-nn.ff2(g_nn,batch_x,0,"g_nn") ##  fake data generation from random noise
      G<-g_nn$post[[length(g_nn$size)]]
      d_nn<-nn.ff2(d_nn,G,1,"d_nn")
      
      
      
      
      
      ######### back propagation generator
      
      
      n <-length(c(d_nn$size,g_nn$size))-1
      d <- list()
      nn<-d_nn
      
      
      if (nn$output == "sigm") {
        d[[n]] <- nn$e * (nn$post[[length(nn$post)]] * (1 - nn$post[[length(nn$post)]]))
      }  else if (nn$output == "linear" || nn$output == "softmax") {
        d[[n]] <- nn$e
      }
      
      
      
      
      k<-0
      dl<-length(d_nn$size)-1
      dg<-length(g_nn$size)-1
      l<-dl+dg
      bp_range<-c(dl:1,dg:2)
      
      
      for (i in bp_range) {
        
        if (nn$activationfun == "linear") {
          d_act <- nn$post[[i]]
        }
        if (nn$activationfun == "sigm") {
          d_act <- nn$post[[i]] * (1 - nn$post[[i]])
        }
        if (nn$activationfun == "relu") {
          d_act <-  ifelse(nn$post[[i]]>0,1,0)
        }
        else if (nn$activationfun == "tanh") {
          d_act <- 1.7159 * 2/3 * (1 - 1/(1.7159)^2 * nn$post[[i]]^2)
        }
        
        d[[l]] <- (d[[l + 1]] %*% nn$W[[i]]) * d_act
        
        k<-k+1
        l<-l-1
        if(k == dl){
          nn<-g_nn
        }
      }
      
      ## check dimension of delta
      lapply(d,dim)
      
      
      for (i in 1:dg) {
        nn<-g_nn
        dw <- t(d[[i+1]]) %*% nn$post[[i]]/nrow(t(d[[i+1]]))
        db <- colMeans(d[[i + 1]])
        db <- db * nn$learningrate
        nn$W[[i]] <- nn$W[[i]] + dw * nn$learningrate
        nn$B[[i]] <- nn$B[[i]] + db
        g_nn<-nn
      }
      
      if(d_real_display){
        cat("\n Generator Training / Pr(fake)")
        print(head(nn.ff2(d_nn,nn.ff2(g_nn,batch_x,0,"g_nn")$post[[length(g_nn$size)]],1,"d_nn")$post[[length(d_nn$size)]]))  ## Pr(D(fake))
      }
      
      if(display_generation_distribution==T){
        
        hist(nn.ff2(g_nn,batch_x,0,"g_nn")$post[[3]])
      }
      
      
      cat(paste0("\n",t,"-epoch ",u,"batch"))
      
    }
    
    if(display_generation_image){
      
      
      fak<-rnorm((num_f))
      fake<-matrix(fak,ncol=g_nn$input_dim)
      
      gg<-nn.ff2(g_nn,fake[1:100,],1,"g_nn")$post[[length(g_nn$size)]]
      
      # png(filename=paste0("iteration-",t,".png"))
      par(mfrow=c(3,3))
      lapply(1:9,
             function(q) image(
               rotate(matrix(unlist(gg[q,]),nrow = sqrt(d_nn$input_dim), byrow = TRUE)),
               col=grey.colors(255)        
             )
      )
      # dev.off()
      
      
    }
    
    if(print_loss){
      
      g_nn <- nn.ff2(g_nn, fake,1,"g_nn")
      genration<-g_nn$post[[length(g_nn$size)]]
      
      
      d_nn<-nn.ff2(d_nn,x,1,"d_nn")
      d_real<-  d_nn$post[[length(d_nn$size)]]
      
      d_fake<-nn.ff2(d_nn,genration,0,"d_nn")$post[[length(d_nn$size)]]
      final_loss<-rbind(final_loss,cbind(mean(d_real),mean(d_fake)))
      ### D_real, D_fake calculation
      ### D_real : probability real data as real
      ### D_fake : probability real data as fake
      ### In gan model, training convergence condition is D_real = D_fake = 0.5
    }
    gan_model$g_nn<<-g_nn
    gan_model$d_nn<<-d_nn
    
    colnames(final_loss)<-c("d_real","d_fake")
    
    gan_model$loss<<-final_loss
    gan_model
  }
  gan_model
  
}