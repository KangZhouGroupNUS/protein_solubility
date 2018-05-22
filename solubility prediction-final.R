#prepare dataset for solubility prediction

#find 90 proteins without gene name in sequence files
library("readxl")
traningdata <- read_excel("all_data.xlsx")
traningdata[ ,1:3] <- NULL
traningdata[ ,2:3] <- NULL
traningdata[ ,3:23] <- NULL
x<-na.omit(traningdata)
sequence <- read_excel("seqinfo.xlsx")
colnames(x)[1] <- "gene"

library(MASS)
merged.train <- merge(x = x, y = sequence,
                      by.x = ("gene"),
                      by.y = ("gene")
)

merged.train[merged.train=='NA'] <- NA
new_DF <- merged.train[is.na(merged.train$sequence),]
new_DF[,2]<-NULL
write.csv(new_DF,"missing_sequence.csv")


#add missing sequence
library("readxl")
traningdata <- read_excel("all_data.xlsx")
traningdata[ ,1:3] <- NULL
traningdata[ ,2:3] <- NULL
traningdata[ ,3:23] <- NULL
x<-na.omit(traningdata)
sequence <- read_excel("seqinfo.xlsx")
colnames(x)[1] <- "gene"

no_sequence=read.csv("no_sequence1.csv")
no_sequence$sequence <- as.character(no_sequence$sequence)
myfunction <- function(y){
  w=gsub(" ", "", y, fixed = TRUE)
  return(gsub("\n", "", w, fixed = TRUE))
}
for (i in seq(1,519,1)) { 
  no_sequence[i,2]<-myfunction(no_sequence[i,2])
}
class(no_sequence$sequence)
write.csv(no_sequence,"no_sequence5.csv")
data <- read.csv("no_sequence5.csv",na.strings=c(""," ","NA"),header=TRUE)
summary(data)
new_DF3 <- data[is.na(data$sequence),]
new_DF4 <- data[!(is.na(data$sequence)),]
new_DF5 <- new_DF4[!(new_DF4$gene=="yjgX"),]
new_DF6 <- new_DF5[!(new_DF5$gene=="gatR"),]

library(MASS)
merged.train3 <- merge(x = x, y = new_DF6,
                      by.x = ("gene"),
                      by.y = ("gene")
)
merged.train3[,3]<-NULL
total <- rbind(new_DF, merged.train3)
#Total proteins 3152=3173-21(17+4)
write.csv(total,"solubility1.csv")

total[,1]<-NULL
total[ ,1][total[ ,1]>100]<-100
total[ ,1]<-total[ ,1]/100
colnames(total)[1] <- "solubility"
total[,3]<-total[,1]
total[,1]<-NULL
x<-na.omit(total)
x[ ,1]<-as.character(x[ ,1])
length(x[ ,1])
y=x[1:3152,1]
y = t(sapply(y, extractAAC))
y = y[(sapply(y, protcheck))]
length(y)
y<-as.data.frame(y)
a=x[,1]
b=y[,1]
ab=a[!(a %in% b)]
dismatch<-as.data.frame(ab)

total=read.csv("solubility1.csv")
total[,1]<-NULL
total <- total[!(total$gene=="ybfH"),]
total <- total[!(total$gene=="yhdW"),]
total <- total[!(total$gene=="ymgH"),]
total <- total[!(total$gene=="yohG"),]
#total proteins 3148(4 proteins include undertemined amino acid X)
write.csv(total,"solubility2.csv")

total=read.csv("solubility2.csv")
total[,1]<-NULL
total[,1]<-NULL
total[ ,1][total[ ,1]>100]<-100
total[ ,1]<-total[ ,1]/100
colnames(total)[1] <- "solubility"
total[,3]<-total[,1]
total[,1]<-NULL
x<-na.omit(total)
x[ ,1]<-as.character(x[ ,1])
length(x[ ,1])
library("protr")
x[ ,1] = x[ ,1][(sapply(x[ ,1], protcheck))]
length(x[ ,1])
x1 = t(sapply(x[ ,1], extractAAC))
x[ ,3:22] <- x1[ ,1:20]
x[ ,1]<-NULL
x[ ,22]<-x[ ,1]
x[ ,1]<-NULL
colnames(x)[21] <- "solubility"
write.csv(x,"cleandata1.csv")


#cut-off point optimization
clean7=read.csv("cleandata1.csv")
clean7[ ,1]<-NULL
library('caret')
set.seed(5)
TrainingDataIndex <- createDataPartition(clean7[ ,21], p=0.9, list = FALSE)
train7 <- clean7[TrainingDataIndex,]
test7<- clean7[-TrainingDataIndex,]
testdata<-test7 
test7[ ,21]<-NULL
library("e1071")
svm_model <- svm(solubility ~ ., train7, method="eps-regression", cost=1.5, gamma=0.026, epsilon = 0.198, cross=10)
pred <- predict(svm_model, newdata = test7)
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
write.csv(d,"accuracy_cutoff1.csv")

##Plot for accuracy vs cutoff points using SVM
plot(x = d$x,y = d$overall.accuracy,
     xlab = "cutoff points",
     ylab = "accuracy",
     xlim = c(0.3,0.7),
     ylim = c(0.4,0.8),		 
     main = "Correlation between accuracy and cutoff points"
)
# when cut off point is 0.36 or 0.39, accuracy is 0.7628
#
yourData=clean7
yourData<-yourData[sample(nrow(yourData)),]
folds <- cut(seq(1,nrow(yourData)),breaks=10,labels=FALSE)


## SVM tuning
clean7=read.csv("cleandata1.csv")
clean7[ ,1]<-NULL
library('caret')
set.seed(5)
TrainingDataIndex <- createDataPartition(clean7[ ,21], p=0.9, list = FALSE)
train7 <- clean7[TrainingDataIndex,]
test7<- clean7[-TrainingDataIndex,]

d = NULL
for (x in seq(0.02,2,0.02)) {
  test18=test7
  testdata<-test18 
  test18[ ,21]<-NULL
  library("e1071")
  svm_model <- svm(solubility ~ ., train7, method="eps-regression", cost=x, cross=10)
  pred <- predict(svm_model, newdata = test18)
  prediction<-data.frame(pred)
  rsq <- function(x, y) summary(lm(y~x))$r.squared
  R2=rsq(testdata$solubility, prediction[,1])
  
  testdata$solubility[testdata$solubility >= 0.39] <-1
  testdata$solubility[testdata$solubility < 0.39] <- 0
  prediction$pred[prediction$pred >= 0.39] <-1
  prediction$pred[prediction$pred < 0.39] <- 0
  library('caret')
  cm<-confusionMatrix(prediction$pred, reference = testdata$solubility)
  overall <- cm$overall
  overall.accuracy <- overall['Accuracy'] 
  d = rbind(d, data.frame(x,R2,overall.accuracy))
}
write.csv(d,"cost2.csv")

plot(x = d$x,y = d$overall.accuracy,
     xlab = "cost",
     ylab = "SVM-accuracy",
     xlim = c(0.02,2),
     ylim = c(0.717,0.758),		 
     main = "Correlation between SVM accuracy and cost"
)
d1 = NULL
for (y in seq(0.001,0.1,0.001)) {
  test18=test7
  testdata<-test18 
  test18[ ,21]<-NULL
  library("e1071")
  svm_model <- svm(solubility ~ ., train7, method="eps-regression", cost=1.26, gamma=y, cross=10)
  pred <- predict(svm_model, newdata = test18)
  prediction<-data.frame(pred)
  rsq <- function(x, y) summary(lm(y~x))$r.squared
  R2=rsq(testdata$solubility, prediction[,1])
  
  testdata$solubility[testdata$solubility >= 0.39] <-1
  testdata$solubility[testdata$solubility < 0.39] <- 0
  prediction$pred[prediction$pred >= 0.39] <-1
  prediction$pred[prediction$pred < 0.39] <- 0
  library('caret')
  cm<-confusionMatrix(prediction$pred, reference = testdata$solubility)
  overall <- cm$overall
  overall.accuracy <- overall['Accuracy'] 
  d1 = rbind(d1, data.frame(y,R2,overall.accuracy))
}
write.csv(d1,"gamma2.csv")
plot(x = d1$y,y = d1$overall.accuracy,
     xlab = "gamma",
     ylab = "SVM-accuracy",
     xlim = c(0.001,0.1),
     ylim = c(0.692,0.758),		 
     main = "Correlation between SVM accuracy and gamma"
)
d2 = NULL
for (z in seq(0.001,0.3,0.004)) {
  test18=test7
  testdata<-test18 
  test18[ ,21]<-NULL
  library("e1071")
  svm_model <- svm(solubility ~ ., train7, method="eps-regression", cost=1.26, gamma=0.029, epsilon = z, cross=10)
  pred <- predict(svm_model, newdata = test18)
  prediction<-data.frame(pred)
  rsq <- function(x, y) summary(lm(y~x))$r.squared
  R2=rsq(testdata$solubility, prediction[,1])
  
  testdata$solubility[testdata$solubility >= 0.39] <-1
  testdata$solubility[testdata$solubility < 0.39] <- 0
  prediction$pred[prediction$pred >= 0.39] <-1
  prediction$pred[prediction$pred < 0.39] <- 0
  library('caret')
  cm<-confusionMatrix(prediction$pred, reference = testdata$solubility)
  overall <- cm$overall
  overall.accuracy <- overall['Accuracy'] 
  d2 = rbind(d2, data.frame(z,R2,overall.accuracy))
}
write.csv(d2,"epsilon2.csv")
plot(x = d2$z,y = d2$overall.accuracy,
     xlab = "epsilon",
     ylab = "SVM-accuracy",
     xlim = c(0.001,0.3),
     ylim = c(0.746,0.759),		 
     main = "Correlation between SVM accuracy and epsilon"
)

#binary accuracy 0.7001
clean7=read.csv("cleandata1.csv")
clean7[ ,1]<-NULL
clean7$solubility[clean7$solubility >= 0.39] <-1
clean7$solubility[clean7$solubility < 0.39] <-0

library('caret')
set.seed(5)
TrainingDataIndex <- createDataPartition(clean7[ ,21], p=0.9, list = FALSE)
train7 <- clean7[TrainingDataIndex,]
test7<- clean7[-TrainingDataIndex,]
testdata<-test7 
test7[ ,21]<-NULL
library("e1071")
svm_model <- svm(solubility ~ ., train7, method="eps-regression", cost=1.5, gamma=0.026, epsilon = 0.198, cross=10)
pred <- predict(svm_model, newdata = test7)
prediction1<-data.frame(pred)
prediction1$pred[prediction1$pred >= 0.39] <-1
prediction1$pred[prediction1$pred < 0.39] <- 0
library("ROSE")
accuracy.meas(testdata$solubility, prediction1[,1])
roc.curve(testdata$solubility, prediction1[,1], plotit = F)
library(EvaluationMeasures)
EvaluationMeasures.MCC(testdata$solubility, prediction1[,1])
library('caret')
confusionMatrix(prediction1[,1], reference = testdata$solubility)

# use 0.3 and 0.7 to produce binary solubility and 0.5 as cut-off points  accuracy 0.8117 
clean7=read.csv("cleandata1.csv")
clean7[ ,1]<-NULL
clean7 <- clean7[!(clean7$solubility>0.3 & clean7$solubility<0.7),]
clean7$solubility[clean7$solubility >= 0.7] <-1
clean7$solubility[clean7$solubility <= 0.3] <-0
soluable<-clean7$solubility[clean7$solubility == 1] 
soluable<-as.data.frame(soluable)
library('caret')
set.seed(5)
TrainingDataIndex <- createDataPartition(clean7[ ,21], p=0.9, list = FALSE)
train7 <- clean7[TrainingDataIndex,]
test7<- clean7[-TrainingDataIndex,]
testdata<-test7 
test7[ ,21]<-NULL
library("e1071")
svm_model <- svm(solubility ~ ., train7, method="eps-regression", cost=1.5, gamma=0.026, epsilon = 0.198, cross=10)
pred <- predict(svm_model, newdata = test7)
prediction1<-data.frame(pred)
prediction1$pred[prediction1$pred >= 0.5] <-1
prediction1$pred[prediction1$pred < 0.5] <- 0
library("ROSE")
accuracy.meas(testdata$solubility, prediction1[,1])
roc.curve(testdata$solubility, prediction1[,1], plotit = F)
library(EvaluationMeasures)
EvaluationMeasures.MCC(testdata$solubility, prediction1[,1])
library('caret')
confusionMatrix(prediction1[,1], reference = testdata$solubility)

#different descriptors for solubility and normailization
x=read.csv("solubility2.csv")
x[,1]<-NULL
#there are some descriptors have limits for lenght of sequence, so shorter one need to be removed.
x <- x[!(x$gene=="kdpF"),]
write.csv(x,"solubility3.csv")
x=read.csv("solubility3.csv")
x[,1:2]<-NULL
x[,3]<-x[,1]
x[,1]<-NULL
x[ ,1]<-as.character(x[ ,1])
library("protr")
# calculate extractAPAAC descriptors
x1 = t(sapply(x[ ,1], extractAPAAC))
x[ ,3:82] <- x1[ ,1:80]
x[,1]<-NULL
x[ ,82]<-x[ ,1]
x[ ,1]<-NULL
colnames(x)[81] <- "solubility"
write.csv(x,"cleandata2.csv")
# calculate extractMoreauBroto descriptors
x=read.csv("solubility3.csv")
x[,1:2]<-NULL
x[,3]<-x[,1]
x[,1]<-NULL
x[ ,1]<-as.character(x[ ,1])

x1 = t(sapply(x[ ,1], extractMoreauBroto))
x[ ,3:242] <- x1[ ,1:240]
x[,1]<-NULL
x[ ,242]<-x[ ,1]
x[ ,1]<-NULL
colnames(x)[241] <- "solubility"
write.csv(x,"cleandata3.csv")
# calculate extractMoran descriptors
x=read.csv("solubility3.csv")
x[,1:2]<-NULL
x[,3]<-x[,1]
x[,1]<-NULL
x[ ,1]<-as.character(x[ ,1]) 

x1 = t(sapply(x[ ,1], extractMoran))
x[ ,3:242] <- x1[ ,1:240]
x[,1]<-NULL
x[ ,242]<-x[ ,1]
x[ ,1]<-NULL
colnames(x)[241] <- "solubility"
write.csv(x,"cleandata4.csv")
# calculate extractGeary descriptors
x=read.csv("solubility3.csv")
x[,1:2]<-NULL
x[,3]<-x[,1]
x[,1]<-NULL
x[ ,1]<-as.character(x[ ,1])

x1 = t(sapply(x[ ,1], extractGeary))
x[ ,3:242] <- x1[ ,1:240]
x[,1]<-NULL
x[ ,242]<-x[ ,1]
x[ ,1]<-NULL
colnames(x)[241] <- "solubility"
write.csv(x,"cleandata5.csv")
# calculate three descriptors
x=read.csv("solubility3.csv")
x[,1:2]<-NULL
x[,3]<-x[,1]
x[,1]<-NULL
x[ ,1]<-as.character(x[ ,1])

x1 = t(sapply(x[ ,1], extractCTDC))
x2 = t(sapply(x[ ,1], extractCTDT))
x3 = t(sapply(x[ ,1], extractCTDD))
x[ ,3:23] <- x1[ ,1:21]
x[ ,24:44] <- x2[ ,1:21]
x[ ,45:149] <- x3[ ,1:105]
x[,1]<-NULL
x[ ,149]<-x[ ,1]
x[ ,1]<-NULL
colnames(x)[148] <- "solubility"
write.csv(x,"cleandata6.csv")
# calculate extractCTriad descriptors
x1 = t(sapply(x[ ,1], extractCTriad))
x[ ,3:345] <- x1[ ,1:343]
x[,1]<-NULL
x[ ,345]<-x[ ,1]
x[ ,1]<-NULL
colnames(x)[344] <- "solubility"
write.csv(x,"cleandata7.csv")
# calculate extractSOCN descriptors
x1 = t(sapply(x[ ,1], extractSOCN))
x[ ,3:62] <- x1[ ,1:60]
x[,1]<-NULL
x[ ,62]<-x[ ,1]
x[ ,1]<-NULL
colnames(x)[61] <- "solubility"
write.csv(x,"cleandata8.csv")
# calculate extractQSO descriptors
x1 = t(sapply(x[ ,1], extractQSO))
x[ ,3:102] <- x1[ ,1:100]
x[,1]<-NULL
x[ ,102]<-x[ ,1]
x[ ,1]<-NULL
colnames(x)[101] <- "solubility"
write.csv(x,"cleandata9.csv")
# calculate extractPAAC descriptors
x1 = t(sapply(x[ ,1], extractPAAC))
x[ ,3:52] <- x1[ ,1:50]
x[,1]<-NULL
x[ ,52]<-x[ ,1]
x[ ,1]<-NULL
colnames(x)[51] <- "solubility"
write.csv(x,"cleandata10.csv")

#normalisation
x=read.csv("cleandata10.csv")
x[,1]<-NULL
x[ ,51][x[ ,51]>100]<-100
x[ ,51]<-x[ ,51]/100
colnames(x) <- sub("^X","",colnames(x))
x[,sub("Value","Normed",colnames(x))] <-(apply(x,2,function(x){(x-min(x)) / diff(range(x))}))
summary(x)
write.csv(x,"cleandataN10.csv")
x=read.csv("cleandata1.csv")
x[ ,1]<-NULL
#remove gene with short length manually
write.csv(x,"cleandataN1.csv")

#test accuracy and r square for different descrptors with svm
SVM1 = NULL
for (i in seq(1,10,1)) {
  n=paste("cleandataN",i,".csv",sep = "", collapse = "")
  clean7=read.csv(n)
  clean7[,1]<-NULL
  library('caret')
  set.seed(5)
  TrainingDataIndex <- createDataPartition(clean7$solubility, p=0.9, list = FALSE)
  train7 <- clean7[TrainingDataIndex,]
  test7<- clean7[-TrainingDataIndex,]
  testdata<-test7 
  test7$solubility<-NULL
  library("e1071")
  svm_model <- svm(solubility ~ ., train7, method="eps-regression", cost=1.5, gamma=0.026, epsilon = 0.198, cross=10)
  pred <- predict(svm_model, newdata = test7)
  prediction<-data.frame(pred)
  rsq <- function(x, y) summary(lm(y~x))$r.squared
  R2=rsq(testdata$solubility, prediction[,1])
  
  testdata$solubility[testdata$solubility >= 0.39] <-1
  testdata$solubility[testdata$solubility < 0.39] <- 0
  prediction$pred[prediction$pred >= 0.39] <-1
  prediction$pred[prediction$pred < 0.39] <- 0
  library('caret')
  cm<-confusionMatrix(prediction$pred, reference = testdata$solubility)
  overall <- cm$overall
  overall.accuracy <- overall['Accuracy'] 
  SVM1 = rbind(SVM1, data.frame(i,R2,overall.accuracy))
}
write.csv(SVM1,"SVM1.csv")
#ccombine cleandataN10 and cleandataN1  R2 0.3718951  Accuracy 0.7156
clean7=read.csv("cleandataN1.csv")
x=read.csv("cleandataN10.csv")
x[,1]<-NULL
clean7[,1]<-NULL
clean7[,21]<-NULL
clean7[,21:71]<-x[,1:51]
library('caret')
set.seed(5)
TrainingDataIndex <- createDataPartition(clean7$solubility, p=0.9, list = FALSE)
train7 <- clean7[TrainingDataIndex,]
test7<- clean7[-TrainingDataIndex,]
testdata<-test7 
test7$solubility<-NULL
library("e1071")
svm_model <- svm(solubility ~ ., train7, method="eps-regression", cost=1.5, gamma=0.026, epsilon = 0.198, cross=10)
pred <- predict(svm_model, newdata = test7)
prediction<-data.frame(pred)
rsq <- function(x, y) summary(lm(y~x))$r.squared
R2=rsq(testdata$solubility, prediction[,1])

testdata$solubility[testdata$solubility >= 0.39] <-1
testdata$solubility[testdata$solubility < 0.39] <- 0
prediction$pred[prediction$pred >= 0.39] <-1
prediction$pred[prediction$pred < 0.39] <- 0
library('caret')
cm<-confusionMatrix(prediction[,1], reference = testdata$solubility)
overall <- cm$overall
overall.accuracy <- overall['Accuracy'] 
#tune SVM for cleandataN10
clean7=read.csv("cleandataN10.csv")
clean7[ ,1]<-NULL
library('caret')
set.seed(5)
TrainingDataIndex <- createDataPartition(clean7$solubility, p=0.9, list = FALSE)
train7 <- clean7[TrainingDataIndex,]
test7<- clean7[-TrainingDataIndex,]

d6 = NULL
for (x in seq(0.01,2,0.02)) {
  test18=test7
  testdata<-test18 
  test18$solubility<-NULL
  library("e1071")
  svm_model <- svm(solubility ~ ., train7, method="eps-regression", cost=x, gamma=0.026, epsilon = 0.198, cross=10)
  pred <- predict(svm_model, newdata = test18)
  prediction<-data.frame(pred)
  rsq <- function(x, y) summary(lm(y~x))$r.squared
  R2=rsq(testdata$solubility, prediction[,1])
  
  testdata$solubility[testdata$solubility >= 0.39] <-1
  testdata$solubility[testdata$solubility < 0.39] <- 0
  prediction$pred[prediction$pred >= 0.39] <-1
  prediction$pred[prediction$pred < 0.39] <- 0
  library('caret')
  cm<-confusionMatrix(prediction$pred, reference = testdata$solubility)
  overall <- cm$overall
  overall.accuracy <- overall['Accuracy'] 
  d6 = rbind(d6, data.frame(x,R2,overall.accuracy))
}

d7 = NULL
for (y in seq(0.001,0.1,0.001)) {
  test18=test7
  testdata<-test18 
  test18$solubility<-NULL
  library("e1071")
  svm_model <- svm(solubility ~ ., train7, method="eps-regression", cost=1.5, gamma=y, epsilon = 0.198, cross=10)
  pred <- predict(svm_model, newdata = test18)
  prediction<-data.frame(pred)
  rsq <- function(x, y) summary(lm(y~x))$r.squared
  R2=rsq(testdata$solubility, prediction[,1])
  
  testdata$solubility[testdata$solubility >= 0.39] <-1
  testdata$solubility[testdata$solubility < 0.39] <- 0
  prediction$pred[prediction$pred >= 0.39] <-1
  prediction$pred[prediction$pred < 0.39] <- 0
  library('caret')
  cm<-confusionMatrix(prediction$pred, reference = testdata$solubility)
  overall <- cm$overall
  overall.accuracy <- overall['Accuracy'] 
  d7 = rbind(d7, data.frame(y,R2,overall.accuracy))
}

d8 = NULL
for (z in seq(0.001,0.3,0.004)) {
  test18=test7
  testdata<-test18 
  test18$solubility<-NULL
  library("e1071")
  svm_model <- svm(solubility ~ ., train7, method="eps-regression", cost=1.5, gamma=0.026, epsilon = z, cross=10)
  pred <- predict(svm_model, newdata = test18)
  prediction<-data.frame(pred)
  rsq <- function(x, y) summary(lm(y~x))$r.squared
  R2=rsq(testdata$solubility, prediction[,1])
  
  testdata$solubility[testdata$solubility >= 0.39] <-1
  testdata$solubility[testdata$solubility < 0.39] <- 0
  prediction$pred[prediction$pred >= 0.39] <-1
  prediction$pred[prediction$pred < 0.39] <- 0
  library('caret')
  cm<-confusionMatrix(prediction$pred, reference = testdata$solubility)
  overall <- cm$overall
  overall.accuracy <- overall['Accuracy'] 
  d8 = rbind(d8, data.frame(z,R2,overall.accuracy))
}
#max R2 0.4047999 in d7,max accuracy 0.7385204 in d7

#other machine learning methods, 10 fold cross validation
clean7=read.csv("cleandata1.csv")
clean7[ ,1]<-NULL
yourData=clean7
#Randomly shuffle the data
yourData<-yourData[sample(nrow(yourData)),]

#Create 10 equally size folds
folds <- cut(seq(1,nrow(yourData)),breaks=10,labels=FALSE)

#Perform 10 fold cross validation
d = NULL #accuracy
dd=NULL  #R square
for (x in seq(1,10,1)) {
  #Segement your data by fold using the which() function 
  testIndexes <- which(folds==x,arr.ind=TRUE)
  test7 <- yourData[testIndexes, ]
  train7 <- yourData[-testIndexes, ]
  testdata<-test7 
  test7[ ,21]<-NULL
  #Logistical Regression
  #input <- train7
  #solubility.data = glm(formula = solubility ~ ., data = input, family = binomial)
  #pred <- predict(solubility.data, newdata = test7)
  #library(party)
  #output.tree <- ctree(solubility ~ ., data=train7)
  #pred <- predict(output.tree, newdata = test7)
  #naive bayes
  #train7$solubility<-as.factor(train7$solubility)
  #library("e1071")
  #model <- naiveBayes(solubility ~ ., data=train7)
  #pred <- predict(model, newdata = test7)
  #prediction<-data.frame(pred)
  #transfer factor into numeric but do not change value
  #prediction$pred<-as.numeric(paste(prediction$pred))
  cf <- cforest(solubility ~ ., data = train7)
  pred <- predict(cf, newdata = test7, type = "response")
  #library("e1071")
  #svm_model <- svm(solubility ~ ., train7, method="eps-regression", cost=1.5, gamma=0.026, epsilon = 0.198)
  #pred <- predict(svm_model, newdata = test7)
  prediction<-data.frame(pred)
  rsq <- function(x, y) summary(lm(y~x))$r.squared
  a=rsq(testdata$solubility, prediction[,1])
  dd = rbind(dd, data.frame(x,a))
  testdata$solubility[testdata$solubility >= 0.39] <-1
  testdata$solubility[testdata$solubility < 0.39] <- 0
  prediction[,1][prediction[,1] >= 0.39] <-1
  prediction[,1][prediction[,1] < 0.39] <- 0
  
  library('caret')
  cm<-confusionMatrix(prediction[,1], reference = testdata$solubility)
  overall <- cm$overall
  overall.accuracy <- overall['Accuracy'] 
  d = rbind(d, data.frame(x,overall.accuracy))
} 
average_d=sum(d[1:10,2])/10 
average_dd=sum(dd[1:10,2])/10 
#binary
yourData$solubility[yourData$solubility >= 0.39] <-1
yourData$solubility[yourData$solubility < 0.39] <- 0
yourData<-yourData[sample(nrow(yourData)),]
folds <- cut(seq(1,nrow(yourData)),breaks=10,labels=FALSE)
#Perform 10 fold cross validation
d = NULL #accuracy
for (x in seq(1,10,1)) {
  testIndexes <- which(folds==x,arr.ind=TRUE)
  test7 <- yourData[testIndexes, ]
  train7 <- yourData[-testIndexes, ]
  #Use the test and train data partitions however you desire...
  testdata<-test7 
  test7[ ,21]<-NULL
  #Logistical Regression
  #input <- train7
  #solubility.data = glm(formula = solubility ~ ., data = input, family = binomial)
  #pred <- predict(solubility.data, newdata = test7)
  #decsion tree
  #library(party)
  #output.tree <- ctree(solubility ~ ., data=train7)
  #pred <- predict(output.tree, newdata = test7)
  #naive bayes
  #train7$solubility<-as.factor(train7$solubility)
  #library("e1071")
  #model <- naiveBayes(solubility ~ ., data=train7)
  #pred <- predict(model, newdata = test7)
  #prediction<-data.frame(pred)
  #transfer factor into numeric but do not change value
  #prediction$pred<-as.numeric(paste(prediction$pred))
  cf <- cforest(solubility ~ ., data = train7)
  pred <- predict(cf, newdata = test7, type = "response")
  prediction<-data.frame(pred)
  prediction[,1][prediction[,1] >= 0.39] <-1
  prediction[,1][prediction[,1] < 0.39] <- 0
  library('caret')
  cm<-confusionMatrix(prediction[,1], reference = testdata$solubility)
  overall <- cm$overall
  overall.accuracy <- overall['Accuracy'] 
  d = rbind(d, data.frame(x,overall.accuracy))
} 
average_d=sum(d[1:10,2])/10 

##optimize the parameters of SVM--kernel:figure 19-radial, figure 20-linear, 21-sigmoid
## figure 22-RBF(default), figure 23-polynomial
clean7=read.csv("cleandata1.csv")
clean7[ ,1]<-NULL
##use 90% as training data
library('caret')
set.seed(5)
TrainingDataIndex <- createDataPartition(clean7[ ,21], p=0.9, list = FALSE)
train7 <- clean7[TrainingDataIndex,]
test7<- clean7[-TrainingDataIndex,]
testdata<-test7 
test7[ ,21]<-NULL
library("e1071")
svm_model <- svm(solubility ~ ., train7, method="eps-regression",cross=10)
pred <- predict(svm_model, newdata = test7)
prediction<-data.frame(pred)
rsq <- function(x, y) summary(lm(y~x))$r.squared
R2=rsq(testdata$solubility, prediction[,1])
print(R2)
testdata$solubility[testdata$solubility >= 0.39] <-1
testdata$solubility[testdata$solubility < 0.39] <- 0
prediction<-data.frame(pred)
prediction$pred[prediction$pred >= 0.39] <-1
prediction$pred[prediction$pred < 0.39] <- 0
library("ROSE")
accuracy.meas(testdata$solubility, prediction$pred)
roc.curve(testdata$solubility, prediction$pred, plotit = F)
library(EvaluationMeasures)
EvaluationMeasures.MCC(testdata$solubility, prediction$pred)
library('caret')
confusionMatrix(prediction$pred, reference = testdata$solubility)

#produce input data train21-23 into python and leave 25% as test data
clean7=read.csv("cleandata1.csv")
clean7[,1]<-NULL
library('caret')
set.seed(7)
TrainingDataIndex <- createDataPartition(clean7[ ,21], p=0.75, list = FALSE)
train7 <- clean7[TrainingDataIndex,]
test7<- clean7[-TrainingDataIndex,]
write.csv(train7,"train23.csv")
write.csv(test7,"test23.csv")
#different GAN for SVM, use hold out test data to test

test19=read.csv("test23.csv")
test19[,1]<-NULL
#change column name for test 19
for (i in seq(1,20,1)) {
  q=i-1
  colnames(test19)[i] <- paste("X", q, sep="")
}
write.csv(test19,"test23.csv")

test19=read.csv("test23.csv")
test19[,1]<-NULL
test=test19

WCGAN1 = NULL
for (i in seq(1,5,1)) {
  n=paste("WCGAN_",i,"00_0408.csv",sep = "", collapse = "")
  clean7=read.csv(n)
  clean7[,1]<-NULL
  colnames(clean7)[21] <- "solubility"
  test18=test
  testdata<-test18 
  test18[ ,21]<-NULL
  library("e1071")
  svm_model <- svm(solubility ~ ., clean7, method="eps-regression", cost=1.26, gamma=0.029, epsilon = 0.057, cross=10)
  pred <- predict(svm_model, newdata = test18)
  prediction<-data.frame(pred)
  rsq <- function(x, y) summary(lm(y~x))$r.squared
  R2=rsq(testdata$solubility, prediction[,1])
  
  testdata$solubility[testdata$solubility >= 0.39] <-1
  testdata$solubility[testdata$solubility < 0.39] <- 0
  prediction$pred[prediction$pred >= 0.39] <-1
  prediction$pred[prediction$pred < 0.39] <- 0
  library('caret')
  cm<-confusionMatrix(prediction$pred, reference = testdata$solubility)
  overall <- cm$overall
  overall.accuracy <- overall['Accuracy'] 
  WCGAN1 = rbind(WCGAN1, data.frame(i,R2,overall.accuracy))
}
write.csv(GAN1,"GAN0408_1.csv")
GAN1=read.csv("WCGAN0408_1.csv")
write.csv(CGAN1,"CGAN0408_1.csv")
write.csv(WGAN1,"WGAN0408_1.csv")
write.csv(WCGAN1,"WCGAN0408_1.csv")

#5000 iterations for GAN
clean7=read.csv("cleandata1.csv")
clean7[ ,1]<-NULL

library('caret')
set.seed(5)
TrainingDataIndex <- createDataPartition(clean7[ ,21], p=0.9, list = FALSE)
train7 <- clean7[TrainingDataIndex,]
test7<- clean7[-TrainingDataIndex,]
testdata<-test7 
test7[ ,21]<-NULL
library("e1071")
svm_model <- svm(solubility ~ ., train7, method="eps-regression", cost=1.26, gamma=0.029, epsilon = 0.057, cross=10)
pred <- predict(svm_model, newdata = test7)
prediction<-data.frame(pred)
rsq <- function(x, y) summary(lm(y~x))$r.squared
R2=rsq(testdata$solubility, prediction[,1])

testdata$solubility[testdata$solubility >= 0.39] <-1
testdata$solubility[testdata$solubility < 0.39] <- 0
prediction$pred[prediction$pred >= 0.39] <-1
prediction$pred[prediction$pred < 0.39] <- 0
library('caret')
cm<-confusionMatrix(prediction$pred, reference = testdata$solubility)
overall <- cm$overall
overall.accuracy <- overall['Accuracy'] 


GAN6 = NULL
for (i in seq(1,50,1)) {
  n=paste("GAN_",i,"00_0408_2.csv",sep = "", collapse = "")
  clean7=read.csv(n)
  clean7[,1]<-NULL
  colnames(clean7)[21] <- "solubility"
  test18=test
  testdata<-test18 
  test18[ ,21]<-NULL
  library("e1071")
  svm_model <- svm(solubility ~ ., clean7, method="eps-regression", cost=1.26, gamma=0.029, epsilon = 0.057, cross=10)
  pred <- predict(svm_model, newdata = test18)
  prediction<-data.frame(pred)
  rsq <- function(x, y) summary(lm(y~x))$r.squared
  R2=rsq(testdata$solubility, prediction[,1])
  
  testdata$solubility[testdata$solubility >= 0.39] <-1
  testdata$solubility[testdata$solubility < 0.39] <- 0
  prediction$pred[prediction$pred >= 0.39] <-1
  prediction$pred[prediction$pred < 0.39] <- 0
  library('caret')
  cm<-confusionMatrix(prediction$pred, reference = testdata$solubility)
  overall <- cm$overall
  overall.accuracy <- overall['Accuracy'] 
  GAN6 = rbind(GAN6, data.frame(i,R2,overall.accuracy))
}
write.csv(GAN4,"GAN0408_2.csv")
write.csv(GAN5,"GAN0408_3.csv")
write.csv(GAN6,"GAN0408_4.csv")

# yield dataset(use 3147 items for all yield data)
library("readxl")
traningdata <- read_excel("all_data.xlsx")
traningdata[ ,1:3] <- NULL
traningdata[ ,2:3] <- NULL
traningdata[ ,3] <- NULL
traningdata[ ,4:22] <- NULL
x<-na.omit(traningdata)
sequence <- read_excel("seqinfo.xlsx")
colnames(x)[1] <- "gene"

no_sequence=read.csv("no_sequence1.csv")
no_sequence$sequence <- as.character(no_sequence$sequence)
myfunction <- function(y){
  w=gsub(" ", "", y, fixed = TRUE)
  return(gsub("\n", "", w, fixed = TRUE))
}
for (i in seq(1,519,1)) { 
  no_sequence[i,2]<-myfunction(no_sequence[i,2])
}
class(no_sequence$sequence)
write.csv(no_sequence,"no_sequence5.csv")
data <- read.csv("no_sequence5.csv",na.strings=c(""," ","NA"),header=TRUE)
summary(data)
new_DF3 <- data[is.na(data$sequence),]
new_DF4 <- data[!(is.na(data$sequence)),]
new_DF5 <- new_DF4[!(new_DF4$gene=="yjgX"),]
new_DF6 <- new_DF5[!(new_DF5$gene=="gatR"),]

library(MASS)
merged.train3 <- merge(x = x, y = new_DF6,
                       by.x = ("gene"),
                       by.y = ("gene")
)

merge=read.csv("merge.csv")
merge$solubility <- merge[ ,8]
merge$yield <- merge[ ,10]
merge[ ,1] <- NULL
merge[ ,2:28] <- NULL
## romove the rows including NA
x<-na.omit(merge)
merged.train3$solubility<-merged.train3[,2]
merged.train3$yield<-merged.train3[,3]
merged.train3[,2:4]<-NULL
totalyield <- rbind(x, merged.train3)
#Total proteins 3152=3173-21(17+4)
total=totalyield
total[ ,3][total[ ,3]>100]<-100
total[ ,3]<-total[ ,3]/100

total <- total[!(total$gene=="ybfH"),]
total <- total[!(total$gene=="yhdW"),]
total <- total[!(total$gene=="ymgH"),]
total <- total[!(total$gene=="yohG"),]
write.csv(total,"solubility4.csv")
#total proteins 3148(4 proteins include undertemined amino acid X)
x=read.csv("solubility4.csv")
x[,1]<-NULL
x <- x[!(x$gene=="kdpF"),]
write.csv(x,"solubility5.csv")


rr=read.csv("solubility5.csv")
#plot the distribution of yield and solubility
plot(rr[,1],rr[,4], main="solubility", 
     xlab="index ", ylab="solubility ")

plot(rr[,1],rr[,5], main="yield", 
     xlab="index ", ylab="yield ")

#remove outliers  remove 99 points 1.5 times IQR
clean19=read.csv("solubility5.csv")
x=clean19$yield
quantile(x)
IQR(x)
qnt = quantile(x, probs=c(.25, .75))
iqt = 1.5 * IQR(x)
y = x 
y[x < (qnt[1] - iqt)] = NA
y[x > (qnt[2] + iqt)] = NA
y[complete.cases(y)]
summary(y)

ww=read.csv("cleandataN1.csv")
ww[,1]<-NULL
ww$solubility<-NULL
ww$yield<-y
x<-na.omit(ww)

colnames(x) <- sub("^X","",colnames(x))
x[,sub("Value","Normed",colnames(x))] <-(apply(x,2,function(x){(x-min(x)) / diff(range(x))}))
summary(x)
write.csv(x,"cleanyieldN1.csv")

#try SVM
clean7=read.csv("cleanyieldN1.csv")
clean7[ ,1]<-NULL
##use 90% as training data
library('caret')
set.seed(5)
TrainingDataIndex <- createDataPartition(clean7[ ,21], p=0.9, list = FALSE)
train7 <- clean7[TrainingDataIndex,]
test7<- clean7[-TrainingDataIndex,]
testdata<-test7 
test7[ ,21]<-NULL

library("e1071")
svm_model <- svm(yield ~ ., train7, method="eps-regression", cross=10)
pred <- predict(svm_model, newdata = test7)

d = NULL
for (x in seq(0.1,0.9,0.01)) {
  testdata$yield[testdata$yield >= x] <-1
  testdata$yield[testdata$yield < x] <- 0
  prediction<-data.frame(pred)
  prediction$pred[prediction$pred >= x] <-1
  prediction$pred[prediction$pred < x] <- 0
  library("ROSE")
  accuracy.meas(testdata$yield, prediction$pred)
  roc.curve(testdata$yield, prediction$pred, plotit = F)
  library(EvaluationMeasures)
  EvaluationMeasures.MCC(testdata$yield, prediction$pred)
  library('caret')
  cm<-confusionMatrix(prediction$pred, reference = testdata$yield)
  overall <- cm$overall
  overall.accuracy <- overall['Accuracy'] 
  print(x)
  print(overall.accuracy)
  d = rbind(d, data.frame(x,overall.accuracy))
}
write.csv(d,"accuracy_cutoff2.csv")

#try different models
clean7=read.csv("cleanyieldN1.csv")
clean7[ ,1]<-NULL
yourData=clean7
#Randomly shuffle the data
yourData<-yourData[sample(nrow(yourData)),]

#Create 10 equally size folds
folds <- cut(seq(1,nrow(yourData)),breaks=10,labels=FALSE)

#Perform 10 fold cross validation
d = NULL #accuracy
dd=NULL  #R square
for (x in seq(1,10,1)) {
  #Segement your data by fold using the which() function 
  testIndexes <- which(folds==x,arr.ind=TRUE)
  test7 <- yourData[testIndexes, ]
  train7 <- yourData[-testIndexes, ]
  testdata<-test7 
  test7[ ,21]<-NULL
  #Logistical Regression
  #input <- train7
  #solubility.data = glm(formula = yield ~ ., data = input, family = binomial)
  #pred <- predict(solubility.data, newdata = test7)
  #library(party)
  #output.tree <- ctree(yield ~ ., data=train7)
  #pred <- predict(output.tree, newdata = test7)
  #naive bayes
  #train7$yield<-as.factor(train7$yield)
  #library("e1071")
  #model <- naiveBayes(yield ~ ., data=train7)
  #pred <- predict(model, newdata = test7)
  #prediction<-data.frame(pred)
  #transfer factor into numeric but do not change value
  #prediction$pred<-as.numeric(paste(prediction$pred))
  cf <- cforest(yield ~ ., data = train7)
  pred <- predict(cf, newdata = test7, type = "response")
  #library("e1071")
  #svm_model <- svm(solubility ~ ., train7, method="eps-regression", cost=1.5, gamma=0.026, epsilon = 0.198)
  #pred <- predict(svm_model, newdata = test7)
  prediction<-data.frame(pred)
  rsq <- function(x, y) summary(lm(y~x))$r.squared
  a=rsq(testdata$yield, prediction[,1])
  dd = rbind(dd, data.frame(x,a))
  testdata$yield[testdata$yield >= 0.2706] <-1
  testdata$yield[testdata$yield < 0.2706] <- 0
  prediction[,1][prediction[,1] >= 0.2706] <-1
  prediction[,1][prediction[,1] < 0.2706] <- 0
  
  library('caret')
  cm<-confusionMatrix(prediction[,1], reference = testdata$yield)
  overall <- cm$overall
  overall.accuracy <- overall['Accuracy'] 
  d = rbind(d, data.frame(x,overall.accuracy))
} 
average_d=sum(d[1:10,2])/10 
average_dd=sum(dd[1:10,2])/10 

#tune SVM
clean7=read.csv("cleanyieldN1.csv")
clean7[ ,1]<-NULL
library('caret')
set.seed(5)
TrainingDataIndex <- createDataPartition(clean7[ ,21], p=0.9, list = FALSE)
train7 <- clean7[TrainingDataIndex,]
test7<- clean7[-TrainingDataIndex,]

#different kernels
test18=test7
testdata<-test18 
test18[ ,21]<-NULL
library("e1071")
svm_model <- svm(yield ~ ., train7, method="eps-regression", cross=10)
#"polynomial","linear",
pred <- predict(svm_model, newdata = test18)
prediction<-data.frame(pred)
rsq <- function(x, y) summary(lm(y~x))$r.squared
R2=rsq(testdata$yield, prediction[,1])
print(R2)

d = NULL
for (x in seq(0.02,2,0.02)) {
  test18=test7
  testdata<-test18 
  test18[ ,21]<-NULL
  library("e1071")
  svm_model <- svm(yield ~ ., train7, method="eps-regression", cost=x, cross=10)
  pred <- predict(svm_model, newdata = test18)
  prediction<-data.frame(pred)
  rsq <- function(x, y) summary(lm(y~x))$r.squared
  R2=rsq(testdata$yield, prediction[,1])
  
  testdata$yield[testdata$yield >= 0.2706] <-1
  testdata$yield[testdata$yield < 0.2706] <- 0
  prediction[,1][prediction[,1] >= 0.2706] <-1
  prediction[,1][prediction[,1] < 0.2706] <- 0
  library('caret')
  cm<-confusionMatrix(prediction$pred, reference = testdata$yield)
  overall <- cm$overall
  overall.accuracy <- overall['Accuracy'] 
  d = rbind(d, data.frame(x,R2,overall.accuracy))
}
write.csv(d,"cost3.csv")

plot(x = d$x,y = d$R2,
     xlab = "cost",
     ylab = "SVM-R Square",
     xlim = c(0.02,2),
     ylim = c(min(d$R2),max(d$R2)),		 
     main = "Correlation between SVM R Square and cost"
)
d1 = NULL
for (y in seq(0.001,0.1,0.001)) {
  test18=test7
  testdata<-test18 
  test18[ ,21]<-NULL
  library("e1071")
  svm_model <- svm(yield ~ ., train7, method="eps-regression", cost=0.2, gamma=y, cross=10)
  pred <- predict(svm_model, newdata = test18)
  prediction<-data.frame(pred)
  rsq <- function(x, y) summary(lm(y~x))$r.squared
  R2=rsq(testdata$yield, prediction[,1])
  
  testdata$yield[testdata$yield >= 0.2706] <-1
  testdata$yield[testdata$yield < 0.2706] <- 0
  prediction[,1][prediction[,1] >= 0.2706] <-1
  prediction[,1][prediction[,1] < 0.2706] <- 0
  library('caret')
  cm<-confusionMatrix(prediction$pred, reference = testdata$yield)
  overall <- cm$overall
  overall.accuracy <- overall['Accuracy'] 
  d1 = rbind(d1, data.frame(y,R2,overall.accuracy))
}
write.csv(d1,"gamma3.csv")
plot(x = d1$y,y = d1$R2,
     xlab = "gamma",
     ylab = "SVM-R Square",
     xlim = c(0.001,0.1),
     ylim = c(min(d1$R2),max(d1$R2)),	 
     main = "Correlation between SVM R Square and gamma"
)
d2 = NULL
for (z in seq(0.001,0.3,0.004)) {
  test18=test7
  testdata<-test18 
  test18[ ,21]<-NULL
  library("e1071")
  svm_model <- svm(yield ~ ., train7, method="eps-regression", cost=0.2, gamma=0.027, epsilon = z, cross=10)
  pred <- predict(svm_model, newdata = test18)
  prediction<-data.frame(pred)
  rsq <- function(x, y) summary(lm(y~x))$r.squared
  R2=rsq(testdata$yield, prediction[,1])
  
  testdata$yield[testdata$yield >= 0.2706] <-1
  testdata$yield[testdata$yield < 0.2706] <- 0
  prediction[,1][prediction[,1] >= 0.2706] <-1
  prediction[,1][prediction[,1] < 0.2706] <- 0
  library('caret')
  cm<-confusionMatrix(prediction$pred, reference = testdata$yield)
  overall <- cm$overall
  overall.accuracy <- overall['Accuracy'] 
  d2 = rbind(d2, data.frame(z,R2,overall.accuracy))
}
write.csv(d2,"epsilon3.csv")
plot(x = d2$z,y = d2$R2,
     xlab = "epsilon",
     ylab = "SVM-R Square",
     xlim = c(0.001,0.3),
     ylim = c(min(d2$R2),max(d2$R2)),	 
     main = "Correlation between SVM R Square and epsilon"
)

#different discriptors
clean19=read.csv("solubility5.csv")
x=clean19$yield
quantile(x)
IQR(x)
qnt = quantile(x, probs=c(.25, .75))
iqt = 1.5 * IQR(x)
y = x 
y[x < (qnt[1] - iqt)] = NA
y[x > (qnt[2] + iqt)] = NA
y[complete.cases(y)]
summary(y)

for (i in seq(2,10,1)) {
  n=paste("cleandataN",i,".csv",sep = "", collapse = "")
  ww=read.csv(n)
  ww[,1]<-NULL
  ww$solubility<-NULL
  ww$yield<-y
  x<-na.omit(ww)
  colnames(x) <- sub("^X","",colnames(x))
  x[,sub("Value","Normed",colnames(x))] <-(apply(x,2,function(x){(x-min(x)) / diff(range(x))}))
  m=paste("cleanyieldN",i,".csv",sep = "", collapse = "")
  write.csv(x,m)
}

SVM1 = NULL
for (i in seq(1,10,1)) {
  n=paste("cleanyieldN",i,".csv",sep = "", collapse = "")
  clean7=read.csv(n)
  clean7[,1]<-NULL
  library('caret')
  set.seed(5)
  TrainingDataIndex <- createDataPartition(clean7$yield, p=0.9, list = FALSE)
  train7 <- clean7[TrainingDataIndex,]
  test7<- clean7[-TrainingDataIndex,]
  testdata<-test7 
  test7$yield<-NULL
  library("e1071")
  svm_model <- svm(yield ~ ., train7, method="eps-regression", cost=0.2, gamma=0.027, epsilon = 0.241, cross=10)
  pred <- predict(svm_model, newdata = test7)
  prediction<-data.frame(pred)
  rsq <- function(x, y) summary(lm(y~x))$r.squared
  R2=rsq(testdata$yield, prediction[,1])
  
  testdata$yield[testdata$yield >= 0.2706] <-1
  testdata$yield[testdata$yield < 0.2706] <- 0
  prediction[,1][prediction[,1] >= 0.2706] <-1
  prediction[,1][prediction[,1] < 0.2706] <- 0
  library('caret')
  cm<-confusionMatrix(prediction$pred, reference = testdata$yield)
  overall <- cm$overall
  overall.accuracy <- overall['Accuracy'] 
  SVM1 = rbind(SVM1, data.frame(i,R2,overall.accuracy))
}
write.csv(SVM1,"SVM2.csv") 

#tune SVM for another descriptor
clean7=read.csv("cleanyieldN2.csv")
clean7[ ,1]<-NULL
library('caret')
set.seed(5)
TrainingDataIndex <- createDataPartition(clean7$yield, p=0.9, list = FALSE)
train7 <- clean7[TrainingDataIndex,]
test7<- clean7[-TrainingDataIndex,]

d = NULL
for (x in seq(0.02,2,0.02)) {
  test18=test7
  testdata<-test18 
  test18$yield<-NULL
  library("e1071")
  svm_model <- svm(yield ~ ., train7, method="eps-regression", cost=x, cross=10)
  pred <- predict(svm_model, newdata = test18)
  prediction<-data.frame(pred)
  rsq <- function(x, y) summary(lm(y~x))$r.squared
  R2=rsq(testdata$yield, prediction[,1])
  
  testdata$yield[testdata$yield >= 0.2706] <-1
  testdata$yield[testdata$yield < 0.2706] <- 0
  prediction[,1][prediction[,1] >= 0.2706] <-1
  prediction[,1][prediction[,1] < 0.2706] <- 0
  library('caret')
  cm<-confusionMatrix(prediction$pred, reference = testdata$yield)
  overall <- cm$overall
  overall.accuracy <- overall['Accuracy'] 
  d = rbind(d, data.frame(x,R2,overall.accuracy))
}
write.csv(d,"cost4.csv")

plot(x = d$x,y = d$R2,
     xlab = "cost",
     ylab = "SVM-R Square",
     xlim = c(0.02,2),
     ylim = c(min(d$R2),max(d$R2)),		 
     main = "Correlation between SVM R Square and cost"
)
d1 = NULL
for (y in seq(0.001,0.1,0.001)) {
  test18=test7
  testdata<-test18 
  test18$yield<-NULL
  library("e1071")
  svm_model <- svm(yield ~ ., train7, method="eps-regression", cost=0.38, gamma=y, cross=10)
  pred <- predict(svm_model, newdata = test18)
  prediction<-data.frame(pred)
  rsq <- function(x, y) summary(lm(y~x))$r.squared
  R2=rsq(testdata$yield, prediction[,1])
  
  testdata$yield[testdata$yield >= 0.2706] <-1
  testdata$yield[testdata$yield < 0.2706] <- 0
  prediction[,1][prediction[,1] >= 0.2706] <-1
  prediction[,1][prediction[,1] < 0.2706] <- 0
  library('caret')
  cm<-confusionMatrix(prediction$pred, reference = testdata$yield)
  overall <- cm$overall
  overall.accuracy <- overall['Accuracy'] 
  d1 = rbind(d1, data.frame(y,R2,overall.accuracy))
}
write.csv(d1,"gamma4.csv")
plot(x = d1$y,y = d1$R2,
     xlab = "gamma",
     ylab = "SVM-R Square",
     xlim = c(0.001,0.1),
     ylim = c(min(d1$R2),max(d1$R2)),	 
     main = "Correlation between SVM R Square and gamma"
)
d2 = NULL
for (z in seq(0.001,0.3,0.004)) {
  test18=test7
  testdata<-test18 
  test18$yield<-NULL
  library("e1071")
  svm_model <- svm(yield ~ ., train7, method="eps-regression", cost=0.2, gamma=0.013, epsilon = z, cross=10)
  pred <- predict(svm_model, newdata = test18)
  prediction<-data.frame(pred)
  rsq <- function(x, y) summary(lm(y~x))$r.squared
  R2=rsq(testdata$yield, prediction[,1])
  
  testdata$yield[testdata$yield >= 0.2706] <-1
  testdata$yield[testdata$yield < 0.2706] <- 0
  prediction[,1][prediction[,1] >= 0.2706] <-1
  prediction[,1][prediction[,1] < 0.2706] <- 0
  library('caret')
  cm<-confusionMatrix(prediction$pred, reference = testdata$yield)
  overall <- cm$overall
  overall.accuracy <- overall['Accuracy'] 
  d2 = rbind(d2, data.frame(z,R2,overall.accuracy))
}
write.csv(d2,"epsilon4.csv")
#epsilon=0.001 is best
plot(x = d2$z,y = d2$R2,
     xlab = "epsilon",
     ylab = "SVM-R Square",
     xlim = c(0.001,0.3),
     ylim = c(min(d2$R2),max(d2$R2)),	 
     main = "Correlation between SVM R Square and epsilon"
)

##remove outliers by 3 times sd
clean19=read.csv("solubility5.csv")
x=clean19$yield
x<-data.frame(x)
findOutlier <- function(data, cutoff = 3) {
  ## Calculate the sd
  sds <- apply(x, 2, sd, na.rm = TRUE)
  ## Identify the cells with value greater than cutoff * sd (column wise)
  result <- mapply(function(d, s) {
    which(d > cutoff * s)
  }, data, sds)
  result
}
outliers <- findOutlier(x)
outliernumber <-as.vector(outliers[,1])

for (i in outliers[,1]) {x[i,1]<-NA}
summary(x)
ww=read.csv("cleandataN1.csv")
ww[,1]<-NULL
ww$solubility<-NULL
ww$yield<-x[,1]
y<-na.omit(ww)

#remove 164 outliers
x<-y
colnames(x) <- sub("^X","",colnames(x))
x[,sub("Value","Normed",colnames(x))] <-(apply(x,2,function(x){(x-min(x)) / diff(range(x))}))
summary(x)
write.csv(x,"cleanyieldN1_3sd.csv")

clean7=read.csv("solubility5.csv")
clean7=read.csv("cleanyieldN1_3sd.csv")
clean7[ ,1]<-NULL
library('caret')
set.seed(5)
TrainingDataIndex <- createDataPartition(clean7[ ,21], p=0.9, list = FALSE)
train7 <- clean7[TrainingDataIndex,]
test7<- clean7[-TrainingDataIndex,]


test18=test7
testdata<-test18 
test18[ ,21]<-NULL
library("e1071")
svm_model <- svm(yield ~ ., train7, method="eps-regression", cross=10)

pred <- predict(svm_model, newdata = test18)
prediction<-data.frame(pred)
rsq <- function(x, y) summary(lm(y~x))$r.squared
R2=rsq(testdata$yield, prediction[,1])
print(R2)

