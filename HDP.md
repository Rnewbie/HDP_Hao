---
title: "Large-scale structure-activity relationship of host defense peptides"
author: Saw Simeon, Hao Li, Theeraphon Piacham, Virapong Prachayasittikul and Chanin
  Nantasanamat
date: "October 13, 2015"
output: html_document
---

###HDPs having their respective activity were merge with negative data set to build seperate models


```
## $cancer
##           Training set Internal validation External set
## ACC_Mean       99.2613             97.5583      97.6046
## ACC_SD          0.1433              0.3102       0.6682
## Sens_Mean      99.3600             98.3807      98.3492
## Sens_SD         0.1644              0.2511       0.6887
## Spec_Mean      98.8989             94.5317      94.9041
## Spec_SD         0.5444              0.9232       2.1142
## MCC_Mean        0.9780              0.9274       0.9287
## MCC_SD          0.0043              0.0092       0.0201
## 
## $fungus
##           Training set Internal validation External set
## ACC_Mean       99.1317             96.9126      97.1210
## ACC_SD          0.1578              0.2887       0.5901
## Sens_Mean      99.3237             96.8298      96.9850
## Sens_SD         0.1842              0.3971       0.9161
## Spec_Mean      98.9841             96.9786      97.2420
## Spec_SD         0.2723              0.3417       0.8210
## MCC_Mean        0.9824              0.9373       0.9416
## MCC_SD          0.0032              0.0059       0.0120
## 
## $bacteria
##           Training set Internal validation External set
## ACC_Mean       99.1115             96.9422      97.0605
## ACC_SD          0.1624              0.2981       0.6091
## Sens_Mean      99.3654             96.9297      97.1218
## Sens_SD         0.1644              0.4589       1.0143
## Spec_Mean      98.9161             96.9535      97.0302
## Spec_SD         0.2763              0.2843       0.8232
## MCC_Mean        0.9820              0.9379       0.9404
## MCC_SD          0.0033              0.0061       0.0124
## 
## $virus
##           Training set Internal validation External set
## ACC_Mean       99.1646             96.0039      95.9189
## ACC_SD          0.1988              0.4324       0.8257
## Sens_Mean      98.4089             91.5197      91.3986
## Sens_SD         0.4843              1.0566       2.3605
## Spec_Mean      99.3926             97.3520      97.2890
## Spec_SD         0.2455              0.3853       0.9133
## MCC_Mean        0.9765              0.8876       0.8848
## MCC_SD          0.0056              0.0123       0.0238
```


###HDPs having anti--bacteria, anti--cancer, anti--fungus and virus-- were combined into one to represent HDPs while the negative set was used as non--HDPs.###


```
##           Training set Internal validation External set
## ACC_Mean       99.3873             98.1257      98.1922
## ACC_SD          0.0763              0.1237       0.2867
## Sens_Mean      99.4387             98.5458      98.5944
## Sens_SD         0.0878              0.1001       0.2640
## Spec_Mean      99.1291             95.9758      96.1536
## Spec_SD         0.1609              0.5045       1.1747
## MCC_Mean        0.9781              0.9325       0.9350
## MCC_SD          0.0027              0.0045       0.0104
```


###HDPs having their respective activity were combined together to build a union predictive model###

Results for training 


```
##           Bacteria  Cancer  Fungus   Virus   Overall
## ACC_Mean   96.0307 96.0505 96.0500 96.0401 96.042825
## ACC_SD      0.2174  0.1845  0.2185  0.2011  0.205375
## Sens_Mean  96.1914 96.2048 96.1826 96.1909 96.192425
## Sens_SD     0.2364  0.2113  0.2100  0.2265  0.221050
## Spec_Mean  93.3186 93.4566 93.7839 93.5015 93.515150
## Spec_SD     1.5234  1.5337  1.5985  1.5726  1.557050
## MCC_Mean    0.7225  0.7241  0.7239  0.7234  0.723475
## MCC_SD      0.0169  0.0145  0.0170  0.0158  0.016050
```

Results for CV


```
##           Bacteria  Cancer  Fungus   Virus   Overall
## ACC_Mean   92.9778 93.0262 93.0751 93.0796 93.039675
## ACC_SD      0.2731  0.2657  0.2426  0.2499  0.257825
## Sens_Mean  97.9807 98.0357 98.0257 98.0608 98.025725
## Sens_SD     0.2337  0.2109  0.2006  0.2191  0.216075
## Spec_Mean  42.9199 42.8544 43.5194 43.2379 43.132900
## Spec_SD     1.8937  2.0311  1.8305  1.9767  1.933000
## MCC_Mean    0.5057  0.5078  0.5130  0.5124  0.509725
## MCC_SD      0.0201  0.0202  0.0181  0.0189  0.019325
```

Results for Testing


```
##           Bacteria  Cancer  Fungus   Virus   Overall
## ACC_Mean   93.3841 93.3183 93.3104 93.2103 93.305775
## ACC_SD      0.5439  0.5368  0.5752  0.5747  0.557650
## Sens_Mean  94.6636 94.7143 94.5986 94.5303 94.626700
## Sens_SD     0.4644  0.3942  0.3781  0.4320  0.417175
## Spec_Mean  71.5641 70.7074 71.3498 70.4909 71.028050
## Spec_SD     5.1636  6.2726  6.3225  6.0771  5.958950
## MCC_Mean    0.5263  0.5260  0.5206  0.5114  0.521075
## MCC_SD      0.0452  0.0397  0.0418  0.0449  0.042900
```


###Function and Scriptrs for Analysis###
#####HDPs having their respective activity were merge with negative data set to build seperate models#####
```
### read FASTA
library(protr)
cancer <- readFASTA("cancer.fasta")
fungus <- readFASTA("fungus.fasta")
bacteria <- readFASTA("bacteria.fasta")
virus <- readFASTA("virus.fasta")
negative <- readFASTA("negative_Protein.fasta")

### removed wired protein
cancer <- cancer[(sapply(cancer, protcheck))]
fungus <- fungus[(sapply(fungus, protcheck))]
bacteria <- bacteria[(sapply(bacteria, protcheck))]
virus <- virus[(sapply(virus, protcheck))]
negative <- negative[(sapply(negative, protcheck))]


### Descriptors generation
cancer_des <- t(sapply(cancer, extractAAC))
fungus_des <- t(sapply(fungus, extractAAC))
bacteria_des <- t(sapply(bacteria, extractAAC))
virus_des <- t(sapply(virus, extractAAC))
negative_des <- t(sapply(negative, extractAAC))
#### label the labels 
cancer_des <- as.data.frame(cancer_des)
cancer_des$Label <- "Cancer"
fungus_des <- as.data.frame(fungus_des)
fungus_des$Label <- "Fungus"
bacteria_des <- as.data.frame(bacteria_des)
bacteria_des$Label <- "Bacteria"
virus_des <- as.data.frame(virus_des)
virus_des$Label <- "Virus"
negative_des <- as.data.frame(negative_des)
negative_des$Label <- "Negative"

cancer <- rbind(cancer_des, negative_des)
cancer$Label <- as.factor(cancer$Label)
fungus <- rbind(fungus_des, negative_des)
fungus$Label <- as.factor(fungus$Label)
bacteria <- rbind(fungus_des, negative_des)
bacteria$Label <- as.factor(bacteria$Label)
virus <- rbind(virus_des, negative_des)
virus$Label <- as.factor(virus$Label)



input <- list(cancer = cancer, fungus = fungus, bacteria = bacteria, 
              virus = virus)



#### training results using J48
J48_training <- function(x) {
  
  results <- list(100)
  for (i in 1:100) {
    in_train <- caret::createDataPartition(x$Label, p = 0.80, list = FALSE)
    train <- x[in_train, ]
    test <- x[-in_train, ]
    model_train <- RWeka::J48(Label~., data = train)
    summary <- summary(model_train)
    confusionmatrix <- summary$confusionMatrix
    results[[i]] <- as.numeric(confusionmatrix)
  }
  return(results)
}

mean_and_sd <- function(x) {
  c(round(mean(x, na.rm = TRUE), digits = 4),
    round(sd(x, na.rm = TRUE), digits = 4))
}

J48_train <- function(x) {
  ok <- J48_training(x)
  results <- data.frame(ok)
  data <- data.frame(results)
  m = ncol(data)
  ACC  <- matrix(nrow = m, ncol = 1)
  SENS  <- matrix(nrow = m, ncol = 1)
  SPEC  <-matrix(nrow = m, ncol = 1)
  MCC <- matrix(nrow = m, ncol = 1)
  
  for(i in 1:m){ 
    ACC[i,1]  = (data[1,i]+data[4,i])/(data[1,i]+data[2,i]+data[3,i]+data[4,i])*100
    SENS[i,1]  =  (data[4,i])/(data[3,i]+data[4,i])*100
    SPEC[i,1]  = (data[1,i]/(data[1,i]+data[2,i]))*100
    MCC1      = (data[1,i]*data[4,i]) - (data[2,i]*data[3,i])
    MCC2      =  (data[4,i]+data[2,i])*(data[4,i]+data[3,i])
    MCC3      =  (data[1,i]+data[2,i])*(data[1,i]+data[3,i])
    MCC4  =  sqrt(MCC2)*sqrt(MCC3)
    
    
    MCC[i,1]  = MCC1/MCC4
  }
  results_ACC <- mean_and_sd(ACC)
  results_SENS <- mean_and_sd(SENS)
  results_SPEC <- mean_and_sd(SPEC)
  results_MCC <- mean_and_sd(MCC)
  results_all <- (data.frame(c(results_ACC, results_SENS, results_SPEC, results_MCC)))
  rownames(results_all) <- c("ACC_Mean", "ACC_SD", "Sens_Mean", "Sens_SD", "Spec_Mean", "Spec_SD",
                             "MCC_Mean", "MCC_SD")
  return(results_all)
}

#### 10-fold results using J48
J48_10fold <- function(x) {
  results <- list(100)
  for (i in 1:100) {
    in_train <- caret::createDataPartition(x$Label, p = 0.80, list = FALSE)
    train <- x[in_train, ]
    test <- x[-in_train, ]
    model_train <- RWeka::J48(Label~., data = train)
    eval_j48 <- RWeka::evaluate_Weka_classifier(model_train, numFolds = 10, complexity = FALSE, seed = 1, class = TRUE)
    confusionmatrix <- eval_j48$confusionMatrix
    results[[i]] <- as.numeric(confusionmatrix)
  }
  return(results)
}

J48_cross_validation <- function(x) {
  ok <- J48_10fold(x)
  results <- data.frame(ok)
  data <- data.frame(results)
  m = ncol(data)
  ACC  <- matrix(nrow = m, ncol = 1)
  SENS  <- matrix(nrow = m, ncol = 1)
  SPEC  <-matrix(nrow = m, ncol = 1)
  MCC <- matrix(nrow = m, ncol = 1)
  
  for(i in 1:m){ 
    ACC[i,1]  = (data[1,i]+data[4,i])/(data[1,i]+data[2,i]+data[3,i]+data[4,i])*100
    SENS[i,1]  =  (data[4,i])/(data[3,i]+data[4,i])*100
    SPEC[i,1]  = (data[1,i]/(data[1,i]+data[2,i]))*100
    MCC1      = (data[1,i]*data[4,i]) - (data[2,i]*data[3,i])
    MCC2      =  (data[4,i]+data[2,i])*(data[4,i]+data[3,i])
    MCC3      =  (data[1,i]+data[2,i])*(data[1,i]+data[3,i])
    MCC4  =  sqrt(MCC2)*sqrt(MCC3)
    
    
    MCC[i,1]  = MCC1/MCC4
  }
  results_ACC <- mean_and_sd(ACC)
  results_SENS <- mean_and_sd(SENS)
  results_SPEC <- mean_and_sd(SPEC)
  results_MCC <- mean_and_sd(MCC)
  results_all <- (data.frame(c(results_ACC, results_SENS, results_SPEC, results_MCC)))
  rownames(results_all) <- c("ACC_Mean", "ACC_SD", "Sens_Mean", "Sens_SD", "Spec_Mean", "Spec_SD",
                             "MCC_Mean", "MCC_SD")
  return(results_all)
}

### testing results using J48
J48_testing <- function(x) {
  results <- list(100)
  for (i in 1:100) {
    in_train <- caret::createDataPartition(x$Label, p = 0.80, list = FALSE)
    train <- x[in_train, ]
    test <- x[-in_train, ]
    model_train <- RWeka::J48(Label~., data = train)
    eval_external <- RWeka::evaluate_Weka_classifier(model_train, newdata = test, numFolds = 0, complexity = FALSE, seed = 1, class = TRUE)
    confusionmatrix <- eval_external$confusionMatrix
    results[[i]] <- as.numeric(confusionmatrix)
  }
  return(results)
}


J48_external <- function(x) {
  ok <- J48_testing(x)
  results <- data.frame(ok)
  data <- data.frame(results)
  m = ncol(data)
  ACC  <- matrix(nrow = m, ncol = 1)
  SENS  <- matrix(nrow = m, ncol = 1)
  SPEC  <-matrix(nrow = m, ncol = 1)
  MCC <- matrix(nrow = m, ncol = 1)
  
  for(i in 1:m){ 
    ACC[i,1]  = (data[1,i]+data[4,i])/(data[1,i]+data[2,i]+data[3,i]+data[4,i])*100
    SENS[i,1]  =  (data[4,i])/(data[3,i]+data[4,i])*100
    SPEC[i,1]  = (data[1,i]/(data[1,i]+data[2,i]))*100
    MCC1      = (data[1,i]*data[4,i]) - (data[2,i]*data[3,i])
    MCC2      =  (data[4,i]+data[2,i])*(data[4,i]+data[3,i])
    MCC3      =  (data[1,i]+data[2,i])*(data[1,i]+data[3,i])
    MCC4  =  sqrt(MCC2)*sqrt(MCC3)
    
    
    MCC[i,1]  = MCC1/MCC4
  }
  results_ACC <- mean_and_sd(ACC)
  results_SENS <- mean_and_sd(SENS)
  results_SPEC <- mean_and_sd(SPEC)
  results_MCC <- mean_and_sd(MCC)
  results_all <- (data.frame(c(results_ACC, results_SENS, results_SPEC, results_MCC)))
  rownames(results_all) <- c("ACC_Mean", "ACC_SD", "Sens_Mean", "Sens_SD", "Spec_Mean", "Spec_SD",
                             "MCC_Mean", "MCC_SD")
  return(results_all)
}

results_J48 <- function(x) {
  training <- J48_train(x)
  cross_validation <- J48_cross_validation(x)
  testing <- J48_external(x)
  results <- data.frame(Training = training, Cross_Validation = cross_validation, Testing = testing)
  colnames(results) <- c("Training set", "Internal validation", "External set")
  return(results)
}

suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(doSNOW))
cl <- makeCluster(8)
registerDoSNOW(cl)
clusterExport(cl = cl, ls())

result_J48_performance <- parLapply(cl = cl, input, function(x) {
  models <- suppressWarnings(results_J48(x))
  return(models)
})
print(result_J48_performance)
stopCluster(cl)

```

#####HDPs having anti--bacteria, anti--cancer, anti--fungus and virus-- were combined into one to represent HDPs while the negative set was used as non--HDPs#####

```
### read FASTA
library(protr)
cancer <- readFASTA("cancer.fasta")
fungus <- readFASTA("fungus.fasta")
bacteria <- readFASTA("bacteria.fasta")
virus <- readFASTA("virus.fasta")
negative <- protr::readFASTA("negative_Protein.fasta")

### removed wired protein
cancer <- cancer[(sapply(cancer, protcheck))]
fungus <- fungus[(sapply(fungus, protcheck))]
bacteria <- bacteria[(sapply(bacteria, protcheck))]
virus <- virus[(sapply(virus, protcheck))]
negative <- negative[(sapply(negative, protcheck))]


### Descriptors generation
cancer_des <- t(sapply(cancer, extractAAC))
fungus_des <- t(sapply(fungus, extractAAC))
bacteria_des <- t(sapply(bacteria, extractAAC))
virus_des <- t(sapply(virus, extractAAC))
negative_des <- t(sapply(negative, extractAAC))
#### label the labels 
cancer_des <- as.data.frame(cancer_des)
fungus_des <- as.data.frame(fungus_des)
bacteria_des <- as.data.frame(bacteria_des)
virus_des <- as.data.frame(virus_des)
hdp <- rbind(cancer_des, fungus_des, bacteria_des, virus_des)
hdp$Label <- "Positive"
hdp <- as.data.frame(hdp)
negative_des <- as.data.frame(negative_des)
negative_des$Label <- "Negative"

data <- rbind(hdp, negative_des)
data$Label <- as.factor(data$Label)


results_J48 <- function(x) {
  training <- J48_train(x)
  cross_validation <- J48_cross_validation(x)
  testing <- J48_external(x)
  results <- data.frame(Training = training, Cross_Validation = cross_validation, Testing = testing)
  colnames(results) <- c("Training set", "Internal validation", "External set")
  return(results)
}

result_J48_performance <- suppressMessages(results_J48(data))
print(result_J48_performance)

```

#####HDPs having their respective activity were combined together to build a union predictive model#####

```
library(protr)
cancer <- readFASTA("cancer.fasta")
fungus <- readFASTA("fungus.fasta")
bacteria <- readFASTA("bacteria.fasta")
virus <- readFASTA("virus.fasta")
negative <- protr::readFASTA("negative_Protein.fasta")

### removed wired protein
cancer <- cancer[(sapply(cancer, protcheck))]
fungus <- fungus[(sapply(fungus, protcheck))]
bacteria <- bacteria[(sapply(bacteria, protcheck))]
virus <- virus[(sapply(virus, protcheck))]
negative <- negative[(sapply(negative, protcheck))]

### function for amino acid composition
composition <- function(x) {
  library(protr)
  c(extractAAC(x))
}
### generation the Descriptors
cancer_des <- t(sapply(cancer, extractAAC))
fungus_des <- t(sapply(fungus, extractAAC))
bacteria_des <- t(sapply(bacteria, extractAAC))
virus_des <- t(sapply(virus, extractAAC))
negative_des <- t(sapply(negative, extractAAC))
#### label the labels 
cancer_des <- as.data.frame(cancer_des)
cancer_des$Label <- "Cancer"
fungus_des <- as.data.frame(fungus_des)
fungus_des$Label <- "Fungus"
bacteria_des <- as.data.frame(bacteria_des)
bacteria_des$Label <- "Bacteria"
virus_des <- as.data.frame(virus_des)
virus_des$Label <- "Virus"
combine_data <- rbind(cancer_des, fungus_des,
                      bacteria_des, virus_des)
combine_data$Label <- as.factor(combine_data$Label)

#### training results using J48
J48_training <- function(x, Label){
  if (Label == "Bacteria") {
    suppressPackageStartupMessages(library(parallel))
    suppressPackageStartupMessages(library(doSNOW))
    cl <- makeCluster(8)
    registerDoSNOW(cl)
    
    ok <- list(100)
    ok <- foreach(i = 1:100) %dopar% { 
      in_train <- caret::createDataPartition(x$Label, p = 0.80, list = FALSE)
      train <- x[in_train, ]
      test <- x[-in_train, ]
      rm(in_train)
      rm(test)
      model_train <- RWeka::J48(Label~., data = train)
      actual <- train$Label
      prediction <- predict(model_train, train)
      rm(train)
      rm(model_train)
      results <- caret::confusionMatrix(prediction, actual)
      results <- results$table
      results <- table(prediction, actual)
      rm(prediction)
      rm(actual)
      results <- as.numeric(results)
      ok[[i]] <- cbind(results[[1]], (results[[5]] + results[[9]] + results[[13]]),
                        (results[[2]] + results[[3]] + results[[4]]), (results[[6]] + results[[11]] + results[[16]]))
      Cancer <- cbind(results[[6]], (results[[2]] + results[[10]] + results[[14]]), 
                      (results[[5]] + results[[7]] + results[[8]]), (results[[1]] + results[[11]] + results[[16]]))
      Fungus <- cbind(results[[11]], (results[[3]] + results[[7]] + results[[15]]),
                      (results[[9]] + results[[10]] + results[[12]]), (results[[1]] + results[[6]] + results[[16]]))
      Virus <- cbind(results[[16]], (results[[4]] + results[[8]] + results[[12]]), 
                     (results[[13]] + results[[14]] + results[[15]]), (results[[1]] + results[[6]] + results[[11]]))
    }
  }  else if (Label == "Cancer") {
    cl <- makeCluster(8)
    registerDoSNOW(cl)
    
    ok <- list(100)
    ok <- foreach(i = 1:100) %dopar% { 
      in_train <- caret::createDataPartition(x$Label, p = 0.80, list = FALSE)
      train <- x[in_train, ]
      test <- x[-in_train, ]
      rm(in_train)
      rm(test)
      model_train <- RWeka::J48(Label~., data = train)
      actual <- train$Label
      prediction <- predict(model_train, train)
      rm(train)
      rm(model_train)
      results <- caret::confusionMatrix(prediction, actual)
      results <- results$table
      results <- table(prediction, actual)
      results <- as.numeric(results)
      rm(prediction)
      rm(actual)
      Bacteria <- cbind(results[[1]], (results[[5]] + results[[9]] + results[[13]]),
                        (results[[2]] + results[[3]] + results[[4]]), (results[[6]] + results[[11]] + results[[16]]))
      ok[[i]] <- cbind(results[[6]], (results[[2]] + results[[10]] + results[[14]]), 
                      (results[[5]] + results[[7]] + results[[8]]), (results[[1]] + results[[11]] + results[[16]]))
      Fungus <- cbind(results[[11]], (results[[3]] + results[[7]] + results[[15]]),
                      (results[[9]] + results[[10]] + results[[12]]), (results[[1]] + results[[6]] + results[[16]]))
      Virus <- cbind(results[[16]], (results[[4]] + results[[8]] + results[[12]]), 
                     (results[[13]] + results[[14]] + results[[15]]), (results[[1]] + results[[6]] + results[[11]]))    
      } 
  }  else if (Label == "Fungus") {
    cl <- makeCluster(8)
    registerDoSNOW(cl)
    
    ok <- list(100)
    ok <- foreach(i = 1:100) %dopar% { 
      in_train <- caret::createDataPartition(x$Label, p = 0.80, list = FALSE)
      train <- x[in_train, ]
      test <- x[-in_train, ]
      rm(test)
      rm(in_train)
      model_train <- RWeka::J48(Label~., data = train)
      actual <- train$Label
      prediction <- predict(model_train, train)
      rm(model_train)
      rm(train)
      results <- caret::confusionMatrix(prediction, actual)
      results <- results$table
      results <- table(prediction, actual)
      results <- as.numeric(results)
      rm(actual)
      rm(prediction)
      
      Bacteria <- cbind(results[[1]], (results[[5]] + results[[9]] + results[[13]]),
                        (results[[2]] + results[[3]] + results[[4]]), (results[[6]] + results[[11]] + results[[16]]))
      Cancer <- cbind(results[[6]], (results[[2]] + results[[10]] + results[[14]]), 
                      (results[[5]] + results[[7]] + results[[8]]), (results[[1]] + results[[11]] + results[[16]]))
      ok[[i]] <- cbind(results[[11]], (results[[3]] + results[[7]] + results[[15]]),
                      (results[[9]] + results[[10]] + results[[12]]), (results[[1]] + results[[6]] + results[[16]]))
      Virus <- cbind(results[[16]], (results[[4]] + results[[8]] + results[[12]]), 
                     (results[[13]] + results[[14]] + results[[15]]), (results[[1]] + results[[6]] + results[[11]]))    
    } 
  }  else if (Label == "Virus") {
    cl <- makeCluster(8)
    registerDoSNOW(cl)
    
    ok <- list(100)
    ok <- foreach(i = 1:100) %dopar% { 
      in_train <- caret::createDataPartition(x$Label, p = 0.80, list = FALSE)
      train <- x[in_train, ]
      test <- x[-in_train, ]
      rm(in_train)
      rm(test)
      model_train <- RWeka::J48(Label~., data = train)
      actual <- train$Label
      prediction <- predict(model_train, train)
      rm(model_train)
      rm(train)
      results <- caret::confusionMatrix(prediction, actual)
      results <- results$table
      results <- table(prediction, actual)
      results <- as.numeric(results)
      rm(prediction)
      rm(actual)
      Bacteria <- cbind(results[[1]], (results[[5]] + results[[9]] + results[[13]]),
                        (results[[2]] + results[[3]] + results[[4]]), (results[[6]] + results[[11]] + results[[16]]))
      Cancer <- cbind(results[[6]], (results[[2]] + results[[10]] + results[[14]]), 
                      (results[[5]] + results[[7]] + results[[8]]), (results[[1]] + results[[11]] + results[[16]]))
      Fungus <- cbind(results[[11]], (results[[3]] + results[[7]] + results[[15]]),
                      (results[[9]] + results[[10]] + results[[12]]), (results[[1]] + results[[6]] + results[[16]]))
      ok[[i]] <- cbind(results[[16]], (results[[4]] + results[[8]] + results[[12]]), 
                     (results[[13]] + results[[14]] + results[[15]]), (results[[1]] + results[[6]] + results[[11]]))    
    }
    return(ok)
    stopCluster(cl)
  } }

mean_and_sd <- function(x) {
  c(round(mean(x, na.rm = TRUE), digits = 4),
    round(sd(x, na.rm = TRUE), digits = 4))
}


results_training_Bacteria <- function(x) {
  yes <- J48_training(x, Label = "Bacteria")
  great <- data.frame(yes)
  rm(yes)
  TP <- seq(from = 1, to = 400, by = 4)
  FN <- seq(from = 2, to = 400, by = 4)
  FP <- seq(from = 3, to = 400, by = 4)
  TN <- seq(from = 4, to = 400, by = 4)
  results <- mapply(c, great[TP], great[FN], great[FP], great[TN])
  rm(great)
  data <- data.frame(results)
  rm(results)
  m = ncol(data)
  ACC  <- matrix(nrow = m, ncol = 1)
  SENS  <- matrix(nrow = m, ncol = 1)
  SPEC  <-matrix(nrow = m, ncol = 1)
  MCC <- matrix(nrow = m, ncol = 1)
  
  for(i in 1:m){ 
    ACC[i,1]  = (data[1,i]+data[4,i])/(data[1,i]+data[2,i]+data[3,i]+data[4,i])*100
    SENS[i,1]  =  (data[4,i])/(data[3,i]+data[4,i])*100
    SPEC[i,1]  = (data[1,i]/(data[1,i]+data[2,i]))*100
    MCC1      = (data[1,i]*data[4,i]) - (data[2,i]*data[3,i])
    MCC2      =  (data[4,i]+data[2,i])*(data[4,i]+data[3,i])
    MCC3      =  (data[1,i]+data[2,i])*(data[1,i]+data[3,i])
    MCC4  =  sqrt(MCC2)*sqrt(MCC3)
    
    
    MCC[i,1]  = MCC1/MCC4
  }
  rm(TP)
  rm(FP)
  rm(TN)
  rm(FN)
  results_ACC <- mean_and_sd(ACC)
  results_SENS <- mean_and_sd(SENS)
  results_SPEC <- mean_and_sd(SPEC)
  results_MCC <- mean_and_sd(MCC)
  rm(ACC)
  rm(SENS)
  rm(SPEC)
  rm(MCC)
  results_all <- (data.frame(c(results_ACC, results_SENS, results_SPEC, results_MCC)))
  rownames(results_all) <- c("ACC_Mean", "ACC_SD", "Sens_Mean", "Sens_SD", "Spec_Mean", "Spec_SD",
                             "MCC_Mean", "MCC_SD")
  rm(results_ACC)
  rm(results_SENS)
  rm(results_SPEC)
  rm(results_MCC)
  return(results_all)
}

results_training_Cancer <- function(x) {
  yes <- J48_training(x, Label = "Cancer")
  great <- data.frame(yes)
  rm(yes)
  TP <- seq(from = 1, to = 400, by = 4)
  FN <- seq(from = 2, to = 400, by = 4)
  FP <- seq(from = 3, to = 400, by = 4)
  TN <- seq(from = 4, to = 400, by = 4)
  results <- mapply(c, great[TP], great[FN], great[FP], great[TN])
  rm(great)
  data <- data.frame(results)
  rm(results)
  m = ncol(data)
  ACC  <- matrix(nrow = m, ncol = 1)
  SENS  <- matrix(nrow = m, ncol = 1)
  SPEC  <-matrix(nrow = m, ncol = 1)
  MCC <- matrix(nrow = m, ncol = 1)
  
  for(i in 1:m){ 
    ACC[i,1]  = (data[1,i]+data[4,i])/(data[1,i]+data[2,i]+data[3,i]+data[4,i])*100
    SENS[i,1]  =  (data[4,i])/(data[3,i]+data[4,i])*100
    SPEC[i,1]  = (data[1,i]/(data[1,i]+data[2,i]))*100
    MCC1      = (data[1,i]*data[4,i]) - (data[2,i]*data[3,i])
    MCC2      =  (data[4,i]+data[2,i])*(data[4,i]+data[3,i])
    MCC3      =  (data[1,i]+data[2,i])*(data[1,i]+data[3,i])
    MCC4  =  sqrt(MCC2)*sqrt(MCC3)
    
    
    MCC[i,1]  = MCC1/MCC4
  }
  rm(TP)
  rm(FP)
  rm(TN)
  rm(FN)
  results_ACC <- mean_and_sd(ACC)
  results_SENS <- mean_and_sd(SENS)
  results_SPEC <- mean_and_sd(SPEC)
  results_MCC <- mean_and_sd(MCC)
  rm(ACC)
  rm(SENS)
  rm(SPEC)
  rm(MCC)
  results_all <- (data.frame(c(results_ACC, results_SENS, results_SPEC, results_MCC)))
  rownames(results_all) <- c("ACC_Mean", "ACC_SD", "Sens_Mean", "Sens_SD", "Spec_Mean", "Spec_SD",
                             "MCC_Mean", "MCC_SD")
  rm(results_ACC)
  rm(results_SENS)
  rm(results_SPEC)
  rm(results_MCC)
  return(results_all)
}



results_training_Fungus <- function(x) {
  yes <- J48_training(x, Label = "Fungus")
  great <- data.frame(yes)
  rm(yes)
  TP <- seq(from = 1, to = 400, by = 4)
  FN <- seq(from = 2, to = 400, by = 4)
  FP <- seq(from = 3, to = 400, by = 4)
  TN <- seq(from = 4, to = 400, by = 4)
  results <- mapply(c, great[TP], great[FN], great[FP], great[TN])
  rm(great)
  data <- data.frame(results)
  rm(results)
  m = ncol(data)
  ACC  <- matrix(nrow = m, ncol = 1)
  SENS  <- matrix(nrow = m, ncol = 1)
  SPEC  <-matrix(nrow = m, ncol = 1)
  MCC <- matrix(nrow = m, ncol = 1)
  
  for(i in 1:m){ 
    ACC[i,1]  = (data[1,i]+data[4,i])/(data[1,i]+data[2,i]+data[3,i]+data[4,i])*100
    SENS[i,1]  =  (data[4,i])/(data[3,i]+data[4,i])*100
    SPEC[i,1]  = (data[1,i]/(data[1,i]+data[2,i]))*100
    MCC1      = (data[1,i]*data[4,i]) - (data[2,i]*data[3,i])
    MCC2      =  (data[4,i]+data[2,i])*(data[4,i]+data[3,i])
    MCC3      =  (data[1,i]+data[2,i])*(data[1,i]+data[3,i])
    MCC4  =  sqrt(MCC2)*sqrt(MCC3)
    
    
    MCC[i,1]  = MCC1/MCC4
  }
  rm(TP)
  rm(FP)
  rm(TN)
  rm(FN)
  results_ACC <- mean_and_sd(ACC)
  results_SENS <- mean_and_sd(SENS)
  results_SPEC <- mean_and_sd(SPEC)
  results_MCC <- mean_and_sd(MCC)
  rm(ACC)
  rm(SENS)
  rm(SPEC)
  rm(MCC)
  results_all <- (data.frame(c(results_ACC, results_SENS, results_SPEC, results_MCC)))
  rownames(results_all) <- c("ACC_Mean", "ACC_SD", "Sens_Mean", "Sens_SD", "Spec_Mean", "Spec_SD",
                             "MCC_Mean", "MCC_SD")
  rm(results_ACC)
  rm(results_SENS)
  rm(results_SPEC)
  rm(results_MCC)
  return(results_all)
}


results_training_Virus <- function(x) {
  yes <- J48_training(x, Label = "Virus")
  great <- data.frame(yes)
  rm(yes)
  TP <- seq(from = 1, to = 400, by = 4)
  FN <- seq(from = 2, to = 400, by = 4)
  FP <- seq(from = 3, to = 400, by = 4)
  TN <- seq(from = 4, to = 400, by = 4)
  results <- mapply(c, great[TP], great[FN], great[FP], great[TN])
  rm(great)
  data <- data.frame(results)
  rm(results)
  m = ncol(data)
  ACC  <- matrix(nrow = m, ncol = 1)
  SENS  <- matrix(nrow = m, ncol = 1)
  SPEC  <-matrix(nrow = m, ncol = 1)
  MCC <- matrix(nrow = m, ncol = 1)
  
  for(i in 1:m){ 
    ACC[i,1]  = (data[1,i]+data[4,i])/(data[1,i]+data[2,i]+data[3,i]+data[4,i])*100
    SENS[i,1]  =  (data[4,i])/(data[3,i]+data[4,i])*100
    SPEC[i,1]  = (data[1,i]/(data[1,i]+data[2,i]))*100
    MCC1      = (data[1,i]*data[4,i]) - (data[2,i]*data[3,i])
    MCC2      =  (data[4,i]+data[2,i])*(data[4,i]+data[3,i])
    MCC3      =  (data[1,i]+data[2,i])*(data[1,i]+data[3,i])
    MCC4  =  sqrt(MCC2)*sqrt(MCC3)
    
    
    MCC[i,1]  = MCC1/MCC4
  }
  rm(TP)
  rm(FP)
  rm(TN)
  rm(FN)
  results_ACC <- mean_and_sd(ACC)
  results_SENS <- mean_and_sd(SENS)
  results_SPEC <- mean_and_sd(SPEC)
  results_MCC <- mean_and_sd(MCC)
  rm(ACC)
  rm(SENS)
  rm(SPEC)
  rm(MCC)
  results_all <- (data.frame(c(results_ACC, results_SENS, results_SPEC, results_MCC)))
  rownames(results_all) <- c("ACC_Mean", "ACC_SD", "Sens_Mean", "Sens_SD", "Spec_Mean", "Spec_SD",
                             "MCC_Mean", "MCC_SD")
  rm(results_ACC)
  rm(results_SENS)
  rm(results_SPEC)
  rm(results_MCC)
  return(results_all)
}


J48_training_all <- function(x) {
  bacteria <- results_training_Bacteria(x)
  cancer <- results_training_Cancer(x)
  fungus <- results_training_Fungus(x)
  virus <- results_training_Virus(x)
  results_all <- cbind(bacteria, cancer, fungus, virus)
  rm(bacteria)
  rm(cancer)
  rm(fungus)
  rm(virus)
  total <- apply(results_all, 1, mean)
  results_all_mean <- cbind(results_all, total)
  rm(results_all)
  rm(total)
  colnames(results_all_mean) <- c("Bacteria", "Cancer", "Fungus", "Virus", "Overall")
  return(results_all_mean)
}


### 10 fold Cross Validations

J48_10_CV <- function(x, Label){
  if (Label == "Bacteria") {
    suppressPackageStartupMessages(library(parallel))
    suppressPackageStartupMessages(library(doSNOW))
    cl <- makeCluster(8)
    registerDoSNOW(cl)
    
    ok <- list(100)
    ok <- foreach(i = 1:100) %dopar% { 
      in_train <- caret::createDataPartition(x$Label, p = 0.80, list = FALSE)
      train <- x[in_train, ]
      test <- x[-in_train, ]
      rm(in_train)
      rm(test)
      model_train <- RWeka::J48(Label~., data = train)
      results <- RWeka::evaluate_Weka_classifier(model_train, newata = NULL, numFolds = 10, complexity = FALSE)
      rm(train)
      confusionMatrix <- results$confusionMatrix
      rm(model_train)
      results <- as.numeric(confusionMatrix)
      rm(confusionMatrix)
      ok[[i]] <- cbind(results[[1]], (results[[5]] + results[[9]] + results[[13]]),
                       (results[[2]] + results[[3]] + results[[4]]), (results[[6]] + results[[11]] + results[[16]]))
      Cancer <- cbind(results[[6]], (results[[2]] + results[[10]] + results[[14]]), 
                      (results[[5]] + results[[7]] + results[[8]]), (results[[1]] + results[[11]] + results[[16]]))
      Fungus <- cbind(results[[11]], (results[[3]] + results[[7]] + results[[15]]),
                      (results[[9]] + results[[10]] + results[[12]]), (results[[1]] + results[[6]] + results[[16]]))
      Virus <- cbind(results[[16]], (results[[4]] + results[[8]] + results[[12]]), 
                     (results[[13]] + results[[14]] + results[[15]]), (results[[1]] + results[[6]] + results[[11]]))
    }
  }  else if (Label == "Cancer") {
    cl <- makeCluster(8)
    registerDoSNOW(cl)
    
    ok <- list(100)
    ok <- foreach(i = 1:100) %dopar% { 
      in_train <- caret::createDataPartition(x$Label, p = 0.80, list = FALSE)
      train <- x[in_train, ]
      test <- x[-in_train, ]
      rm(in_train)
      rm(test)
      model_train <- RWeka::J48(Label~., data = train)
      results <- RWeka::evaluate_Weka_classifier(model_train, newata = NULL, numFolds = 10, complexity = FALSE)
      rm(train)
      confusionMatrix <- results$confusionMatrix
      rm(model_train)
      results <- as.numeric(confusionMatrix)
      rm(confusionMatrix)
      Bacteria <- cbind(results[[1]], (results[[5]] + results[[9]] + results[[13]]),
                        (results[[2]] + results[[3]] + results[[4]]), (results[[6]] + results[[11]] + results[[16]]))
      ok[[i]] <- cbind(results[[6]], (results[[2]] + results[[10]] + results[[14]]), 
                       (results[[5]] + results[[7]] + results[[8]]), (results[[1]] + results[[11]] + results[[16]]))
      Fungus <- cbind(results[[11]], (results[[3]] + results[[7]] + results[[15]]),
                      (results[[9]] + results[[10]] + results[[12]]), (results[[1]] + results[[6]] + results[[16]]))
      Virus <- cbind(results[[16]], (results[[4]] + results[[8]] + results[[12]]), 
                     (results[[13]] + results[[14]] + results[[15]]), (results[[1]] + results[[6]] + results[[11]]))    
    } 
  }  else if (Label == "Fungus") {
    cl <- makeCluster(8)
    registerDoSNOW(cl)
    
    ok <- list(100)
    ok <- foreach(i = 1:100) %dopar% { 
      in_train <- caret::createDataPartition(x$Label, p = 0.80, list = FALSE)
      train <- x[in_train, ]
      test <- x[-in_train, ]
      rm(in_train)
      rm(test)
      model_train <- RWeka::J48(Label~., data = train)
      results <- RWeka::evaluate_Weka_classifier(model_train, newata = NULL, numFolds = 10, complexity = FALSE)
      rm(train)
      rm(model_train)
      confusionMatrix <- results$confusionMatrix
      rm(results)
      results <- as.numeric(confusionMatrix)
      rm(confusionMatrix)
      Bacteria <- cbind(results[[1]], (results[[5]] + results[[9]] + results[[13]]),
                        (results[[2]] + results[[3]] + results[[4]]), (results[[6]] + results[[11]] + results[[16]]))
      Cancer <- cbind(results[[6]], (results[[2]] + results[[10]] + results[[14]]), 
                      (results[[5]] + results[[7]] + results[[8]]), (results[[1]] + results[[11]] + results[[16]]))
      ok[[i]] <- cbind(results[[11]], (results[[3]] + results[[7]] + results[[15]]),
                       (results[[9]] + results[[10]] + results[[12]]), (results[[1]] + results[[6]] + results[[16]]))
      Virus <- cbind(results[[16]], (results[[4]] + results[[8]] + results[[12]]), 
                     (results[[13]] + results[[14]] + results[[15]]), (results[[1]] + results[[6]] + results[[11]]))    
    } 
  }  else if (Label == "Virus") {
    cl <- makeCluster(8)
    registerDoSNOW(cl)
    
    ok <- list(100)
    ok <- foreach(i = 1:100) %dopar% { 
      in_train <- caret::createDataPartition(x$Label, p = 0.80, list = FALSE)
      train <- x[in_train, ]
      test <- x[-in_train, ]
      rm(in_train)
      rm(test)
      model_train <- RWeka::J48(Label~., data = train)
      results <- RWeka::evaluate_Weka_classifier(model_train, newata = NULL, numFolds = 10, complexity = FALSE)
      rm(model_train)
      rm(train)
      confusionMatrix <- results$confusionMatrix
      results <- as.numeric(confusionMatrix)
      rm(confusionMatrix)
      Bacteria <- cbind(results[[1]], (results[[5]] + results[[9]] + results[[13]]),
                        (results[[2]] + results[[3]] + results[[4]]), (results[[6]] + results[[11]] + results[[16]]))
      Cancer <- cbind(results[[6]], (results[[2]] + results[[10]] + results[[14]]), 
                      (results[[5]] + results[[7]] + results[[8]]), (results[[1]] + results[[11]] + results[[16]]))
      Fungus <- cbind(results[[11]], (results[[3]] + results[[7]] + results[[15]]),
                      (results[[9]] + results[[10]] + results[[12]]), (results[[1]] + results[[6]] + results[[16]]))
      ok[[i]] <- cbind(results[[16]], (results[[4]] + results[[8]] + results[[12]]), 
                       (results[[13]] + results[[14]] + results[[15]]), (results[[1]] + results[[6]] + results[[11]]))    
    }
    return(ok)
    stopCluster(cl)
  } }

mean_and_sd <- function(x) {
  c(round(mean(x, na.rm = TRUE), digits = 4),
    round(sd(x, na.rm = TRUE), digits = 4))
}


results_CV_Bacteria <- function(x) {
  yes <- J48_10_CV(x, Label = "Bacteria")
  great <- data.frame(yes)
  rm(yes)
  TP <- seq(from = 1, to = 400, by = 4)
  FN <- seq(from = 2, to = 400, by = 4)
  FP <- seq(from = 3, to = 400, by = 4)
  TN <- seq(from = 4, to = 400, by = 4)
  results <- mapply(c, great[TP], great[FN], great[FP], great[TN])
  data <- data.frame(results)
  rm(great)
  rm(results)
  m = ncol(data)
  ACC  <- matrix(nrow = m, ncol = 1)
  SENS  <- matrix(nrow = m, ncol = 1)
  SPEC  <-matrix(nrow = m, ncol = 1)
  MCC <- matrix(nrow = m, ncol = 1)
  
  for(i in 1:m){ 
    ACC[i,1]  = (data[1,i]+data[4,i])/(data[1,i]+data[2,i]+data[3,i]+data[4,i])*100
    SENS[i,1]  =  (data[4,i])/(data[3,i]+data[4,i])*100
    SPEC[i,1]  = (data[1,i]/(data[1,i]+data[2,i]))*100
    MCC1      = (data[1,i]*data[4,i]) - (data[2,i]*data[3,i])
    MCC2      =  (data[4,i]+data[2,i])*(data[4,i]+data[3,i])
    MCC3      =  (data[1,i]+data[2,i])*(data[1,i]+data[3,i])
    MCC4  =  sqrt(MCC2)*sqrt(MCC3)
    
    
    MCC[i,1]  = MCC1/MCC4
  }
  rm(TP)
  rm(FP)
  rm(TN)
  rm(FN)
  results_ACC <- mean_and_sd(ACC)
  results_SENS <- mean_and_sd(SENS)
  results_SPEC <- mean_and_sd(SPEC)
  results_MCC <- mean_and_sd(MCC)
  rm(data)
  rm(ACC)
  rm(SENS)
  rm(SPEC)
  rm(MCC)
  results_all <- (data.frame(c(results_ACC, results_SENS, results_SPEC, results_MCC)))
  rownames(results_all) <- c("ACC_Mean", "ACC_SD", "Sens_Mean", "Sens_SD", "Spec_Mean", "Spec_SD",
                             "MCC_Mean", "MCC_SD")
  rm(results_ACC)
  rm(results_SENS)
  rm(results_SPEC)
  rm(results_MCC)
  return(results_all)
}

results_CV_Cancer <- function(x) {
  yes <- J48_10_CV(x, Label = "Cancer")
  great <- data.frame(yes)
  rm(yes)
  TP <- seq(from = 1, to = 400, by = 4)
  FN <- seq(from = 2, to = 400, by = 4)
  FP <- seq(from = 3, to = 400, by = 4)
  TN <- seq(from = 4, to = 400, by = 4)
  results <- mapply(c, great[TP], great[FN], great[FP], great[TN])
  data <- data.frame(results)
  rm(great)
  rm(results)
  m = ncol(data)
  ACC  <- matrix(nrow = m, ncol = 1)
  SENS  <- matrix(nrow = m, ncol = 1)
  SPEC  <-matrix(nrow = m, ncol = 1)
  MCC <- matrix(nrow = m, ncol = 1)
  
  for(i in 1:m){ 
    ACC[i,1]  = (data[1,i]+data[4,i])/(data[1,i]+data[2,i]+data[3,i]+data[4,i])*100
    SENS[i,1]  =  (data[4,i])/(data[3,i]+data[4,i])*100
    SPEC[i,1]  = (data[1,i]/(data[1,i]+data[2,i]))*100
    MCC1      = (data[1,i]*data[4,i]) - (data[2,i]*data[3,i])
    MCC2      =  (data[4,i]+data[2,i])*(data[4,i]+data[3,i])
    MCC3      =  (data[1,i]+data[2,i])*(data[1,i]+data[3,i])
    MCC4  =  sqrt(MCC2)*sqrt(MCC3)
    
    
    MCC[i,1]  = MCC1/MCC4
  }
  rm(TP)
  rm(FP)
  rm(TN)
  rm(FN)
  results_ACC <- mean_and_sd(ACC)
  results_SENS <- mean_and_sd(SENS)
  results_SPEC <- mean_and_sd(SPEC)
  results_MCC <- mean_and_sd(MCC)
  rm(data)
  rm(ACC)
  rm(SENS)
  rm(SPEC)
  rm(MCC)
  results_all <- (data.frame(c(results_ACC, results_SENS, results_SPEC, results_MCC)))
  rownames(results_all) <- c("ACC_Mean", "ACC_SD", "Sens_Mean", "Sens_SD", "Spec_Mean", "Spec_SD",
                             "MCC_Mean", "MCC_SD")
  rm(results_ACC)
  rm(results_SENS)
  rm(results_SPEC)
  rm(results_MCC)
  return(results_all)
}


results_CV_Fungus <- function(x) {
  yes <- J48_10_CV(x, Label = "Fungus")
  great <- data.frame(yes)
  rm(yes)
  TP <- seq(from = 1, to = 400, by = 4)
  FN <- seq(from = 2, to = 400, by = 4)
  FP <- seq(from = 3, to = 400, by = 4)
  TN <- seq(from = 4, to = 400, by = 4)
  results <- mapply(c, great[TP], great[FN], great[FP], great[TN])
  data <- data.frame(results)
  rm(great)
  rm(results)
  m = ncol(data)
  ACC  <- matrix(nrow = m, ncol = 1)
  SENS  <- matrix(nrow = m, ncol = 1)
  SPEC  <-matrix(nrow = m, ncol = 1)
  MCC <- matrix(nrow = m, ncol = 1)
  
  for(i in 1:m){ 
    ACC[i,1]  = (data[1,i]+data[4,i])/(data[1,i]+data[2,i]+data[3,i]+data[4,i])*100
    SENS[i,1]  =  (data[4,i])/(data[3,i]+data[4,i])*100
    SPEC[i,1]  = (data[1,i]/(data[1,i]+data[2,i]))*100
    MCC1      = (data[1,i]*data[4,i]) - (data[2,i]*data[3,i])
    MCC2      =  (data[4,i]+data[2,i])*(data[4,i]+data[3,i])
    MCC3      =  (data[1,i]+data[2,i])*(data[1,i]+data[3,i])
    MCC4  =  sqrt(MCC2)*sqrt(MCC3)
    
    
    MCC[i,1]  = MCC1/MCC4
  }
  rm(TP)
  rm(FP)
  rm(TN)
  rm(FN)
  results_ACC <- mean_and_sd(ACC)
  results_SENS <- mean_and_sd(SENS)
  results_SPEC <- mean_and_sd(SPEC)
  results_MCC <- mean_and_sd(MCC)
  rm(data)
  rm(ACC)
  rm(SENS)
  rm(SPEC)
  rm(MCC)
  results_all <- (data.frame(c(results_ACC, results_SENS, results_SPEC, results_MCC)))
  rownames(results_all) <- c("ACC_Mean", "ACC_SD", "Sens_Mean", "Sens_SD", "Spec_Mean", "Spec_SD",
                             "MCC_Mean", "MCC_SD")
  rm(results_ACC)
  rm(results_SENS)
  rm(results_SPEC)
  rm(results_MCC)
  return(results_all)
}

results_CV_Virus <- function(x) {
  yes <- J48_10_CV(x, Label = "Virus")
  great <- data.frame(yes)
  rm(yes)
  TP <- seq(from = 1, to = 400, by = 4)
  FN <- seq(from = 2, to = 400, by = 4)
  FP <- seq(from = 3, to = 400, by = 4)
  TN <- seq(from = 4, to = 400, by = 4)
  results <- mapply(c, great[TP], great[FN], great[FP], great[TN])
  data <- data.frame(results)
  rm(great)
  rm(results)
  m = ncol(data)
  ACC  <- matrix(nrow = m, ncol = 1)
  SENS  <- matrix(nrow = m, ncol = 1)
  SPEC  <-matrix(nrow = m, ncol = 1)
  MCC <- matrix(nrow = m, ncol = 1)
  
  for(i in 1:m){ 
    ACC[i,1]  = (data[1,i]+data[4,i])/(data[1,i]+data[2,i]+data[3,i]+data[4,i])*100
    SENS[i,1]  =  (data[4,i])/(data[3,i]+data[4,i])*100
    SPEC[i,1]  = (data[1,i]/(data[1,i]+data[2,i]))*100
    MCC1      = (data[1,i]*data[4,i]) - (data[2,i]*data[3,i])
    MCC2      =  (data[4,i]+data[2,i])*(data[4,i]+data[3,i])
    MCC3      =  (data[1,i]+data[2,i])*(data[1,i]+data[3,i])
    MCC4  =  sqrt(MCC2)*sqrt(MCC3)
    
    
    MCC[i,1]  = MCC1/MCC4
  }
  rm(TP)
  rm(FP)
  rm(TN)
  rm(FN)
  results_ACC <- mean_and_sd(ACC)
  results_SENS <- mean_and_sd(SENS)
  results_SPEC <- mean_and_sd(SPEC)
  results_MCC <- mean_and_sd(MCC)
  rm(data)
  rm(ACC)
  rm(SENS)
  rm(SPEC)
  rm(MCC)
  results_all <- (data.frame(c(results_ACC, results_SENS, results_SPEC, results_MCC)))
  rownames(results_all) <- c("ACC_Mean", "ACC_SD", "Sens_Mean", "Sens_SD", "Spec_Mean", "Spec_SD",
                             "MCC_Mean", "MCC_SD")
  rm(results_ACC)
  rm(results_SENS)
  rm(results_SPEC)
  rm(results_MCC)
  return(results_all)
}

J48_CV_all <- function(x) {
  bacteria <- results_CV_Bacteria(x)
  cancer <- results_CV_Cancer(x)
  fungus <- results_CV_Fungus(x)
  virus <- results_CV_Virus(x)
  results_all <- cbind(bacteria, cancer, fungus, virus)
  rm(bacteria)
  rm(cancer)
  rm(fungus)
  rm(virus)
  total <- apply(results_all, 1, mean)
  results_all_mean <- cbind(results_all, total)
  rm(results_all)
  rm(total)
  colnames(results_all_mean) <- c("Bacteria", "Cancer", "Fungus", "Virus", "Overall")
  return(results_all_mean)
}

#### function for testing

J48_testing <- function(x, Label){
  if (Label == "Bacteria") {
    suppressPackageStartupMessages(library(parallel))
    suppressPackageStartupMessages(library(doSNOW))
    cl <- makeCluster(8)
    registerDoSNOW(cl)
    
    ok <- list(100)
    ok <- foreach(i = 1:100) %dopar% { 
      in_train <- caret::createDataPartition(x$Label, p = 0.80, list = FALSE)
      train <- x[in_train, ]
      test <- x[-in_train, ]
      rm(in_train)
      model_train <- RWeka::J48(Label~., data = train)
      rm(train)
      actual <- test$Label
      prediction <- predict(model_train, test)
      rm(model_train)
      rm(test)
      results <- caret::confusionMatrix(prediction, actual)
      results <- results$table
      results <- table(prediction, actual)
      rm(prediction)
      rm(actual)
      results <- as.numeric(results)
      ok[[i]] <- cbind(results[[1]], (results[[5]] + results[[9]] + results[[13]]),
                       (results[[2]] + results[[3]] + results[[4]]), (results[[6]] + results[[11]] + results[[16]]))
      Cancer <- cbind(results[[6]], (results[[2]] + results[[10]] + results[[14]]), 
                      (results[[5]] + results[[7]] + results[[8]]), (results[[1]] + results[[11]] + results[[16]]))
      Fungus <- cbind(results[[11]], (results[[3]] + results[[7]] + results[[15]]),
                      (results[[9]] + results[[10]] + results[[12]]), (results[[1]] + results[[6]] + results[[16]]))
      Virus <- cbind(results[[16]], (results[[4]] + results[[8]] + results[[12]]), 
                     (results[[13]] + results[[14]] + results[[15]]), (results[[1]] + results[[6]] + results[[11]]))
    }
  }  else if (Label == "Cancer") {
    cl <- makeCluster(8)
    registerDoSNOW(cl)
    
    ok <- list(100)
    ok <- foreach(i = 1:100) %dopar% { 
      in_train <- caret::createDataPartition(x$Label, p = 0.80, list = FALSE)
      train <- x[in_train, ]
      test <- x[-in_train, ]
      rm(in_train)
      model_train <- RWeka::J48(Label~., data = train)
      rm(train)
      actual <- test$Label
      prediction <- predict(model_train, test)
      rm(model_train)
      rm(test)
      results <- caret::confusionMatrix(prediction, actual)
      results <- results$table
      results <- table(prediction, actual)
      results <- as.numeric(results)
      rm(prediction)
      rm(actual)
      Bacteria <- cbind(results[[1]], (results[[5]] + results[[9]] + results[[13]]),
                        (results[[2]] + results[[3]] + results[[4]]), (results[[6]] + results[[11]] + results[[16]]))
      ok[[i]] <- cbind(results[[6]], (results[[2]] + results[[10]] + results[[14]]), 
                       (results[[5]] + results[[7]] + results[[8]]), (results[[1]] + results[[11]] + results[[16]]))
      Fungus <- cbind(results[[11]], (results[[3]] + results[[7]] + results[[15]]),
                      (results[[9]] + results[[10]] + results[[12]]), (results[[1]] + results[[6]] + results[[16]]))
      Virus <- cbind(results[[16]], (results[[4]] + results[[8]] + results[[12]]), 
                     (results[[13]] + results[[14]] + results[[15]]), (results[[1]] + results[[6]] + results[[11]]))    
    } 
  }  else if (Label == "Fungus") {
    cl <- makeCluster(8)
    registerDoSNOW(cl)
    
    ok <- list(100)
    ok <- foreach(i = 1:100) %dopar% { 
      in_train <- caret::createDataPartition(x$Label, p = 0.80, list = FALSE)
      train <- x[in_train, ]
      test <- x[-in_train, ]
      rm(in_train)
      model_train <- RWeka::J48(Label~., data = train)
      rm(train)
      actual <- test$Label
      prediction <- predict(model_train, test)
      rm(model_train)
      rm(test)
      results <- caret::confusionMatrix(prediction, actual)
      results <- results$table
      results <- table(prediction, actual)
      results <- as.numeric(results)
      rm(prediction)
      rm(actual)
      Bacteria <- cbind(results[[1]], (results[[5]] + results[[9]] + results[[13]]),
                        (results[[2]] + results[[3]] + results[[4]]), (results[[6]] + results[[11]] + results[[16]]))
      Cancer <- cbind(results[[6]], (results[[2]] + results[[10]] + results[[14]]), 
                      (results[[5]] + results[[7]] + results[[8]]), (results[[1]] + results[[11]] + results[[16]]))
      ok[[i]] <- cbind(results[[11]], (results[[3]] + results[[7]] + results[[15]]),
                       (results[[9]] + results[[10]] + results[[12]]), (results[[1]] + results[[6]] + results[[16]]))
      Virus <- cbind(results[[16]], (results[[4]] + results[[8]] + results[[12]]), 
                     (results[[13]] + results[[14]] + results[[15]]), (results[[1]] + results[[6]] + results[[11]]))    
    } 
  }  else if (Label == "Virus") {
    cl <- makeCluster(8)
    registerDoSNOW(cl)
    
    ok <- list(100)
    ok <- foreach(i = 1:100) %dopar% { 
      in_train <- caret::createDataPartition(x$Label, p = 0.80, list = FALSE)
      train <- x[in_train, ]
      test <- x[-in_train, ]
      rm(in_train)
      model_train <- RWeka::J48(Label~., data = train)
      rm(train)
      actual <- test$Label
      prediction <- predict(model_train, test)
      rm(model_train)
      rm(test)
      results <- caret::confusionMatrix(prediction, actual)
      results <- results$table
      results <- table(prediction, actual)
      results <- as.numeric(results)
      rm(prediction)
      rm(actual)
      Bacteria <- cbind(results[[1]], (results[[5]] + results[[9]] + results[[13]]),
                        (results[[2]] + results[[3]] + results[[4]]), (results[[6]] + results[[11]] + results[[16]]))
      Cancer <- cbind(results[[6]], (results[[2]] + results[[10]] + results[[14]]), 
                      (results[[5]] + results[[7]] + results[[8]]), (results[[1]] + results[[11]] + results[[16]]))
      Fungus <- cbind(results[[11]], (results[[3]] + results[[7]] + results[[15]]),
                      (results[[9]] + results[[10]] + results[[12]]), (results[[1]] + results[[6]] + results[[16]]))
      ok[[i]] <- cbind(results[[16]], (results[[4]] + results[[8]] + results[[12]]), 
                       (results[[13]] + results[[14]] + results[[15]]), (results[[1]] + results[[6]] + results[[11]]))    
    }
    return(ok)
    stopCluster(cl)
  } }

mean_and_sd <- function(x) {
  c(round(mean(x, na.rm = TRUE), digits = 4),
    round(sd(x, na.rm = TRUE), digits = 4))
}


results_testing_Bacteria <- function(x) {
  yes <- J48_testing(x, Label = "Bacteria")
  great <- data.frame(yes)
  rm(yes)
  TP <- seq(from = 1, to = 400, by = 4)
  FN <- seq(from = 2, to = 400, by = 4)
  FP <- seq(from = 3, to = 400, by = 4)
  TN <- seq(from = 4, to = 400, by = 4)
  results <- mapply(c, great[TP], great[FN], great[FP], great[TN])
  data <- data.frame(results)
  rm(great)
  rm(results)
  m = ncol(data)
  ACC  <- matrix(nrow = m, ncol = 1)
  SENS  <- matrix(nrow = m, ncol = 1)
  SPEC  <-matrix(nrow = m, ncol = 1)
  MCC <- matrix(nrow = m, ncol = 1)
  
  for(i in 1:m){ 
    ACC[i,1]  = (data[1,i]+data[4,i])/(data[1,i]+data[2,i]+data[3,i]+data[4,i])*100
    SENS[i,1]  =  (data[4,i])/(data[3,i]+data[4,i])*100
    SPEC[i,1]  = (data[1,i]/(data[1,i]+data[2,i]))*100
    MCC1      = (data[1,i]*data[4,i]) - (data[2,i]*data[3,i])
    MCC2      =  (data[4,i]+data[2,i])*(data[4,i]+data[3,i])
    MCC3      =  (data[1,i]+data[2,i])*(data[1,i]+data[3,i])
    MCC4  =  sqrt(MCC2)*sqrt(MCC3)
    
    
    MCC[i,1]  = MCC1/MCC4
  }
  rm(TP)
  rm(FP)
  rm(TN)
  rm(FN)
  results_ACC <- mean_and_sd(ACC)
  results_SENS <- mean_and_sd(SENS)
  results_SPEC <- mean_and_sd(SPEC)
  results_MCC <- mean_and_sd(MCC)
  rm(ACC)
  rm(SENS)
  rm(SPEC)
  rm(MCC)
  rm(data)
  results_all <- (data.frame(c(results_ACC, results_SENS, results_SPEC, results_MCC)))
  rownames(results_all) <- c("ACC_Mean", "ACC_SD", "Sens_Mean", "Sens_SD", "Spec_Mean", "Spec_SD",
                             "MCC_Mean", "MCC_SD")
  rm(results_ACC)
  rm(results_SENS)
  rm(results_SPEC)
  rm(results_MCC)
  return(results_all)
}

results_testing_Cancer <- function(x) {
  yes <- J48_testing(x, Label = "Cancer")
  great <- data.frame(yes)
  rm(yes)
  TP <- seq(from = 1, to = 400, by = 4)
  FN <- seq(from = 2, to = 400, by = 4)
  FP <- seq(from = 3, to = 400, by = 4)
  TN <- seq(from = 4, to = 400, by = 4)
  results <- mapply(c, great[TP], great[FN], great[FP], great[TN])
  data <- data.frame(results)
  rm(great)
  rm(results)
  m = ncol(data)
  ACC  <- matrix(nrow = m, ncol = 1)
  SENS  <- matrix(nrow = m, ncol = 1)
  SPEC  <-matrix(nrow = m, ncol = 1)
  MCC <- matrix(nrow = m, ncol = 1)
  
  for(i in 1:m){ 
    ACC[i,1]  = (data[1,i]+data[4,i])/(data[1,i]+data[2,i]+data[3,i]+data[4,i])*100
    SENS[i,1]  =  (data[4,i])/(data[3,i]+data[4,i])*100
    SPEC[i,1]  = (data[1,i]/(data[1,i]+data[2,i]))*100
    MCC1      = (data[1,i]*data[4,i]) - (data[2,i]*data[3,i])
    MCC2      =  (data[4,i]+data[2,i])*(data[4,i]+data[3,i])
    MCC3      =  (data[1,i]+data[2,i])*(data[1,i]+data[3,i])
    MCC4  =  sqrt(MCC2)*sqrt(MCC3)
    
    
    MCC[i,1]  = MCC1/MCC4
  }
  rm(TP)
  rm(FP)
  rm(TN)
  rm(FN)
  results_ACC <- mean_and_sd(ACC)
  results_SENS <- mean_and_sd(SENS)
  results_SPEC <- mean_and_sd(SPEC)
  results_MCC <- mean_and_sd(MCC)
  rm(ACC)
  rm(SENS)
  rm(SPEC)
  rm(MCC)
  rm(data)
  results_all <- (data.frame(c(results_ACC, results_SENS, results_SPEC, results_MCC)))
  rownames(results_all) <- c("ACC_Mean", "ACC_SD", "Sens_Mean", "Sens_SD", "Spec_Mean", "Spec_SD",
                             "MCC_Mean", "MCC_SD")
  rm(results_ACC)
  rm(results_SENS)
  rm(results_SPEC)
  rm(results_MCC)
  return(results_all)
}


results_testing_Fungus <- function(x) {
  yes <- J48_testing(x, Label = "Fungus")
  great <- data.frame(yes)
  rm(yes)
  TP <- seq(from = 1, to = 400, by = 4)
  FN <- seq(from = 2, to = 400, by = 4)
  FP <- seq(from = 3, to = 400, by = 4)
  TN <- seq(from = 4, to = 400, by = 4)
  results <- mapply(c, great[TP], great[FN], great[FP], great[TN])
  data <- data.frame(results)
  rm(great)
  rm(results)
  m = ncol(data)
  ACC  <- matrix(nrow = m, ncol = 1)
  SENS  <- matrix(nrow = m, ncol = 1)
  SPEC  <-matrix(nrow = m, ncol = 1)
  MCC <- matrix(nrow = m, ncol = 1)
  
  for(i in 1:m){ 
    ACC[i,1]  = (data[1,i]+data[4,i])/(data[1,i]+data[2,i]+data[3,i]+data[4,i])*100
    SENS[i,1]  =  (data[4,i])/(data[3,i]+data[4,i])*100
    SPEC[i,1]  = (data[1,i]/(data[1,i]+data[2,i]))*100
    MCC1      = (data[1,i]*data[4,i]) - (data[2,i]*data[3,i])
    MCC2      =  (data[4,i]+data[2,i])*(data[4,i]+data[3,i])
    MCC3      =  (data[1,i]+data[2,i])*(data[1,i]+data[3,i])
    MCC4  =  sqrt(MCC2)*sqrt(MCC3)
    
    
    MCC[i,1]  = MCC1/MCC4
  }
  rm(TP)
  rm(FP)
  rm(TN)
  rm(FN)
  results_ACC <- mean_and_sd(ACC)
  results_SENS <- mean_and_sd(SENS)
  results_SPEC <- mean_and_sd(SPEC)
  results_MCC <- mean_and_sd(MCC)
  rm(ACC)
  rm(SENS)
  rm(SPEC)
  rm(MCC)
  rm(data)
  results_all <- (data.frame(c(results_ACC, results_SENS, results_SPEC, results_MCC)))
  rownames(results_all) <- c("ACC_Mean", "ACC_SD", "Sens_Mean", "Sens_SD", "Spec_Mean", "Spec_SD",
                             "MCC_Mean", "MCC_SD")
  rm(results_ACC)
  rm(results_SENS)
  rm(results_SPEC)
  rm(results_MCC)
  return(results_all)
}

results_testing_Virus <- function(x) {
  yes <- J48_testing(x, Label = "Virus")
  great <- data.frame(yes)
  rm(yes)
  TP <- seq(from = 1, to = 400, by = 4)
  FN <- seq(from = 2, to = 400, by = 4)
  FP <- seq(from = 3, to = 400, by = 4)
  TN <- seq(from = 4, to = 400, by = 4)
  results <- mapply(c, great[TP], great[FN], great[FP], great[TN])
  data <- data.frame(results)
  rm(great)
  rm(results)
  m = ncol(data)
  ACC  <- matrix(nrow = m, ncol = 1)
  SENS  <- matrix(nrow = m, ncol = 1)
  SPEC  <-matrix(nrow = m, ncol = 1)
  MCC <- matrix(nrow = m, ncol = 1)
  
  for(i in 1:m){ 
    ACC[i,1]  = (data[1,i]+data[4,i])/(data[1,i]+data[2,i]+data[3,i]+data[4,i])*100
    SENS[i,1]  =  (data[4,i])/(data[3,i]+data[4,i])*100
    SPEC[i,1]  = (data[1,i]/(data[1,i]+data[2,i]))*100
    MCC1      = (data[1,i]*data[4,i]) - (data[2,i]*data[3,i])
    MCC2      =  (data[4,i]+data[2,i])*(data[4,i]+data[3,i])
    MCC3      =  (data[1,i]+data[2,i])*(data[1,i]+data[3,i])
    MCC4  =  sqrt(MCC2)*sqrt(MCC3)
    
    
    MCC[i,1]  = MCC1/MCC4
  }
  rm(TP)
  rm(FP)
  rm(TN)
  rm(FN)
  results_ACC <- mean_and_sd(ACC)
  results_SENS <- mean_and_sd(SENS)
  results_SPEC <- mean_and_sd(SPEC)
  results_MCC <- mean_and_sd(MCC)
  rm(ACC)
  rm(SENS)
  rm(SPEC)
  rm(MCC)
  rm(data)
  results_all <- (data.frame(c(results_ACC, results_SENS, results_SPEC, results_MCC)))
  rownames(results_all) <- c("ACC_Mean", "ACC_SD", "Sens_Mean", "Sens_SD", "Spec_Mean", "Spec_SD",
                             "MCC_Mean", "MCC_SD")
  rm(results_ACC)
  rm(results_SENS)
  rm(results_SPEC)
  rm(results_MCC)
  return(results_all)
}

J48_testing_all <- function(x) {
  bacteria <- results_testing_Bacteria(x)
  cancer <- results_testing_Cancer(x)
  fungus <- results_testing_Fungus(x)
  virus <- results_testing_Virus(x)
  results_all <- cbind(bacteria, cancer, fungus, virus)
  rm(bacteria)
  rm(cancer)
  rm(fungus)
  rm(virus)
  total <- apply(results_all, 1, mean)
  results_all_mean <- cbind(results_all, total)
  rm(results_all)
  rm(total)
  colnames(results_all_mean) <- c("Bacteria", "Cancer", "Fungus", "Virus", "Overall")
  return(results_all_mean)
}
###Results for training 

training_results <- suppressWarnings(J48_training_all(combine_data))
print(training_results)

####Results for CV

CV_results <- suppressWarnings(J48_CV_all(combine_data))
print(CV_results)

####Results for Testing
{r, echo=FALSE, cache = TRUE}
testing_results <- suppressWarnings(J48_testing_all(combine_data))
print(testing_results)
```




