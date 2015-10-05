### read FASTA
library(protr)
cancer <- readFASTA("cancer.fasta")
fungus <- readFASTA("fungus.fasta")
bacteria <- readFASTA("bacteria.fasta")
virus <- readFASTA("virus.fasta")
mammal <- readFASTA("mammal.fasta")
negative <- protr::readFASTA("negative_Protein.fasta")

### removed wired protein
cancer <- cancer[(sapply(cancer, protcheck))]
fungus <- fungus[(sapply(fungus, protcheck))]
bacteria <- bacteria[(sapply(bacteria, protcheck))]
virus <- virus[(sapply(virus, protcheck))]
mammal <- mammal[(sapply(mammal, protcheck))]
negative <- negative[(sapply(negative, protcheck))]

### function for amino acid composition, dipeptide composition, tripeptide composition
composition <- function(x) {
  library(protr)
  c(extractAAC(x))
}
### generation the Descriptors
cancer_des <- t(sapply(cancer, extractAAC))
fungus_des <- t(sapply(fungus, extractAAC))
bacteria_des <- t(sapply(bacteria, extractAAC))
virus_des <- t(sapply(virus, extractAAC))
mammal_des <- t(sapply(mammal, extractAAC))
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
#mammal_des <- as.data.frame(mammal_des)
#mammal_des$Label <- "Mammal"
#negative_des <- as.data.frame(negative_des)
#negative_des$Label <- "Negative"
combine_data <- rbind(cancer_des, fungus_des,
                      bacteria_des, virus_des)
combine_data$Label <- as.factor(combine_data$Label)



#### training results using J48
J48_training <- function(x, Label){
  if (Label == "Bacteria") {
    library(parallel)
    library(doSNOW)
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
      model_train <- RWeka::J48(Label~., data = train)
      actual <- train$Label
      prediction <- predict(model_train, train)
      results <- caret::confusionMatrix(prediction, actual)
      results <- results$table
      results <- table(prediction, actual)
      results <- as.numeric(results)
      rm(model_train)
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
      model_train <- RWeka::J48(Label~., data = train)
      actual <- train$Label
      prediction <- predict(model_train, train)
      results <- caret::confusionMatrix(prediction, actual)
      results <- results$table
      results <- table(prediction, actual)
      results <- as.numeric(results)
      rm(model_train)
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
  TP <- seq(from = 1, to = 400, by = 4)
  FN <- seq(from = 2, to = 400, by = 4)
  FP <- seq(from = 3, to = 400, by = 4)
  TN <- seq(from = 4, to = 400, by = 4)
  results <- mapply(c, great[TP], great[FN], great[FP], great[TN])
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
  rm(yes)
  rm(great)
  rm(data)
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
  return(results_all)
}

results_training_Cancer <- function(x) {
  yes <- J48_training(x, Label = "Cancer")
  great <- data.frame(yes)
  TP <- seq(from = 1, to = 400, by = 4)
  FN <- seq(from = 2, to = 400, by = 4)
  FP <- seq(from = 3, to = 400, by = 4)
  TN <- seq(from = 4, to = 400, by = 4)
  results <- mapply(c, great[TP], great[FN], great[FP], great[TN])
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
  rm(yes)
  rm(great)
  rm(results)
  results_ACC <- mean_and_sd(ACC)
  results_SENS <- mean_and_sd(SENS)
  results_SPEC <- mean_and_sd(SPEC)
  results_MCC <- mean_and_sd(MCC)
  results_all <- (data.frame(c(results_ACC, results_SENS, results_SPEC, results_MCC)))
  rownames(results_all) <- c("ACC_Mean", "ACC_SD", "Sens_Mean", "Sens_SD", "Spec_Mean", "Spec_SD",
                             "MCC_Mean", "MCC_SD")
  return(results_all)
}


results_training_Fungus <- function(x) {
  yes <- J48_training(x, Label = "Fungus")
  great <- data.frame(yes)
  TP <- seq(from = 1, to = 400, by = 4)
  FN <- seq(from = 2, to = 400, by = 4)
  FP <- seq(from = 3, to = 400, by = 4)
  TN <- seq(from = 4, to = 400, by = 4)
  results <- mapply(c, great[TP], great[FN], great[FP], great[TN])
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
  rm(yes)
  rm(great)
  rm(data)
  results_ACC <- mean_and_sd(ACC)
  results_SENS <- mean_and_sd(SENS)
  results_SPEC <- mean_and_sd(SPEC)
  results_MCC <- mean_and_sd(MCC)
  results_all <- (data.frame(c(results_ACC, results_SENS, results_SPEC, results_MCC)))
  rownames(results_all) <- c("ACC_Mean", "ACC_SD", "Sens_Mean", "Sens_SD", "Spec_Mean", "Spec_SD",
                             "MCC_Mean", "MCC_SD")
  return(results_all)
}

results_training_Virus <- function(x) {
  yes <- J48_training(x, Label = "Virus")
  great <- data.frame(yes)
  TP <- seq(from = 1, to = 400, by = 4)
  FN <- seq(from = 2, to = 400, by = 4)
  FP <- seq(from = 3, to = 400, by = 4)
  TN <- seq(from = 4, to = 400, by = 4)
  results <- mapply(c, great[TP], great[FN], great[FP], great[TN])
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
  rm(yes)
  results_ACC <- mean_and_sd(ACC)
  results_SENS <- mean_and_sd(SENS)
  results_SPEC <- mean_and_sd(SPEC)
  results_MCC <- mean_and_sd(MCC)
  results_all <- (data.frame(c(results_ACC, results_SENS, results_SPEC, results_MCC)))
  rownames(results_all) <- c("ACC_Mean", "ACC_SD", "Sens_Mean", "Sens_SD", "Spec_Mean", "Spec_SD",
                             "MCC_Mean", "MCC_SD")
  return(results_all)
}

J48_training_all <- function(x) {
  bacteria <- results_training_Bacteria(x)
  cancer <- results_training_Cancer(x)
  fungus <- results_training_Fungus(x)
  virus <- results_training_Virus(x)
  results_all <- cbind(bacteria, cancer, fungus, virus)
  total <- apply(results_all, 1, mean)
  results_all_mean <- cbind(results_all, total)
  colnames(results_all_mean) <- c("Bacteria", "Cancer", "Fungus", "Virus", "Overall")
  return(results_all_mean)
}


### 10 fold Cross Validations

J48_10_CV <- function(x, Label){
  if (Label == "Bacteria") {
    library(parallel)
    library(doSNOW)
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
      confusionMatrix <- results$confusionMatrix
      rm(model_train)
      rm(results)
      results <- as.numeric(confusionMatrix)
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
      confusionMatrix <- results$confusionMatrix
      rm(model_train)
      rm(results)
      results <- as.numeric(confusionMatrix)
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
      model_train <- RWeka::J48(Label~., data = train)
      results <- RWeka::evaluate_Weka_classifier(model_train, newata = NULL, numFolds = 10, complexity = FALSE)
      confusionMatrix <- results$confusionMatrix
      rm(model_train)
      rm(results)
      results <- as.numeric(confusionMatrix)
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
      model_train <- RWeka::J48(Label~., data = train)
      results <- RWeka::evaluate_Weka_classifier(model_train, newata = NULL, numFolds = 10, complexity = FALSE)
      confusionMatrix <- results$confusionMatrix
      rm(model_train)
      rm(results)
      results <- as.numeric(confusionMatrix)
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
  TP <- seq(from = 1, to = 400, by = 4)
  FN <- seq(from = 2, to = 400, by = 4)
  FP <- seq(from = 3, to = 400, by = 4)
  TN <- seq(from = 4, to = 400, by = 4)
  results <- mapply(c, great[TP], great[FN], great[FP], great[TN])
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
  rm(yes)
  rm(great)
  rm(data)
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
  return(results_all)
}

results_CV_Cancer <- function(x) {
  yes <- J48_10_CV(x, Label = "Cancer")
  great <- data.frame(yes)
  TP <- seq(from = 1, to = 400, by = 4)
  FN <- seq(from = 2, to = 400, by = 4)
  FP <- seq(from = 3, to = 400, by = 4)
  TN <- seq(from = 4, to = 400, by = 4)
  results <- mapply(c, great[TP], great[FN], great[FP], great[TN])
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
  rm(yes)
  rm(great)
  rm(results)
  results_ACC <- mean_and_sd(ACC)
  results_SENS <- mean_and_sd(SENS)
  results_SPEC <- mean_and_sd(SPEC)
  results_MCC <- mean_and_sd(MCC)
  results_all <- (data.frame(c(results_ACC, results_SENS, results_SPEC, results_MCC)))
  rownames(results_all) <- c("ACC_Mean", "ACC_SD", "Sens_Mean", "Sens_SD", "Spec_Mean", "Spec_SD",
                             "MCC_Mean", "MCC_SD")
  return(results_all)
}


results_CV_Fungus <- function(x) {
  yes <- J48_10_CV(x, Label = "Fungus")
  great <- data.frame(yes)
  TP <- seq(from = 1, to = 400, by = 4)
  FN <- seq(from = 2, to = 400, by = 4)
  FP <- seq(from = 3, to = 400, by = 4)
  TN <- seq(from = 4, to = 400, by = 4)
  results <- mapply(c, great[TP], great[FN], great[FP], great[TN])
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
  rm(yes)
  rm(great)
  rm(data)
  results_ACC <- mean_and_sd(ACC)
  results_SENS <- mean_and_sd(SENS)
  results_SPEC <- mean_and_sd(SPEC)
  results_MCC <- mean_and_sd(MCC)
  results_all <- (data.frame(c(results_ACC, results_SENS, results_SPEC, results_MCC)))
  rownames(results_all) <- c("ACC_Mean", "ACC_SD", "Sens_Mean", "Sens_SD", "Spec_Mean", "Spec_SD",
                             "MCC_Mean", "MCC_SD")
  return(results_all)
}

results_CV_Virus <- function(x) {
  yes <- J48_10_CV(x, Label = "Virus")
  great <- data.frame(yes)
  TP <- seq(from = 1, to = 400, by = 4)
  FN <- seq(from = 2, to = 400, by = 4)
  FP <- seq(from = 3, to = 400, by = 4)
  TN <- seq(from = 4, to = 400, by = 4)
  results <- mapply(c, great[TP], great[FN], great[FP], great[TN])
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
  rm(yes)
  results_ACC <- mean_and_sd(ACC)
  results_SENS <- mean_and_sd(SENS)
  results_SPEC <- mean_and_sd(SPEC)
  results_MCC <- mean_and_sd(MCC)
  results_all <- (data.frame(c(results_ACC, results_SENS, results_SPEC, results_MCC)))
  rownames(results_all) <- c("ACC_Mean", "ACC_SD", "Sens_Mean", "Sens_SD", "Spec_Mean", "Spec_SD",
                             "MCC_Mean", "MCC_SD")
  return(results_all)
}

J48_CV_all <- function(x) {
  bacteria <- results_CV_Bacteria(x)
  cancer <- results_CV_Cancer(x)
  fungus <- results_CV_Fungus(x)
  virus <- results_CV_Virus(x)
  results_all <- cbind(bacteria, cancer, fungus, virus)
  total <- apply(results_all, 1, mean)
  results_all_mean <- cbind(results_all, total)
  colnames(results_all_mean) <- c("Bacteria", "Cancer", "Fungus", "Virus", "Overall")
  return(results_all_mean)
}

#### function for testing

J48_testing <- function(x, Label){
  if (Label == "Bacteria") {
    library(parallel)
    library(doSNOW)
    cl <- makeCluster(8)
    registerDoSNOW(cl)
    
    ok <- list(100)
    ok <- foreach(i = 1:100) %dopar% { 
      in_train <- caret::createDataPartition(x$Label, p = 0.80, list = FALSE)
      train <- x[in_train, ]
      test <- x[-in_train, ]
      rm(in_train)
      model_train <- RWeka::J48(Label~., data = train)
      actual <- test$Label
      prediction <- predict(model_train, test)
      rm(train)
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
      actual <- test$Label
      prediction <- predict(model_train, test)
      rm(train)
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
      model_train <- RWeka::J48(Label~., data = train)
      actual <- test$Label
      prediction <- predict(model_train, test)
      results <- caret::confusionMatrix(prediction, actual)
      results <- results$table
      results <- table(prediction, actual)
      results <- as.numeric(results)
      rm(model_train)
      rm(train)
      rm(test)
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
      model_train <- RWeka::J48(Label~., data = train)
      actual <- train$Label
      prediction <- predict(model_train, train)
      results <- caret::confusionMatrix(prediction, actual)
      results <- results$table
      results <- table(prediction, actual)
      results <- as.numeric(results)
      rm(model_train)
      rm(in_train)
      rm(train)
      rm(test)
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
  TP <- seq(from = 1, to = 400, by = 4)
  FN <- seq(from = 2, to = 400, by = 4)
  FP <- seq(from = 3, to = 400, by = 4)
  TN <- seq(from = 4, to = 400, by = 4)
  results <- mapply(c, great[TP], great[FN], great[FP], great[TN])
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
  rm(yes)
  rm(great)
  rm(data)
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
  return(results_all)
}

results_testing_Cancer <- function(x) {
  yes <- J48_testing(x, Label = "Cancer")
  great <- data.frame(yes)
  TP <- seq(from = 1, to = 400, by = 4)
  FN <- seq(from = 2, to = 400, by = 4)
  FP <- seq(from = 3, to = 400, by = 4)
  TN <- seq(from = 4, to = 400, by = 4)
  results <- mapply(c, great[TP], great[FN], great[FP], great[TN])
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
  rm(yes)
  rm(great)
  rm(results)
  results_ACC <- mean_and_sd(ACC)
  results_SENS <- mean_and_sd(SENS)
  results_SPEC <- mean_and_sd(SPEC)
  results_MCC <- mean_and_sd(MCC)
  results_all <- (data.frame(c(results_ACC, results_SENS, results_SPEC, results_MCC)))
  rownames(results_all) <- c("ACC_Mean", "ACC_SD", "Sens_Mean", "Sens_SD", "Spec_Mean", "Spec_SD",
                             "MCC_Mean", "MCC_SD")
  return(results_all)
}


results_testing_Fungus <- function(x) {
  yes <- J48_testing(x, Label = "Fungus")
  great <- data.frame(yes)
  TP <- seq(from = 1, to = 400, by = 4)
  FN <- seq(from = 2, to = 400, by = 4)
  FP <- seq(from = 3, to = 400, by = 4)
  TN <- seq(from = 4, to = 400, by = 4)
  results <- mapply(c, great[TP], great[FN], great[FP], great[TN])
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
  rm(yes)
  rm(great)
  rm(data)
  results_ACC <- mean_and_sd(ACC)
  results_SENS <- mean_and_sd(SENS)
  results_SPEC <- mean_and_sd(SPEC)
  results_MCC <- mean_and_sd(MCC)
  results_all <- (data.frame(c(results_ACC, results_SENS, results_SPEC, results_MCC)))
  rownames(results_all) <- c("ACC_Mean", "ACC_SD", "Sens_Mean", "Sens_SD", "Spec_Mean", "Spec_SD",
                             "MCC_Mean", "MCC_SD")
  return(results_all)
}

results_testing_Virus <- function(x) {
  yes <- J48_testing(x, Label = "Virus")
  great <- data.frame(yes)
  TP <- seq(from = 1, to = 400, by = 4)
  FN <- seq(from = 2, to = 400, by = 4)
  FP <- seq(from = 3, to = 400, by = 4)
  TN <- seq(from = 4, to = 400, by = 4)
  results <- mapply(c, great[TP], great[FN], great[FP], great[TN])
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
  rm(yes)
  results_ACC <- mean_and_sd(ACC)
  results_SENS <- mean_and_sd(SENS)
  results_SPEC <- mean_and_sd(SPEC)
  results_MCC <- mean_and_sd(MCC)
  results_all <- (data.frame(c(results_ACC, results_SENS, results_SPEC, results_MCC)))
  rownames(results_all) <- c("ACC_Mean", "ACC_SD", "Sens_Mean", "Sens_SD", "Spec_Mean", "Spec_SD",
                             "MCC_Mean", "MCC_SD")
  return(results_all)
}

J48_testing_all <- function(x) {
  bacteria <- results_testing_Bacteria(x)
  cancer <- results_testing_Cancer(x)
  fungus <- results_testing_Fungus(x)
  virus <- results_testing_Virus(x)
  results_all <- cbind(bacteria, cancer, fungus, virus)
  total <- apply(results_all, 1, mean)
  results_all_mean <- cbind(results_all, total)
  colnames(results_all_mean) <- c("Bacteria", "Cancer", "Fungus", "Virus", "Overall")
  return(results_all_mean)
}

#### Results for training, 10-fold and testing
#J48_testing_all(combine_data)
J48_training_all(combine_data)
#J48_CV_all(combine_data)

