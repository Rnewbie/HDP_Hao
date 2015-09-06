library(Biostrings)
library(RCurl)

getwd()
#setwd("/Users/sawsimeon/Downloads")
#data <- read.csv("HDP_DATA.csv", header = TRUE)
link <- "https://raw.githubusercontent.com/Rnewbie/HDP_Hao/master/HDP_DATA.csv"
csv <- getURL(link)
data <- read.csv(text = csv, header = TRUE)
df <- apply(data, 2, unique)
ok <- do.call(rbind.data.frame, df)
data <- t(ok)
data <- data.frame(data)
### substracting
cancer <- data$cancer
fungus <- data$fungus
gram_positive <- data$gram_positive
gram_negative <- data$gram_negative
virus <- data$virus
mammal <- data$mammal
#setwd("/Volumes/SAW SIMEON/HDP")


### writing fasta files
cancer <- AAStringSet(cancer)
names_cancer <- rep("cancer", length(cancer))
fungus <- AAStringSet(fungus)
gram_positive <- AAStringSet(gram_positive)
gram_negative <- AAStringSet(gram_negative)
virus <- AAStringSet(virus)
mammal <- AAStringSet(mammal)
###
writeXStringSet(cancer, file = "cancer.fasta", width = 80)
writeXStringSet(fungus, file = "fungus.fasta", width = 80)
writeXStringSet(gram_positive, file = "gram_positive.fasta", width = 80)
writeXStringSet(gram_negative, file = "gram_negative.fasta", width = 80)
writeXStringSet(virus, file = "virus.fasta", width = 80)
writeXStringSet(mammal, file = "mammal.fasta", width = 80)
### read FASta
library(protr)
cancer <- readFASTA("cancer.fasta")
fungus <- readFASTA("fungus.fasta")
gram_positive <- readFASTA("gram_positive.fasta")
gram_negative <- readFASTA("gram_negative.fasta")
virus <- readFASTA("virus.fasta")
mammal <- readFASTA("mammal.fasta")
### removed wired protein
cancer <- cancer[(sapply(cancer, protcheck))]
fungus <- fungus[(sapply(fungus, protcheck))]
gram_positive <- gram_positive[(sapply(gram_positive, protcheck))]
gram_negative <- gram_negative[(sapply(gram_negative, protcheck))]
virus <- virus[(sapply(virus, protcheck))]
mammal <- mammal[(sapply(mammal, protcheck))]
### function for amino acid composition, dipeptide composition, tripeptide composition
composition <- function(x) {
  library(protr)
  c(extractAAC(x), extractDC(x), extractTC(x))
}
### generation the Descriptors
cancer_des <- t(sapply(cancer, composition))
fungus_des <- t(sapply(fungus, composition))
gram_positive_des <- t(sapply(gram_positive, composition))
gram_negative_des <- t(sapply(gram_negative, composition))
virus_des <- t(sapply(virus, composition))
mammal_des <- t(sapply(mammal, composition))
#### label the labels 
cancer_des$Label <- "cancer"
fungus_des$Label <- "fugus"
gram_positive_des$Label <- "gram_positive"
gram_negative$Label <- "gram_negative"
virus_des$Label <- "virus"
big_data <- rbind(cancer_des, fungus_des, gram_positive_des, gram_negative_des, virus_des, mammal_des)
#### training results using J48
J48_training <- function(x) {
  
  results <- list(100)
  for (i in 1:100) {
    in_trian <- createDataPartition(x$Label, p = 0.80, list = FALSE)
    train <- x[in_train, ]
    test <- x[-in_train, ]
    model_train <- J48(Label~., data = train)
    summary <- summary(model_train)
    confusionmatrix <- summary$confusionMatrix
    results[[i]] <- as.numeric(confusionmatrix)
  }
  return(results)
}
#### 10-fold results using J48
J48_10fold <- function(x) {
  results <- list(100)
  for (i in 1:100) {
    in_trian <- createDataPartition(x$Label, p = 0.80, list = FALSE)
    train <- x[in_train, ]
    test <- x[-in_train, ]
    model_train <- J48(Label~., data = train)
    eval_j48 <- evaluate_Weka_classifier(model_train, numFolds = 10, complexity = FALSE, seed = 1, class = TRUE)
    confusionmatrix <- eval_j48$confusionMatrix
    results[[i]] <- as.numeric(confusionmatrix)
  }
  return(results)
}

### testing results using J48
J48_testing <- function(x) {
  results <- list(100)
  for (i in 1:100) {
    in_trian <- createDataPartition(x$Label, p = 0.80, list = FALSE)
    train <- x[in_train, ]
    test <- x[-in_train, ]
    model_train <- J48(Label~., data = train)
    eval_external <- evaluate_Weka_classifier(model_train, newdata = test, numFolds = 0, complexity = FALSE, seed = 1, class = TRUE)
    confusionmatrix <- eval_external$confusionMatrix
    results[[i]] <- as.numeric(confusionmatrix)
  }
  return(results)
}














