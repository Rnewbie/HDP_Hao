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

### Physiochemical Properties


library(protr)
cancer <- readFASTA("cancer.fasta")
fungus <- readFASTA("fungus.fasta")
gram_positive <- readFASTA("gram_positive.fasta")
gram_negative <- readFASTA("gram_negative.fasta")
virus <- readFASTA("virus.fasta")
mammal <- readFASTA("mammal.fasta")

cancer <- cancer[(sapply(cancer, protcheck))]
fungus <- fungus[(sapply(fungus, protcheck))]
gram_positive <- gram_positive[(sapply(gram_positive, protcheck))]
gram_negative <- gram_negative[(sapply(gram_negative, protcheck))]
virus <- virus[(sapply(virus, protcheck))]
mammal <- mammal[(sapply(mammal, protcheck))]
### function for PCP from AAindex
pcp_des <- function(x) {
  library(protr)
  PCP <- t(sapply(x, extractMoreauBroto, props = AAindex$AccNo, nlag = 1L))
  PCP <- PCP[, colSums(is.na(PCP)) != nrow(PCP)]
  return(PCP)
}
### generation the Descriptors
cancer_PCP <- pcp_des(cancer)
fungus_PCP <- pcp_des(fungus)
gram_positive_PCP <- pcp_des(gram_positive)
gram_negative_PCP <- pcp_des(gram_negative)
virus_PCP <- pcp_des(virus)
mammal_PCP <- pcp_des(mammal)

#### label the labels 
cancer_PCP <- as.data.frame(cancer_PCP)
Label <- c("cancer")
cancer_PCP <- cbind(Label, cancer_PCP)

fungus_PCP <- as.data.frame(fungus_PCP)
Label <- c("fungus")
fungus_PCP <- cbind(Label, fungus_PCP)

gram_positive_PCP <- as.data.frame(gram_positive_PCP)
Label <- c("gram_positive")
gram_positive_PCP <- cbind(Label, gram_positive_PCP)

gram_negative_PCP <- as.data.frame(gram_negative_PCP)
Label <- c("gram_negative")
gram_negative_PCP <- cbind(Label, gram_negative_PCP)

virus_PCP <- as.data.frame(virus_PCP)
Label <- c("virus")
virus_PCP <- cbind(Label, virus_PCP)

mammal_PCP <- as.data.frame(mammal_PCP)
Label <- c("mammal")
mammal_PCP <- cbind(Label, mammal_PCP)

PCP_data <- rbind(cancer_PCP, fungus_PCP, gram_positive_PCP, gram_negative_PCP, virus_PCP, mammal_PCP)
#### REAdy

ctdDes <- function(x) {
  c(extractCTDC(x), ### composition
    extractCTDT(x), ### transition
    extractCTDD(x)) ### distribution
}

cancer_CTD <- t(sapply(cancer, ctdDes))
fungus_CTD <- t(sapply(fungus, ctdDes))
gram_positive_CTD <- t(sapply(gram_positive, ctdDes))
gram_negative_CTD <- t(sapply(gram_negative, ctdDes))
virus_CTD <- t(sapply(virus, ctdDes))
mammal_CTD <- t(sapply(mammal, ctdDes))


cancer_CTD <- as.data.frame(cancer_CTD)
Label <- c("cancer")
cancer_CTD <- cbind(Label, cancer_CTD)

fungus_CTD <- as.data.frame(fungus_CTD)
Label <- c("fungus")
fungus_CTD <- cbind(Label, fungus_CTD)

gram_positive_CTD <- as.data.frame(gram_positive_CTD)
Label <- c("gram_positive")
gram_positive_CTD <- cbind(Label, gram_positive_CTD)

gram_negative_CTD <- as.data.frame(gram_negative_CTD)
Label <- c("gram_negative")
gram_negative_CTD <- cbind(Label, gram_negative_CTD)

virus_CTD <- as.data.frame(virus_CTD)
Label <- c("virus")
virus_CTD <- cbind(Label, virus_CTD)

mammal_CTD <- as.data.frame(mammal_CTD)
Label <- c("mammal")
mammal_CTD <- cbind(Label, mammal_CTD)

CTD_data <- rbind(cancer_CTD, fungus_CTD, gram_positive_CTD, gram_negative_CTD, virus_CTD, mammal_CTD)

input <- list(composition = big_data,
              AAindex = PCP_data,
              CTD = CTD_data)
results_train <- lapply(input, function(x) {
  model <- J48_training(x)
  return(model)
})

















extractCTDD = function (x) {
  
  if (protcheck(x) == FALSE) stop('x has unrecognized amino acid type')
  
  group1 = list(hydrophobicity  = c('R', 'K', 'E', 'D', 'Q', 'N'),
                normwaalsvolume = c('G', 'A', 'S', 'T', 'P', 'D', 'C'),
                polarity        = c('L', 'I', 'F', 'W', 'C', 'M', 'V', 'Y'),
                polarizability  = c('G', 'A', 'S', 'D', 'T'),
                charge          = c('K', 'R'),
                secondarystruct = c('E', 'A', 'L', 'M', 'Q', 'K', 'R', 'H'),
                solventaccess   = c('A', 'L', 'F', 'C', 'G', 'I', 'V', 'W'))
  
  group2 = list(hydrophobicity  = c('G', 'A', 'S', 'T', 'P', 'H', 'Y'),
                normwaalsvolume = c('N', 'V', 'E', 'Q', 'I', 'L'),
                polarity        = c('P', 'A', 'T', 'G', 'S'),
                polarizability  = c('C', 'P', 'N', 'V', 'E', 'Q', 'I', 'L'),
                charge          = c('A', 'N', 'C', 'Q', 'G', 'H', 'I', 'L', 
                                    'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'),
                secondarystruct = c('V', 'I', 'Y', 'C', 'W', 'F', 'T'),
                solventaccess   = c('R', 'K', 'Q', 'E', 'N', 'D'))
  
  group3 = list(hydrophobicity  = c('C', 'L', 'V', 'I', 'M', 'F', 'W'),
                normwaalsvolume = c('M', 'H', 'K', 'F', 'R', 'Y', 'W'),
                polarity        = c('H', 'Q', 'R', 'K', 'N', 'E', 'D'),
                polarizability  = c('K', 'M', 'H', 'F', 'R', 'Y', 'W'),
                charge          = c('D', 'E'),
                secondarystruct = c('G', 'N', 'P', 'S', 'D'),
                solventaccess   = c('M', 'S', 'P', 'T', 'H', 'Y'))
  
  xSplitted = strsplit(x, split = '')[[1]]
  n  = nchar(x)
  
  G = vector('list', 7)
  for (i in 1:7) G[[i]] = rep(NA, n)
  
  # Get groups for each property & each amino acid
  
  for (i in 1:7) {
    try(G[[i]][which(xSplitted %in% group1[[i]])] <- 'G1')
    try(G[[i]][which(xSplitted %in% group2[[i]])] <- 'G2')
    try(G[[i]][which(xSplitted %in% group3[[i]])] <- 'G3')
  }
  
  # Compute Distribution
  
  D = vector('list', 7)
  for (i in 1:7) D[[i]] = matrix(ncol = 5, nrow = 3)
  
  for (i in 1:7) {
    inds = which(G[[i]] == 'G1')
    quartiles <- floor(length(inds) * c(0.25, 0.5, 0.75))
    D[[i]][1, ] = ifelse(length(inds)>0,
                         (inds[c(1, ifelse(quartiles>0, quartiles, 1), length(inds))])*100/n,
                         0)
    inds = which(G[[i]] == 'G2')
    quartiles <- floor(length(inds) * c(0.25, 0.5, 0.75))
    D[[i]][2, ] = ifelse(length(inds)>0,
                         (inds[c(1, ifelse(quartiles>0, quartiles, 1), length(inds))])*100/n,
                         0)
    inds = which(G[[i]] == 'G3')
    quartiles <- floor(length(inds) * c(0.25, 0.5, 0.75))
    D[[i]][3, ] = ifelse(length(inds)>0,
                         (inds[c(1, ifelse(quartiles>0, quartiles, 1), length(inds))])*100/n,
                         0)
  }
  
  D = do.call(rbind, D)
  D = as.vector(t(D))
  
  names(D) = paste(rep(paste('prop', 1:7, sep = ''), each = 15),
                   rep(rep(c('.G1', '.G2', '.G3'), each = 5), times = 7),
                   rep(paste('.residue', c('0', '25', '50', '75', '100'), 
                             sep = ''), times = 21), sep = '')
  
  return(D)
  
}

