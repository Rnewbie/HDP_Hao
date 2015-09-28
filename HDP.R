library(Biostrings)
library(RCurl)
library(RWeka)
library(caret)

getwd()
#setwd("/Users/sawsimeon/Downloads")
#data <- read.csv("HDP_DATA.csv", header = TRUE)
link <- "https://raw.githubusercontent.com/Rnewbie/HDP_Hao/master/HDP_DATA.csv"
csv <- getURL(link)
data <- read.csv(text = csv, header = TRUE)
### substracting
cancer <- data$cancer
cancer <- as.data.frame(cancer)
fungus <- data$fungus
fungus <- as.data.frame(fungus)
bacteria <- data$bacteria
bacteria <- as.data.frame(bacteria)
virus <- data$virus
virus <- as.data.frame(virus)
mammal <- data$mammal
mammal <- as.data.frame(mammal)
#### remove duplicated 
index_cancer <- which(duplicated(cancer))
cancer <- cancer[-index_cancer, ]
index_fungus <- which(duplicated(fungus))
fungus <- fungus[-index_fungus, ]
index_bacteria <- which(duplicated(bacteria))
bacteria <- bacteria[-index_bacteria, ]
index_virus <- which(duplicated(virus))
virus <- virus[-index_virus, ]
index_mammal <- which(duplicated(mammal))
mammal <- mammal[-index_mammal, ]


setwd("/Volumes/SAW SIMEON 1/HDP_Hao")

### getting Negative Data set
library(RCurl)
library(seqinr)
library(protr)

link <- "http://cs.gmu.edu/~ashehu/sites/default/files/tools/EFC-FCBF_AMPRecognition/data/fasta_files/Xiao_GBMR4_Training_DECOYS.fasta"
x <- getURL(link)
write.fasta(sequence = x, names(x), 
            nbchar = 80,  file.out = "negative_DNA.fasta")

DNA <- readFASTA("negative_DNA.fasta")
DNA <- t(sapply(DNA, s2c))
Protein <- t(sapply(DNA, seqinr::translate))
Protein <- lapply(Protein, str_replace_all,"[[:punct:]]", "*")
write.fasta(sequence=Protein, names(Protein), 
            nbchar = 80, file.out = "negative_Protein.fasta")



### writing fasta files
cancer <- AAStringSet(cancer)
names_cancer <- rep("cancer", length(cancer))
fungus <- AAStringSet(fungus)
bacteria <- AAStringSet(bacteria)
virus <- AAStringSet(virus)
mammal <- AAStringSet(mammal)
###
writeXStringSet(cancer, file = "cancer.fasta", width = 80)
writeXStringSet(fungus, file = "fungus.fasta", width = 80)
writeXStringSet(bacteria, file = "bacteria.fasta", width = 80)
writeXStringSet(virus, file = "virus.fasta", width = 80)
writeXStringSet(mammal, file = "mammal.fasta", width = 80)
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
mammal_des <- as.data.frame(mammal_des)
mammal_des$Label <- "Mammal"
negative_des <- as.data.frame(negative_des)
negative_des$Label <- "Negative"

cancer <- rbind(cancer_des, negative_des)
fungus <- rbind(fungus_des, negative_des)
bacteria <- rbind(fungus_des, negative_des)
virus <- rbind(virus_des, negative_des)
mammal <- rbind(mammal_des, negative_des)

input <- list(cancer = cancer, funus = fungus, bacteira = bacteria, 
              virus = virus, mammal = mammal)



#### training results using J48
J48_training <- function(x) {
  
  results <- list(100)
  for (i in 1:100) {
    in_trian <- caret::createDataPartition(x$Label, p = 0.80, list = FALSE)
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
  return(data.frame(c(results_ACC, results_SENS, results_SPEC, results_MCC)))
}

#### 10-fold results using J48
J48_10fold <- function(x) {
  results <- list(100)
  for (i in 1:100) {
    in_trian <- caret::createDataPartition(x$Label, p = 0.80, list = FALSE)
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
  return(data.frame(c(results_ACC, results_SENS, results_SPEC, results_MCC)))
}

### testing results using J48
J48_testing <- function(x) {
  results <- list(100)
  for (i in 1:100) {
    in_trian <- caret::createDataPartition(x$Label, p = 0.80, list = FALSE)
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
  return(data.frame(c(results_ACC, results_SENS, results_SPEC, results_MCC)))
}

results_J48 <- function(x) {
  c(J48_train(x), J48_cross_validation(x), J48_external(x))
}

library(parallel)
library(doSNOW)
cl <- makeCluster(8)
registerDoSNOW(cl)
clusterExport(cl = cl, ls())

result_J48_performance <- parLapply(cl = cl, input, function(x) {
  models <- results_J48(x)
  return(models)
})

stopCluster(cl)

















### Physiochemical Properties


library(protr)
cancer <- readFASTA("cancer.fasta")[[1]]
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


num_cores <- detectCores()
cl <- makeCluster(num_cores)
cancer_CTD <- t(parSapply(cancer, ctdDes))


stopCluster(cl)




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

### Paran analysis to obtain paren analysis  and PCA
### df shoudld have the descriptors for the
df <- load("composition.RData")
df = df[, -nearZeroVar(df, uniquecut = 1)]
library(paran)
paran(df, iterations = 5000)
pca <- prcomp(df, retx=TRUE,scale.=TRUE)
summary(pca)


### Correlation Descriptors
extractMoran = function (x, props = c('CIDH920105', 'BHAR880101',
                                      'CHAM820101', 'CHAM820102',
                                      'CHOC760101', 'BIGC670101',
                                      'CHAM810101', 'DAYM780201'), 
                         nlag = 30L, customprops = NULL) {
  
  if (protcheck(x) == FALSE) stop('x has unrecognized amino acid type')
  
  if (nchar(x) <= nlag){
    warning("extractMoran: The length of the sequence is <= nlag; NA's will result.")
  }
  
  
  # 1. Compute Pr values for each type of property
  
  AAidx = read.csv(system.file('sysdata/AAidx.csv', package = 'protr'), header = TRUE)
  
  if (!is.null(customprops)) AAidx = rbind(AAidx, customprops)
  
  aaidx = AAidx[, -1]
  row.names(aaidx) = AAidx[, 1]
  n = length(props)
  pmean = rowMeans(aaidx[props, ])
  psd   = apply(aaidx[props, ], 1, sd) * sqrt((20 - 1)/20) # sd() uses (n-1)
  
  Pr = data.frame(matrix(ncol = 20, nrow = n))
  for (i in 1:n) Pr[i, ] = (aaidx[props[i], ] - pmean[i])/psd[i]
  
  # 2. Replace character with numbers, also applies to less than 20 AA occured
  
  xSplitted = strsplit(x, split = '')[[1]]
  AADict = c('A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 
             'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V')
  names(Pr) = AADict
  
  P = vector('list', n)
  for (i in 1:n) P[[i]] = xSplitted
  
  for (i in 1:n) {
    for (j in AADict) {
      try(P[[i]][which(P[[i]] == j)] <- Pr[i, j], silent = TRUE)
    }
  }
  
  P = lapply(P, as.numeric)
  
  # 3. Compute Moran Autocorrelation Descriptor
  
  Moran = vector('list', n)
  N  = length(xSplitted)
  
  Pbar = sapply(P, mean)
  
  for (i in 1:n) {
    for (j in 1:nlag) {
      Moran[[i]][j] = ifelse(N-j > 0,
                             (N/(N - j)) * ((sum((P[[i]][1:(N - j)] - Pbar[i]) * (P[[i]][(1:(N - j)) + j] - Pbar[i])))/(sum((P[[i]] - Pbar[i])^2))),
                             NA)
    }
  }
  
  Moran = unlist(Moran)
  
  names(Moran) = as.vector(t(outer(props, 
                                   paste('.lag', 1:nlag, sep = ''), 
                                   paste, sep = '')))
  
  return(Moran)
  
}


extractMoreauBroto = function (x, props = c('CIDH920105', 'BHAR880101',
                                            'CHAM820101', 'CHAM820102',
                                            'CHOC760101', 'BIGC670101',
                                            'CHAM810101', 'DAYM780201'), 
                               nlag = 30L, customprops = NULL) {
  
  if (protcheck(x) == FALSE) stop('x has unrecognized amino acid type')
  
  if (nchar(x) <= nlag){
    warning("extractMoreauBroto: The length of the sequence is <= nlag; NA's will result.")
  }
  
  # 1. Compute Pr values for each type of property
  
  AAidx = read.csv(system.file('sysdata/AAidx.csv', package = 'protr'), header = TRUE)
  
  if (!is.null(customprops)) AAidx = rbind(AAidx, customprops)
  
  aaidx = AAidx[, -1]
  row.names(aaidx) = AAidx[, 1]
  n = length(props)
  pmean = rowMeans(aaidx[props, ])
  psd   = apply(aaidx[props, ], 1, sd) * sqrt((20 - 1)/20) # sd() uses (n-1)
  
  Pr = data.frame(matrix(ncol = 20, nrow = n))
  for (i in 1:n) Pr[i, ] = (aaidx[props[i], ] - pmean[i])/psd[i]
  
  # 2. Replace character with numbers, also applies to less than 20 AA occured
  
  xSplitted = strsplit(x, split = '')[[1]]
  AADict = c('A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 
             'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V')
  names(Pr) = AADict
  
  P = vector('list', n)
  for (i in 1:n) P[[i]] = xSplitted
  
  for (i in 1:n) {
    for (j in AADict) {
      try(P[[i]][which(P[[i]] == j)] <- Pr[i, j], silent = TRUE)
    }
  }
  
  P = lapply(P, as.numeric)
  
  # 3. Compute Moreau-Broto Autocorrelation Descriptor
  
  MB = vector('list', n)
  N  = length(xSplitted)
  
  for (i in 1:n) {
    for (j in 1:nlag) {
      MB[[i]][j] = ifelse(N-j > 0,
                          sum(P[[i]][1:(N - j)] * P[[i]][(1:(N - j)) + j])/(N - j),
                          NA)
    }
  }
  
  MB = unlist(MB)
  
  names(MB) = as.vector(t(outer(props, 
                                paste('.lag', 1:nlag, sep = ''), 
                                paste, sep = '')))
  
  return(MB)
  
}

extractGeary = function (x, props = c('CIDH920105', 'BHAR880101',
                                      'CHAM820101', 'CHAM820102',
                                      'CHOC760101', 'BIGC670101',
                                      'CHAM810101', 'DAYM780201'), 
                         nlag = 30L, customprops = NULL) {
  
  if (protcheck(x) == FALSE) stop('x has unrecognized amino acid type')
  
  if (nchar(x) <= nlag){
    warning("extractGeary: The length of the sequence is <= nlag; NA's will result.")
  }
  
  # 1. Compute Pr values for each type of property
  
  AAidx = read.csv(system.file('sysdata/AAidx.csv', package = 'protr'), header = TRUE)
  
  if (!is.null(customprops)) AAidx = rbind(AAidx, customprops)
  
  aaidx = AAidx[, -1]
  row.names(aaidx) = AAidx[, 1]
  n = length(props)
  pmean = rowMeans(aaidx[props, ])
  psd   = apply(aaidx[props, ], 1, sd) * sqrt((20 - 1)/20) # sd() uses (n-1)
  
  Pr = data.frame(matrix(ncol = 20, nrow = n))
  for (i in 1:n) Pr[i, ] = (aaidx[props[i], ] - pmean[i])/psd[i]
  
  # 2. Replace character with numbers, also applies to less than 20 AA occured
  
  xSplitted = strsplit(x, split = '')[[1]]
  AADict = c('A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 
             'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V')
  names(Pr) = AADict
  
  P = vector('list', n)
  for (i in 1:n) P[[i]] = xSplitted
  
  for (i in 1:n) {
    for (j in AADict) {
      try(P[[i]][which(P[[i]] == j)] <- Pr[i, j], silent = TRUE)
    }
  }
  
  P = lapply(P, as.numeric)
  
  # 3. Compute Geary Autocorrelation Descriptor
  
  Geary = vector('list', n)
  N  = length(xSplitted)
  
  Pbar = sapply(P, mean)
  
  for (i in 1:n) {
    for (j in 1:nlag) {
      Geary[[i]][j] = ifelse(N-j > 0,
                             ((N - 1)/(2 * (N - j))) * ((sum((P[[i]][1:(N - j)] - P[[i]][(1:(N - j)) + j])^2))/(sum((P[[i]] - Pbar[i])^2))),
                             NA)
    }
  }
  
  Geary = unlist(Geary)
  
  names(Geary) = as.vector(t(outer(props, 
                                   paste('.lag', 1:nlag, sep = ''), 
                                   paste, sep = '')))
  
  return(Geary)
  
}



randomforest_crossvalidation <- function(x) {
  results <- list(100)
  for (i in 1:100) {
    in_trian <- createDataPartition(x$Label, p = 0.80, list = FALSE)
    myData <- x[in_train, ]
    test <- x[-in_train, ]
    k = 10
    index <- sample(1:k, nrow(myData), replace = TRUE)
    model_train <- J48(Label~., data = train)
    eval_j48 <- evaluate_Weka_classifier(model_train, numFolds = 10, complexity = FALSE, seed = 1, class = TRUE)
    confusionmatrix <- eval_j48$confusionMatrix
    results[[i]] <- as.numeric(confusionmatrix)
  }
  return(results)
}


k = 10
index <- sample(1:k, nrow(myData), replace = TRUE)
folds <- 1:k
myRes <- data.frame()
for (j in 1:k)
  training <- subset(myData, index %in% folds[-j])
testing <- subset(myData, index %in% c(j))
ctrl <- caret::trainControl(method = "repeatedcv", number = 10, repeats = 1)
tune <- train(pIC50 ~., data = training, method = "pls", tunLength = 10, trControl = ctrl)
model <- plsr(pIC5~., data = training, ncomp = tune$bestTune[[1]])
prediction <- predict(model, testing, ncomp = tune$bestTune[[1]])




