### getting Negative Data set
library(RCurl)
library(seqinr)
library(protr)

link <- "http://cs.gmu.edu/~ashehu/sites/default/files/tools/EFC-FCBF_AMPRecognition/data/fasta_files/Xiao_GBMR4_Training_DECOYS.fasta"
x <- getURL(link)
write.fasta(sequence = x, names(x), 
            nbchar = 80, , file.out = "negative_DNA.fasta")

DNA <- readFASTA("negative_DNA.fasta")
DNA <- t(sapply(DNA, s2c))
Protein <- t(sapply(DNA, seqinr::translate))
Protein <- lapply(Protein, str_replace_all,"[[:punct:]]", "*")
write.fasta(sequence=Protein, names(Protein), 
                                    nbchar = 80, file.out = "negative_Protein.fasta")

negative <- protr::readFASTA("negative_Protein.fasta")
negative <- negative[(sapply(negative, protcheck))]

training <- c(99.38, 0.07, 99.44, 0.08, 99.08, 0.21, 0.98, 0.00)
cross_validata <- c(98.12, 0.13, 98.55, 0.10)




