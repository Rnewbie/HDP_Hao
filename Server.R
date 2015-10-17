library(shiny)
library(seqinr)
library(protr)
library(caret)
library(randomForest)
library(shinyjs)
library(dplyr)
library(conformal)


fit <- get(load("fit.RData"))
#data <- readRDS("data.Rda")
#bacteria <- subset(data, Label == "Bacteria")
#cancer <- subset(data, Label == "Cancer")
#fungus <- subset(data, Label == "Fungus")
#virus <- subset(data, Label == "Virus")
#bacteria <- dplyr::sample_n(bacteria, size = 500, replace = TRUE)
#cancer <- dplyr::sample_n(cancer, size = 500, replace = TRUE)
#fungus <- dplyr::sample_n(fungus, size = 500, replace = TRUE)
#virus <- dplyr::sample_n(virus, size = 500, replace = TRUE)
#wtf <- dplyr::sample_n(data, size = 4)
#data <- rbind(bacteria, cancer, fungus, virus)
example <- get(load("comformal.RData"))
#Label <- data$Label
#train <- data[, 1:20]
#data <- readRDS("data.Rda")
#showClass("ConformalClassification")
#trControl <- trainControl(method = "cv", number = 5, savePredictions = TRUE)
#model <- train(train, Label, data = data, method = "rf", trControl = trControl, predict.all = TRUE)
#save(model, file = "model.RData")
#save(example, file = "example.RData")
#save(model, file = "fit.RData")
x <- readFASTA("text.fasta")
x <- x[(sapply(x, protcheck))]
ACC <- t(sapply(x, extractAAC))
test <- ACC
#Prediction <- predict(fit, wtf)
#Prediction <- as.data.frame(Prediction)
#Protein <- cbind(Name = rownames(test, test))
#results <- cbind(Protein, Prediction)
#results <- data.frame(results, row.names=NULL)
#example$CalculateCVScores(model = fit)
example$CalculatePValues(new.data = test)
p_values <- example$p.values$P.values
comformal_predictor <- ConformalClassification$new()
comformal_predictor$CalculateCVScores(model = fit)
comformal_predictor$CalculatePValues(new.data = test)
comformal_predictor$p.values
save(comformal_predictor, file = "comformal.RData")
#p_values <- data.frame(p_values)
#p_values <- head(p_values, n = nrow(results))
#results_all <- cbind(results, p_values)
#print(results_all)

#example <- ConformalClassification$new()
#suppressMessages(example$CalculateCVScores(model = model))

#save(example, file = "example.RData")

shinyServer(function(input, output, session) {
  
  
  observe({
    FASTADATA <- ''
    fastaexample <- '>HDP-Bacteria
KVLKAAAKAALNAVLVGANA
'
    if(input$addlink>0) {
      isolate({
        FASTADATA <- fastaexample
        updateTextInput(session, inputId = "Sequence", value = FASTADATA)
      })
    }
  })
  
  observe({
    emptyDATA <- ""
    if(input$clearbutton>0) {
      isolate({
        updateTextInput(session, inputId = "Sequence", value = emptyDATA) 
        is.null(datasetInput())
        
      })
    }
  })
  

  datasetInput <- reactive({
    
    inFile <- input$file1 
    inTextbox <- input$Sequence
    
    
    if (is.null(inTextbox)) {
      return("Please insert/upload sequence in FASTA format")
    } else {
      if (is.null(inFile)) {
        x <- inTextbox
        write.fasta(sequence = x, names = names(x),
                    nbchar = 80, , file.out = "text.fasta")
        x <- readFASTA("text.fasta")
        x <- x[(sapply(x, protcheck))]
        ACC <- t(sapply(x, extractAAC))
        test <- data.frame(ACC)
        Prediction <- predict(fit, test)
        Prediction <- as.data.frame(Prediction)
        Protein <- cbind(Name = rownames(test, test))
        results <- cbind(Protein, Prediction)
        results <- data.frame(results, row.names=NULL)
        suppressMessages(example$CalculatePValues(new.data = test))
        p_values <- example$p.values$P.values
        p_values <- data.frame(p_values)
        p_values <- head(p_values, n = nrow(results))
        results_all <- cbind(results, p_values)
        print(results_all)
      } 
      else {     
        x <- readFASTA(inFile$datapath)
        x <- x[(sapply(x, protcheck))]
        ACC <- t(sapply(x, extractAAC))
        test <- data.frame(ACC)
        Prediction <- predict(fit, test)
        Prediction <- as.data.frame(Prediction)
        Protein <- cbind(Name = rownames(test, test))
        results <- cbind(Protein, Prediction)
        results <- data.frame(results, row.names=NULL)
        suppressMessages(example$CalculatePValues(new.data = test))
        p_values <- example$p.values$P.values
        p_values <- data.frame(p_values)
        p_values <- head(p_values, n = nrow(results))
        results_all <- cbind(results, p_values)
        print(results_all)
        
      }
    }
    
    
  })
  
  
  
  output$contents <- renderPrint({
    if (input$submitbutton>0) { 
    isolate(datasetInput()) 
    } else {
      if (input$clearbutton>0) {
      isolate(is.null(datasetInput()))
    } else {
      return("Please insert/upload sequence in FASTA format")
    }
    }
  })
  
  output$downloadData <- downloadHandler(
    filename = function() { paste('Predicted_Results', '.csv', sep='') },
    content = function(file) {
      write.csv(datasetInput(), file, row.names=FALSE)
    })
  
  

  })

