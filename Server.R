library(shiny)
library(seqinr)
library(protr)
library(caret)
library(randomForest)
library(shinyjs)
library(dplyr)
library(conformal)
library(RWeka)

combine_data <- readRDS("data.Rds")

fit <- J48(Label~., data = combine_data)

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
        print(results)
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
        print(results)
        
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

