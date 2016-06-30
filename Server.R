library(shiny)
library(seqinr)
library(protr)
library(caret)
library(randomForest)
library(shinyjs)
library(dplyr)
library(conformal)

cancer <- readFASTA("cancer.fasta")
fungus <- readFASTA("fungus.fasta")
bacteria <- readFASTA("bacteria.fasta")
virus <- readFASTA("virus.fasta")
negative <- readFASTA("negative_Protein.fasta")

## removed wired protein
cancer <- cancer[(sapply(cancer, protcheck))]
fungus <- fungus[(sapply(fungus, protcheck))]
bacteira <- bacteira[(sapply(bacteria, protcheck))]
virus <- virus[(sapply(virus, protcheck))]
negative <- negative[(sapply(negative, protcheck))]

## function for aminoa acid composition
composition <- function(x) {
  library(protr)
  c(extractAAC(x))
}

### Generate the descriptors

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
        suppressMessages(example$CalculatePValues(new.data = test))
        p_values <- example$p.values$P.values
        p_values <- data.frame(p_values)
        p_values <- head(p_values, n = nrow(results))
        results_all <- cbind(results, p_values)
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

