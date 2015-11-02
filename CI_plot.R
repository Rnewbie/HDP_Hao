library(caret)
library(conformal)
library(randomForest)
showClass("ConformalClassification")
data <- readRDS("data.Rda")
in_train <- createDataPartition(data$Label, p = 0.8, list = FALSE)
train <- data[in_train, ]
test <- data[-in_train, ]
in_proper <- createDataPartition(train$Label, p = 0.7, list = FALSE)
proper <- train[in_proper, ]
cal <- train[-in_proper, ]
trControl <- trainControl(method = "cv", number = 5, savePredictions = TRUE)
set.seed(3)
### use test set only to save time
model <- train(Label~., data = cal, method = "rf",  trControl = trControl, predict.all = TRUE)

### confidence 0.80
example <- ConformalClassification$new()
example$CalculateCVScores(model = model)
example$CalculatePValues(new.data = test)
results_confidence_0.8 <- example$ClassPredictions$aggregate
results_prediction_0.8 <- as.data.frame(results_confidence_0.8)
actual <- test$Label
data_confidence_0.8 <- cbind(results_prediction_0.8, actual)
data_confidence_0.8 <- example$p.values$Significance_p.values
write.csv(data_confidence_0.8, file  = "data_confidence_0.8.csv", row.names = FALSE)

set.seed(333)
library(reshape2)
data <- read.csv("data_confidence_0.8.csv", header = TRUE)
data <- as.data.frame(data)
data <- cbind(sample = rownames(data), data)
data_significant <- example$p.values$Significance_p.values
data_significant <- as.data.frame(data_significant)
bacteria <- data_significant$Bacteria
cancer <- data_significant$Cancer
fungus <- data_significant$Fungus
virus <- data_significant$Virus
count <- rbind(bacteria, cancer, fungus, virus)
ok <- apply(data_significant[, c("Bacteria", "Cancer", "Fungus", "Virus")], MARGIN = 1, FUN = sum)
bacteria_sum <- apply(data_significant[ c("Bacteria")], MARGIN = 2, FUN = sum)
cancer_sum <- apply(data_significant[ c("Cancer")], MARGIN = 2, FUN = sum)
fungus_sum <- apply(data_significant[ c("Fungus")], MARGIN = 2, FUN = sum)
virus_sum <- apply(data_significant[ c("Virus")], MARGIN = 2, FUN = sum)
data_boss <- rbind(bacteria_sum, cancer_sum, fungus_sum, virus_sum)
ok <- as.data.frame(ok)
ok <- ok$ok
data.m <- melt(data, id.vars = "sample")
data.m <- cbind(data.m, ok, count)
data.m$sample <- as.factor(data.m$sample)
data.m$ok <- as.factor(data.m$ok)
#levels(ok) <- ordered(ok, levels = c("0", "1", "2", "3", "4"))
#ok <- as.data.frame(ok)
### plot the stacked bar plot
#data <- data.m[, c("sample", "variable", "value", "ok")]
#data_2 <- data.frame(sample = factor(1), variable = "Null", value = 0, ok = "Zero")
#data <- rbind(data, data_2)
data <- data.m
data$ok <- as.factor(ok)
data <- data.frame(data)
p.80  <- ggplot(data) + geom_bar(aes(x = variable, y = value, fill = ok, order = desc(ok)), stat = "identity") + 
  geom_bar(aes(x = ok, y = value, fill = variable, desc(variable)), stat = "identity") +
  guides(fill = FALSE) + ggtitle("C") + scale_y_continuous(limits = c(0, 5000), breaks = c(0, 1000, 2000, 3000, 4000, 5000)) +
  theme(
    plot.title = element_text(size = 30, face = "bold", vjust = -1.9, hjust = -.3),
    panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size = 1),
    legend.position = ("none"),
    axis.text.y = element_text(color = "black", size = 20),
    axis.text.x = element_text(color = "black", size = 20),
    axis.title.y = element_text(color = "black", size = 30),
    axis.text.y = element_text(color = "black", size = 20)) + xlab(" ") + ylab("Count")
p.80


### confidence 0.85
library(caret)
library(conformal)
library(randomForest)
showClass("ConformalClassification")
data <- readRDS("data.Rda")
in_train <- createDataPartition(data$Label, p = 0.8, list = FALSE)
train <- data[in_train, ]
test <- data[-in_train, ]
in_proper <- createDataPartition(train$Label, p = 0.7, list = FALSE)
proper <- train[in_proper, ]
cal <- train[-in_proper, ]
trControl <- trainControl(method = "cv", number = 5, savePredictions = TRUE)
set.seed(3)
### use test set only to save time
model <- train(Label~., data = cal, method = "rf",  trControl = trControl, predict.all = TRUE)

example_2 <- ConformalClassification$new()
example_2$initialize(confi = 0.85)
example_2$CalculateCVScores(model = model)
example_2$CalculatePValues(new.data = test)
results_confidence_0.85 <- example_2$ClassPredictions$aggregate
results_prediction_0.85 <- as.data.frame(results_confidence_0.85)
actual <- test$Label
data_confidence_0.85 <- cbind(results_prediction_0.85, actual)
data_confidence_0.85 <- example_2$p.values$Significance_p.values

write.csv(data_confidence_0.85, file = "data_confidence_0.85.csv", row.names = FALSE)
set.seed(333)
library(reshape2)
data <- read.csv("data_confidence_0.85.csv", header = TRUE)
data <- as.data.frame(data)
data <- cbind(sample = rownames(data), data)
data_significant <- example_2$p.values$Significance_p.values
data_significant <- as.data.frame(data_significant)
bacteria <- data_significant$Bacteria
cancer <- data_significant$Cancer
fungus <- data_significant$Fungus
virus <- data_significant$Virus
count <- rbind(bacteria, cancer, fungus, virus)

#bacteria <- nrow(subset(data_significant, Bacteria == 1))
#cancer <- nrow(subset(data_significant, Cancer == 1))
#fungus <- nrow(subset(data_significant, Fungus ==1))
#virus <- nrow(subset(data_significant, Virus ==1))

ok <- apply(data_significant[, c("Bacteria", "Cancer", "Fungus", "Virus")], MARGIN = 1, FUN = sum)
ok <- as.data.frame(ok)
ok <- ok$ok
data.m <- melt(data, id.vars = "sample")
data.m <- cbind(data.m, ok, count)
data.m$sample <- as.factor(data.m$sample)
data.m$ok <- as.factor(data.m$ok)
#levels(ok) <- ordered(ok, levels = c("0", "1", "2", "3", "4"))
#ok <- as.data.frame(ok)
### plot the stacked bar plot
#data <- data.m[, c("sample", "variable", "value", "ok")]
#data_2 <- data.frame(sample = factor(1), variable = "Null", value = 0, ok = "Zero")
#data <- rbind(data, data_2)
data <- data.m
data$ok <- as.factor(ok)
data <- data.frame(data)
p.85  <- ggplot(data) +  geom_bar(aes(x = variable, y = value, fill = ok, order = desc(ok)), stat = "identity") + 
  geom_bar(aes(x = ok, y = value, fill = variable, desc(variable)), stat = "identity") +
  guides(fill = FALSE) + ggtitle("D") + scale_y_continuous(limits = c(0, 5000), breaks = c(0, 1000, 2000, 3000, 4000, 5000)) +
  theme(
    plot.title = element_text(size = 30, face = "bold", vjust = -1.9, hjust = -.3),
    panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size = 1),
    legend.position = ("none"),
    axis.text.y = element_text(color = "black", size = 20),
    axis.text.x = element_text(color = "black", size = 20),
    axis.title.y = element_text(color = "black", size = 30),
    axis.text.y = element_text(color = "black", size = 20)) + xlab(" ") + ylab("Count")
p.85


### confidence 0.90
library(caret)
library(conformal)
library(randomForest)
showClass("ConformalClassification")
data <- readRDS("data.Rda")
in_train <- createDataPartition(data$Label, p = 0.8, list = FALSE)
train <- data[in_train, ]
test <- data[-in_train, ]
in_proper <- createDataPartition(train$Label, p = 0.7, list = FALSE)
proper <- train[in_proper, ]
cal <- train[-in_proper, ]
trControl <- trainControl(method = "cv", number = 5, savePredictions = TRUE)
set.seed(3)
### use test set only to save time
model <- train(Label~., data = cal, method = "rf",  trControl = trControl, predict.all = TRUE)

example_3 <- ConformalClassification$new()
example_3$initialize(confi = 0.90)
example_3$CalculateCVScores(model = model)
example_3$CalculatePValues(new.data = test)
data_confidence_0.9 <- example_3$p.values$Significance_p.values
write.csv(data_confidence_0.9, file = "data_confidence_0.9.csv", row.names = FALSE)


set.seed(333)
library(reshape2)
data <- read.csv("data_confidence_0.9.csv", header = TRUE)
data <- as.data.frame(data)
data <- cbind(sample = rownames(data), data)
data_significant <- example_3$p.values$Significance_p.values
data_significant <- as.data.frame(data_significant)
bacteria <- data_significant$Bacteria
cancer <- data_significant$Cancer
fungus <- data_significant$Fungus
virus <- data_significant$Virus
count <- rbind(bacteria, cancer, fungus, virus)

#bacteria <- nrow(subset(data_significant, Bacteria == 1))
#cancer <- nrow(subset(data_significant, Cancer == 1))
#fungus <- nrow(subset(data_significant, Fungus ==1))
#virus <- nrow(subset(data_significant, Virus ==1))

ok <- apply(data_significant[, c("Bacteria", "Cancer", "Fungus", "Virus")], MARGIN = 1, FUN = sum)
ok <- as.data.frame(ok)
ok <- ok$ok
data.m <- melt(data, id.vars = "sample")
data.m <- cbind(data.m, ok, count)
data.m$sample <- as.factor(data.m$sample)
data.m$ok <- as.factor(data.m$ok)
#levels(ok) <- ordered(ok, levels = c("0", "1", "2", "3", "4"))
#ok <- as.data.frame(ok)
### plot the stacked bar plot
#data <- data.m[, c("sample", "variable", "value", "ok")]
#data_2 <- data.frame(sample = factor(1), variable = "Null", value = 0, ok = "Zero")
#data <- rbind(data, data_2)
data <- data.m
data$ok <- as.factor(ok)
data <- data.frame(data)
p.90  <- ggplot(data)  + geom_bar(aes(x = variable, y = value, fill = ok, order = desc(ok)), stat = "identity") + 
  geom_bar(aes(x = ok, y = value, fill = variable, desc(variable)), stat = "identity") +
  guides(fill = FALSE) + ggtitle("E") + scale_y_continuous(limits = c(0, 5000), breaks = c(0, 1000, 2000, 3000, 4000, 5000)) +
  theme(
    plot.title = element_text(size = 30, face = "bold", vjust = -1.9, hjust = -.3),
    panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size = 1),
    legend.position = ("none"),
    axis.text.y = element_text(color = "black", size = 20),
    axis.text.x = element_text(color = "black", size = 20),
    axis.title.y = element_text(color = "black", size = 30),
    axis.text.y = element_text(color = "black", size = 20)) + xlab(" ") + ylab("Count")
p.90

### p values 0.7
library(caret)
library(conformal)
library(randomForest)
showClass("ConformalClassification")
data <- readRDS("data.Rda")
in_train <- createDataPartition(data$Label, p = 0.8, list = FALSE)
train <- data[in_train, ]
test <- data[-in_train, ]
in_proper <- createDataPartition(train$Label, p = 0.7, list = FALSE)
proper <- train[in_proper, ]
cal <- train[-in_proper, ]
trControl <- trainControl(method = "cv", number = 5, savePredictions = TRUE)
set.seed(3)
### use test set only to save time
model <- train(Label~., data = cal, method = "rf",  trControl = trControl, predict.all = TRUE)

example_4 <- ConformalClassification$new()
example_4$initialize(confi = 0.70)
example_4$CalculateCVScores(model = model)
example_4$CalculatePValues(new.data = test)
data_confidence_0.7 <- example_4$p.values$Significance_p.values

write.csv(data_confidence_0.7, file = "data_confidence_0.7.csv", row.names = FALSE)


set.seed(333)
library(reshape2)
data <- read.csv("data_confidence_0.7.csv", header = TRUE)
data <- as.data.frame(data)
data <- cbind(sample = rownames(data), data)
data_significant <- example_4$p.values$Significance_p.values
data_significant <- as.data.frame(data_significant)
bacteria <- data_significant$Bacteria
cancer <- data_significant$Cancer
fungus <- data_significant$Fungus
virus <- data_significant$Virus
count <- rbind(bacteria, cancer, fungus, virus)

ok <- apply(data_significant[, c("Bacteria", "Cancer", "Fungus", "Virus")], MARGIN = 1, FUN = sum)
ok <- as.data.frame(ok)
ok <- ok$ok
data.m <- melt(data, id.vars = "sample")
data.m <- cbind(data.m, ok, count)
data.m$sample <- as.factor(data.m$sample)
data.m$ok <- as.factor(data.m$ok)
#levels(ok) <- ordered(ok, levels = c("0", "1", "2", "3", "4"))
#ok <- as.data.frame(ok)
### plot the stacked bar plot
#data <- data.m[, c("sample", "variable", "value", "ok")]
#data_2 <- data.frame(sample = factor(1), variable = "Null", value = 0, ok = "Zero")
#data <- rbind(data, data_2)
data <- data.m
data$ok <- as.factor(ok)
data <- data.frame(data)
p.70  <- ggplot(data)  + geom_bar(aes(x = variable, y = value, fill = ok, order = desc(ok)), stat = "identity") + 
  geom_bar(aes(x = ok, y = value, fill = variable, desc(variable)), stat = "identity") +
  guides(fill = FALSE) + ggtitle("A") + scale_y_continuous(limits = c(0, 5000), breaks = c(0, 1000, 2000, 3000, 4000, 5000)) +
  theme(
    plot.title = element_text(size = 30, face = "bold", vjust = -1.9, hjust = -.3),
    panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size = 1),
    legend.position = ("none"),
    axis.text.y = element_text(color = "black", size = 20),
    axis.text.x = element_text(color = "black", size = 20),
    axis.title.y = element_text(color = "black", size = 30),
    axis.text.y = element_text(color = "black", size = 20)) + xlab(" ") + ylab("Count")
p.70


### p values 7.5

library(caret)
library(conformal)
library(randomForest)
showClass("ConformalClassification")
data <- readRDS("data.Rda")
in_train <- createDataPartition(data$Label, p = 0.8, list = FALSE)
train <- data[in_train, ]
test <- data[-in_train, ]
in_proper <- createDataPartition(train$Label, p = 0.7, list = FALSE)
proper <- train[in_proper, ]
cal <- train[-in_proper, ]
trControl <- trainControl(method = "cv", number = 5, savePredictions = TRUE)
set.seed(3)
### use test set only to save time
model <- train(Label~., data = cal, method = "rf",  trControl = trControl, predict.all = TRUE)


example_5 <- ConformalClassification$new()
example_5$initialize(confi = 0.75)
example_5$CalculateCVScores(model = model)
example_5$CalculatePValues(new.data = test)
data_confidence_0.75 <- example_5$p.values$Significance_p.values

write.csv(data_confidence_0.75, file = "data_confidence_0.75.csv", row.names = FALSE)


set.seed(333)
library(reshape2)
data <- read.csv("data_confidence_0.75.csv", header = TRUE)
data <- as.data.frame(data)
data <- cbind(sample = rownames(data), data)
data_significant <- example_5$p.values$Significance_p.values
data_significant <- as.data.frame(data_significant)
bacteria <- data_significant$Bacteria
cancer <- data_significant$Cancer
fungus <- data_significant$Fungus
virus <- data_significant$Virus
count <- rbind(bacteria, cancer, fungus, virus)

ok <- apply(data_significant[, c("Bacteria", "Cancer", "Fungus", "Virus")], MARGIN = 1, FUN = sum)
ok <- as.data.frame(ok)
ok <- ok$ok
data.m <- melt(data, id.vars = "sample")
data.m <- cbind(data.m, ok, count)
data.m$sample <- as.factor(data.m$sample)
data.m$ok <- as.factor(data.m$ok)
#levels(ok) <- ordered(ok, levels = c("0", "1", "2", "3", "4"))
#ok <- as.data.frame(ok)
### plot the stacked bar plot
#data <- data.m[, c("sample", "variable", "value", "ok")]
#data_2 <- data.frame(sample = factor(1), variable = "Null", value = 0, ok = "Zero")
#data <- rbind(data, data_2)
data <- data.m
data$ok <- as.factor(ok)
data <- data.frame(data)
p.75  <- ggplot(data) +  geom_bar(aes(x = variable, y = value, fill = ok, order = desc(ok)), stat = "identity") + 
  geom_bar(aes(x = ok, y = value, fill = variable, desc(variable)), stat = "identity") +
  guides(fill = FALSE) + ggtitle("B") + scale_y_continuous(limits = c(0, 5000), breaks = c(0, 1000, 2000, 3000, 4000, 5000)) +
  theme(
    plot.title = element_text(size = 30, face = "bold", vjust = -1.9, hjust = -.3),
    panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size = 1),
    legend.position = ("none"),
    axis.text.y = element_text(color = "black", size = 20),
    axis.text.x = element_text(color = "black", size = 20),
    axis.title.y = element_text(color = "black", size = 30),
    axis.text.y = element_text(color = "black", size = 20)) + xlab(" ") + ylab("Count")
p.75
#####

conformal_plot <- plot_grid(p.70, p.75, p.80, p.85, p.90, labels = c("", "", "", "", ""), nrow = 5)
print(conformal_plot)
conformal_plot_2 <- plot_grid(p.70, p.75, p.80, p.85, p.90, labels = c("A", "B", "C", "D", "E"), nrow = 5, label_size = 30)
print(conformal_plot_2)
