library("caret")
library("ggplot2")
library("mixOmics")
library("pracma")

nan_norm <- function(x) {
  return(x / sqrt(sum(x[!is.na(x)]^2)))
}
accuracy_score <- function(y_true, y_pred) {
  return(sum(y_true == y_pred) / length(y_true))
}

phospho <- read.table("data/phospho.txt.gz", sep = ",", header = 1, row.names = 1)
global <- read.table("data/global.txt.gz", sep = ",", header = 1, row.names = 1)
race_labels <- read.table("data/race_label.txt.gz", sep = ",", header = 1, row.names = 1)

race_labels = subset(race_labels, race_labels$Race_label %in% list("Black", "White"))
phospho = phospho[rownames(race_labels),]
global = global[rownames(race_labels),]
phospho = t(apply(phospho, 1, nan_norm))
global = t(apply(global, 1, nan_norm))

race_labels = as.vector(race_labels$Race_label)
race_labels = as.integer(factor(race_labels)) -  1

folds <- createFolds(
  race_labels,
  k = 10, 
  list = TRUE,
  returnTrain = TRUE
)
ranks = c(1:10)
accuracies = numeric(length(ranks))

for (rank in ranks) {
  predicted = integer(length(race_labels))
  for (train_index in folds) {
    train_data <- list(
      "phospho" = phospho[train_index,],
      "global" = global[train_index,]
    )
    train_labels = race_labels[train_index]
    test_data <- list(
      "phospho" = phospho[-train_index,],
      "global" = global[-train_index,]
    )
    diablo <- block.splsda(
      train_data, 
      train_labels, 
      all.outputs = FALSE, 
      ncomp = rank
    )
    train_components <- data.frame(
      predict(diablo, train_data)$WeightedPredict[,1,]
    )
    colnames(train_components) = c(1:rank)
    train_components[,"race_label"] = race_labels[train_index]
    lr <- glm(race_label ~ ., data = train_components, family = binomial("logit"))
    
    test_components <- data.frame(
      predict(diablo, test_data)$WeightedPredict[,1,]
    )
    colnames(test_components) = c(1:rank)
    predicted[-train_index] <- ifelse(
      predict(lr, test_components) >= 0.5,
      1,
      0
    )
  }
  accuracies[rank] = accuracy_score(race_labels, predicted)
}

plot_df = data.frame(
  Ranks = ranks,
  Accuracy = accuracies
)
plot <- ggplot(data = plot_df, aes(x = Ranks, y = Accuracy))
plot = plot + geom_line()
plot = plot + scale_x_continuous(breaks = 1:10)
plot = plot + scale_y_continuous(breaks = linspace(0.5, 1, 6))
plot = plot + ylim(0.5, 1)
plot = plot + xlab("Number of Components")

# 3 components is optimal for race prediction accuracy

data = list(
  "phospho" = phospho,
  "global" = global
)
diablo <- block.splsda(
  data, 
  race_labels, 
  all.outputs = TRUE, 
  ncomp = 3
)
components = data.frame(
  predict(diablo, data)$WeightedPredict[,1,]
)
colnames(components) = c(1:3)
components[,"race_label"] = race_labels

lr <- glm(race_label ~ ., data = components, family = binomial("logit"))
comp_order <- order(abs(lr$coefficients[2:4]), decreasing = TRUE)

diablo$loadings$Y = diablo$loadings$Y[,comp_order]
diablo$loadings$phospho = diablo$loadings$phospho[,comp_order]
diablo$loadings$global = diablo$loadings$global[,comp_order]
diablo$variates$phospho = diablo$variates$phospho[,comp_order]
diablo$variates$global = diablo$variates$global[,comp_order]

colnames(diablo$loadings$Y) = c("comp1", "comp2", "comp3")
colnames(diablo$loadings$phospho) = c("comp1", "comp2", "comp3")
colnames(diablo$loadings$global) = c("comp1", "comp2", "comp3")
colnames(diablo$variates$phospho) = c("comp1", "comp2", "comp3")
colnames(diablo$variates$global) = c("comp1", "comp2", "comp3")

write.csv(diablo$loadings$Y, "data/diablo_y_loadings.csv")
write.csv(diablo$loadings$phospho, "data/diablo_phospho_loadings.csv")
write.csv(diablo$loadings$global, "data/diablo_global_loadings.csv")
write.csv(diablo$variates$phospho, "data/diablo_phospho_scores.csv")
write.csv(diablo$variates$global, "data/diablo_global_scores.csv")
