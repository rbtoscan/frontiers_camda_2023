"""
boruta_results_processing.R

This script processes the results of isolate prediction using the Boruta feature selection method.
It calculates the mean probabilities for predictions, identifies the top predictions, and prints them.

Usage: source('boruta_results_processing.R')
"""

# Load required libraries
library(randomForest)
library(Boruta)
library(pROC)
library(MDFS)
library(mltools)
library(stringr)

# Define the function to process the Boruta results
boruta_wojtek <- function(file) {
  # Set the main directory and results directory
  main_dir <- "data_directory_path/scripts"
  results_dir <- "data_directory_path/results"
  setwd(main_dir)
  
  # Read input data files
  input <- read.csv(file)
  iso.file <- "data_directory_path/data/isolates.csv"
  iso <- read.csv(iso.file)
  
  input_file <- basename(file)
  method <- 'boruta'
  method_dir <- paste0(results_dir, '/', str_replace(input_file, ".csv", ""), "/", method)
  
  # Create necessary directories
  if (!dir.exists(results_dir)) dir.create(results_dir)
  if (!dir.exists(paste0(results_dir, '/', str_replace(input_file, ".csv", "")))) dir.create(paste0(results_dir, '/', str_replace(input_file, ".csv", "")))
  if (!dir.exists(method_dir)) dir.create(method_dir)
  
  # Define output file paths
  output.file <- paste0(method_dir, "/result_auc_", input_file)
  output.txt.file <- paste0(method_dir, "/important_variables_", str_replace(input_file, "csv", "txt"))
  prediction.output.file <- paste0(method_dir, "/prediction_", str_replace(input_file, "csv", "txt"))
  
  # Read and preprocess the data
  city.data.full <- read.csv(file = file)
  y <- city.data.full$city
  city.data.full$city <- NULL
  row.names(isolates) <- isolates$X
  isolates <- isolates[,-1]
  
  city.data.full <- city.data.full[,-1]
  mask <- names(city.data.full) %in% names(isolates)
  city.data <- city.data.full[,mask]
  x <- city.data
  
  # 1 v all analysis
  unique.dec <- unique(y)
  i <- 1
  rf.models.1all.list <- list()
  auc.1all.list <- list()
  one.all.list <- list()
  prediction.1all.list <- list()
  prediction.1all.list.prob <- list()
  important.vars.1vall <- list()
  
  for (city in unique.dec) {
    mask <- y %in% city
    decision <- ifelse(mask, 1, 0)
    print(city)
    print.city <- paste(city, "versus all")
    one.all.list[i] <- print.city
    
    boruta <- Boruta(x = x, y = decision, ntree = 1000)
    important.variables <- names(x[boruta$finalDecision == "Confirmed"])
    important.vars.1vall[[i]] <- important.variables
    
    if (length(important.variables) > 1) {
      model <- randomForest(x = x[, important.variables], y = as.factor(decision), ntree = 1000)
    } else {
      model <- randomForest(x = x, y = as.factor(decision), ntree = 1000)
    }
    
    rf.models.1all.list[[i]] <- model
    auc.1all.list[[i]] <- auc_roc(model$votes[,2], decision)
    pred <- predict(model, newdata = isolates)
    prediction.1all.list[[i]] <- pred
    pred.prob <- predict(model, newdata = isolates, type = "prob")[,2]
    prediction.1all.list.prob[[i]] <- pred.prob
    
    i <- i + 1
  }
  
  # 1 v 1 analysis
  cities_combination <- combn(unique.dec, 2)
  rf.models.1v1.list <- list()
  auc.1v1.list <- list()
  prediction.1v1.list <- list()
  prediction.1v1.list.prob <- list()
  pair.list <- list()
  
  for (i in 1:dim(cities_combination)[2]) {
    city1 <- cities_combination[1, i]
    city2 <- cities_combination[2, i]
    pair <- paste('+', city1, city2, '+')
    print(pair)
    pair.list[[i]] <- pair
    
    mask1 <- y %in% city1
    mask2 <- y %in% city2
    maskd <- mask1 | mask2
    decision <- ifelse(mask1, 1, 2)
    data <- cbind(decision, x[maskd, ])
    
    d <- data[, -1]
    boruta <- Boruta(x = d, y = data$decision, ntree = 1000)
    important.variables <- names(d[boruta$finalDecision == "Confirmed"])
    important.vars.1v1[[i]] <- important.variables
    
    if (length(important.variables) > 1) {
      model <- randomForest(x = data[, important.variables], y = as.factor(data$decision), ntree = 1000)
    } else {
      model <- randomForest(x = data[, -1], y = as.factor(data$decision), ntree = 1000)
    }
    
    pred <- predict(model, newdata = isolates)
    prediction.1v1.list[[i]] <- pred
    pred.prob <- predict(model, newdata = isolates, type = "prob")[,2]
    prediction.1v1.list.prob[[i]] <- pred.prob
    
    rf.models.1v1.list[[i]] <- model
    auc.1v1.list[[i]] <- auc_roc(model$votes[,2], data[,1] - 1)
  }
  
  # Save results
  result.names <- c(unlist(one.all.list), unlist(pair.list))
  result.all <- c(unlist(auc.1all.list), unlist(auc.1v1.list))
  result.all <- as.data.frame(result.all)
  row.names(result.all) <- result.names
  write.csv(result.all, output.file)
  
  var.list <- c(important.vars.1vall, important.vars.1v1)
  names(var.list) <- result.names
  sink(output.txt.file)
  print(var.list)
  sink()
  
  # Create prediction data frames
  prediction.1all.data.frame <- as.data.frame(matrix(unlist(prediction.1all.list), nrow = length(prediction.1all.list[[1]]), byrow = FALSE))
  colnames(prediction.1all.data.frame) <- unlist(one.all.list)
  row.names(prediction.1all.data.frame) <- row.names(isolates)
  
  prediction.1v1.data.frame <- as.data.frame(matrix(unlist(prediction.1v1.list), nrow = length(prediction.1v1.list[[1]]), byrow = FALSE))
  colnames(prediction.1v1.data.frame) <- unlist(pair.list)
  row.names(prediction.1v1.data.frame) <- row.names(isolates)
  
  prediction.frame <- cbind(prediction.1all.data.frame, prediction.1v1.data.frame)
  
  prob.prediction.1all.data.frame <- as.data.frame(matrix(unlist(prediction.1all.list.prob), nrow = length(prediction.1all.list.prob[[1]]), byrow = FALSE))
  colnames(prob.prediction.1all.data.frame) <- unlist(one.all.list)
  row.names(prob.prediction.1all.data.frame) <- row.names(isolates)
  
  prob.prediction.1v1.data.frame <- as.data.frame(matrix(unlist(prediction.1v1.list.prob), nrow = length(prediction.1v1.list.prob[[1]]), byrow = FALSE))
  colnames(prob.prediction.1v1.data.frame) <- unlist(pair.list)
  row.names(prob.prediction.1v1.data.frame) <- row.names(isolates)
  
  prob.prediction.frame <- cbind(prob.prediction.1all.data.frame, prob.prediction.1v1.data.frame)
  
  write.csv(prediction.frame, prediction.output)
  write.csv(prob.prediction.frame, prob.prediction.output)
}

# Example usage
# boruta_wojtek('data_directory_path/data/training.csv')