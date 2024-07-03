"""
This script predicts isolate origins using the Boruta feature selection method and Random Forest for classification. 
It processes input data to create binary and multiclass classification models for predicting city labels based on isolate data.

Usage: source('isolate_prediction_with_boruta.R')
"""

library(randomForest)
library(Boruta)
library(pROC)
library(MDFS)
library(mltools)
library(stringr)

setwd('data_directory_path')

file <- 'tools_iso_filt_city/bowtie_city_iso_filt.csv'
iso.file <- 'data/isolates.csv'

isolates_prediction_boruta <- function(file) {
  
  # Read input data
  city.data.full <- read.csv(file)
  y <- as.factor(city.data.full$city)
  city.data.full$city <- NULL
  isolates <- read.csv(iso.file)
  
  # Removing decision and sample information
  row.names(isolates) <- isolates$X 
  isolates <- isolates[,-1]
  tmp.y <- city.data.full$X
  tmp.yy <- strsplit(tmp.y, split = "_")
  
  city.data.full <- city.data.full[,-1]
  names.isolates <- names(isolates)
  names.city.data <- names(city.data.full)
  mask <- names.city.data %in% names.isolates 
  city.data <- city.data.full[,mask]
  x <- city.data
  
  # Information about classes
  summary(as.factor(y))
  
  # Prepare directories
  input_file <- basename(file)
  method <- 'boruta'
  results_dir <- 'results_isolates_prediction'
  method_dir <- paste0(results_dir, '/', str_replace(input_file, ".csv", ""), "/", method)
  
  if (!dir.exists(results_dir)) dir.create(results_dir)
  if (!dir.exists(paste0(results_dir, '/', str_replace(input_file, ".csv", "")))) dir.create(paste0(results_dir, '/', str_replace(input_file, ".csv", "")))
  if (!dir.exists(method_dir)) dir.create(method_dir)
  
  output.file <- paste0(method_dir, "/result_auc_", input_file)
  output.txt.file <- paste0(method_dir, "/important_variables_", str_replace(input_file, "csv", "txt"))
  prediction.output.file <- paste0(method_dir, "/prediction_", str_replace(input_file, "csv", "txt"))
  
  ################# One-vs-All Classification ##################
  
  unique.dec <- unique(y)
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
    one.all.list[[i]] <- print.city
    
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
    prediction.1all.list[[i]] <- predict(model, newdata = isolates)
    prediction.1all.list.prob[[i]] <- predict(model, newdata = isolates, type = "prob")[,2]
    i <- i + 1
  }
  
  ############## One-vs-One Classification #####################
  
  cities_combination <- combn(unique.dec, 2)
  rf.models.1v1.list <- list()
  auc.1v1.list <- list()
  prediction.1v1.list <- list()
  prediction.1v1.list.prob <- list()
  pair.list <- list()
  
  for (i in 1:ncol(cities_combination)) {
    city1 <- cities_combination[1, i]
    city2 <- cities_combination[2, i]
    pair <- paste('+', city1, city2, '+')
    pair.list[[i]] <- pair
    
    mask1 <- y %in% city1
    mask2 <- y %in% city2
    maskd <- mask1 + mask2
    decision <- ifelse(mask1, 1, 2)
    data <- cbind(decision, x)[maskd == 1, ]
    d <- data[,-1]
    
    boruta <- Boruta(x = d, y = data$decision, ntree = 1000)
    important.variables <- names(d[boruta$finalDecision == "Confirmed"])
    important.vars.1v1[[i]] <- important.variables
    
    if (length(important.variables) > 1) {
      model <- randomForest(x = d[, important.variables], y = as.factor(data$decision), ntree = 1000)
    } else {
      model <- randomForest(x = d, y = as.factor(data$decision), ntree = 1000)
    }
    
    rf.models.1v1.list[[i]] <- model
    auc.1v1.list[[i]] <- auc_roc(model$votes[,2], data[,1] - 1)
    prediction.1v1.list[[i]] <- predict(model, newdata = isolates)
    prediction.1v1.list.prob[[i]] <- predict(model, newdata = isolates, type = "prob")[,2]
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
  print("end")
  
  # Create prediction data frames
  prediction.1all.data.frame <- matrix(nrow = dim(isolates)[1], ncol = 6)
  for (i in 1:length(prediction.1all.list)) {
    prediction.1all.data.frame[,i] <- prediction.1all.list[[i]]
  }
  prediction.1all.data.frame <- as.data.frame(prediction.1all.data.frame)
  names(prediction.1all.data.frame) <- unlist(one.all.list)
  row.names(prediction.1all.data.frame) <- row.names(isolates)
  
  prediction.1v1.data.frame <- matrix(nrow = dim(isolates)[1], ncol = 15)
  for (i in 1:length(prediction.1v1.list)) {
    prediction.1v1.data.frame[,i] <- prediction.1v1.list[[i]]
  }
  prediction.1v1.data.frame <- as.data.frame(prediction.1v1.data.frame)
  names(prediction.1v1.data.frame) <- unlist(pair.list)
  row.names(prediction.1v1.data.frame) <- row.names(isolates)
  
  prediction.frame <- cbind(prediction.1all.data.frame, prediction.1v1.data.frame)
  
  prob.prediction.1all.data.frame <- matrix(nrow = dim(isolates)[1], ncol = 6)
  for (i in 1:length(prediction.1all.list)) {
    prob.prediction.1all.data.frame[,i] <- prediction.1all.list.prob[[i]]
  }
  prob.prediction.1all.data.frame <- as.data.frame(prob.prediction.1all.data.frame)
  names(prob.prediction.1all.data.frame) <- unlist(one.all.list)
  row.names(prob.prediction.1all.data.frame) <- row.names(isolates)
  
  prob.prediction.1v1.data.frame <- matrix(nrow = dim(isolates)[1], ncol = 15)
  for (i in 1:length(prediction.1v1.list)) {
    prob.prediction.1v1.data.frame[,i] <- prediction.1v1.list.prob[[i]]
  }
  prob.prediction.1v1.data.frame <- as.data.frame(prob.prediction.1v1.data.frame)
  names(prob.prediction.1v1.data.frame) <- unlist(pair.list)
  row.names(prob.prediction.1v1.data.frame) <- row.names(isolates)
  
  prob.prediction.frame <- cbind(prob.prediction.1all.data.frame, prob.prediction.1v1.data.frame)
  
  prediction.output <- paste0(method_dir, "/prediction_result_boruta_", input_file)
  prob.prediction.output <- paste0(method_dir, "/prob_prediction_result_boruta_", input_file)
  
  write.csv(prediction.frame, prediction.output)
  write.csv(prob.prediction.frame, prob.prediction.output)
}
