"""
This script performs out-of-bag (OOB) modeling using Random Forest and MDFS for feature selection. 
It processes input data to create binary and multiclass classification models for predicting city labels 
based on antimicrobial resistance (AMR) gene data. The results include important variables and AUC scores.

Usage: source('oob_modeling_with_mdfs.R')
"""

library(randomForest)
library(Boruta)
library(pROC)
library(MDFS)
library(mltools)
library(stringr)

# File with data
file <- "amrpp_full.csv"

output.file <- paste0("result_auc_", gsub(".csv", "", file), ".csv")
output.txt.file <- paste0("important_variables_", gsub(".csv", ".txt", file))

# Reading data
city.data <- read.csv(file = file)
x <- city.data[,-(1:2)]  # Removing decision and sample information
y <- as.factor(city.data$city)

# Information about classes
summary(as.factor(y))

# Initialize lists for storing results
unique.dec <- unique(y)
rf.models.1all.list <- list()
auc.1all.list <- list()
one.all.list <- list()
important.vars.1vall <- list()
important.vars.1v1 <- list()

# Loop for binary classification (one-vs-all)
for (i in 1:length(unique.dec)) {
  city <- unique.dec[i]
  mask <- y %in% city
  decision <- ifelse(mask, 1, 0)
  print(city)
  print.city <- paste(city, "versus all")
  one.all.list[i] <- print.city
  
  mdfs <- MDFS(data = x, decision = decision)
  important.variables <- names(x[,mdfs$relevant.variables])
  important.vars.1vall[[i]] <- important.variables
  
  model <- randomForest(x = x[, important.variables, drop = FALSE], y = as.factor(decision), ntree = 1000)
  rf.models.1all.list[[i]] <- model
  print(model$confusion)
  print(auc_roc(model$votes[,2], decision))
  auc.1all.list[[i]] <- auc_roc(model$votes[,2], decision)
  
  print(i)
  print("--------#######-------")
}

# Generate combinations of city pairs for binary classification (one-vs-one)
cities_combination <- combn(unique.dec, 2)
rf.models.1v1.list <- list()
auc.1v1.list <- list()
pair.list <- list()

for (i in 1:ncol(cities_combination)) {
  city1 <- cities_combination[1, i]
  city2 <- cities_combination[2, i]
  pair <- paste('+', city1, city2, '+')
  print(pair)
  pair.list[[i]] <- pair
  
  mask1 <- y %in% city1
  mask2 <- y %in% city2
  maskd <- mask1 + mask2
  decision <- ifelse(mask1, 1, 2)
  data <- cbind(decision, x)[maskd == 1, ]
  d <- data[,-1]
  
  mdfs <- MDFS(data = d, decision = data$decision)
  important.variables <- names(d[,mdfs$relevant.variables])
  important.vars.1v1[[i]] <- important.variables
  
  model <- randomForest(x = data[, important.variables, drop = FALSE], y = as.factor(data$decision), ntree = 1000)
  rf.models.1v1.list[[i]] <- model
  print(model$confusion)
  print(auc_roc(model$votes[,2], data$decision - 1))
  auc.1v1.list[[i]] <- auc_roc(model$votes[,2], data$decision - 1)
  
  print(i)
  print("-----#######------")
}

# Saving results
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