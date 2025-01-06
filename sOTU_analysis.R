##### from OTUs #####
Sys.setenv(LANG = "en")

# Load necessary libraries
library(phyloseq)
library(qiime2R)
library(caret)
library(randomForest)
library(rpart)
library(glmnet)
library(kernlab)
library(ggplot2)
library(pROC)
getwd()

# Load phyloseq object named physeq
physeq <-qza_to_phyloseq(
  features="merged-table.qza",
  tree="tree.qza",
  taxonomy="taxonomy.qza",
  metadata="metadata.tsv"
)


# Assuming  phyloseq file is already loaded and object is named 'physeq'
# Ensure  metadata is included in the phyloseq object

# Extract OTU table and metadata from phyloseq object
otu_table <- as.matrix(otu_table(physeq))
metadata <- as.data.frame(sample_data(physeq))
otu_table <- t(otu_table)


######## only run in order to normalize data ########
ot1 = as.data.frame(otu_table(physeq))
ottransposed = as.data.frame(otu_table)  

# Calculate total abundance per sample
total_abundance_per_sample <- rowSums(otu_table)

# Convert absolute abundances to relative abundances (percentages)
relative_abundances <- sweep(otu_table, 1, total_abundance_per_sample, "/")

# Calculate sums of relative abundances per sample
sums_per_sample <- rowSums(relative_abundances)

# Check if all sums are approximately equal to 1
all_close_to_one <- all(abs(sums_per_sample - 1) < 1e-6)

# Print the results
print(all_close_to_one)

# Check the dimensions and structure of relative_abundances
dim(relative_abundances)
head(relative_abundances)
otu_table = as.matrix(relative_abundances)

  ########end of normalization########


#######   transform to log10 if needed    ######
small_constant <- 1e-6
abundance_matrix <- otu_table + small_constant
log10otu_table = log10(abundance_matrix)
otu_table = log10otu_table

# Define the outcome variable
outcome_var <- "diabetes"  # Replace with  actual metadata column name for diabetes diagnosis
metadata$outcome <- ifelse(metadata[, outcome_var] == "Diagnosed by doctor", "positive", "negative")
colnames(metadata) <- c("diabetes_type", "diabetes", "outcome")




# Initialize an empty list to store results
all_results <- list()

# Initialize a data frame to store aggregated variable importance scores

agg_var_imp <- data.frame(OTU = colnames(otu_table), 
                          "Decision Tree" = 0, 
                          "Elastic Net" = 0, 
                          "Random Forest" = 0, 
                          "SVM Radial Kernel" = 0)



####### Automatic loop  ###########


# Repeat the process 50 times
for (i in 1:50) {
  # Split data into training (70%) and testing (30%) datasets
  set.seed(i)  
shuffle_index <- sample(nrow(otu_table))  # create a random shuffle of row indices
otu_table_shuffled <- otu_table[shuffle_index, ]
metadata_shuffled <- metadata[shuffle_index, ]
otu_table = as.matrix(otu_table_shuffled)
metadata = as.data.frame(metadata_shuffled)


train_index <- createDataPartition(metadata_shuffled$outcome, p = 0.7, list = FALSE)
train_data <- otu_table[train_index, ]
test_data <- otu_table[-train_index, ]
dim(train_data) ####### Check the dimensions of the data correspond to the actual subsample dimensions  ###########
dim(test_data)

train_outcome <- metadata$outcome[train_index]
test_outcome <- metadata$outcome[-train_index]
train_outcome <- as.factor(train_outcome)
test_outcome = as.factor(test_outcome)
train_outcome
test_outcome
length(train_outcome)


train_data = as.data.frame(train_data)
train_data

# Calculate OTU-wise variance and select top 500 OTUs, this can be changed to select different number of OTUs

otu_variance <- apply(train_data, 2, var)
top_otus <- names(sort(otu_variance, decreasing = TRUE)[1:500])
train_data_subset <- train_data[,top_otus]
test_data_subset <- test_data[,top_otus]
train_data_subset = as.matrix(train_data_subset)
test_data_subset = as.matrix(test_data_subset)
dim(train_data_subset)
dim(test_data_subset)
td = as.data.frame(train_data_subset)

# Define training control using 5-fold cross-validation with 5 repeats
train_control <- trainControl(method = "repeatedcv", number = 5, repeats = 5, classProbs = TRUE, summaryFunction = twoClassSummary)

# Define hyperparameter tuning grids for each model
tune_grids <- list(
  "Decision Tree" = expand.grid(cp = seq(0.001, 0.1, length = 10)),
  "Elastic Net" = expand.grid(alpha = seq(0, 1, length = 10), lambda = seq(0.0001, 1, length = 10)),
  "Random Forest" = expand.grid(mtry = seq(1, ncol(train_data_subset), length = 10)),
  "SVM Radial Kernel" = expand.grid(C = 2 ^ seq(-5, 2, length = 10), sigma = 2 ^ seq(-5, 2, length = 10))
)

# Define the models with hyperparameter tuning
models <- list(
  "Decision Tree" = train(train_data_subset, train_outcome, method = "rpart", trControl = train_control, metric = "ROC", tuneGrid = tune_grids[["Decision Tree"]]),
  "Elastic Net" = train(train_data_subset, train_outcome, method = "glmnet", trControl = train_control, metric = "ROC", tuneGrid = tune_grids[["Elastic Net"]]),
  "Random Forest" = train(train_data_subset, train_outcome, method = "rf", trControl = train_control, metric = "ROC", tuneGrid = tune_grids[["Random Forest"]]),
  "SVM Radial Kernel" = train(train_data_subset, train_outcome, method = "svmRadial", trControl = train_control, metric = "ROC", tuneGrid = tune_grids[["SVM Radial Kernel"]])
)


# Evaluate models on test set
results <- lapply(models, function(model) {
  prob_predictions <- predict(model, newdata = test_data_subset, type = "prob")[, 2]
  roc_curve <- roc(test_outcome, prob_predictions)
  auc <- auc(roc_curve)
  predictions <- ifelse(prob_predictions > 0.7, "positive", "negative") ####### select the desired threshold  ###########
  confusion_matrix <- confusionMatrix(as.factor(predictions), as.factor(test_outcome))
  sensitivity <- confusion_matrix$byClass["Sensitivity"]
  specificity <- confusion_matrix$byClass["Specificity"]
  precision <- confusion_matrix$byClass["Pos Pred Value"]
  accuracy <- confusion_matrix$overall["Accuracy"]
  f1 <- confusion_matrix$byClass["F1"]
  
  return(c(AUC = auc, Sensitivity = sensitivity, Specificity = specificity, Precision = precision, Accuracy = accuracy, F1 = f1))
})

# Store results in a list
all_results[[i]] <- data.frame(Model = names(models), do.call(rbind, results))

# Assuming var_imp_scores is a list of data frames with variable importance scores
var_imp_scores <- lapply(models, varImp)




# Aggregate variable importance scores across models
for (model_name in names(var_imp_scores)) {
  if (model_name == "SVM Radial Kernel") {
    next  # Skip to the next iteration if model_name is "SVM Radial Kernel"
  }
  var_imp_df <- var_imp_scores[[model_name]]$importance
  # Match OTU names and aggregate scores
  agg_var_imp[match(var_imp_df$OTU, agg_var_imp$OTU), model_name] <- agg_var_imp[match(var_imp_df$OTU, agg_var_imp$OTU), model_name] + var_imp_df$Overall
}

}



########## START FROM HERE WHEN THE LOOP FINISHES to extract metrics, the files produced differ depending on the files needed for upstream analyses #####################


# Combine results from all iterations into a single data frame
combined_results <- do.call(rbind, all_results)

# Print combined results
print(combined_results)

# Write combined results to CSV
write.csv(combined_results, "status_500.csv", row.names = FALSE)


rf_scores  = var_imp_scores[["Random Forest"]]
en_scores  = var_imp_scores[["Elastic Net"]]
dt_scores = var_imp_scores[["Decision Tree"]]

rf_scores  = rf_scores[["importance"]]
en_scores = en_scores[["importance"]]
dt_scores = dt_scores[["importance"]]

rf_scores <- tibble::rownames_to_column(rf_scores, "OTU")
en_scores <- tibble::rownames_to_column(en_scores, "OTU")
dt_scores <- tibble::rownames_to_column(dt_scores, "OTU")

write.csv(rf_scores, "status_500_rf_scores.csv", row.names = FALSE)
write.csv(rf_scores, "status_500_en_scores.csv", row.names = FALSE)
write.csv(rf_scores, "status_500_dt_scores.csv", row.names = FALSE)


# Rename second columns
colnames(en_scores)[2] <- "Elastic_Net"
colnames(dt_scores)[2] <- "Decision_Tree"
colnames(rf_scores)[2] <- "Random_Forest"



# Merge data frames by OTU, including only second columns
merged_data <- merge(en_scores[, c("OTU", "Elastic_Net")],
                     dt_scores[, c("OTU", "Decision_Tree")],
                     by = "OTU", all = TRUE)
merged_data <- merge(merged_data,
                     rf_scores[, c("OTU", "Random_Forest")],
                     by = "OTU", all = TRUE)

write.csv(merged_data, "status_importance_scores_list.csv", row.names = FALSE)





#add taxonomy
file_path <- "taxonomyotus.tsv"
taxonomyotus <- read.delim(file_path, stringsAsFactors = FALSE)



# Ensure OTU column is named appropriately

colnames(taxonomyotus) <- c("OTU", "taxa")

# Merge the top OTUs with their taxonomy
Importancescores_with_taxonomy <- merge(merged_data, taxonomyotus, by.x = "OTU", by.y = "OTU")

# Print the top OTUs with taxonomy
print(Importancescores_with_taxonomy)

# Write the top OTUs with taxonomy to CSV
write.csv(Importancescores_with_taxonomy, "status_500_Importancescores_with_taxonomy.csv", row.names = FALSE)



# Average the variable importance scores over all iterations
agg_var_imp[,-1] <- agg_var_imp[,-1] / 50

# Print the aggregated variable importance scores
print(agg_var_imp)

# Write aggregated variable importance scores to CSV
write.csv(agg_var_imp, "status_500_aggregated_var_imp.csv", row.names = FALSE)



top_n <- 50  # Set this to the desired number of top OTUs to compile

# Average the variable importance scores over all iterations
agg_var_imp[,-1] <- agg_var_imp[,-1] / 50

# Calculate the average importance score across all models
agg_var_imp$AverageImportance <- rowMeans(agg_var_imp[,-1])

# Select the top N OTUs based on the average importance score
top_otus_importance <- agg_var_imp[order(agg_var_imp$AverageImportance, decreasing = TRUE), ][1:top_n, ]

# Print the top N hyper-important OTUs
print(top_otus_importance)

# Write the top N hyper-important OTUs to CSV
write.csv(top_otus_importance, "status_500_30pc_top_hyper_important_otus.csv", row.names = FALSE)



file_path <- "taxonomyotus.tsv"
taxonomyotus <- read.delim(file_path, stringsAsFactors = FALSE)



# Ensure OTU column is named appropriately

colnames(taxonomyotus) <- c("OTU", "taxa")

# Merge the top OTUs with their taxonomy
top_otus_with_taxonomy <- merge(top_otus_importance, taxonomyotus, by.x = "OTU", by.y = "OTU")

# Print the top OTUs with taxonomy
print(top_otus_with_taxonomy)

# Write the top OTUs with taxonomy to CSV
write.csv(top_otus_with_taxonomy, "status_500_30pc_top_hyper_important_otus_with_taxonomy.csv", row.names = FALSE)



# Extract the best lambda value from the final Elastic Net model 
#repeat with other models

best_lambda <- models[["Elastic Net"]]$bestTune$lambda

# Extract the coefficients using the best lambda value
elastic_net_coefs <- as.data.frame(as.matrix(coef(models[["Elastic Net"]]$finalModel, s = best_lambda)))

# Add a column for OTU names
elastic_net_coefs$OTU <- rownames(elastic_net_coefs)

# Rename columns for clarity
colnames(elastic_net_coefs) <- c("Coefficient", "OTU")

# Print the coefficients
print(elastic_net_coefs)
filtered_elastic_net_coefs500 <- elastic_net_coefs[elastic_net_coefs$OTU %in% top_otus_with_taxonomy$OTU, ]
filtered_elastic_net_coefs_with_taxonomy <- elastic_net_coefs[elastic_net_coefs$OTU %in% top_otus_with_taxonomy$OTU, ]
write.csv(filtered_elastic_net_coefs500, "status_500_filtered_elastic_net_coefs500.csv", row.names = FALSE)

# Merge the dataframes by the OTU column
filtered_elastic_net_coefs_with_taxonomy <- merge(filtered_elastic_net_coefs500, top_otus_with_taxonomy, by = "OTU")

# Write the merged dataframe to a CSV file
write.csv(filtered_elastic_net_coefs_with_taxonomy, "status_500_filtered_elastic_net_coefs_with_taxonomy.csv", row.names = FALSE)




