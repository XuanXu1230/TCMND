# Load required libraries
library(dplyr)
library(data.table)
library(parallel)

# Load data
load("~/AD/PD/GEOdata/deg0.1/GSE12649-DEG0.1step2output.Rdata")

# Initialize results storage
enrichment_results <- list()

# Process combined_data to merge genes by ingredient
processed_data <- combined_data %>%
  group_by(INGREDIENT_ENG, HERB_CHN, HERB_ENG, HERB_PY, PUBCHEM_CID) %>%
  summarise(
    Combined_Targets = paste(unique(RELATED_TARGET), collapse = ", "),
    .groups = "drop"
  )

# Convert to data.table for faster operations
module_summary <- as.data.table(module_summary)
processed_data <- as.data.table(processed_data)

# Define total number of genes
total_genes <- dim(exp_sub)[1]

# Function to process a single module
process_module <- function(module_name) {
  # Get genes in the module
  module_genes <- unlist(strsplit(module_summary[Module == module_name, ]$Genes, split = ";"))
  m <- length(module_genes)  # Number of genes in the module
  n <- total_genes - m       # Number of genes not in the module
  
  # Perform enrichment analysis for each ingredient
  results <- processed_data[, {
    target_genes <- unlist(strsplit(Combined_Targets, split = ",\\s*"))
    overlap_genes <- intersect(module_genes, target_genes)
    k <- length(overlap_genes)
    N <- length(target_genes)
    
    # Hypergeometric test for enrichment
    p_value <- phyper(k - 1, m, n, N, lower.tail = FALSE)
    
    # Return results
    list(
      Module = module_name,
      INGREDIENT_CHN = HERB_CHN,
      INGREDIENT_ENG = INGREDIENT_ENG,
      Related_Herbs = HERB_ENG,
      Related_Herbs_pinyin = HERB_PY,
      Overlap_Genes = paste(overlap_genes, collapse = ", "),
      P_Value = p_value,
      PUBCHEM_CID = PUBCHEM_CID
    )
  }, by = 1:nrow(processed_data)]
  
  return(results)
}

# Parallel processing of all modules
cl <- makeCluster(10)  # Create a cluster with 10 cores
clusterExport(cl, c("module_summary", "processed_data", "exp_sub", "total_genes", "phyper"))
clusterEvalQ(cl, library(data.table))

enrichment_results <- parLapply(cl, module_summary$Module, process_module)
stopCluster(cl)  # Stop the cluster

# Combine results from all modules
final_results <- rbindlist(enrichment_results)

# Sort results by P-value
final_results <- final_results[order(P_Value)]

# Filter significant results (P-value <= 0.05)
top_ingredient <- final_results[P_Value <= 0.05]

# Calculate module importance sum
top_ingredient$Module_Importance_Sum <- sapply(1:nrow(top_ingredient), function(i) {
  module <- top_ingredient$Module[i]
  overlap_genes <- unlist(strsplit(top_ingredient$Overlap_Genes[i], ",\\s*"))
  
  # Sum importance values for the module
  module_summary[Module == module, sum(Total_Importance, na.rm = TRUE)]
})

# Sort by module importance sum
top_ingredient <- top_ingredient[order(-Module_Importance_Sum)]

# Join with additional data and process
ingredient_result2 <- top_ingredient %>%
  left_join(all_data, by = c("Overlap_Genes" = "Gene", "Module" = "Module")) %>%
  group_by(Module, INGREDIENT_ENG, INGREDIENT_CHN, PUBCHEM_CID, Related_Herbs_pinyin, Overlap_Genes, P_Value, Module_Importance_Sum, Importance) %>%
  summarise(Total_Importance = sum(Importance, na.rm = TRUE), .groups = 'drop') %>%
  arrange(desc(Module_Importance_Sum), desc(Total_Importance)) %>%
  separate_rows(Related_Herbs_pinyin, sep = ";") %>%
  mutate(Related_Herbs_pinyin = trimws(Related_Herbs_pinyin)) %>%
  mutate(Module = gsub(".csv", "", Module)) %>%
  distinct()

# Summarize ingredient scores
Ingredient_result_sum2 <- ingredient_result2 %>%
  group_by(INGREDIENT_ENG) %>%
  summarise(INGREDIENT_Total_Importance = sum(Total_Importance, na.rm = TRUE), .groups = 'drop') %>%
  arrange(desc(INGREDIENT_Total_Importance)) %>%
  distinct()

# Handle NA values
ingredient_result2 <- ingredient_result2 %>%
  mutate(Importance = ifelse(is.na(Importance), 0, Importance))

# Define weights for scoring
w1 <- 1/3
w2 <- 1/3
w3 <- 1/3

# Normalize and calculate scores
ingredient_result2_mut <- ingredient_result2 %>%
  mutate(
    normalized_P_Value = (P_Value - min(P_Value, na.rm = TRUE)) / 
      (max(P_Value, na.rm = TRUE) - min(P_Value, na.rm = TRUE)),
    normalized_Module_Importance_Sum = (Module_Importance_Sum - min(Module_Importance_Sum, na.rm = TRUE)) / 
      (max(Module_Importance_Sum, na.rm = TRUE) - min(Module_Importance_Sum, na.rm = TRUE)),
    normalized_Importance = ifelse(
      max(Importance, na.rm = TRUE) == min(Importance, na.rm = TRUE),
      0,
      (Importance - min(Importance, na.rm = TRUE)) / 
        (max(Importance, na.rm = TRUE) - min(Importance, na.rm = TRUE))
    ),
    Score = w1 * normalized_P_Value + w2 * normalized_Module_Importance_Sum + w3 * normalized_Importance
  )

# Summarize ingredient scores
setwd("~/AD/PD")
ingredient_scores <- ingredient_result2_mut %>%
  group_by(INGREDIENT_CHN) %>%
  summarise(
    Sum_Score = sum(Score, na.rm = TRUE),
    Related_INGREDIENT = paste(unique(INGREDIENT_ENG), collapse = "; "),
    Related_INGREDIENT_CID = paste(unique(PUBCHEM_CID), collapse = "; "),
    Related_Herbs_pinyin = paste(unique(Related_Herbs_pinyin), collapse = "; "),
    Overlap_Genes = paste(unique(Overlap_Genes), collapse = "; "),
    .groups = "drop"
  ) %>%
  arrange(desc(Sum_Score))

# Save results
write.csv(ingredient_scores, "ingredient_scores.csv", row.names = FALSE)

# Print results
print(ingredient_scores)