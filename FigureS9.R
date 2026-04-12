library(numbat)
library(SingleCellExperiment)
library(dplyr)
library(data.table)

# ==========================================
# Numbat Workflow Execution
# ==========================================
# This script runs the Numbat pipeline for CNV inference on specific patient batches.
# It loads raw expression matrices and allele counts, filters for common cells, 
# and executes the Numbat algorithm.

# NOTE: We have already executed this workflow and generated the result files. 
# You can proceed directly to the benchmarking section to evaluate the results.

# Load the processed SingleCellExperiment object
sce <- readRDS("./numbat_results/NKTCL_NK.rds")

# Define the list of patient IDs to process
ids <- c("P03_0", "P04_0", "P09_0", "P10_0")

# ==========================================
# Main Loop: Run Numbat per Patient
# ==========================================
for (id in ids){
  cat(sprintf("Running Numbat for %s ...\n", id))
  
  # 1. Extract Raw Expression Matrix
  # Subset the SCE object by batch ID and retrieve the 'Raw' assay
  expr_mat <- assays(sce[, colData(sce)$batch == id])$Raw
  expr_cell_names <- colnames(expr_mat)
  
  # 2. Load Allele Counts Data
  # Construct path to the allele counts file for the current batch
  df_allele_path <- paste0("./numbat_results/NKTCL_EBV_", id, "_allele_counts.tsv.gz")
  df_allele <- fread(df_allele_path)
  
  # 3. Harmonize Cell Names
  # Append batch ID suffix to allele count cell names to match SCE format
  df_allele$cell <- paste0(df_allele$cell, "-", id)
  
  # 4. Filter Allele Data
  # Retain only cells that are present in the expression matrix
  df_allele_filtered <- df_allele[df_allele$cell %in% expr_cell_names, ]
  
  # Log filtering statistics
  cat(sprintf("Batch %s finished: loaded %d rows → remained %d rows after filtering\n", 
              id, nrow(df_allele), nrow(df_allele_filtered)))
  
  # 5. Execute Numbat Pipeline
  numbat_out <- numbat::run_numbat(
    count_mat = expr_mat,           # Input: Raw count matrix
    lambdas_ref = ref_hca,          # Input: Reference lambda values (e.g., from HCA)
    df_allele = df_allele_filtered, # Input: Filtered allele counts data frame
    gtf = gtf_hg38,                 # Input: Gene annotation file (GTF)
    out_dir = paste0("./numbat_results/", id), # Output: Directory for results
    genome = "hg38",                # Parameter: Reference genome version
    ncores = 24,                    # Parameter: Number of CPU cores to use
    plot = TRUE,                    # Parameter: Generate diagnostic plots
    verbose = TRUE                  # Parameter: Print detailed logs
  )
}

# ==========================================
# Numbat CNV Inference & Benchmarking
# ==========================================
# This script compares Numbat-predicted CNV status with pre-annotated labels.
# It calculates per-sample ARI and global classification metrics (Accuracy, F1).

library(ggplot2)
library(numbat)
library(dplyr)
library(glue)
library(data.table)
library(ggtree)
library(scater)
library(stringr)
library(tidygraph)
library(patchwork)
library(purrr)      # For list operations
library(mclust)     # For ARI calculation (adjustedRandIndex)

# ==========================================
# Initialization
# ==========================================
# Initialize storage lists
plist <- list()            # Store plots
all_pred_list <- list()    # Store prediction results for each sample (for merging)
ari_results <- data.frame() # Store ARI results

# Define ground truth labels based on pre-annotation
# Classify as "Tumor" if annotation contains "Malig NK", otherwise "Normal"
sce$cnv_group <- ifelse(
  str_detect(as.character(sce$anno.cnv), "Malig NK"), 
  "Tumor", 
  "Normal"
)

# Define patient list to process
patients <- c("P03_0", "P04_0", "P09_0", "P10_0")

# ==========================================
# Loop Through Patients
# ==========================================
for (patient in patients) {
  cat(sprintf("\n=== Processing Sample: %s ===\n", patient))
  
  # 1. Load Numbat Object
  nb_path <- file.path('./numbat_results', patient)
  if (!dir.exists(nb_path)) { next }
  nb <- Numbat$new(out_dir = nb_path)
  
  # 2. Cell Name Suffix Completion Logic
  # Ensure cell names in Numbat match the format in the SCE object
  raw_target_cells <- nb$clone_post %>% pull(cell)
  suffix <- paste0("-", patient)
  needs_suffix <- !endsWith(raw_target_cells, suffix)
  target_cells <- ifelse(needs_suffix, paste0(raw_target_cells, suffix), raw_target_cells)
  
  # 3. Find Intersection
  # Keep only cells present in both Numbat results and the SCE object
  valid_cells <- intersect(target_cells, colnames(sce))
  if (length(valid_cells) == 0) { warning("No matching cells found"); next }
  
  # 4. Subset SCE Object
  sce_sub <- sce[, valid_cells]
  
  # 5. Assign Numbat Predicted Results
  clone_df <- nb$clone_post %>%
    mutate(cell_corrected = ifelse(!endsWith(cell, suffix), paste0(cell, suffix), cell))
  clone_map <- clone_df %>% select(cell_corrected, compartment_opt) %>% tibble::deframe()
  sce_sub$cnv <- as.character(clone_map[valid_cells])
  
  # --- Plotting ---
  # Visualize Numbat CNV status on UMAP
  p <- plotReducedDim(sce_sub, dimred = "UMAP", colour_by = "cnv") +
    labs(title = paste0(patient, " - Numbat CNV status")) + theme(plot.title = element_text(hjust = 0.5))
  
  ggsave(file.path("./numbat_results", paste0(patient, "_umap.pdf")), p, width = 5, height = 5)
  plist[[patient]] <- p
  
  # --- Calculate ARI ---
  # Compare predicted CNV status (pred) vs. ground truth (cnv_group)
  true_labels <- sce_sub$cnv_group 
  pred_labels <- sce_sub$cnv
  
  if (length(unique(true_labels)) > 1 && length(unique(pred_labels)) > 1) {
    ari_val <- mclust::adjustedRandIndex(true_labels, pred_labels)
    ari_results <- rbind(ari_results, data.frame(Patient = patient, ARI = ari_val, N_Cells = length(valid_cells)))
    cat(sprintf("  -> ARI: %.4f\n", ari_val))
  }
  
  # Collect data for global calculation
  all_pred_list[[patient]] <- data.frame(
    cell = valid_cells,
    patient = patient,
    pred = pred_labels,
    true = true_labels,
    stringsAsFactors = FALSE
  )
}

# ==========================================
# Global Metrics Calculation
# ==========================================
cat("\n=== Calculating Global Consistency (All Samples Combined) ===\n")

# Merge data from all samples
all_data <- bind_rows(all_pred_list)

library(caret)
library(dplyr)

# --- Data Preprocessing ---
# CARET requires factors with consistent level ordering
# We define "Tumor" as the Positive class
df <- all_data %>%
  mutate(
    # Standardize naming (handle case sensitivity)
    pred_clean = ifelse(tolower(pred) %in% c("tumor", "malignant"), "Tumor", "Normal"),
    true_clean = ifelse(tolower(true) %in% c("tumor", "malignant", "malig nk"), "Tumor", "Normal"),
    
    # Convert to factor with fixed order: Normal (Negative), Tumor (Positive)
    pred_fac = factor(pred_clean, levels = c("Normal", "Tumor")),
    true_fac = factor(true_clean, levels = c("Normal", "Tumor"))
  )

# --- Calculate Metrics ---
# data = Prediction, reference = Truth, positive = Class of interest
cm <- confusionMatrix(data = df$pred_fac, reference = df$true_fac, positive = "Tumor")

# Extract Results
acc_val <- cm$overall['Accuracy']
f1_val <- 2 * cm$byClass['Sensitivity'] * cm$byClass['Pos Pred Value'] / 
  (cm$byClass['Sensitivity'] + cm$byClass['Pos Pred Value'])

cat("=== Evaluation Results ===\n")
cat(sprintf("Accuracy (ACC): %.4f\n", acc_val))
cat(sprintf("F1-Score:       %.4f\n", f1_val))
cat(sprintf("Precision:      %.4f\n", cm$byClass['Pos Pred Value']))
cat(sprintf("Recall:         %.4f\n", cm$byClass['Sensitivity']))

# Print Confusion Matrix
print(cm$table)

