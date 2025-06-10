#' Predict Metabolite Expression from Genotype Data
#'
#' @description Predicts metabolite expression using genotype data and a pre-trained model.
#'
#' @param geno_file Path to genotype file (PLINK --recode A format).
#' @param metabo_code Metabolite code (e.g., "p23400").
#' @param output_file Output file name for predictions (default: "predictions.csv").
#' @param impute_missing Logical, whether to impute missing SNPs (default: FALSE).
#' @return A data frame with sample IDs and predicted expression values.
#' @export
#' @examples
#' # Interactive R usage:
#' predict_metabo("/path/to/genotype_test.raw", "p23400", "predictions_test.csv")
#' # With imputation:
#' predict_metabo("/path/to/genotype_test.raw", "p23400", "predictions_test.csv", impute_missing = TRUE)
#' metabo_file <- data.table::fread(system.file("metabo_id.csv", package = "MetaboRidgePredict"), data.table = FALSE)
#' # Command-line usage via Rscript:
#' # Rscript run_predict.R --input=genotype_test.raw --id=p23400 --out=predictions_test.csv --impute=TRUE
#' # Rscript run_predict.R --list
#' @note
#' This function can also be run from the command line using the provided script run_predict.R:
#' \preformatted{
#' Usage: Rscript run_predict.R [options]
#' Options:
#'   --input=<path>    Path to the genotype file in PLINK --recode A format (.raw).
#'   --id=<code>       Metabolite code (e.g., "p23400").
#'   --out=<path>      Output file path for predictions (default: "predictions.csv").
#'   --impute          Whether to impute missing SNPs (TRUE/FALSE, default: FALSE).
#'   --list            List all available metabolite codes and their corresponding model files.
#'   --help            Display this help message and exit.
#' }
#' For more details, see the package documentation or help.txt.
predict_metabo <- function(geno_file = NULL, metabo_code = NULL, output_file = "predictions.csv", impute_missing = FALSE) {
  library(data.table)
  library(DBI)
  library(RSQLite)

  if (is.null(geno_file) || is.null(metabo_code)) {
    stop("ERROR: Missing genotype file or metabolite code.")
  }

  if (!file.exists(geno_file)) {
    stop("ERROR: Genotype file not found: ", geno_file)
  }

  db_path <- system.file("imputation.db", package = "MetaboRidgePredict")
  con <- dbConnect(RSQLite::SQLite(), db_path)

  weights_db_path <- system.file("weights.db", package = "MetaboRidgePredict")
  weights_con <- dbConnect(RSQLite::SQLite(), weights_db_path)

  available_tables <- dbListTables(weights_con)
  if (!(metabo_code %in% available_tables)) {
    dbDisconnect(con)
    stop("ERROR: The provided metabolite code is not found in the database as a table.")
  }

  model_query <- sprintf("SELECT rsid, weight FROM '%s'", metabo_code)
  model_data <- dbGetQuery(weights_con, model_query)

  intercept_query <- sprintf("SELECT intercept FROM model_intercepts WHERE model_name = '%s'", metabo_code)
  intercept_value <- dbGetQuery(weights_con, intercept_query)
  dbDisconnect(weights_con)

  intercept <- if (nrow(intercept_value) == 0) 0 else as.numeric(intercept_value$intercept[1])

  geno_data <- fread(geno_file, data.table = FALSE)
  sample_ids <- geno_data$IID
  geno_matrix <- as.matrix(geno_data[, -(1:6)])
  colnames(geno_matrix) <- gsub("_.*", "", colnames(geno_matrix))  
  geno_matrix <- apply(geno_matrix, 2, as.numeric)

  model_snps <- model_data$rsid
  model_weights <- model_data$weight
  common_snps <- intersect(colnames(geno_matrix), model_snps)

  snp_ratio <- length(common_snps) / length(model_snps)
  cat(sprintf("WARNING: Found %.2f%% of required SNPs are available.\n", snp_ratio * 100))

  if (impute_missing && snp_ratio < 1) {
    cat("Imputing missing SNPs...\n")
    imputation_query <- sprintf("SELECT snp, fill_value FROM '%s'", metabo_code)
    fill_data <- dbGetQuery(con, imputation_query)
    missing_snps <- setdiff(model_snps, common_snps)

    for (snp in missing_snps) {
      if (snp %in% fill_data$snp) {
        fill_value <- fill_data$fill_value[fill_data$snp == snp]
        geno_matrix <- cbind(geno_matrix, fill_value)
        colnames(geno_matrix)[ncol(geno_matrix)] <- snp
      }
    }
    common_snps <- intersect(colnames(geno_matrix), model_snps)
  } else if (!impute_missing) {
    cat("Continuing with available SNPs only (no imputation).\n")
  }

  dbDisconnect(con)

  if (length(common_snps) == 0) {
    stop("ERROR: No overlapping SNPs found between genotype data and model.")
  }

  snp_weights <- model_weights[match(common_snps, model_snps)]
  pred_values <- as.vector(intercept + geno_matrix[, common_snps] %*% snp_weights)

  trans_file <- system.file("skewness_transformation_method.csv", package = "MetaboRidgePredict")
  if (!file.exists(trans_file)) {
    warning("Transformation file not found, skipping reverse transformation.")
    reverse_method <- "none"
  } else {
    trans_info <- fread(trans_file, data.table = FALSE)
    reverse_method <- trans_info$transformation[trans_info$id == metabo_code]
    if (length(reverse_method) == 0) {
      reverse_method <- "none"
    }
  }

  if (reverse_method == "log10") {
    cat("Applying reverse transformation: 10^x (log10 undo)\n")
    pred_values <- 10 ^ pred_values
  } else if (reverse_method == "square") {
    cat("Applying reverse transformation: sqrt(x) (square undo)\n")
    pred_values <- sqrt(pred_values)
  } else {
    cat("No reverse transformation applied.\n")
  }

  result <- data.frame(ID = sample_ids, Predicted_Value = pred_values)
  write.csv(result, output_file, row.names = FALSE)
  cat("Prediction results saved to:", output_file, "\n")

  return(result)
}
