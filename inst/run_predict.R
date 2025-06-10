library(MetaboRidgePredict)

args <- commandArgs(trailingOnly = TRUE)

if ("--help" %in% args) {
  help_file <- system.file("help.txt", package = "MetaboRidgePredict")
  if (file.exists(help_file)) {
    cat("Usage: Rscript run_predict.R [options]\n\n")
    cat(readLines(help_file), sep = "\n")
  } else {
    cat("Help file not found. Basic usage:\n")
    cat("Rscript run_predict.R --input=<path> --id=<code> --out=<path> [--impute=TRUE|FALSE]\n")
    cat("Options:\n")
    cat("  --input=<path>    Path to the genotype file in PLINK --recode A format (.raw).\n")
    cat("  --id=<code>       Metabolite code (e.g., 'p23400').\n")
    cat("  --out=<path>      Output file path for predictions (default: 'predictions.csv').\n")
    cat("  --impute=TRUE|FALSE Whether to impute missing SNPs (default: FALSE).\n")
    cat("  --list            List all available metabolite codes and their corresponding model files.\n")
    cat("  --help            Display this help message and exit.\n")
  }
  quit(status = 0)
}

metabo_file <- system.file("metabo_id.csv", package = "MetaboRidgePredict")
metabo_codes <- data.table::fread(metabo_file, data.table = FALSE)

if ("--list" %in% args) {
  cat("Listing all available Metabolite Codes and Phenotype Names:\n")
  for (i in 1:nrow(metabo_codes)) {
    cat(metabo_codes$Metabolite_Code[i], ":", metabo_codes$Metabolite_Name[i], "\n")
  }
  
  quit(status = 0)
}


input <- gsub("--input=", "", args[grep("--input=", args)])
id <- gsub("--id=", "", args[grep("--id=", args)])
out <- gsub("--out=", "", args[grep("--out=", args)])


if (length(out) == 0) {
  out <- "predictions.csv"
}


impute <- gsub("--impute=", "", args[grep("--impute=", args)])
if (length(impute) == 0) {
  impute <- FALSE
} else {
  impute <- tolower(impute) == "true"  
}

if (length(input) == 0 || length(id) == 0) {
  stop("Missing --input or --id. Use --help for usage.")
}

if (!(id %in% metabo_codes$Metabolite_Code)) {
  stop("Invalid Metabolite Code. Please select a valid code from the list above.")
}

predict_metabo(
  geno_file = input, 
  metabo_code = id, 
  output_file = out, 
  impute_missing = impute
)