# Genotype-Based Prediction of Metabolite Levels and Ratios
### Predict metabolite expression values from genotype data using pre-trained ridge regression models.

<pre>
--input= (Required)
Path to the genotype file in PLINK `--recode A` format. This file should contain sample IDs and SNP data, where each row represents an individual and columns include metadata (`FID`, `IID`, etc.) followed by SNP genotypes (0, 1, 2, or NA).

--id= (Required)
Metabolite code (e.g., `p23400`) corresponding to a pre-trained model. This code must match an entry in `metabo_codes.csv` and a model file named `.rds` (e.g., `p23400.rds`) in the `models/` directory.

--out= (Optional)
Output file path for predictions (default: `predictions.csv`). The result will be a CSV file with two columns: `ID` and `Predicted_Value`.

--impute= (Optional)
Whether to impute missing SNPs (`TRUE` or `FALSE`, default: `FALSE`).

--list (Optional)
List all available metabolite codes and their corresponding model files.

--help
Display this help message and exit.
</pre>


## Examples

<pre>
# Predict using genotype input and a specific metabolite model
Rscript run_predict.R --input=test.raw --id=p23400 --out=predictions_test.csv

# List available metabolite models
Rscript run_predict.R --list
</pre>
## Notes
<pre>
The genotype file must be in PLINK --recode A format.

Predictions use the lambda.1se value from the pre-trained ridge regression model.

Ensure the MetaboRidgePredict package is installed before running.
</pre>
