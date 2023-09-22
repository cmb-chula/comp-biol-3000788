## Load sleuth library
library('sleuth')

## Set working directory (so that files will be read from/written to this location)
setwd('PATH')

## Load metadata table
s2c <- read.table(file.path('METADATA_FILE'), header = TRUE, stringsAsFactors=FALSE)

## Preprocess data into sleuth format
#### We use bootstrapping data from kallisto here
so <- sleuth_prep(s2c, extra_bootstrap_summary = TRUE, read_bootstrap_tpm = TRUE)

## Fitting of the alternative hypothesis = condition-specific expression
so <- sleuth_fit(so, ~condition, 'full')

## Fitting of the null hypothesis = no difference across condition
so <- sleuth_fit(so, ~1, 'reduced')

## Perform Likelihood-Ratio-Test (LRT)
so <- sleuth_lrt(so, 'reduced', 'full')

## Generate result as table sorted by q-value (similar to FDR)
sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)

## Extract significant transcripts at 5% FDR
sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05)

## View the top 20 transcripts
head(sleuth_significant, 20)

## View expression level of specific transcripts
plot_bootstrap(so, 'GENE_ID', units = 'est_counts', color_by = 'condition')
plot_bootstrap(so, 'GENE_ID', units = 'tpm', color_by = 'condition')

## Output differential expression test result to file
write.table(sleuth_table, file = 'sleuth_differential_expression.txt', sep = '\t', quote = FALSE)

## Output transcript expression level to file
kallisto_tpm <- sleuth_to_matrix(so, 'obs_norm', 'tpm')
write.table(kallisto_tpm, file = 'kallisto_expression.txt', sep = '\t', quote = FALSE)