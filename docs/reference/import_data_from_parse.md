# Import data from Parse (experimental)

`import_data_from_parse()` imports data from the output of the Parse
count matrix generation program. See [Chapter 1 of the
manual](https://timothy-barry.github.io/sceptre-book/import-data.html#import-from-the-parse-program-experimental)
for more information about this function.

## Usage

``` r
import_data_from_parse(
  gene_mat_fp,
  grna_mat_fp,
  all_genes_fp,
  all_grnas_fp,
  moi,
  grna_target_data_frame,
  extra_covariates = data.frame()
)
```

## Arguments

- gene_mat_fp:

  file path to the gene `count_matrix.mtx` file.

- grna_mat_fp:

  file path to the gRNA `count_matrix.mtx` file.

- all_genes_fp:

  file path to the `all_genes.csv` file containing the gene IDs.

- all_grnas_fp:

  file path to the `all_guides.csv` file containing the gRNA IDs. The
  gRNA IDs are assumed to be in the second column (i.e., the "gene_name"
  column) of this file.

- moi:

  a string indicating the MOI of the dataset, either "low" or "high".

- grna_target_data_frame:

  a data frame containing columns `grna_id` and `grna_target` mapping
  each individual gRNA to its target. Non-targeting gRNAs should be
  assigned a label of "non-targeting". Optionally,
  `grna_target_data_frame` can contain columns `chr`, `start`, and
  `end`, giving the chromosome, start coordinate, and end coordiante,
  respectively, of each gRNA. Additionally, `grna_target_data_frame` can
  contain the column `vector_id` specifying the vector to which a given
  gRNA belongs.

- extra_covariates:

  (optional) a data frame containing extra covariates (e.g., batch,
  biological replicate) beyond those that `sceptre` can compute.

## Value

an initialized `sceptre_object`

## Details

`import_data_from_parse()` is experimental, and the API of this function
is subject to change. We expect the API to solidify as we learn more
about the Parse platform and the structure of the Parse count matrix
generation program output.

The example data shipped with the package (used in the examples below)
are a small subset of Parse Biosciences' "Technical Performance of
CRISPR Detect in Cell Lines" demo dataset, used under the Creative
Commons Attribution 4.0 International (CC BY 4.0) license, with
attribution to Parse Biosciences. See
<https://www.parsebiosciences.com/datasets/technical-performance-of-crispr-detect-in-cell-lines/>.

## Examples

``` r
directory <- paste0(
  system.file("extdata", package = "sceptre"),
  "/parse_example/"
)
gene_mat_fp <- paste0(directory, "gene_mat.mtx")
grna_mat_fp <- paste0(directory, "grna_mat.mtx")
all_genes_fp <- paste0(directory, "all_genes.csv")
all_grnas_fp <- paste0(directory, "all_grnas.csv")
grna_target_data_frame <- data.frame(
  grna_id = c("guide_A", "guide_B", "guide_C"),
  grna_target = c("target-A", "target-B", "non-targeting")
)
sceptre_object <- import_data_from_parse(
  gene_mat_fp = gene_mat_fp,
  grna_mat_fp = grna_mat_fp,
  all_genes_fp = all_genes_fp,
  all_grnas_fp = all_grnas_fp,
  moi = "low",
  grna_target_data_frame = grna_target_data_frame
)
```
