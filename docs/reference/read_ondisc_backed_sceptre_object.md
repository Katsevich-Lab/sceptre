# Write or read an `ondisc`-backed `sceptre_object`

`write_ondisc_backed_sceptre_object()` and
`read_ondisc_backed_sceptre_object()` enable the writing and reading of
`ondisc`-backed `sceptre_object`s, respectively. First,
`write_ondisc_backed_sceptre_object()` writes an `ondisc`-backed
`sceptre_object` to disk, creating a file `sceptre_object.rds` in the
specified directory. Next, `read_ondisc_backed_sceptre_object()` reads
and initializes a `sceptre_object` from a `sceptre_object.rds` file,
`response.odm` file, and `grna.odm` file stored on disk.

## Usage

``` r
read_ondisc_backed_sceptre_object(
  sceptre_object_fp,
  response_odm_file_fp,
  grna_odm_file_fp
)

write_ondisc_backed_sceptre_object(sceptre_object, directory_to_write)
```

## Arguments

- sceptre_object_fp:

  file path to a `sceptre_object.rds` file

- response_odm_file_fp:

  file path to the backing `.odm` file for the response modality

- grna_odm_file_fp:

  file path to the backing `.odm` file for the gRNA modality

- sceptre_object:

  a `sceptre_object`

- directory_to_write:

  the directory in which to write the `sceptre_object.rds` file

## Value

`write_ondisc_backed_sceptre_object()` returns NULL, and
`read_ondisc_backed_sceptre_object()` returns an `ondisc`-backed
`sceptre_object`

## Examples

``` r
data(lowmoi_example_data)
# 1. create ondisc-backed sceptre_object
sceptre_object <- import_data(
  response_matrix = lowmoi_example_data$response_matrix,
  grna_matrix = lowmoi_example_data$grna_matrix,
  grna_target_data_frame = lowmoi_example_data$grna_target_data_frame,
  extra_covariates = lowmoi_example_data$extra_covariates,
  moi = "low",
  use_ondisc = TRUE,
  directory_to_write = tempdir()
)

# 2. write
write_ondisc_backed_sceptre_object(
  sceptre_object = sceptre_object,
  directory_to_write = tempdir()
)

# 3. read
rm(sceptre_object)
sceptre_object <- read_ondisc_backed_sceptre_object(
  sceptre_object_fp = paste0(tempdir(), "/sceptre_object.rds"),
  response_odm_file_fp = paste0(tempdir(), "/response.odm"),
  grna_odm_file_fp = paste0(tempdir(), "/grna.odm")
)
```
