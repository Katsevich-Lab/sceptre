# Gene position data frames

`gene_position_data_frame_grch38` and `gene_position_data_frame_grch37`
contain the coordinate and transcription start site position of each
gene relative to reference genome GRCh38 and GRCh37, respectively. Both
`gene_position_data_frame_grch38` and `gene_position_data_frame_grch37`
were constructed from reference genomes available on the 10x Genomics
website. The GRCh38 reference genome has been used by 10x Cell Ranger
since 2020.

## Usage

``` r
data(gene_position_data_frame_grch38)

data(gene_position_data_frame_grch37)
```

## Format

An object of class `data.table` (inherits from `data.frame`) with 36572
rows and 4 columns.

An object of class `data.table` (inherits from `data.frame`) with 57773
rows and 4 columns.

## Examples

``` r
data(gene_position_data_frame_grch38)
data(gene_position_data_frame_grch37)
head(gene_position_data_frame_grch38)
#>        response_id response_name    chr position
#>             <char>        <char> <fctr>    <int>
#> 1: ENSG00000243485   MIR1302-2HG   chr1    29554
#> 2: ENSG00000237613       FAM138A   chr1    36081
#> 3: ENSG00000186092         OR4F5   chr1    65419
#> 4: ENSG00000239945    AL627309.3   chr1    91105
#> 5: ENSG00000238009    AL627309.1   chr1   133723
#> 6: ENSG00000239906    AL627309.2   chr1   140339
head(gene_position_data_frame_grch37)
#>        response_id response_name    chr position
#>             <char>        <char> <fctr>    <int>
#> 1: ENSG00000223972       DDX11L1   chr1    11869
#> 2: ENSG00000227232        WASH7P   chr1    29806
#> 3: ENSG00000243485    MIR1302-10   chr1    29554
#> 4: ENSG00000237613       FAM138A   chr1    36081
#> 5: ENSG00000268020        OR4G4P   chr1    52473
#> 6: ENSG00000240361       OR4G11P   chr1    62948
```
