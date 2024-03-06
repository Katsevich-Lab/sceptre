cp ~/research_code/sceptre-book/sceptre.qmd ~/research_code/sceptre/vignettes/sceptre.Rmd

sed -i -e '/::: callout-note/,/:::/d' sceptre.Rmd
sed -i -e 's/# The whole game {#sec-whole_game .unnumbered}//g' sceptre.Rmd
sed -i -e 's|\[@gasperini2019\]|\[(Gasperini, 2019)\](https://pubmed.ncbi.nlm.nih.gov/30612741/)|g' sceptre.Rmd
sed -i -e 's|\[@replogle2022\]|\[(Replogle, 2022)\](https://pubmed.ncbi.nlm.nih.gov/35688146/)|g' sceptre.Rmd
sed -i -e 's/{#[^}]*}//g' sceptre.Rmd
sed -i -e 's|\[`sceptre` Nextflow pipeline\](at-scale.qmd)|\[`sceptre` Nextflow pipeline\](https://timothy-barry.github.io/sceptre-book/)|g' sceptre.Rmd
sed -i -e 's/[^.]*@sec-[^.]*\.//g' sceptre.Rmd
sed -i -e 's|this book|this vignette|g' sceptre.Rmd
sed -i -e 's|`sceptre` may not exhibit good calibration initially, which is OK.|`sceptre` may not exhibit good calibration initially, which is OK. See the \[book\](https://timothy-barry.github.io/sceptre-book/) for strategies for improving calibration.|g' sceptre.Rmd
sed -i -e 's/parallel = TRUE/parallel = TRUE, n_processors = 2/g' sceptre.Rmd
sed -i -e '/## Further reading/,$d' sceptre.Rmd
echo "## Further reading
We encourage readers interested in learning more to consult the [sceptre manual](https://timothy-barry.github.io/sceptre-book/).
\`\`\`{r}
sessioninfo::session_info()
\`\`\`" >> sceptre.Rmd
sed -i '' '1i\
--- \
title: "Getting started with sceptre" \
output: rmarkdown::html_vignette \
vignette: > \
  %\\VignetteIndexEntry{Getting started with sceptre} \
  %\\VignetteEngine{knitr::rmarkdown} \
  %\\VignetteEncoding{UTF-8} \
--- \
' sceptre.Rmd

rm sceptre.Rmd-e
