cp ~/research_code/sceptre-book/sceptre.qmd ~/research_code/sceptre/vignettes/sceptre.Rmd

sed -i -e '/::: callout-note/,/:::/d' sceptre.Rmd
sed -i -e 's/# The whole game {#sec-whole_game .unnumbered}//g' sceptre.Rmd

rm sceptre.Rmd-e
