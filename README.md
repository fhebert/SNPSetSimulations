The `SNPSetSimulations` R package provides functions to generate GWAS related data. It aims namely at simulating genotypic
profiles under specific dependence structures and genering binary phenotypes under logistic models. To install and use the 
package and read its vignette, run the following commands:

```{r,eval=FALSE}
install.packages("devtools")
devtools::install_github("fhebert/SNPSetSimulations",build_opts = c("--no-resave-data", "--no-manual"))
library(SNPSetSimulations)
vignette("SNPSetSimulations")
```