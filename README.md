The R package `SNPSetSimulations` provides functions to generate GWAS related data. It aims namely at simulating genotypic
profiles under specific dependence structures and genering binary phenotypes under logistic models. To install and use the
package and read its vignette, run the following commands:

```{r,eval=FALSE}
install.packages("devtools")
devtools::install_github("fhebert/SNPSetSimulations",build_opts=NULL)
library(SNPSetSimulations)
vignette("SNPSetSimulations")
```

This package was used in simulation studies conducted in some research papers, namely:
- Hébert, F., Causeur, D., & Emily, M. (2022). Omnibus testing approach for gene‐based gene‐gene interaction. Statistics in Medicine, 41(15), 2854-2878.
- Wang, S., Qian, G., & Hopper, J. (2023). Integrated logistic ridge regression and random forest for phenotype-genotype association analysis in categorical genomic data containing non-ignorable missing values. Applied Mathematical Modelling, 123, 1-22.
- John M, Lencz T. Potential application of elastic nets for shared polygenicity detection with adapted threshold selection. The International Journal of Biostatistics. 2023 Nov;19(2):417-438. DOI: 10.1515/ijb-2020-0108. PMID: 36327464; PMCID: PMC10154439.
- Qian, G., Sun, PY. Association rule mining for genome-wide association studies through Gibbs sampling. Int J Data Sci Anal (2023).

If you use this package in your research, please cite it as follows:

Hébert F., Emily M., Causeur D. (2019). SNPSetSimulations: Simulation of genotypic profiles and binary phenotypes for GWAS. R package version 1.0.