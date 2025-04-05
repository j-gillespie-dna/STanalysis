# Future Directions and Conclusion



## Deconvolution

As this step is optional, it will not be covered here in detail. However, we did want to mention two R packages, RCTD and CARD that are useful for deconvolution. RCTD uses expression profiles from a reference dataset and supervised learning to determine cell type proportions in a spot. CARD uses a gene expression reference and spatial correlation to determine cell type at each spot across a tissue. Either package requires our Seurat object with spatial data and an annotated reference RNA-seq dataset appropriate for the tissue under study.

As we have RNA-seq data, analyses typically performed with these data can also be done here. Differential gene expression (DEG) is often performed. We can examine DEG per cell group defined in the clustering step or we can also compare expression in one cluster versus the rest of the cell population in the tissue. Seurat has methods to perform this analysis. Gene set enrichment analysis (GSEA) or a pathway analysis are also possible areas for further study.

## Conclusion

Spatial transcriptomics adds another level to RNA analysis, combining gene expression and location information. While we present a simple workflow with mostly default options, there are many nuances to ST analysis. For example, many software has the option to pre-process data before clustering yet there has been evidence that pre-processing can greatly affect analysis outcome18. Similarly, normalization can affect deconvolution and it is suggested only raw spatial data should be used. For SVGs and CCC, the list of results returned and their associated statistics greatly depends on the software and database choice respectively. It has been demonstrated in multiple studies that analysis performance is highly dependent on the dataset/software pairing and further algorithm refinement is necessary.

Here we focus on one platform and a handful of software for analysis but there are many options available. Visium is only one sequencing platform with other popular options being Slide-seq, MERFISH, seqFISH, Visium HD, and the emerging 10X Genomics Xenium. Each platform has its own pros and cons and offers slightly different analysis options. There are also many software options for each analysis step mentioned here. Again, each has distinct advantages and strengths depending on the dataset and sequencing platform. In the future, we would like to expand this protocol to make it applicable to a wider variety of platforms and analysis options.



``` r
sessionInfo()
```

```
## R version 4.4.3 (2025-02-28)
## Platform: x86_64-pc-linux-gnu
## Running under: Linux Mint 21
## 
## Matrix products: default
## BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.10.0 
## LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.10.0
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## time zone: America/New_York
## tzcode source: system (glibc)
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## loaded via a namespace (and not attached):
##  [1] digest_0.6.37     R6_2.6.1          bookdown_0.42     fastmap_1.2.0    
##  [5] xfun_0.51         cachem_1.1.0      knitr_1.49        htmltools_0.5.8.1
##  [9] rmarkdown_2.29    lifecycle_1.0.4   cli_3.6.4         sass_0.4.9       
## [13] jquerylib_0.1.4   compiler_4.4.3    rstudioapi_0.17.1 tools_4.4.3      
## [17] evaluate_1.0.3    bslib_0.9.0       yaml_2.3.10       jsonlite_1.9.0   
## [21] rlang_1.1.5
```
