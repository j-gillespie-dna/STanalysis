# Materials



## Computer Infrastructure

As R will be the main working environment any operating system capable of running R is sufficient. ST datasets can be quite large (10’s of GB) and saving data at several stages in the process may be desired so 10-20GB of space may be desired. A multi-core processor for faster analysis may be beneficial at some steps, though not required. Finally, a fairly large amount of RAM may be required (>16GB). High Performance Computing services are available allowing users to analyze their data on a cluster, although that procedure is not detailed here.

## Software

R is the coding environment in which the analysis will take place. R can be run via the command-line or used in the RStudio GUI. Packages can then be installed from within R. Version numbers used in this protocol are provided but these specific versions are not mandatory.

[R (v4.4.2)](https://www.r-project.org/)
[RStudio (RStudio 2024.12.0+467 "Kousa Dogwood")](https://posit.co/downloads/)
[Seurat5 (v5.2.0)](https://satijalab.org/seurat/)
[SPARK-X6](https://xzhoulab.github.io/SPARK/)
[CellChat7 (v2.1.0)](https://github.com/jinworks/CellChat)

## Collect necessary files

### Space Ranger output

After sequencing is completed, 10X Visium data will be provided as either base call files (.BCL) or FASTQ files, either of which can be used for input to 10X Genomics’ Space Ranger software. The output of Space Ranger is the input for our first analysis method so we will start there. Space Ranger provides a variety of files containing data and various metrics on the samples and sequencing run. We are interested in one file and one folder. We need the filtered feature-barcode matrices: hdf5 which holds the barcodes, features, and RNA count data. We will also need the spatial folder which contains data pertaining to the image that accompanies the sample. For a full picture of output file structure, see the [10X Genomics website](https://www.10xgenomics.com/support/software/space-ranger/latest/analysis/outputs/output-overview).

### Reference data

In order to mark cell types, we need to compare gene expression to that of already determined cells. This is frequently done with a reference dataset. While the software mentioned here has some reference sets built in, there are others available for download. These can be used after the clustering step for cell-type annotation and in deconvolution. Some annotated datasets are available from the [HuBMAP Consortium](https://azimuth.hubmapconsortium.org/).

### Example datasets

10X Genomics has freely available many datasets produced with their technology for users to download and analyze. These are useful for practice and learning. A mouse anterior brain dataset will be used in this chapter and is available [here](https://www.10xgenomics.com/datasets/mouse-brain-serial-section-2-sagittal-anterior-1-standard). NOTE: If you click this link from an institution's network you will be asked to sign up for an account. However, if you click this link from a home network, you won't have to provide an email to download a free dataset.

**To follow along with this tutorial**
After clicking the link, scroll down on the page and click on "Output and supplemental files". Download the "Feature / barcode matrix HDF5 (filtered)" and place it in the `datasets` folder. Also download the "Spatial imagine data", a compressed folder named "V1_Mouse_Brain_Sagittal_Anterior_Section_2_spatial.tar.gz" into the `dtasets` folder. Extracting this should give you a folder named `spatial` with various files inside. Leave this as is for now.


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
