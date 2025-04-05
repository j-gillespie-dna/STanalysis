---
title: "A Gentle Introduction to Spatial Transcriptomics Analysis"
author: "Jessica Gillespie"
date: "2025-04-04"
site: bookdown::bookdown_site
documentclass: book
link-citations: yes
output: html_document
---



# Welcome!
Future home of a simple, (hopefully) easy to understand spatial transcriptomics analysis pipeline for those who have never before done ST analysis.

This started as part of my masters thesis but has been continued as a way to provide some guidance to all of those trying to teach themselves a new NGS analysis techniques. The first version of this protocol will soon be published (pending link).

## Who is this NOT for
This protocol is not for anyone with extensive ST analysis experience. We will not be discussing intricate details or highly technical aspects beyond general notes to new users or "tips and tricks".

Conversely, this protocol is also not for someone completely new to next-generation sequencing (NGS), bioinformatics (BMI) analysis, or R. We will not be covering in detail datatypes, file handling, basic biological concepts involved with ST samples, or ST concepts prior to completed sequencing. There are places to find such information and repeating them here would be redundant and beyond the scope of this project. Instead, two good resources are [W3Schools](https://www.w3schools.com/r/default.asp) which will teach you the very basics of R and [DataQuest](https://www.dataquest.io/blog/learn-r-for-data-science/) which has more basics and also some simple projects where you can practice coding.

## So who is this actually for then?
The idea of this protocol is to run someone through their "first time" doing spatial transcriptomics analysis. Maybe you want to see if this is a new experiment type you would like to pursue in your lab. Maybe you just want to dabble in a new method to expand your BMI repertoire. If you want to take a simple data set and get a small hands-on experience with ST analysis, this is for you.

## Want to contribute?
The world of ST analysis is ever changing and constantly expanding. There is no guarantee this page or project will contain the most up to date information. This project is also currently managed by a single person in their spare time. If you would like to make updates, corrections, report problems, or suggest improvements, please feel free to do so using the issue tracker or discussion boards here on github or (if you are comfortable with git and coding) contributing directly to the repository.

There are several plans for the expansion of this protocol and they will be implemented as time allows. However, I hope you find this project helpful and informative in its current state.

## About this actual site
These pages are all written in R using a package called [Bookdown](https://bookdown.org/). This package allows you to write R scripts that can then be converted into html files. These files can then be uploaded to GitHub and published via GitPages. The R files are available for download on GitHub if you want a "behind the scenes" view of these pages or you want to run them without copy/pasting code from these pages.


At the bottom of every chapter, you will see output from `sessionInfo()`. This is a command in R that prints a list of current R environment infomration as well as packages that were loaded at the time of creating this page. While it is not necessary if you are creating simple webpages like this one, it is provided as part of the protocol as it can be useful in troubleshooting.


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
