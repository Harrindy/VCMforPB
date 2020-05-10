# Varying-coefficient regression analysis for pooled biomonitoring

This repository contains R programs for the article, “Varying-coefficient regression analysis for pooled biomonitorin,” by Dewei Wang, Xichen Mou, and Yan Liu. This article has been submitted for publication.

The required R packges are Rcpp and RcppArmadillo. You can install these two R packages by

        install.packages(c("RcppArmadillo","Rcpp"))
        
The necessary files are "VCM.cpp" and "Functions.R". Please download these two files into the directory of your R or Rstudio, and use the following codes to source both of them.

        source("Functions.R")
        
The "VCM.cpp" file will be sourced in "Function.R". Now you are ready to fit our varying-coefficient model.
