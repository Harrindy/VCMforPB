# Varying-coefficient regression analysis for pooled biomonitoring

This repository contains R programs for the article, “Varying-coefficient regression analysis for pooled biomonitorin,” by Dewei Wang, Xichen Mou, and Yan Liu. This article has been submitted for publication.

The required R packges are Rcpp and RcppArmadillo. You can install these two R packages by

        install.packages(c("RcppArmadillo","Rcpp"))
        
The necessary files are "VCM.cpp" and "Functions.R". Please download these two files into the directory of your R or Rstudio, and use the following codes to source both of them.

        source("Functions.R")
        
The "VCM.cpp" file is sourced in "Function.R" automatically. Now you are ready to fit our varying-coefficient model.

## When indiviudal-level biomonitoring is used

        Individual.fit=IndT(u_grid, Y, U, X, W, c(a,b))
        # u_grid is a grid of u-values
        # Y: a size-J vector of individual-level measurements
        # U: a size-J vector of indivdiual-level U-covariates
        # X: a J*(p+1) matrix of indiviudal-level X-covariates
        # W: a size-J vector of individual-level sampling weights
        # (a,b): the cross-validation will search h within the interval (a,b)
        
        #The output Individual.fit is a list of the selected bandwidth and the estimates
        Individual.fit$h # a single value
        Individual.fit$fit # a (p+1)*length(u_grid) matrix
        
        # You can plot the estimates by 
        par(mfrow=c(1,nrow(Individual.fit$fit)))
        for(k in 1:nrow(Individual.fit$fit))
        {
                plot(u_grid,Individual.fit$fit[k,],type="l")
        } 
             
## When randomly pooled biomonitoring is used

        Random.fit=RT(u_grid,Z,U,X,POOLID,W,c(a,b))
        # u_grid is a grid of u-values
        # Z: a size-N vector of pooled-level measurements for each individual,
        #       if two indiviudals are in the same pool, the corresponding Z-value are the same
        # U: a size-N vector of indivdiual-level U-covariates
        # X: a N*(p+1) matrix of indiviudal-level X-covariates
        # PoolID: a size-N vector of individuals' pool ID,
        #       if two individuals are in the same pool, their pool IDs are the same
        # W: a size-N vector of individual-level sampling weights
        # (a,b): the cross-validation will search h within the interval (a,b)
        
        #The output Random.fit is a list of the selected bandwidth and the estimates
        Random.fit$h # a single value
        Random.fit$fit # a (p+1)*length(u_grid) matrix
        
        # You can plot the estimates by 
        par(mfrow=c(1,nrow(Random.fit$fit)))
        for(k in 1:nrow(Random.fit$fit))
        {
                plot(u_grid,Random.fit$fit[k,],type="l")
        }
        
## When homogeneously pooled biomonitoring is used

        Homogenous.fit=HT(u_grid,Z,U,X,POOLID,W,c(a,b))
        # u_grid is a grid of u-values
        # Z: a size-N vector of pooled-level measurements for each individual,
        #       if two indiviudals are in the same pool, the corresponding Z-value are the same
        # U: a size-N vector of indivdiual-level U-covariates
        # X: a N*(p+1) matrix of indiviudal-level X-covariates
        # PoolID: a size-N vector of individuals' pool ID,
        #       if two individuals are in the same pool, their pool IDs are the same
        # W: a size-N vector of individual-level sampling weights
        # (a,b): the cross-validation will search h within the interval (a,b)
        
        #The output Homogeneous.fit is a list of the selected bandwidth and the estimates
        Homogeneous.fit$h # a single value
        Homogeneous.fit$fit # a (p+1)*length(u_grid) matrix
        
        # You can plot the estimates by 
        par(mfrow=c(1,nrow(Homogeneous.fit$fit)))
        for(k in 1:nrow(Homogeneous.fit$fit))
        {
                plot(u_grid,Homogeneous.fit$fit[k,],type="l)
        }
