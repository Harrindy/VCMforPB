# Varying-coefficient regression analysis for pooled biomonitoring

This repository contains R programs for the article, “Varying-coefficient regression analysis for pooled biomonitorin,” by Dewei Wang, Xichen Mou, and Yan Liu. This article has been submitted for publication.

The required R packges are Rcpp and RcppArmadillo. You can install these two R packages by

        install.packages(c("RcppArmadillo","Rcpp"))
        
The necessary files are "VCM.cpp" and "Functions.R". Please download the two files to the directory of your R/Rstudio. Use the following codes to source both of them.

        library(Rcpp)
        library(RcppArmadillo)
        sourceCpp("VCM.cpp")
        source("Functions.R")
        
Now you are ready to fit our varying-coefficient model.

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
        par(mfrow=c(ceiling(nrow(Individual.fit$fit)/2),2))
        par(mar=c(4,4,.1,.1))
        for(k in 1:nrow(Individual.fit$fit))
        {
                plot(u_grid,Individual.fit$fit[k,],type="l",xlab="u",ylab=expression(beta(u)))
        }
### An example
Download the "individual.cvs" file to the directory. Run the following programs.

        individual.data=read.csv("individual.csv")
        u_grid=seq(-1.5,1.5,length=400)
        Y=individual.data$Y
        U=individual.data$U
        X=cbind(individual.data$X1,individual.data$X2,individual.data$X3,individual.data$X4)
        W=individual.data$W
        Individual.fit=IndT(u_grid,Y,U,X,W,c(0.01,1))
        Individual.fit$h
        par(mfrow=c(ceiling(nrow(Individual.fit$fit)/2),2))
        par(mar=c(4,4,.1,.1))
        for(k in 1:nrow(Individual.fit$fit))
        {
                plot(u_grid,Individual.fit$fit[k,],type="l",xlab="u",ylab=expression(beta(u)))
        }

Output are

        >  Individual.fit$h
        [1] 0.5753412
        
![Optional Text](../master/individual.png)
     
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
        par(mfrow=c(ceiling(nrow(Random.fit$fit)/2),2))
        par(mar=c(4,4,.1,.1))
        for(k in 1:nrow(Random.fit$fit))
        {
                plot(u_grid,Random.fit$fit[k,],type="l",xlab="u",ylab=expression(beta(u)))
        }
        
### An example
Download the "random.cvs" file to the directory. Run the following programs.

        random.data=read.csv("random.csv")
        u_grid=seq(-1.5,1.5,length=400)
        Z=random.data$Z
        U=random.data$U
        X=cbind(random.data$X1,random.data$X2,random.data$X3,random.data$X4)
        W=random.data$W
        PoolID=random.data$poolID
        Random.fit=RT(u_grid,Z,U,X,PoolID,W,c(0.01,1))
        Random.fit$h
        par(mfrow=c(ceiling(nrow(Random.fit$fit)/2),2))
        par(mar=c(4,4,.1,.1))
        for(k in 1:nrow(Random.fit$fit))
        {
                plot(u_grid,Random.fit$fit[k,],type="l",xlab="u",ylab=expression(beta(u)))
        }

Output are

        >  Random.fit$h
        [1] 0.7529637
        
![Optional Text](../master/random.png)
        
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
        par(mfrow=c(ceiling(nrow(Homogeneous.fit$fit)/2),2))
        par(mar=c(4,4,.1,.1))
        for(k in 1:nrow(Homogeneous.fit$fit))
        {
                plot(u_grid,Homogeneous.fit$fit[k,],type="l",xlab="u",ylab=expression(beta(u)))
        }
        
### An example
Download the "homogeneous.cvs" file to the directory. Run the following programs.

        homogeneous.data=read.csv("homogeneous.csv")
        u_grid=seq(-1.5,1.5,length=400)
        Z=homogeneous.data$Z
        U=homogeneous.data$U
        X=cbind(homogeneous.data$X1,homogeneous.data$X2,homogeneous.data$X3,homogeneous.data$X4)
        W=homogeneous.data$W
        PoolID=homogeneous.data$poolID
        Homogeneous.fit=HT(u_grid,Z,U,X,PoolID,W,c(0.01,1))
        Homogeneous.fit$h
        par(mfrow=c(ceiling(nrow(Homogeneous.fit$fit)/2),2))
        par(mar=c(4,4,.1,.1))
        for(k in 1:nrow(Homogeneous.fit$fit))
        {
                plot(u_grid,Homogeneous.fit$fit[k,],type="l",xlab="u",ylab=expression(beta(u)))
        }

Output are

        >  Homogeneous.fit$h
        [1] 0.3141536
        
![Optional Text](../master/homogeneous.png)
