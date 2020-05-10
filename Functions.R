library(Rcpp)
library(RcppArmadillo)
sourceCpp("VCM.cpp")


VaryCurve=function(u,mod=1)
{
  if(mod==2)
  {
    res=rbind(u^2/2,
              exp((-(u>0)*1.44-(u<0)/1.44)*u^2/4),
              exp(-0.72^(-1)*(u+1)^2)+0.7*exp(-0.98^(-1)*(u-1)^2),
              (0.4*sin(pi*(u-0.5)/2.5)+0.8)/exp((u+(u+0.5)^2*I(u>-0.5))/6))
  }
  if(mod==1)
  {
    res=rbind(u^3/8,
              sin(pi*u/2),
              u*(1-u)/4,
              exp(-u^2))
  }
  return(res)
}

# 
# IndT=function(u_grid,Y,U,X,W_sample,W_cv,start_ind,end_ind,kernel_type=2)
# {
#   CV.object<-function(h) #objective function of cross-validation selection of h
#   {
#     res=VCM_cv_cpp(X,U,Y,h,W_sample,W_cv,start_ind,end_ind,kernel_type)
#     return(res)
#   }
#   h=optimize(CV.object,interval=c(0.01,1))$minimum 
#   fit_ind=VCM_cpp(X,U,Y,h,u_grid,W_sample,kernel_type)
#   return(list(h=h,fit=fit_ind[1:(ncol(X)),]))
# }

##functions of indT, RT, HT, each consists of 2 versions.
#Ver1. fit the model based on original data and use cross-validation to choose h
#Ber2. Based on ver1, fixed h, bootstrap the data for bs_iter times.
#the return of ver 1 is a list of bandwith h and fit value
#the return of ver 2 is a list containing 1.a vector of h and 2. a list of the fitted values of each iteration (it is a list of matrices)

################
#individual test
################

######original funciton 
IndT=function(u_grid,Y,U,X,W,h_range,kernel_type=2)
{
  W_cv=W
  W_sample=W
  start_ind=seq(1,length(Y)) #starting index of each pool
  end_ind=seq(1,length(Y))   #ending index of each pool
  CV.object<-function(h) #objective function of cross-validation selection of h
  {
    res=VCM_cv_cpp(X,U,Y,h,W_sample,W_cv,start_ind,end_ind,kernel_type)
    return(res)
  }
  h=optimize(CV.object,interval=h_range)$minimum 
  fit=VCM_cpp(X,U,Y,h,u_grid,W_sample,kernel_type)
  return(list(h=h,fit=fit[1:ncol(X),]))
}

######bootstrap function
IndT.bs.fixh=function(u_grid,Y,U,X,W,bs_iter,h,kernel_type=2)
{ 
  original_data=cbind(W,Y,U,X)
  
  #record matrix
  record_h=rep(0,bs_iter)
  record_beta=list()
  record_beta_prime=list()
  for (b in 1:bs_iter){
    new_data=original_data[sample(1:nrow(original_data)),]
    
    W     = new_data[,1]
    Y     = new_data[,2]
    U     = new_data[,3]
    X     = new_data[,4:ncol(original_data)]
    
    W_cv=W
    W_sample=W
    start_ind=seq(1,length(Y)) #starting index of each pool
    end_ind=seq(1,length(Y))   #ending index of each pool
    
    fit=VCM_cpp(X,U,Y,h,u_grid,W_sample,kernel_type)
    record_h[b]=h
    record_beta[[b]]=fit[1:ncol(X),]
    record_beta_prime[[b]]=fit[(ncol(X)+1):(2*ncol(X)),]
    if(b%%25==0){cat("number of iteration=",b, "\n")}
  }
  return(list(h=record_h,fit=record_beta,fit.prime=record_beta_prime))
}

################
#random pooling
################

######original function
RT=function(u_grid,Z,U,X,PoolID,W,h_range,kernel_type=2) #W is the sample weight
{
  
  #start_ind,end_ind, and W_pool
  start_ind=rep(0,length(unique(PoolID)))
  end_ind=rep(0,length(unique(PoolID)))
  W_pool=rep(0,length(U))
  count=0
  for (id in unique(PoolID)){
    count=count+1
    index=which(PoolID==id)
    start_ind[count]=index[1]
    end_ind[count]=index[length(index)]
    W_pool[index]=W[index]/sum(W[index]) #normalize W_pool
  }
  
  #calcuate W_dotj
  W_dotj        =rep(0,length(U))
  count=0
  for (id in unique(PoolID)){
    count=count+1
    index=start_ind[count]:end_ind[count]
    W_dotj[index]=sum(W[index]) #normalize W_pool
  }
  
  mu_star_hat=sum(W_dotj[start_ind]*Z[start_ind]/length(index))/length(start_ind)
  
  #Step1:
  #calculate related input variables from data
  W_cv1=W^{-1}
  W_sample1=W^{-1}
  X1=W*X
  
  R1=rep(0,length(U))
  for (j in 1:length(start_ind)){
    index=start_ind[j]:end_ind[j]
    for (i in 1:length(index)){
      R1[index[i]]=W_dotj[index[i]]*Z[index[i]]-(length(index)-1)*mu_star_hat
    }
  }
  
  #
  CV.object.S1<-function(h) #objective function of cross-validation selection of h
  {
    res=VCM_cv_cpp(X1,U,R1,h,W_sample1,W_cv1,start_ind,end_ind,kernel_type)
    return(res)
  }
  h1=optimize(CV.object.S1,interval=h_range)$minimum 
  fit1_U=VCM_cpp(X1,U,R1,h1,U,W_sample1,kernel_type) #u_grid=U
  
  fit1=VCM_cpp(X1,U,R1,h1,u_grid,W_sample1,kernel_type) #u_grid=U
  
  #Step 2:
  #calculate related input variables from data
  W_cv2=W^{-1}
  W_sample2=W^{-1}
  X2=W*X
  
  R2=rep(0,length(U))
  for (j in 1:length(start_ind)){
    index=start_ind[j]:end_ind[j]
    for (i in 1:length(index)){
      R2[index[i]]=W_dotj[index[i]]*Z[index[i]]-sum(X2[index[-i],]*t(fit1_U[1:ncol(X),index[-i]]))
    }
  }
  
  CV.object.S2<-function(h) #objective function of cross-validation selection of h
  {
    res=VCM_cv_cpp(X2,U,R2,h,W_sample2,W_cv2,start_ind,end_ind,kernel_type)
    return(res)
  }
  h2=optimize(CV.object.S2,interval=h_range)$minimum 
  fit2=VCM_cpp(X2,U,R2,h2,u_grid,W_sample2,kernel_type)
  return(list(h=h1,fit=fit1[1:ncol(X),],h2=h2,fit2=fit2[1:ncol(X),]))
}

######bootstrap function
#bs_num is the number of bootstrap
RT.bs.fixh=function(u_grid,Z,U,X,PoolID,W,bs_iter=10,h1,h2,kernel_type=2) #bs_iter: bootstrap iteration numbers
{
  original_data=cbind(PoolID,W,Z,U,X)
  original_ID_table=table(original_data[,1])
  original_ID_value=as.numeric(levels(as.data.frame(original_ID_table)$Var1))
  original_poolsize=as.numeric(as.data.frame(original_ID_table)$Freq)
  
  #record matrix
  record_h1=rep(0,bs_iter)
  record_h2=rep(0,bs_iter)
  record_beta1=list()
  record_beta2=list()
  record_beta_prime1=list()
  record_beta_prime2=list()
  for (b in 1:bs_iter){
    #bootstrap the pools
    bs_ID_table=table(sort(sample(unique(original_data[,1]),replace=TRUE)))
    bs_ID_value = as.numeric(levels(as.data.frame(bs_ID_table)$Var1))
    bs_ID_freq     = as.numeric(as.data.frame(bs_ID_table)$Freq)
    bs_poolsize=rep(0,length(bs_ID_value))
    for (i in 1:length(bs_ID_value)){
      bs_poolsize[i]=original_poolsize[which(original_ID_value==bs_ID_value[i])]
    }
    new_data=matrix(nrow=sum(bs_poolsize*bs_ID_freq),ncol=ncol(original_data))
    bs_end_ind=cumsum(bs_poolsize*bs_ID_freq)
    bs_start_ind=c(1,bs_end_ind+1)
    bs_start_ind=bs_start_ind[-length(bs_start_ind)]
    for (i in 1:length(bs_ID_value)){
      temp=matrix(nrow=0,ncol=ncol(original_data))
      for (j in 1:bs_ID_freq[i]){
        temp=rbind(temp,original_data[which(original_data[,1]==bs_ID_value[i]),])
      }
      new_data[bs_start_ind[i]:bs_end_ind[i],]=temp
    }
    #reassign poolID
    new_PoolID=c()
    new_PoolID_end_ind=cumsum(bs_ID_freq)
    new_PoolID_start_ind=c(1,new_PoolID_end_ind+1)
    new_PoolID_start_ind=new_PoolID_start_ind[-length(new_PoolID_start_ind)]
    for (i in 1:length(bs_ID_value)){
      new_PoolID=c(new_PoolID,rep(new_PoolID_start_ind[i]:new_PoolID_end_ind[i],each=bs_poolsize[i]))
    }
    new_data[,1]=new_PoolID
    
    PoolID= new_data[,1]
    W     = new_data[,2]
    Z     = new_data[,3]
    U     = new_data[,4]
    X     = new_data[,5:ncol(original_data)]
    
    ####start estimation
    #start_ind,end_ind, and W_pool
    start_ind=rep(0,length(unique(PoolID)))
    end_ind=rep(0,length(unique(PoolID)))
    W_pool=rep(0,length(U))
    count=0
    for (id in unique(PoolID)){
      count=count+1
      index=which(PoolID==id)
      start_ind[count]=index[1]
      end_ind[count]=index[length(index)]
      W_pool[index]=W[index]/sum(W[index]) #normalize W_pool
    }
    
    #calcuate W_dotj
    W_dotj        =rep(0,length(U))
    count=0
    for (id in unique(PoolID)){
      count=count+1
      index=start_ind[count]:end_ind[count]
      W_dotj[index]=sum(W[index]) #normalize W_pool
    }
    
    mu_star_hat=sum(W_dotj[start_ind]*Z[start_ind]/length(index))/length(start_ind)
    
    #Step1:
    #calculate related input variables from data
    W_cv1=W^{-1}
    W_sample1=W^{-1}
    X1=W*X
    
    R1=rep(0,length(U))
    for (j in 1:length(start_ind)){
      index=start_ind[j]:end_ind[j]
      for (i in 1:length(index)){
        R1[index[i]]=W_dotj[index[i]]*Z[index[i]]-(length(index)-1)*mu_star_hat
      }
    }
    
    
    fit1_U=VCM_cpp(X1,U,R1,h1,U,W_sample1,kernel_type) #u_grid=U
    fit1=VCM_cpp(X1,U,R1,h1,u_grid,W_sample1,kernel_type) #u_grid=U
    
    #Step 2:
    #calculate related input variables from data
    W_cv2=W^{-1}
    W_sample2=W^{-1}
    X2=W*X
    
    R2=rep(0,length(U))
    for (j in 1:length(start_ind)){
      index=start_ind[j]:end_ind[j]
      for (i in 1:length(index)){
        R2[index[i]]=W_dotj[index[i]]*Z[index[i]]-sum(X2[index[-i],]*t(fit1_U[1:ncol(X),index[-i]]))
      }
    }
    fit2=VCM_cpp(X2,U,R2,h2,u_grid,W_sample2,kernel_type)
    record_h1[b]=h1
    record_h2[b]=h2
    record_beta1[[b]]=fit1[1:ncol(X),]
    record_beta2[[b]]=fit2[1:ncol(X),]
    record_beta_prime1[[b]]=fit1[(ncol(X)+1):(2*ncol(X)),]
    record_beta_prime2[[b]]=fit2[(ncol(X)+1):(2*ncol(X)),]
    if(b%%25==0){cat("number of iteration=",b, "\n")}
  }
  return(list(h1=record_h1,fit1=record_beta1,fit.prime1=record_beta_prime1,
              h2=record_h2,fit2=record_beta2,fit.prime2=record_beta_prime2))
}

##################
HT<-function(u_grid,Z,U,X,PoolID,W,h_range,kernel_type=2) #W is the sample weight; h_range is the range of h in cross validation
{
  #calculate related input variables from data
  #start_ind,end_ind, and W_pool
  start_ind=rep(0,length(unique(PoolID)))
  end_ind=rep(0,length(unique(PoolID)))
  W_pool=rep(0,length(U))
  W_cv=rep(0,length(U))
  W_sample=rep(0,length(U))
  
  count=0
  for (id in unique(PoolID)){
    count=count+1
    index=which(PoolID==id)
    start_ind[count]=index[1]
    end_ind[count]=index[length(index)]
    W_pool[index]=W[index]/sum(W[index]) #normalize W_pool
    W_cv[index]=sum(W[index])
    W_sample[index]=sum(W[index])
  }
  
  #pool X 
  X_HP=X
  for (j in 1:ncol(X_HP)){ #the intercept does not need to be pooled, thus starting from
    X_pool_temp=rep(0,length(start_ind))
    X_temp=rep(0,length(U))
    for (i in 1:length(start_ind)){
      index=start_ind[i]:end_ind[i]
      X_pool_temp[i]=sum(W_pool[index]*X_HP[index,j])
      X_temp[index]=X_pool_temp[i]
    }
    X_HP[,j]=X_temp
  }
  
  #choose h
  CV.object<-function(h) 
  {
    res=VCM_cv_cpp(X_HP,U,Z,h,W_sample,W_cv,start_ind,end_ind,kernel_type)
    return(res)
  }
  h=optimize(CV.object,interval=h_range)$minimum
  #fit model
  fit=VCM_cpp(X_HP,U,Z,h,u_grid,W_sample,kernel_type)
  return(list(h=h,fit=fit[1:(ncol(X)),]))
}

######bootstrap function
HT.bs.fixh<-function(u_grid,Z,U,X,PoolID,W,bs_iter,h=0.2,kernel_type=2) #W is the sample weight; h_range is the range of h in cross validation
{
  
  original_data=cbind(PoolID,W,Z,U,X)
  original_ID_table=table(original_data[,1])
  original_ID_value=as.numeric(levels(as.data.frame(original_ID_table)$Var1))
  original_poolsize=as.numeric(as.data.frame(original_ID_table)$Freq)
  
  #record matrix
  record_h=rep(0,bs_iter)
  record_beta=list()
  record_beta_prime=list()
  for (b in 1:bs_iter){
    #bootstrap the pools
    bs_ID_table=table(sort(sample(unique(original_data[,1]),replace=TRUE)))
    bs_ID_value = as.numeric(levels(as.data.frame(bs_ID_table)$Var1))
    bs_ID_freq     = as.numeric(as.data.frame(bs_ID_table)$Freq)
    bs_poolsize=rep(0,length(bs_ID_value))
    for (i in 1:length(bs_ID_value)){
      bs_poolsize[i]=original_poolsize[which(original_ID_value==bs_ID_value[i])]
    }
    new_data=matrix(nrow=sum(bs_poolsize*bs_ID_freq),ncol=ncol(original_data))
    bs_end_ind=cumsum(bs_poolsize*bs_ID_freq)
    bs_start_ind=c(1,bs_end_ind+1)
    bs_start_ind=bs_start_ind[-length(bs_start_ind)]
    for (i in 1:length(bs_ID_value)){
      temp=matrix(nrow=0,ncol=ncol(original_data))
      for (j in 1:bs_ID_freq[i]){
        temp=rbind(temp,original_data[which(original_data[,1]==bs_ID_value[i]),])
      }
      new_data[bs_start_ind[i]:bs_end_ind[i],]=temp
    }
    #reassign poolID
    new_PoolID=c()
    new_PoolID_end_ind=cumsum(bs_ID_freq)
    new_PoolID_start_ind=c(1,new_PoolID_end_ind+1)
    new_PoolID_start_ind=new_PoolID_start_ind[-length(new_PoolID_start_ind)]
    for (i in 1:length(bs_ID_value)){
      new_PoolID=c(new_PoolID,rep(new_PoolID_start_ind[i]:new_PoolID_end_ind[i],each=bs_poolsize[i]))
    }
    new_data[,1]=new_PoolID
    
    PoolID= new_data[,1]
    W     = new_data[,2]
    Z     = new_data[,3]
    U     = new_data[,4]
    X     = new_data[,5:ncol(original_data)]
    
    #begin estimation
    start_ind=rep(0,length(unique(PoolID)))
    end_ind=rep(0,length(unique(PoolID)))
    W_pool=rep(0,length(U))
    W_cv=rep(0,length(U))
    W_sample=rep(0,length(U))
    
    count=0
    for (id in unique(PoolID)){
      count=count+1
      index=which(PoolID==id)
      start_ind[count]=index[1]
      end_ind[count]=index[length(index)]
      W_pool[index]=W[index]/sum(W[index]) #normalize W_pool
      W_cv[index]=sum(W[index])
      W_sample[index]=sum(W[index])
    }
    
    #pool X 
    X_HP=X
    for (j in 1:ncol(X_HP)){ #the intercept does not need to be pooled, thus starting from
      X_pool_temp=rep(0,length(start_ind))
      X_temp=rep(0,length(U))
      for (i in 1:length(start_ind)){
        index=start_ind[i]:end_ind[i]
        X_pool_temp[i]=sum(W_pool[index]*X_HP[index,j])
        X_temp[index]=X_pool_temp[i]
      }
      X_HP[,j]=X_temp
    }
    
    #fit model
    fit=VCM_cpp(X_HP,U,Z,h,u_grid,W_sample,kernel_type)
    record_h[b]=h
    record_beta[[b]]=fit[1:ncol(X),]
    record_beta_prime[[b]]=fit[(ncol(X)+1):(2*ncol(X)),]
    if(b%%25==0){cat("number of iteration=",b, "\n")}
  }
  return(list(h=record_h,fit=record_beta,fit.prime=record_beta_prime))
}

