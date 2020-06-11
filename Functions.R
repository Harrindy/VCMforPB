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
  
  mu_star_hat=sum(W_dotj[start_ind]*Z[start_ind])/length(U)
  #mu_star_hat=sum(W_dotj[start_ind]*Z[start_ind]/length(index))/length(start_ind)
  
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
  
  return(list(h=h1,fit=fit1[1:ncol(X),]))
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

