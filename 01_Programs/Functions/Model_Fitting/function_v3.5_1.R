
set.seed(1234)
suppressPackageStartupMessages({
library(stats)
library(pracma)
library(MASS)
library(geepack)
library(devtools)
})
#library(SMARTutils)
options(digits=6)

CS=function(rho,n){
  return((matrix(rho,n,n)+diag(1-rho,n)))
}

AR1=function(rho,n){
  result=toeplitz(rho^(0:(n-1)))
  return(result)
}

#change from conditional parameters (treatment pathway level mean and variance) to marginal parameters (AI level)
conditional_to_marginal=function(p1,p2,conditional_paras){
  ans=matrix(0,4,2)
  colnames(ans)=c("mu","sigma2")
  ans[1,1]=conditional_paras[1,1]*p1+conditional_paras[2,1]*(1-p1)
  ans[2,1]=conditional_paras[1,1]*p1+conditional_paras[3,1]*(1-p1)
  ans[3,1]=conditional_paras[4,1]*p2+conditional_paras[5,1]*(1-p2)
  ans[4,1]=conditional_paras[4,1]*p2+conditional_paras[6,1]*(1-p2)
  
  ans[1,2]=conditional_paras[1,2]*p1+conditional_paras[2,2]*(1-p1)+p1*(1-p1)*(conditional_paras[1,1]-conditional_paras[2,1])^2
  ans[2,2]=conditional_paras[1,2]*p1+conditional_paras[3,2]*(1-p1)+p1*(1-p1)*(conditional_paras[1,1]-conditional_paras[3,1])^2
  ans[3,2]=conditional_paras[4,2]*p2+conditional_paras[5,2]*(1-p2)+p2*(1-p2)*(conditional_paras[4,1]-conditional_paras[5,1])^2
  ans[4,2]=conditional_paras[4,2]*p2+conditional_paras[6,2]*(1-p2)+p2*(1-p2)*(conditional_paras[4,1]-conditional_paras[6,1])^2
  return(ans)
}


#Given original conditional parameters and desired effect size, adjust the parameters to achieve the aim
#  the adjustment is made by changing the coefficients of covariates, as well as modifying the variance and correlation
#conditional_paras_original is a 6*3 matrix containing the starting mean, variance and ICC for each pathway
#  the conditional mean effect, and the proportion of variance and ICC between pathways will be kept with the adjustment
#If for no covariates, then just not inputing effectsize_eta
#aimed_comparison is for the aimed comparison of Adaptive Interventions (default is comparing AI(1,1) vs AI(-1,-1))
#  format of aimed_comparison: should be a vector length in 4, containing 2 0's and 2 1's
parameter_adjust=function(conditional_paras_original,p1,p2,aimed_comparison=c(1,0,0,1),
                          effectsize_dtr,effectsize_eta=NULL,rhoadjust,sigma2_x=NULL){
  coeff=matrix(c(1,1,1,1,1,1,-1,-1,1,-1,1,-1,1,-1,-1,1),4,4)  #the matrix that changes coefficients to marginal mean
  q=4
  marginal_paras=conditional_to_marginal(p1,p2,conditional_paras_original)
  beta_margin=solve(coeff) %*% as.matrix(marginal_paras[,1],4,1)
  
  if (!(all.equal(aimed_comparison[order(aimed_comparison)],c(0,0,1,1))==TRUE)){
    stop("Error: aimed_comparison should be a length 4 vector with two 0's and two 1's")
  }
  num_1=which(aimed_comparison==1)[1]; num_2=which(aimed_comparison==1)[2]  #which two AIs are being compared
  
  weight_beta=coeff[num_1,]-coeff[num_2,]
  sum_beta=sum(weight_beta*beta_margin)
  
  if (is.null(effectsize_eta)){
    p=0
    eta=NULL
  }else{
    if (sum(effectsize_eta^2)>=1) stop("Error: total effect size of covariates too large, their sum of squares exceeding 1")
    p=length(effectsize_eta)
    if (is.null(sigma2_x)) sigma2_x=rep(1,p)
    if (length(sigma2_x)!=p) stop("Error: different length of effect size and variance of covariates")
    eta=effectsize_eta/effectsize_dtr*sum_beta/sqrt(sigma2_x)
  }
  sum_varx=sum(eta^2*sigma2_x)
  
  coef1=2*((sum_beta/effectsize_dtr)^2-sum_varx)
  
  coeff_new=matrix(0,4,2)
  #the within-pathway part that can be duplicated by varianceadjust
  coeff_new[1,1]=conditional_paras_original[1,2]*p1+conditional_paras_original[2,2]*(1-p1)
  coeff_new[2,1]=conditional_paras_original[1,2]*p1+conditional_paras_original[3,2]*(1-p1)
  coeff_new[3,1]=conditional_paras_original[4,2]*p2+conditional_paras_original[5,2]*(1-p2)
  coeff_new[4,1]=conditional_paras_original[4,2]*p2+conditional_paras_original[6,2]*(1-p2)
  #the between-pathway part that cannot be duplicated
  coeff_new[1,2]=p1*(1-p1)*(conditional_paras_original[1,1]-conditional_paras_original[2,1])^2
  coeff_new[2,2]=p1*(1-p1)*(conditional_paras_original[1,1]-conditional_paras_original[3,1])^2
  coeff_new[3,2]=p2*(1-p2)*(conditional_paras_original[4,1]-conditional_paras_original[5,1])^2
  coeff_new[4,2]=p2*(1-p2)*(conditional_paras_original[4,1]-conditional_paras_original[6,1])^2
  
  a=t(as.matrix(aimed_comparison)) %*% coeff_new
  varianceadjust=(coef1-a[2])/a[1]
  if (varianceadjust<0) stop("Error: total effect size too large, consider reducing effect of covariates or difference between pathway mean effects")
  
  conditional_paras=conditional_paras_original
  conditional_paras[,2]=conditional_paras[,2]*varianceadjust
  conditional_paras[,3]=conditional_paras[,3]*rhoadjust
  marginal_paras=conditional_to_marginal(p1,p2,conditional_paras)
  
  result=list(beta_margin=beta_margin,conditional_paras=conditional_paras,marginal_paras=marginal_paras,eta=eta)
  return(result)
}


#Default setting is N=100, cluster size constant=10, effectsize_DTR=0.5, effectsize_eta=0 (no cluster-level covariates)
#eta_x: an array of coefficients of covariates (if not imputed, assuming there is no cluster-level covariates)
#sigma2_x: an array of variances of each covariate, default is all 1
#detail_x: an array specifying each covariate being cluster-level(0) or individual_level(1), default is cluster level
#p1 and p2 are response rate corresponding to different first-stage treatments
#N is cluster number, Mmin and Mmax are lower and upper bound of cluster size (and the actual size is uniformly drawn between)
generate_SMART=function(conditional_paras=matrix(c(6,8,4,4,4,4,75.4461,61.1113,48.2855,36.9686,27.1606,27.1606,0.20,0.18,0.16,0.14,0.12,0.10),6,3),
                        eta_x=NULL,sigma2_x=NULL,detail_x=NULL,
                        p1=0.25,p2=0.55,N=100,Mmax=10,Mmin=10,seed=1234){
  coeff=matrix(c(1,1,1,1,1,1,-1,-1,1,-1,1,-1,1,-1,-1,1),4,4)
  
  marginal_paras=conditional_to_marginal(p1,p2,conditional_paras)
  beta_margin=rep(0,4)
  beta_margin[1]=(marginal_paras[1,1]+marginal_paras[2,1]+marginal_paras[3,1]+marginal_paras[4,1])/4
  beta_margin[2]=(marginal_paras[1,1]+marginal_paras[2,1]-marginal_paras[3,1]-marginal_paras[4,1])/4
  beta_margin[3]=(marginal_paras[1,1]-marginal_paras[2,1]+marginal_paras[3,1]-marginal_paras[4,1])/4
  beta_margin[4]=(marginal_paras[1,1]-marginal_paras[2,1]-marginal_paras[3,1]+marginal_paras[4,1])/4
  q=4
  
  if (is.null(eta_x)){
    p=0
    beta=beta_margin
  }else{
    p=length(eta_x)
    beta=c(beta_margin,eta_x)
    if (is.null(sigma2_x)) sigma2_x=rep(1,p)
    if (is.null(detail_x)) detail_x=rep(0,p)
  }
  
  #generating the belonging pathway and corresponding AI for each cluster
  set.seed(seed)
  correct=0
  while (correct==0){
    a1=rbinom(N,size=1,prob=0.5)
    r=rep(0,N); r[a1==1]=rbinom(sum(a1),size=1,prob=p1); r[a1==0]=rbinom(N-sum(a1),size=1,prob=p2)
    a2=rep(0,N); a2[r==0]=2*rbinom(N-sum(r),size=1,prob=0.5)-1
    a1=2*a1-1
    label=1+(1-a1)*3/2+(3*a2-1)*a2/2
    correct=as.numeric(length(table(label))==6)
  }
  within=matrix(0,N,4)
  within[,1]=as.numeric((a1==1)&(a2>=0)); within[,2]=as.numeric((a1==1)&(a2<=0))
  within[,3]=as.numeric((a1==-1)&(a2>=0)); within[,4]=as.numeric((a1==-1)&(a2<=0))
  
  
  M=as.integer(runif(N,Mmin,Mmax+1-1e-8))
  
  #Building covariates
  if (p>0){
    x=list()
    for (i in 1:p){
      if (detail_x[i]==0){
        x[[i]]=matrix(rep(rnorm(N,mean=0,sd=sqrt(sigma2_x[i])),Mmax),N,Mmax)
      }else{
        x[[i]]=matrix(rnorm(N*Mmax,mean=0,sd=sqrt(sigma2_x[i])),N,Mmax)
      }
    }
  }
  
  #Building residuals
  epsilon=matrix(0,N,Mmax); mu=matrix(0,N,Mmax)
  for (i in 1:N){
    epsilon[i,]=mvrnorm(1,mu=rep(0,Mmax),Sigma=conditional_paras[label[i],2]*CS(conditional_paras[label[i],3],Mmax))
    mu[i,]=rep(conditional_paras[label[i],1],Mmax)
  }
  
  #merging all data to generate the outcome
  y=mu+epsilon
  if (p>0){
    for (i in 1:p){
      y=y+eta_x[i]*x[[i]]
    }
  }
  
  #changing into long format data
  result=matrix(0,sum(M),6+p)
  colnames(result)=c(1:(6+p))
  colnames(result)[1:6]=c("Y","id","L","A1","A2","R")
  if (p>1) colnames(result)[7:(6+p)]=paste("X",c(1:p),sep="")
  if (p==1) colnames(result)[7]="X"
  result=data.frame(result)
  count=0
  for (i in 1:N){
    for (j in 1:M[i]){
      count=count+1
      result[count,1]=y[i,j]
      result[count,2]=i
      result[count,3]=label[i]
      result[count,4]=ifelse(label[i] %in% c(1,2,3),1,-1)
      result[count,5]=ifelse(label[i] %in% c(1,4),NA,ifelse(label[i] %in% c(2,5),1,-1))
      result[count,6]=as.numeric(is.na(result[count,5]))
      if (p>0){
        for (k in 1:p) result[count,6+k]=x[[k]][i,j]
      }
    }
  }
  if (p>0){
    data=list(Y=result$Y,id=result$id,L=result$L,A1=result$A1,A2=result$A2,R=result$R,X=result[,7:(6+p)])
  }else{
    data=list(Y=result$Y,id=result$id,L=result$L,A1=result$A1,A2=result$A2,R=result$R,X=NULL)
  }
  
  return(list(data=data,marginal_paras=marginal_paras,beta=beta))
}


#changing the formula from response variables (A1,R,A2) to the label of treatment pathway L
response_to_pathway=function(A1,R,A2){
  n=length(A1)
  if (length(R)!=n) stop("Error: Variables array not of same length")
  if (length(A2)!=n) stop("Error: Variables array not of same length")
  A2=ifelse(is.na(A2),0,A2)
  if (!(all.equal(unique(A2[which(R==1)]),0)==TRUE)) stop("Error: Randomizing responders of first-stage treatment")
  label=1+(1-A1)*3/2+(3*A2-1)*A2/2
  return(label)
}

#changing the formula from treatment pathway label L into (A1,R,A2)
pathway_to_response=function(L){
  if (FALSE %in% (unique(L) %in% c(1,2,3,4,5,6))) stop("Error: value of L should be in 1,2,...,6")
  A1=ifelse(L>=4,-1,1)
  R=ifelse(L %in% c(1,4),1,0)
  A2=ifelse(L %in% c(1,4),NA,ifelse(L %in% c(2,5),1,-1))
  result=cbind(A1,R,A2)
  return(result)
}


#Expanding the data to use in geeglm (if we only care about effect estimation and not about inference)
#seperate: Logical. If TRUE then each response cluster is seperated into two clusters, if FALSE the response cluster is expanded to be twice large. Default is FALSE.
data_expand=function(Y,X=NULL,cluster_id,A1,R,A2,seperate=FALSE){
  id=cluster_id
  N=length(unique(id))
  M=as.numeric(table(id))
  id_name=names(table(id))
  
  if (all.equal(sort(unique(A1)),c(-1,1))!=TRUE) stop("A1 must be in -1 and 1")
  if (all.equal(sort(unique(R)),c(0,1))!=TRUE) stop("R must be in 0 and 1")
  if (length(unique(A2))!=3) stop("A2 must be in -1, 1 and NA")
  if (length(sort(unique(A2)))!=2) stop("A2 must be in -1, 1 and NA")
  if (all.equal(sort(unique(A2)),c(-1,1))!=TRUE) stop("A2 must be in -1, 1 and NA")
  if (length(unique(A2[which(R==1)]))>1) stop("The responders cannot be randomized in prototypical SMART")
  if (!is.na(unique(A2[which(R==1)]))) stop("The responders cannot be randomized in prototypical SMART")
  if (length(unique(A2[which(R==0)]))!=2) stop("The non-responders' A2 should be in -1 and 1")
  if (all.equal(sort(unique(A2[which(R==0)])),c(-1,1))!=TRUE) stop("The non-responders' A2 should be in -1 and 1")
  L=ifelse(A1==1,ifelse(R==1,1,ifelse(A2==1,2,3)),ifelse(R==1,4,ifelse(A2==1,5,6)))
  print("Warning:  The analysis using geeglm is only valid when all clusters have same size or when using indepdendence working variance structure")
  
  if (is.null(X)){
    data=cbind(Y,id,L,A1,R,A2)
  }else{
    data=cbind(Y,id,L,A1,R,A2,X)
  }
  data=data.frame(data)
  
  ndata=nrow(data)
  
  if (seperate){
    maxid=max(data$id)
    maxid=10^trunc(log(maxid,base=10)+1)
    for (i in 1:ndata){
      if (is.na(data$A2[i])){
        new_row=data[i,]
        data$A2[i]=1
        new_row$A2=-1
        new_row$id=new_row$id+maxid
        data=rbind(data,new_row)
      }
    }
    data$W=ifelse(data$L %in% c(1,4),2,4)
    return(list(data=data,structure=NULL))
  }
  
  for (i in 1:ndata){
    if (is.na(data$A2[i])){
      new_row=data[i,]
      data$A2[i]=1
      new_row$A2=-1
      new_row$id=new_row$id
      data=rbind(data,new_row)
    }
  }
  data$W=ifelse(data$L %in% c(1,4),2,4)
  data=data[order(data$id),]
  
  if (length(unique(M))==1){
    m=nrow(data)
    wave=rep(0,m)
    wave[1]=1
    for (i in 2:m){
      if (data$id[i]==data$id[i-1]){
        wave[i]=wave[i-1]+1
      }else{
        wave[i]=1
      }
    }
    
    M=M[1]
    zcor1=genZcor(clusz=table(data$id),waves=wave,corstrv=2*M)
    zcor.toep1=matrix(0,nrow(zcor1),1)
    count=0
    for (i in 1:(2*M-1)){
      for (j in (i+1):(2*M)){
        count=count+1
        if ((j<=M)|(i>M)) zcor.toep1[,1]=zcor.toep1[,1]+zcor1[,count]
      }
    }
    return(list(data=data,structure=zcor.toep1))
  }else{
    return(list(data=data,structure=NULL))
  }
}


#Input Format: Y:outcome; X: covariates; id:cluster id; L:treatment pathway
#Tuning Parameters:
#  estimate_weight: 0 for not estimating IPW, 1 for estimating IPW and do the corresponding adjustment (default = 0)
#  variance_structure: 0 for same variance across Adaptive Intervention, 1 for different variance (default = 1)
#  correlation_structure: 0 for independence, 1 for same ICC across AI, 2 for different ICC (default = 2)
#  ICC_lower_thresh: the manual lower bound for the ICC, set to 0 to force ICC to be nonnegative, to -1 for not adjusting. Should be betwen -1 and 1 (default = 0)
#  max_iter: the maximum allowed number of iterations in solving the estimating equation
#  dof_adjustment,use_t,bias_correction are all proposed finite-sample adjustments (default = 1)
#  verbose: 0 for not outputing results, 1 for outputing
#  aimed_comparison: r*4 matrix, each row has 2 0's and 2 1's, corresponding to the proposed comparison between AIs
solve_SMART=function(Y,X=NULL,cluster_id,A1,R,A2,
                     aimed_comparison=matrix(c(1,1,1,0,0,0,1,0,0,1,1,0,0,1,0,1,0,1,0,0,1,0,1,1),6,4),
                     variance_structure=1,correlation_structure=2,ICC_lower_thresh=0,
                     max_iter=10,convergence_thresh=1e-5,alpha=0.05,
                     estimate_weight=FALSE,dof_adjustment=TRUE,use_t=TRUE,bias_correction=TRUE,verbose=2){
  coeff=matrix(c(1,1,1,1,1,1,-1,-1,1,-1,1,-1,1,-1,-1,1),4,4)
  AI_names=c("AI(1,1)","AI(1,-1)","AI(-1,1)","AI(-1,-1)")
  
  #need to check imput validity
  if ((ICC_lower_thresh<(-1))|(ICC_lower_thresh>1)) stop("ICC_lower_thresh should be between -1 and 1")
  
  id=cluster_id
  N=length(unique(id))
  id_name=names(table(id))
  Mmax=max(table(id))
  M=as.numeric(table(id))
  
  if (all.equal(sort(unique(A1)),c(-1,1))!=TRUE) stop("A1 must be in -1 and 1")
  if (all.equal(sort(unique(R)),c(0,1))!=TRUE) stop("R must be in 0 and 1")
  if (length(unique(A2))!=3) stop("A2 must be in -1, 1 and NA")
  if (length(sort(unique(A2)))!=2) stop("A2 must be in -1, 1 and NA")
  if (all.equal(sort(unique(A2)),c(-1,1))!=TRUE) stop("A2 must be in -1, 1 and NA")
  if (length(unique(A2[which(R==1)]))>1) stop("The responders cannot be randomized in prototypical SMART")
  if (!is.na(unique(A2[which(R==1)]))) stop("The responders cannot be randomized in prototypical SMART")
  if (length(unique(A2[which(R==0)]))!=2) stop("The non-responders' A2 should be in -1 and 1")
  if (all.equal(sort(unique(A2[which(R==0)])),c(-1,1))!=TRUE) stop("The non-responders' A2 should be in -1 and 1")
  L=ifelse(A1==1,ifelse(R==1,1,ifelse(A2==1,2,3)),ifelse(R==1,4,ifelse(A2==1,5,6)))
  
  inside_list=list()
  for (i in 1:N) inside_list[[i]]=which(id==id_name[i])
  label=rep(0,N)
  for (i in 1:N){
    ll=L[inside_list[[i]]]
    if (length(unique(ll))!=1){
      stop("id is not nested in treatment pathway")
    }else{
      label[i]=unique(ll)
    }
  }
  if (length(unique(label))<6) stop("not every treatment pathway has observations, postitivity assumption violated")
  a1=ifelse(label>=4,-1,1)
  a2=ifelse(label %in% c(1,4),0,ifelse(label %in% c(2,5),1,-1))
  within=matrix(0,N,4)
  within[,1]=as.numeric((a1==1)&(a2>=0)); within[,2]=as.numeric((a1==1)&(a2<=0))
  within[,3]=as.numeric((a1==-1)&(a2>=0)); within[,4]=as.numeric((a1==-1)&(a2<=0))
  
  q=4
  if (is.null(X)){
    p=0
  }else{
    if (is.null(dim(X))){
      X=matrix(X,nrow=length(X),ncol=1)
      colnames(X)="X"
      p=1
    }else{
      p=ncol(X)
      if (!is.matrix(X)) X=as.matrix(X)
    }
  }
  
  if (estimate_weight==FALSE){
    count=as.numeric(table(label))
    p1_hat=count[1]/sum(count[1:3]); p2_hat=count[4]/sum(count[4:6])
    pa1_hat=1/2; p1a2_hat=1/2; p2a2_hat=1/2
    w=c(2,4,4,2,4,4)
  }else{
    count=as.numeric(table(label))
    pa1_hat=sum(count[1:3])/N
    p1_hat=count[1]/sum(count[1:3]); p2_hat=count[4]/sum(count[4:6])
    p1a2_hat=count[2]/(count[2]+count[3]); p2a2_hat=count[5]/(count[5]+count[6])
    w=c(1/pa1_hat,1/pa1_hat/p1a2_hat,1/pa1_hat/(1-p1a2_hat),1/(1-pa1_hat),1/(1-pa1_hat)/p2a2_hat,1/(1-pa1_hat)/(1-p2a2_hat))
  }
  
  
  #Use iterative method to calculate beta_hat and working covariance matrix (model based V_naive)
  sigma2_hat=rep(1,q); rho_hat=rep(0,q); beta_hat=rep(0,q+p); diff=1
  count_iter=0
  while ((diff>convergence_thresh)&(count_iter<max_iter)){
    count_iter=count_iter+1
    V_naive=matrix(0,q+p,q+p)
    weighty=matrix(0,q+p,1)
    for (i in 1:N){
      yi=as.matrix(Y[inside_list[[i]]],M[i],1)
      for (group in 1:q){
        if (within[i,group]==1){
          Di=matrix(rep(c(coeff[group,]),each=M[i]),M[i],q)
          if (p>0) Di=cbind(Di,matrix(X[inside_list[[i]],],M[i],p))
          Vi=sigma2_hat[group]*CS(rho_hat[group],M[i])
          V_naive=V_naive+w[label[i]]*t(Di)%*%solve(Vi)%*%Di
          weighty=weighty+w[label[i]]*t(Di)%*%solve(Vi)%*%yi
        }
      }
    }
    V_naive=solve(V_naive)
    beta_new=V_naive%*%weighty
    diff=sum(abs(beta_new-beta_hat))
    beta_hat=beta_new
    
    sigma_group=rep(0,q); weight_sigma=rep(0,q)
    rho_group=rep(0,q); weight_rho=rep(0,q)
    
    for (group in 1:q){
      for (i in 1:N){
        if (within[i,group]==1){
          mu_hat=beta_hat[1]*coeff[group,1]+beta_hat[2]*coeff[group,2]+beta_hat[3]*coeff[group,3]+beta_hat[4]*coeff[group,4]
          if (p>1) mu_hat=mu_hat+X[inside_list[[i]],]%*%beta_hat[(q+1):(q+p),1]
          if (p==1) mu_hat=mu_hat+X[inside_list[[i]],]*beta_hat[q+1]
          epsilon_hat=as.matrix(Y[inside_list[[i]]],M[i],1)-mu_hat
          sigma_group[group]=sigma_group[group]+w[label[i]]*sum(epsilon_hat^2)
          weight_sigma[group]=weight_sigma[group]+w[label[i]]*M[i]
          if (correlation_structure>=1){
            rho_group[group]=rho_group[group]+w[label[i]]*(sum(epsilon_hat)^2-sum(epsilon_hat^2))
            weight_rho[group]=weight_rho[group]+w[label[i]]*M[i]*(M[i]-1)}
        }
      }
    }
    if (variance_structure==0){
      sigma2_hat=rep(sum(sigma_group)/sum(weight_sigma),q)}
    else{
      sigma2_hat=sigma_group/weight_sigma}
    weight_rho=weight_rho*sigma2_hat
    if (correlation_structure==0){rho_hat=rep(0,q)}
    
    if (correlation_structure==1){rho_hat=rep(max(ICC_lower_thresh,sum(rho_group)/sum(weight_rho)),q)}
    if (correlation_structure==2){
      rho_hat=rho_group/weight_rho
      for (group in 1:q) rho_hat[group]=max(rho_hat[group],ICC_lower_thresh)
    }
  }
  
  
  if(estimate_weight){
    S_mat=matrix(c(1/p1_hat,-1/(1-p1_hat),-1/(1-p1_hat),0,0,0,
                   0,0,0,1/p2_hat,-1/(1-p2_hat),-1/(1-p2_hat),
                   1/pa1_hat,1/pa1_hat,1/pa1_hat,-1/(1-pa1_hat),-1/(1-pa1_hat),-1/(1-pa1_hat),
                   0,1/p1a2_hat,-1/(1-p1a2_hat),0,0,0,
                   0,0,0,0,1/p2a2_hat,-1/(1-p2a2_hat))
                 ,5,6,byrow=T)
    W_der_mat=matrix(c(0,0,0,0,0,0,
                       0,0,0,0,0,0,
                       -1/(pa1_hat^2),-1/(pa1_hat^2*p1a2_hat),-1/(pa1_hat^2*(1-p1a2_hat)),1/((1-pa1_hat)^2),1/((1-pa1_hat)^2*p2a2_hat),1/((1-pa1_hat)^2*(1-p2a2_hat)),
                       0,-1/(pa1_hat*p1a2_hat^2),1/(pa1_hat*(1-p1a2_hat)^2),0,0,0,
                       0,0,0,0,-1/((1-pa1_hat)*p2a2_hat^2),1/((1-pa1_hat)*(1-p2a2_hat)^2))
                     ,5,6,byrow=T)
    
    B_mat=matrix(0,5,5)
    for (i in 1:N){
      B_mat=B_mat+S_mat[,label[i]]%*%t(S_mat[,label[i]])
    }
  }
  
  
  if (!bias_correction){
    #Calculating Liang&Zeger Robust Sandwich Estimator
    meat_lz=matrix(0,p+q,p+q)
    if (estimate_weight) C_mat=matrix(0,p+q,5)
    for (i in 1:N){
      yi=as.matrix(Y[inside_list[[i]]],M[i],1)
      Ui=matrix(0,p+q,1)
      for (group in 1:q){
        if (within[i,group]==1){
          Di=matrix(rep(c(coeff[group,]),each=M[i]),M[i],q)
          if (p>0) Di=cbind(Di,matrix(X[inside_list[[i]],],M[i],p))
          mu_hat=beta_hat[1]*coeff[group,1]+beta_hat[2]*coeff[group,2]+beta_hat[3]*coeff[group,3]+beta_hat[4]*coeff[group,4]
          if (p>1) mu_hat=mu_hat+X[inside_list[[i]],]%*%beta_hat[(q+1):(q+p),1]
          if (p==1) mu_hat=mu_hat+X[inside_list[[i]],]*beta_hat[q+1]
          Vi=sigma2_hat[group]*CS(rho_hat[group],M[i])
          Ui=Ui+w[label[i]]*t(Di)%*%solve(Vi)%*%(yi-mu_hat)
          if (estimate_weight) C_mat=C_mat+t(Di)%*%solve(Vi)%*%(yi-mu_hat)%*%t(W_der_mat[,label[i]])
        }
      }
      meat_lz=meat_lz+Ui%*%t(Ui)
    }
    
    if (estimate_weight){
      V_lz=V_naive%*%(meat_lz+C_mat%*%solve(B_mat)%*%t(C_mat))%*%V_naive}
    else{
      V_lz=V_naive%*%meat_lz%*%V_naive}
    V_result=V_lz
  }else{
    #Calculating Mancl&Derhoen Bias Corrected Estimator
    meat_md=matrix(0,p+q,p+q)
    if (estimate_weight) C_mat=matrix(0,p+q,5)
    for (i in 1:N){
      yi=as.matrix(Y[inside_list[[i]]],M[i],1)
      Ui=matrix(0,p+q,1)
      for (group in 1:q){
        if (within[i,group]==1){
          Di=matrix(rep(c(coeff[group,]),each=M[i]),M[i],q)
          if (p>0) Di=cbind(Di,matrix(X[inside_list[[i]],],M[i],p))
          mu_hat=beta_hat[1]*coeff[group,1]+beta_hat[2]*coeff[group,2]+beta_hat[3]*coeff[group,3]+beta_hat[4]*coeff[group,4]
          if (p>1)mu_hat=mu_hat+X[inside_list[[i]],]%*%beta_hat[(q+1):(q+p),1]
          if (p==1) mu_hat=mu_hat+X[inside_list[[i]],]*beta_hat[q+1]
          Vi=sigma2_hat[group]*CS(rho_hat[group],M[i])
          Hi=Di%*%V_naive%*%t(Di)%*%solve(Vi)
          Ui=Ui+w[label[i]]*t(Di)%*%solve(Vi)%*%solve(diag(1,M[i],M[i])-Hi)%*%(yi-mu_hat)
          if (estimate_weight) C_mat=C_mat+t(Di)%*%solve(Vi)%*%(yi-mu_hat)%*%t(W_der_mat[,label[i]])
        }
      }
      meat_md=meat_md+Ui%*%t(Ui)
    }
    if (estimate_weight){
      V_md=V_naive%*%(meat_md+C_mat%*%solve(B_mat)%*%t(C_mat))%*%V_naive}
    else{
      V_md=V_naive%*%meat_md%*%V_naive}
    V_result=V_md
  }
  
  if (dof_adjustment) V_result=V_result*N/(N-p-q)
  
  
  
  
  
  #--------------------------------output session------------------------
  #output model structure
  formu="Y ~ a1 + a2 + I(a1 * a2)"
  if (p>0){
    for (i in 1:p) formu=paste(formu," + ",colnames(X)[i],sep="")
  }
  if (verbose>=2){
    cat("Marginal Mean Model: ",formu,"\n")
    cat(N,"Number of Clusters, Minimum Cluster Size =",min(M),", Maximum Cluster Size =",max(M),"\n")
    cat("Algorithm stops after ",count_iter," iterations.\n")
  }
  formu=as.formula(formu)
  
  #construct parameter output
  summary_paras=matrix(0,nrow=p+q,ncol=7)
  colnames(summary_paras)=c("Parameter","Estimate","Std.Err","CI Lower","CI Higher","Z Score","p-value")
  summary_paras=data.frame(summary_paras)
  summary_paras[1:4,1]=c("(Intercept)","a1","a2","I(a1*a2)")
  if (p>0) summary_paras[5:(4+p),1]=colnames(X)
  summary_paras[,2]=beta_hat
  summary_paras[,3]=sqrt(diag(V_result))
  summary_paras[,6]=summary_paras[,2]/summary_paras[,3]
  if (!use_t){
    summary_paras[,4]=beta_hat-qnorm(1-alpha/2)*sqrt(diag(V_result))
    summary_paras[,5]=beta_hat+qnorm(1-alpha/2)*sqrt(diag(V_result))
    summary_paras[,7]=2*(1-pnorm(abs(summary_paras[,6])))
  }else{
    colnames(summary_paras)[6]=c("T Score")
    summary_paras[,4]=beta_hat-qt(1-alpha/2,df=N-p-q)*sqrt(diag(V_result))
    summary_paras[,5]=beta_hat+qt(1-alpha/2,df=N-p-q)*sqrt(diag(V_result))
    summary_paras[,7]=2*(1-pt(abs(summary_paras[,6]),df=N-p-q))
  }
  
  #construct aimed comparison
  if (is.null(aimed_comparison)){
    r=0
  }else{
    if (is.null(dim(aimed_comparison))){
      aimed_comparison=matrix(aimed_comparison,nrow=1,ncol=4)
      r=1
    }else{
      r=nrow(aimed_comparison)
      if (!is.matrix(aimed_comparison)) aimed_comparison=as.matrix(aimed_comparison)
    }
    
    for (rr in 1:r){
      if (all.equal(aimed_comparison[order(aimed_comparison[rr,])],c(0,0,1,1))==FALSE){
        stop("Formula Error: aimed_comparison should be a length 4 vector with two 0's and two 1's")
      }
      num_1=which(aimed_comparison[rr,]==1)[1]; num_2=which(aimed_comparison[rr,]==1)[2]  #which two AIs are being compared
      weight_beta=coeff[num_1,]-coeff[num_2,]
      
      sum_test=sum(weight_beta*beta_hat[1:4])
      var_test=t(matrix(weight_beta))%*%V_result[1:4,1:4]%*%matrix(weight_beta)
      sd_test=sqrt(var_test)
      value_test=sum_test/sd_test
      if (!use_t){
        ci_low=sum_test-qnorm(1-alpha/2)*sd_test
        ci_high=sum_test+qnorm(1-alpha/2)*sd_test
        p_test=2*(1-pnorm(abs(value_test)))
      }else{
        ci_low=sum_test-qt(1-alpha/2,df=N-p-q)*sd_test
        ci_high=sum_test+qt(1-alpha/2,df=N-p-q)*sd_test
        p_test=2*(1-pt(abs(value_test),df=N-p-q))
      }
      addon=c(paste(AI_names[num_1],"-",AI_names[num_2],sep=""),sum_test,sd_test,ci_low,ci_high,value_test,p_test)
      summary_paras=rbind(summary_paras,addon)
    }
  }
  
  #formatting the parameter output
  output_paras=data.frame(cbind(summary_paras,rep(0,nrow(summary_paras))))
  output_paras[,2]=round(as.numeric(output_paras[,2]),digits=5)
  output_paras[,3]=round(as.numeric(output_paras[,3]),digits=5)
  output_paras[,4]=round(as.numeric(output_paras[,4]),digits=5)
  output_paras[,5]=round(as.numeric(output_paras[,5]),digits=5)
  output_paras[,6]=round(as.numeric(output_paras[,6]),digits=2)
  for (i in 1:nrow(summary_paras)){
    t=as.numeric(output_paras[i,7])
    output_paras[i,8]=ifelse(t>0.1," ",ifelse(t>0.05,".",ifelse(t>0.01,"*",ifelse(t>0.001,"**","***"))))
    output_paras[i,7]=ifelse(t<2e-16,"<2e-16",ifelse(t>0.001,as.character(round(t,digits=4)),formatC(t,format="e",digits=2)))
  }
  output_paras=rbind(colnames(output_paras),output_paras)
  output_paras[1,8]=""
  output_paras[1,7]=ifelse(!use_t,"Pr(>|z|)","Pr(>|t|)")
  colnames(output_paras)=NULL
  rownames(output_paras)=NULL
  
  #outputting the model coefficients
  if (verbose>=1){
    cat("Summary of model Coefficients:\n")
    print(output_paras[1:(p+q+1),],row.names=FALSE)
    cat("---\n")
    cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
  }
  
  
  
  marginal_paras=matrix(0,4,3)
  colnames(marginal_paras)=c("Mean","Variance","ICC")
  rownames(marginal_paras)=AI_names
  marginal_paras[,1]=coeff%*%matrix(as.numeric(summary_paras[1:4,2]),4,1)
  marginal_paras[,2]=sigma2_hat
  marginal_paras[,3]=rho_hat
  
  if(verbose>=2){
    cat("\n")
    cat("Working Variance Structure: ")
    if (correlation_structure==0){cat("Independence")}else{cat("Exchangeable")}
    cat("\nDetailed Structure: ")
    if (variance_structure==1){cat("sigma2_{a1,a2} * ")} else{cat("sigma2 * ")}
    if (correlation_structure==0) cat("I_{M[i]}")
    if (correlation_structure==1) cat("EXCH(rho,M[i])")
    if (correlation_structure==2) cat("EXCH(rho_{a1,a2},M[i])")
    cat("\n")
  }
  
  if (verbose==2){
    cat("Marginal Parameters:\n")
    if ((variance_structure==0)&(correlation_structure<=1)){
      if (correlation_structure==0) cat("sigma2 =",marginal_paras[1,2])
      if (correlation_structure==1) cat("sigma2 =",marginal_paras[1,2]," , rho =",marginal_paras[1,3])
    }else{
      if (correlation_structure==0) print(marginal_paras[,2])
      if (correlation_structure>0) print(marginal_paras[,-1])
    }
    cat("\n")
  }
  
  
  if ((verbose>=3)&(p+r>0)){
    cat("\nWarning:\n")
    cat("The following marginal mean and inference are valid only when:\n")
    cat("    The model only contains baseline covariates and their interaction with A1,A2\n")
    cat("    The baseline covaraites should have mean zero\n")
  }
  
  if (verbose>=3){
    cat("\nMarginal Parameters:\n")
    if (correlation_structure==0) print(marginal_paras[,-3]) else print(marginal_paras)
    cat("\n")
  }
  
  if (r>0){
    summary_comparison=summary_paras[(p+q+1):(p+q+r),]
    colnames(summary_comparison)[1]="Comparison"
    rownames(summary_comparison)=c()
    summary_paras=summary_paras[1:(p+q),]
    
    if (verbose>=3){
      output_comparison=output_paras[c(1,(p+q+2):(p+q+r+1)),]
      output_comparison[1,1]="Comparison"
      if (r==1){
        cat("Result of Proposed Comparison between ",AI_names[num_1]," and ",AI_names[num_2],"\n")
      }else{
        cat("Result of Proposed Comparison between Adaptive Interventions:\n")
      }
      print(output_comparison,row.names=FALSE)
      cat("---\n")
      cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
      cat("\n")
    }
  }
  
  
  
  
  rownames(V_result)=summary_paras[1:(p+q),1]
  colnames(V_result)=summary_paras[1:(p+q),1]
  
  if (verbose>=3){
    cat("\n")
    cat("Variance-Covariance matrix of the estimates\n")
    print(V_result)
  }
  
  for (i in 2:7){
    summary_paras[,i]=as.numeric(summary_paras[,i])
  }
  
  if (!is.null(aimed_comparison)){
    for (i in 2:7) {summary_comparison[,i]=as.numeric(summary_comparison[,i])}
    return(list(summary_paras=summary_paras,marginal_paras=marginal_paras,var_estimator=V_result,
                formula=formu,comparison=summary_comparison,N=N))
  }else{
    return(list(summary_paras=summary_paras,marginal_paras=marginal_paras,var_estimator=V_result,
                formula=formu,comparison=NULL,N=N))
  }
}



#aimed test is an r*(p+q+1) matrix, each row corresponding to the combination of coefficients
hypothesis_testing=function(result,aimed_test,alpha=0.05, q=4, use_t=TRUE){
  p=nrow(result$summary_paras)-q
  if (is.null(dim(aimed_test))){
    aimed_test=matrix(aimed_test,nrow=1)
    r=1
  }else{
    r=nrow(aimed_test)
    if (!is.matrix(aimed_test)) aimed_test=as.matrix(aimed_test)
  }
  if (ncol(aimed_test)!=(p+q+1)) stop("Error: aimed_test should have q+p+1 columns")
  
  summary_test=matrix(0,r,7)
  colnames(summary_test)=c("Hypothesis_Estimand","Estimate","Std.Err","CI Lower","CI Higher","Test Score","p-value")
  summary_test=data.frame(summary_test)
  variable_names=result$summary_paras[,1]
  for (rr in 1:r){
    weight_beta=aimed_test[rr,1:(p+q)]
    if (max(as.numeric(weight_beta==0))==0) stop(paste("Error: row ",rr," contains only zero",sep=""))
    
    first_place=which.min(weight_beta==0)
    H0=""
    t=weight_beta[first_place]
    if (t>0) {
      if (t==1){
        H0=variable_names[first_place]
      }else{
        H0=paste(round(t,digits=2)," * ",variable_names[first_place],sep="")
      }
    }
    if (t<0){
      if (t==-1){
        H0=paste("- ",variable_names[first_place],sep="")
      }else{
        H0=paste("- ",round(-t,digits=2)," * ",variable_names[first_place],sep="")
      }
    }
    
    for (i in (first_place+1):(p+q)){
      t=weight_beta[i]
      if (t>0) {
        if (t==1){
          H0=paste(H0," + ",variable_names[i],sep="")
        }else{
          H0=paste(H0," + ",round(t,digits=2)," * ",variable_names[i],sep="")
        }
      }
      if (t<0){
        if (t==-1){
          H0=paste(H0," - ",variable_names[i],sep="")
        }else{
          H0=paste(H0," - ",round(-t,digits=2)," * ",variable_names[i],sep="")
        }
      }
    }
    
    t=aimed_test[rr,p+q+1]
    if (t>0){H0=paste(H0," + ",round(t,digits=2),sep="")}
    if (t<0){H0=paste(H0," - ",round(-t,digits=2),sep="")}
    
    #H0=paste(H0," = 0",sep="")
    summary_test[rr,1]=H0
    
    summary_test[rr,2]=sum(weight_beta*as.numeric(result$summary_paras[,2]))+aimed_test[rr,p+q+1]
    summary_test[rr,3]=sqrt(t(matrix(weight_beta))%*%result$var_estimator%*%matrix(weight_beta))
    summary_test[rr,6]=summary_test[rr,2]/summary_test[rr,3]
    if (!use_t){
      summary_test[rr,4]=summary_test[rr,2]-qnorm(1-alpha/2)*summary_test[rr,3]
      summary_test[rr,5]=summary_test[rr,2]+qnorm(1-alpha/2)*summary_test[rr,3]
      summary_test[rr,7]=2*(1-pnorm(abs(summary_test[rr,4])))
    }else{
      summary_test[rr,4]=summary_test[rr,2]-qt(1-alpha/2,df=result$N-p-q)*summary_test[rr,3]
      summary_test[rr,5]=summary_test[rr,2]+qt(1-alpha/2,df=result$N-p-q)*summary_test[rr,3]
      summary_test[rr,7]=2*(1-pt(abs(summary_test[rr,6]),df=result$N-p-q))
    }
  }
  return(summary_test)
}




#Generating correlation matrix that has multi-layer structure
#  M: number of people within a cluster (number of diagonal blocks)
#  T: number of observations per people (size of diagonal blocks)
#  diagonal_structure: be in "AR1", "Independence, "Exchangeable" or "Unstructured"
#  off_diagonal_structure: be in "Exchangeable" or "Unstructured"
#  para_diagonal: in "AR1" or "Exchangeable" structure, be the parameter; in "Unstructured", be the whole matrix
#  para__off_diagonal: in "Exchangeable" structure, be the parameter; in "Unstructured", be the whole matrix
Corr_Multilayer=function(M,max_T,diagonal_structure="AR1",para_diagonal,
                         off_diagonal_structure="Exchangeable",para_off_diagonal){
  if (diagonal_structure == "Independence"){
    diagonal_structure="Exchangeable"
    para_diagonal=0
  }
  if (off_diagonal_structure == "Independence"){
    off_diagonal_structure="Exchangeable"
    para_off_diagonal=0
  }
  
  if (! (diagonal_structure %in% c("AR1","Exchangeable","Unstructured"))) stop("diagonal_structure should be in AR1, Independence, Exchangeable or Unstructured")
  if (! (off_diagonal_structure %in% c("Exchangeable","Unstructured"))) stop("off_diagonal_structure should be in Independence, Exchangeable or Unstructured")
  
  result_basis=matrix(0,max_T,max_T)
  if (off_diagonal_structure=="Exchangeable"){
    result_basis=matrix(rep(para_off_diagonal,max_T*max_T),max_T,max_T)
  }else{
    result_basis=para_off_diagonal
  }
  
  result=matrix(rep(t(result_basis),M),nrow=max_T*M,ncol=max_T,byrow=TRUE)
  result=matrix(rep(result,M),nrow=max_T*M,ncol=max_T*M)
  
  if (diagonal_structure=="AR1"){
    for (i in 1:M){
      result[((i-1)*max_T+1):(i*max_T),((i-1)*max_T+1):(i*max_T)]=AR1(para_diagonal,max_T)
    }
  }
  if (diagonal_structure=="Exchangeable"){
    for (i in 1:M){
      result[((i-1)*max_T+1):(i*max_T),((i-1)*max_T+1):(i*max_T)]=CS(para_diagonal,max_T)
    }
  }
  if (diagonal_structure=="Unstructured"){
    for (i in 1:M){
      result[((i-1)*max_T+1):(i*max_T),((i-1)*max_T+1):(i*max_T)]=para_diagonal
    }
  }
  return(result)
}
#Corr_Multilayer(4,3,diagonal_structure="AR1",off_diagonal_structure="Exchangeable",para_diagonal=0.2,para_off_diagonal=0.3)



#Solving multi-layer SMART
#N: number of clusters; M: length N array of number of people within a cluster; max_T: number of observation per people; D:number of AIs
#Ind: D*L matrix indicating the value of indicating variables (each column represents an indicating variable, like 1, a1, a2, a1*a2 for prototypical SMART)
#Y: length max(M)*max_T vector, represent individual outcome
#X: a list of length L, each includes a matrix with sum(M)*max(T) rows, each row represent the covariates terms that interacts with the corresponding indicating variables
#mean_covariate: NULL, or a vector of length sum(ncol(X)), indicating the expectation of each covariate when calculating the marginal mean of each AI (should be 1 if the covariate is indicator alone and 0 if mean-zero baseline)
#within: N*D matrix indicating whether a cluster is consistent with an AI
#weight: estimated weight for a given cluster
#var_homo_across_time:  the variance sigma^2(d)_t is homogeneous across time t
#var_homo_across_AI:  the variance sigma^2(d)_t is homogeneous across intervention d
#diagonal_homo_across_AI:  the parameters in the correlation structure is homogeneous across AIs
#correlation_normalize: 0 for the old version, 1 for the new version that forces to use un-homogeneous variance
solve_SMART_Multilayer=function(N,M,max_T,D,
                                Ind,Y,X,
                                #mean_covariate=NULL,
                                within,weight,
                                var_homo_across_time=FALSE,var_homo_across_AI=FALSE,
                                diagonal_structure="AR1",diagonal_homo_across_AI=FALSE,diagonal_ICC_lower_thresh=0,
                                off_diagonal_structure="Exchangeable",off_diagonal_homo_across_AI=FALSE,off_diagonal_ICC_lower_thresh=0,
                                max_iter=10,dof_adjustment=0,use_t=0,bias_correction=0,correlation_normalize=1,verbose=1){
  #need to check imput validity
  if ((diagonal_ICC_lower_thresh<(-1))|(diagonal_ICC_lower_thresh>1)) stop("diagonal_ICC_lower_thresh should be between -1 and 1")
  if ((off_diagonal_ICC_lower_thresh<(-1))|(off_diagonal_ICC_lower_thresh>1)) stop("off_diagonal_ICC_lower_thresh should be between -1 and 1")
  if (! (diagonal_structure %in% c("AR1","Independence","Exchangeable","Unstructured"))) stop("diagonal_structure should be in AR1, Independence, Exchangeable or Unstructured")
  if (! (off_diagonal_structure %in% c("Independence","Exchangeable","Unstructured"))) stop("off_diagonal_structure should be in Independence, Exchangeable or Unstructured")
  
  L=ncol(Ind)
  Ind_names=colnames(Ind)
  p_each=rep(0,L)
  for (l in 1:L){
    if (!is.null(X[[l]])){
      if (is.null(dim(X[[l]]))){
        X[[l]]=matrix(X[[l]],nrow=length(X[[l]]),ncol=1)
        colnames(X[[l]])="1"
        p_each[l]=1
      }else{
        p_each[l]=ncol(X[[l]])
        if (!is.matrix(X[[l]])) X[[l]]=as.matrix(X[[l]])
      }
    }
  }
  p=sum(p_each)
  starting_ind=c(0,cumsum(p_each)[1:(L-1)])+1
  ending_ind=cumsum(p_each)
  
  starting_cluster=c(0,cumsum(M)[1:(N-1)])*max_T+1
  ending_cluster=cumsum(M)*max_T
  
  
  #Use iterative method to calculate beta_hat and working covariance matrix (model based V_naive)
  sigma2_hat=matrix(1,nrow=max_T,ncol=p)
  if (diagonal_structure %in% c("AR1","Independence","Exchangeable")){
    diagonal_rho_hat=rep(0,p)
  }else{
    diagonal_rho_hat=rep(list(diag(1,max_T)),p)
  }
  if (off_diagonal_structure %in% c("Independence","Exchangeable")){
    off_diagonal_rho_hat=rep(0,p)
  }else{
    off_diagonal_rho_hat=rep(list(matrix(0,max_T,max_T)),p)
  }
  
  
  
  beta_hat=rep(0,p); diff=1; count_iter=0
  while ((diff>1e-5)&(count_iter<max_iter)){
    count_iter=count_iter+1
    
    V_naive=matrix(0,p,p)
    weighty=matrix(0,p,1)
    
    for (i in 1:N){
      yi=as.matrix(Y[starting_cluster[i]:ending_cluster[i]],M[i]*max_T,1)
      for (group in 1:D){
        if (within[i,group]==1){
          Di=matrix(0,M[i]*max_T,0)
          for (l in 1:L){
            if (p_each[l]>0) Di=cbind(Di,Ind[group,l]*X[[l]][starting_cluster[i]:ending_cluster[i],])
          }
          
          Vi=diag(rep(sqrt(sigma2_hat[,group]),M[i]))%*%
            Corr_Multilayer(M[i],max_T,diagonal_structure=diagonal_structure,off_diagonal_structure=off_diagonal_structure,para_diagonal=diagonal_rho_hat[[group]],para_off_diagonal=off_diagonal_rho_hat[[group]])%*%
            diag(rep(sqrt(sigma2_hat[,group]),M[i]))
          V_naive=V_naive+weight[i]*t(Di)%*%solve(Vi)%*%Di
          weighty=weighty+weight[i]*t(Di)%*%solve(Vi)%*%yi
        }
      }
    }
    V_naive=solve(V_naive)
    beta_new=V_naive%*%weighty
    diff=sum(abs(beta_new-beta_hat))
    beta_hat=beta_new
    
    sigma_group=matrix(0,max_T,D); weight_sigma=matrix(0,max_T,D)
    diagonal_rho_group=rep(list(matrix(0,max_T,max_T)),D); weight_diagonal_rho=rep(list(matrix(0,max_T,max_T)),D)
    off_diagonal_rho_group=rep(list(matrix(0,max_T,max_T)),D); weight_off_diagonal_rho=rep(list(matrix(0,max_T,max_T)),D)
    
    for (group in 1:D){
      for (i in 1:N){
        if (within[i,group]==1){
          
          mu_hat=matrix(0,M[i]*max_T,1)
          for (l in 1:L){
            if (p_each[l]>1) mu_hat=mu_hat+Ind[group,l]*X[[l]][starting_cluster[i]:ending_cluster[i],1:p_each[l]]%*%beta_hat[starting_ind[l]:ending_ind[l],1]
            if (p_each[l]==1) mu_hat=mu_hat+Ind[group,l]*X[[l]][starting_cluster[i]:ending_cluster[i],1:p_each[l]]*beta_hat[starting_ind[l]]
          }
          
          epsilon_hat=as.matrix(Y[starting_cluster[i]:ending_cluster[i]],M[i]*max_T,1)-mu_hat
          
          for (t in 1:max_T){
            sigma_group[t,group]=sigma_group[t,group]+weight[i]*sum(epsilon_hat[(1:M[i])*max_T-max_T+t]^2)
            weight_sigma[t,group]=weight_sigma[t,group]+weight[i]*M[i]
          }
          
          for (j in 1:max_T){
            for (k in 1:max_T){
              if (j != k){
                diagonal_rho_group[[group]][j,k]=diagonal_rho_group[[group]][j,k]+weight[i]*sum(epsilon_hat[(1:M[i])*max_T-max_T+j]*epsilon_hat[(1:M[i])*max_T-max_T+k])
                weight_diagonal_rho[[group]][j,k]=weight_diagonal_rho[[group]][j,k]+weight[i]*M[i]
              }
            }
          }
          
          for (j in 1:max_T){
            for (k in 1:max_T){
              off_diagonal_rho_group[[group]][j,k]=off_diagonal_rho_group[[group]][j,k]+
                weight[i]*(sum(epsilon_hat[(1:M[i])*max_T-max_T+j])*sum(epsilon_hat[(1:M[i])*max_T-max_T+k])-sum(epsilon_hat[(1:M[i])*max_T-max_T+j]*epsilon_hat[(1:M[i])*max_T-max_T+k]))
              weight_off_diagonal_rho[[group]][j,k]=weight_off_diagonal_rho[[group]][j,k]+weight[i]*M[i]*(M[i]-1)
            }
          }
        }
      }
    }
    
    
    #--------------calculating parameters of variances--------------------
    #using weighted average
    #if (var_homo_across_AI==TRUE){
    #  if (var_homo_across_time==TRUE){
    #    sigma2_hat=matrix(sum(sigma_group)/sum(weight_sigma),max_T,D)
    #  }else{
    #    sigma2_hat=matrix(rep(rowSums(sigma_group)/rowSums(weight_sigma),D),max_T,D)
    #  }
    #}
    #else{
    #  if (var_homo_across_time==TRUE){
    #    sigma2_hat=matrix(rep(colSums(sigma_group)/colSums(weight_sigma),each=max_T),max_T,D)
    #  }else{
    #    sigma2_hat=sigma_group/weight_sigma
    #  }
    #}
    
    
    if (correlation_normalize==0){
      #using direct average
      if (var_homo_across_AI==TRUE){
        if (var_homo_across_time==TRUE){
          sigma2_hat=matrix(mean(rowSums(sigma_group)/rowSums(weight_sigma)),max_T,D)
        }else{
          sigma2_hat=matrix(rep(rowSums(sigma_group)/rowSums(weight_sigma),D),max_T,D)
        }
      }else{
        if (var_homo_across_time==TRUE){
          sigma2_hat=matrix(rep(colMeans(sigma_group/weight_sigma),each=max_T),max_T,D)
        }else{
          sigma2_hat=sigma_group/weight_sigma
        }
      }}else{
        sigma2_hat=sigma_group/weight_sigma
      }
    #----------------------Adjust the correlation using derived variance-----------------
    
    
    for (group in 1:D){
      for (j in 1:max_T){
        for (k in 1:max_T){
          diagonal_rho_group[[group]][j,k]=diagonal_rho_group[[group]][j,k]/sqrt(sigma2_hat[j,group]*sigma2_hat[k,group])
          off_diagonal_rho_group[[group]][j,k]=off_diagonal_rho_group[[group]][j,k]/sqrt(sigma2_hat[j,group]*sigma2_hat[k,group])
        }
      }
    }
    
    if (correlation_normalize==1){
      #using direct average
      if (var_homo_across_AI==TRUE){
        if (var_homo_across_time==TRUE){
          sigma2_hat=matrix(mean(rowSums(sigma_group)/rowSums(weight_sigma)),max_T,D)
        }else{
          sigma2_hat=matrix(rep(rowSums(sigma_group)/rowSums(weight_sigma),D),max_T,D)
        }
      }else{
        if (var_homo_across_time==TRUE){
          sigma2_hat=matrix(rep(colMeans(sigma_group/weight_sigma),each=max_T),max_T,D)
        }else{
          sigma2_hat=sigma_group/weight_sigma
        }
      }
    }
    #-------------calculating parameters of diagonal blocks------------------
    if (diagonal_structure=="Independence"){
      diagonal_rho_hat=rep(0,D)
    }
    
    if (diagonal_structure=="Exchangeable"){
      if (diagonal_homo_across_AI){
        #Using direct Average
        sum_rho=0
        for (group in 1:D){
          sum_rho=sum_rho+sum(diagonal_rho_group[[group]])/sum(weight_diagonal_rho[[group]])
          #Or using direct average on already-truncated data
          #sum_rho=sum_rho+max(sum(diagonal_rho_group[[group]])/sum(weight_diagonal_rho[[group]]),diagonal_ICC_lower_thresh)
        }
        diagonal_rho_hat=rep(max(sum_rho/D,diagonal_ICC_lower_thresh),D)
        
        #using weighted average
        #sum_rho=0; weight_rho=0
        #for (group in 1:D){
        #  sum_rho=sum_rho+sum(diagonal_rho_group[[group]])
        #  weight_rho=weight_rho+sum(weight_diagonal_rho[[group]])
        #}
        #diagonal_rho_hat=rep(max(sum_rho/weight_rho,diagonal_ICC_lower_thresh),D)
        
      }else{
        for (group in 1:D){
          diagonal_rho_hat[group]=max(sum(diagonal_rho_group[[group]])/sum(weight_diagonal_rho[[group]]),diagonal_ICC_lower_thresh)
        }
      }
    }
    
    
    if (diagonal_structure=="AR1"){
      if (diagonal_homo_across_AI){
        
        #Using Direct Average
        for (group in 1:D){
          sum_value=0; sum_weight=0
          for (j in 1:(max_T-1)){
            sum_value=sum_value+diagonal_rho_group[[group]][j,j+1]
            sum_weight=sum_weight+weight_diagonal_rho[[group]][j,j+1]
          }
          diagonal_rho_hat[group]=max(sum_value/sum_weight,diagonal_ICC_lower_thresh)
        }        
        diagonal_rho_hat=rep(mean(diagonal_rho_hat),D)
        
        #Using weighted average
        #sum_value=0; sum_weight=0
        #for (group in 1:D){
        #  for (j in 1:(max_T-1)){
        #    sum_value=sum_value+diagonal_rho_group[[group]][j,j+1]
        #    sum_weight=sum_weight+weight_diagonal_rho[[group]][j,j+1]
        #  }
        #}        
        #diagonal_rho_hat=rep(max(sum_value/sum_weight,diagonal_ICC_lower_thresh),D)
        
      }else{
        for (group in 1:D){
          sum_value=0; sum_weight=0
          for (j in 1:(max_T-1)){
            sum_value=sum_value+diagonal_rho_group[[group]][j,j+1]
            sum_weight=sum_weight+weight_diagonal_rho[[group]][j,j+1]
          }
          diagonal_rho_hat[group]=max(sum_value/sum_weight,diagonal_ICC_lower_thresh)
        }
      }
    }
    
    
    
    if (diagonal_structure=="Unstructured"){
      if (diagonal_homo_across_AI){
        #Using direct average
        for (group in 1:D){
          for (j in 1:max_T){
            for (k in 1:max_T){
              if (j != k){
                diagonal_rho_hat[[group]][j,k]=max(diagonal_rho_group[[group]][j,k]/weight_diagonal_rho[[group]][j,k],diagonal_ICC_lower_thresh)
              }else{
                diagonal_rho_hat[[group]][j,j]=1
              }
            }
          }
        }
        for (j in 1:max_T){
          for (k in 1:max_T){
            entry_sum=0
            for (group in 1:D) entry_sum=entry_sum+diagonal_rho_hat[[group]][j,k]
            for (group in 1:D) diagonal_rho_hat[[group]][j,k]=entry_sum/D
          }
        }
        
      }else{
        for (group in 1:D){
          for (j in 1:max_T){
            for (k in 1:max_T){
              if (j != k){
                diagonal_rho_hat[[group]][j,k]=max(diagonal_rho_group[[group]][j,k]/weight_diagonal_rho[[group]][j,k],diagonal_ICC_lower_thresh)
              }else{
                diagonal_rho_hat[[group]][j,j]=1
              }
            }
          }
        }
      }
    }
    
    
    
    
    
    #-------------------------calculating parameters of off-diagonal blocks---------------------
    if (off_diagonal_structure=="Independence"){
      off_diagonal_rho_hat=rep(0,D)
    }
    
    if (off_diagonal_structure=="Exchangeable"){
      if (off_diagonal_homo_across_AI){
        #Using direct average
        for (group in 1:D){
          off_diagonal_rho_hat[group]=max(sum(off_diagonal_rho_group[[group]])/sum(weight_off_diagonal_rho[[group]]),off_diagonal_ICC_lower_thresh)
        }
        off_diagonal_rho_hat=rep(mean(off_diagonal_rho_hat),D)
        
        #Using weighted average
        #sum_rho=0; weight_rho=0;
        #for (group in 1:D){
        #  sum_rho=sum_rho+sum(off_diagonal_rho_group[[group]])
        #  weight_rho=weight_rho+sum(weight_off_diagonal_rho[[group]])
        #}
        #off_diagonal_rho_hat=rep(max(sum_rho/weight_rho,off_diagonal_ICC_lower_thresh),D)
        
      }else{
        for (group in 1:D){
          off_diagonal_rho_hat[group]=max(sum(off_diagonal_rho_group[[group]])/sum(weight_off_diagonal_rho[[group]]),off_diagonal_ICC_lower_thresh)
        }
      }
    }
    
    if (off_diagonal_structure=="Unstructured"){
      if (off_diagonal_homo_across_AI){
        for (group in 1:D){
          for (j in 1:max_T){
            for (k in 1:max_T){
              off_diagonal_rho_hat[[group]][j,k]=max(off_diagonal_rho_group[[group]][j,k]/weight_off_diagonal_rho[[group]][j,k],off_diagonal_ICC_lower_thresh)
            }
          }
        }
        for (j in 1:max_T){
          for (k in 1:max_T){
            entry_sum=0
            for (group in 1:D) entry_sum=entry_sum+off_diagonal_rho_hat[[group]][j,k]
            for (group in 1:D) off_diagonal_rho_hat[[group]][j,k]=entry_sum/D
          }
        } 
        
      }else{
        for (group in 1:D){
          for (j in 1:max_T){
            for (k in 1:max_T){
              off_diagonal_rho_hat[[group]][j,k]=max(off_diagonal_rho_group[[group]][j,k]/weight_off_diagonal_rho[[group]][j,k],off_diagonal_ICC_lower_thresh)
            }
          }
        }
      }
    }
    
  }
  
  
  rownames(sigma2_hat)=paste("time_",c(1:max_T),sep="")
  colnames(sigma2_hat)=paste("AI(",c(1:D),")",sep="")
  
  
  
  #-------------------------
  #-------------------------
  
  
  #Calculating Sandwich Estimator
  if (!bias_correction){
    meat_lz=matrix(0,p,p)
    for (i in 1:N){
      yi=as.matrix(Y[starting_cluster[i]:ending_cluster[i]],M[i]*max_T,1)
      Ui=matrix(0,p,1)
      for (group in 1:D){
        if (within[i,group]==1){
          Di=matrix(0,M[i]*max_T,0)
          for (l in 1:L){
            if (p_each[l]>0) Di=cbind(Di,Ind[group,l]*X[[l]][starting_cluster[i]:ending_cluster[i],])
          }
          
          mu_hat=matrix(0,M[i]*max_T,1)
          for (l in 1:L){
            if (p_each[l]>1) mu_hat=mu_hat+Ind[group,l]*X[[l]][starting_cluster[i]:ending_cluster[i],1:p_each[l]]%*%beta_hat[starting_ind[l]:ending_ind[l],1]
            if (p_each[l]==1) mu_hat=mu_hat+Ind[group,l]*X[[l]][starting_cluster[i]:ending_cluster[i],1:p_each[l]]*beta_hat[starting_ind[l]]
          }
          
          Vi=diag(rep(sqrt(sigma2_hat[,group]),M[i]))%*%
            Corr_Multilayer(M[i],max_T,diagonal_structure=diagonal_structure,off_diagonal_structure=off_diagonal_structure,para_diagonal=diagonal_rho_hat[[group]],para_off_diagonal=off_diagonal_rho_hat[[group]])%*%
            diag(rep(sqrt(sigma2_hat[,group]),M[i]))
          
          Ui=Ui+weight[i]*t(Di)%*%solve(Vi)%*%(yi-mu_hat)
        }
      }
      meat_lz=meat_lz+Ui%*%t(Ui)
    }
    V_lz=V_naive%*%meat_lz%*%V_naive
    V_result=V_lz
  }else{
    meat_brl=matrix(0,p,p)
    for (i in 1:N){
      yi=as.matrix(Y[starting_cluster[i]:ending_cluster[i]],M[i]*max_T,1)
      Ui=matrix(0,p,1)
      Hi=matrix(0,p,p)
      for (group in 1:D){
        if (within[i,group]==1){
          Di=matrix(0,M[i]*max_T,0)
          for (l in 1:L){
            if (p_each[l]>0) Di=cbind(Di,Ind[group,l]*X[[l]][starting_cluster[i]:ending_cluster[i],])
          }
          
          mu_hat=matrix(0,M[i]*max_T,1)
          for (l in 1:L){
            if (p_each[l]>1) mu_hat=mu_hat+Ind[group,l]*X[[l]][starting_cluster[i]:ending_cluster[i],1:p_each[l]]%*%beta_hat[starting_ind[l]:ending_ind[l],1]
            if (p_each[l]==1) mu_hat=mu_hat+Ind[group,l]*X[[l]][starting_cluster[i]:ending_cluster[i],1:p_each[l]]*beta_hat[starting_ind[l]]
          }
          
          Vi=diag(rep(sqrt(sigma2_hat[,group]),M[i]))%*%
            Corr_Multilayer(M[i],max_T,diagonal_structure=diagonal_structure,off_diagonal_structure=off_diagonal_structure,para_diagonal=diagonal_rho_hat[[group]],para_off_diagonal=off_diagonal_rho_hat[[group]])%*%
            diag(rep(sqrt(sigma2_hat[,group]),M[i]))
          
          Hi=Hi+weight[i]*t(Di)%*%solve(Vi)%*%Di%*%V_naive
          Ui=Ui+weight[i]*t(Di)%*%solve(Vi)%*%(yi-mu_hat)
        }
      }
      Ui=solve(diag(1,p,p)-Hi)%*%Ui
      meat_brl=meat_brl+Ui%*%t(Ui)
    }
    V_brl=V_naive%*%meat_brl%*%V_naive
    V_result=V_brl
    
  }
  
  if (dof_adjustment==1) V_result=V_result*N/(N-p)
  
  
  #------------------------Output-------------------
  
  if (verbose==1){
    cat(N,"Number of Clusters, Minimum Cluster Size =",min(M),", Maximum Cluster Size =",max(M),"\n")
    cat("Algorithm stops after ",count_iter," iterations.\n")
  }
  
  
  summary_paras=matrix(0,nrow=p,ncol=5)
  colnames(summary_paras)=c("Parameter","Estimate","Std.Err","Z Score","p-value")
  summary_paras=data.frame(summary_paras)
  for (l in 1:L){
    if (p_each[l]>0){
      for (i in 1:p_each[l]){
        summary_paras[starting_ind[l]+i-1,1]=paste(Ind_names[l],"*",colnames(X[[l]])[i])
      }
    }
  }
  summary_paras[,2]=beta_hat
  summary_paras[,3]=sqrt(diag(V_result))
  summary_paras[,4]=summary_paras[,2]/summary_paras[,3]
  if (use_t==0){
    summary_paras[,5]=2*(1-pnorm(abs(summary_paras[,4])))
  }else{
    colnames(summary_paras)[4]=c("T Score")
    summary_paras[,5]=2*(1-pt(abs(summary_paras[,4]),df=N-p))
  }
  
  aimed_comparison=NULL
  if (is.null(aimed_comparison)){
    r=0
  }else{
    if (is.null(dim(aimed_comparison))){
      aimed_comparison=matrix(aimed_comparison,nrow=1,ncol=4)
      r=1
    }else{
      r=nrow(aimed_comparison)
      if (!is.matrix(aimed_comparison)) aimed_comparison=as.matrix(aimed_comparison)
    }
    
    for (rr in 1:r){
      if (all.equal(aimed_comparison[order(aimed_comparison[rr,])],c(rep(0,D-2),1,1))==FALSE){
        stop("Formula Error: aimed_comparison should be a length D vector with D-2 0's and two 1's")
      }
      num_1=which(aimed_comparison[rr,]==1)[1]; num_2=which(aimed_comparison[rr,]==1)[2]  #which two AIs are being compared
      weight_beta=coeff[num_1,]-coeff[num_2,]
      
      sum_test=sum(weight_beta*beta_hat[1:q])
      var_test=t(matrix(weight_beta))%*%V_result[1:q,1:q]%*%matrix(weight_beta)
      sd_test=sqrt(var_test)
      value_test=sum_test/sd_test
      if (use_t==0){
        p_test=2*(1-pnorm(abs(value_test)))
      }else{
        p_test=2*(1-pt(abs(value_test),df=N-p-q))
      }
      summary_paras=rbind(summary_paras,c(paste("AI(",as.character(num_1),") - AI(",as.character(num_2),")",sep=""),sum_test,sd_test,value_test,p_test))
    }
  }
  
  
  
  
  if (verbose==1){
    output_paras=data.frame(cbind(summary_paras,rep(0,nrow(summary_paras))))
    output_paras[,2]=round(as.numeric(output_paras[,2]),digits=5)
    output_paras[,3]=round(as.numeric(output_paras[,3]),digits=5)
    output_paras[,4]=round(as.numeric(output_paras[,4]),digits=2)
    for (i in 1:nrow(summary_paras)){
      t=as.numeric(output_paras[i,5])
      output_paras[i,6]=ifelse(t>0.1," ",ifelse(t>0.05,".",ifelse(t>0.01,"*",ifelse(t>0.001,"**","***"))))
      output_paras[i,5]=ifelse(t<2e-16,"<2e-16",ifelse(t>0.001,as.character(round(t,digits=4)),formatC(t,format="e",digits=2)))
    }
    output_paras=rbind(colnames(output_paras),output_paras)
    output_paras[1,6]=""
    output_paras[1,5]=ifelse(use_t==0,"Pr(>|z|)","Pr(>|t|)")
    colnames(output_paras)=NULL
    rownames(output_paras)=NULL
    
    cat("Summary of model Coefficients:\n")
    print(output_paras[1:(p+1),],row.names=FALSE)
    cat("---\n")
    cat("Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1\n")
  }
  
  if (verbose==1){
    cat("Estimated Variance:\n")
    print(sigma2_hat)
    cat("Estimated Correlation for each Adaptive Intervention:\n")
    for (group in 1:D){
      cat("AI(",as.character(group),"), diagonal block\n")
      print(Corr_Multilayer(2,max_T,diagonal_structure=diagonal_structure,off_diagonal_structure=off_diagonal_structure,para_diagonal=diagonal_rho_hat[[group]],para_off_diagonal=off_diagonal_rho_hat[[group]])[1:max_T,1:max_T])
      cat("AI(",as.character(group),"), off-diagonal block\n")
      print(Corr_Multilayer(2,max_T,diagonal_structure=diagonal_structure,off_diagonal_structure=off_diagonal_structure,para_diagonal=diagonal_rho_hat[[group]],para_off_diagonal=off_diagonal_rho_hat[[group]])[1:max_T,(max_T+1):(2*max_T)])
    }
  }
  
  var_paras=list(marginal_var=sigma2_hat,diagonal_structure=diagonal_structure,off_diagonal_structure=off_diagonal_structure,para_diagonal=diagonal_rho_hat,para_off_diagonal=off_diagonal_rho_hat)
  
  #if (r>0){
  #  summary_comparison=summary_paras[(p+q+1):(p+q+r),]
  #  colnames(summary_comparison)[1]="Comparison"
  #  rownames(summary_comparison)=c()
  # summary_paras=summary_paras[1:(p+q),]
  #  if (verbose==1){
  #    output_comparison=output_paras[c(1,(p+q+2):(p+q+r+1)),]
  #    output_comparison[1,1]="Comparison"
  #    if (r==1){
  #      cat("Result of Proposed Comparison between ",AI_names[num_1]," and ",AI_names[num_2],"\n")
  #    }else{
  #      cat("Result of Proposed Comparison between Adaptive Interventions:\n")
  #    }
  #    print(output_comparison,row.names=FALSE)
  #    cat("---\n")
  #    cat("Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1\n")
  #    cat("\n")
  #  }
  #}
  
  if (verbose>0){
    if ((p>0)|(var_homo_across_time)) cat("Warning:\n")
    if (p>0){
      cat("All Inference are valid only when all covariates are mean zero\n")
      cat("All covariates should be baseline only to avoid causal problems\n")
      #warning("All Inference are valid only when all covariates are mean zero",call.=FALSE)
      #warning("All covariates should be baseline only to avoid causal problems",call.=FALSE)
    }
    if (var_homo_across_time){
      cat("Setting working variation homogeneous across time may lead to problems when there are too many timepoints and variance increases with time\n")
    }
  }
  
  rownames(V_result)=summary_paras[1:p,1]
  colnames(V_result)=summary_paras[1:p,1]
  
  if (verbose==1){
    cat("\n")
    cat("Variance-Covariance matrix of the estimates\n")
    print(V_result)
  }
  
  if (!is.null(aimed_comparison)){
    return(list(summary_paras=summary_paras,var_estimator=V_result,
                comparison=summary_comparison,var_paras=var_paras,N=N,max_T=max_T,count_iter=count_iter))
  }else{
    return(list(summary_paras=summary_paras,var_estimator=V_result,
                comparison=NULL,var_paras=var_paras,N=N,max_T=max_T,count_iter=count_iter))
  }
}



#getting working variance/correlation matrix
#result:   the result of solve_SMART_multilayer
#group:  the label of AI
#diagonal:   TRUE (give the diagonal block)  or FALSE (give the off-diagonal block)
#corr_only:  TRUE (give the correlation matrix) or FALSE (give the variance matrix)
get_var=function(result,group=group,diagonal=TRUE,corr_only=TRUE){
  max_T=result$max_T
  sigma2_hat=result$var_paras$marginal_var
  diagonal_structure=result$var_paras$diagonal_structure
  off_diagonal_structure=result$var_paras$off_diagonal_structure
  diagonal_rho_hat=result$var_paras$para_diagonal
  off_diagonal_rho_hat=result$var_paras$para_off_diagonal
  
  if (diagonal==TRUE){
    cor_matrix=Corr_Multilayer(2,max_T,diagonal_structure=diagonal_structure,off_diagonal_structure=off_diagonal_structure,para_diagonal=diagonal_rho_hat[[group]],para_off_diagonal=off_diagonal_rho_hat[[group]])[1:max_T,1:max_T]
  }else{
    cor_matrix=Corr_Multilayer(2,max_T,diagonal_structure=diagonal_structure,off_diagonal_structure=off_diagonal_structure,para_diagonal=diagonal_rho_hat[[group]],para_off_diagonal=off_diagonal_rho_hat[[group]])[1:max_T,(max_T+1):(2*max_T)]
  }
  
  if (corr_only==TRUE){
    return(cor_matrix)
  }else{
    sd_mat=diag(sqrt(sigma2_hat[,group]))
    var_matrix=sd_mat %*% cor_matrix %*% sd_mat
    return(var_matrix)
  }
}



