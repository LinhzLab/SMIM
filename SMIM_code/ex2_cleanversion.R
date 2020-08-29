##### needing to choose different h#####################
##### needing to decide tuning principle ###############


## L2_norm
L2_norm=function(x) return(sqrt(sum(x^2)))

## Epanechnikov kernel
K<-function(u, h){                             
  3/4*(1-(u/h)^2)*ifelse(abs(u/h)<=1, 1, 0) 
}


############# generate data ##########################
######################################################

AR<-function(n,p,rho){
  wn.sd <- sqrt(1 - rho^2)   
  wn.mat <- matrix(rnorm(n*p), nrow=n, ncol=p) * wn.sd  ## wn means white noise
  Xmat <- matrix(0, nrow=n, ncol=p)
  Xmat[,1] <- rnorm(n) 
  for(j in 2:p)
    Xmat[,j] <- rho * Xmat[,j-1] + wn.mat[,j]
  Xmat=matrix(pmax(-1,pmin(1,Xmat)),n,p)   
  return(Xmat)
}

f1=function(x) return(5*x)
f3=function(x) return(0.5^2*(x+2)^2)

gdata=function(n,rho,p){
  Xmat=AR(n,p,rho)
  epsi=rnorm(n)
  X1=c(Xmat%*%beta1) 
  X3=c(Xmat%*%beta3)
  Y=f1(X1)+epsi*sqrt(f3(X3))
  return(list(Xmat=Xmat,Y=Y))
}

######### estimate theta ##########################
#####################################################

## likelihood
lik=function(y,v1,v3){
  v=v3^2+eps
  result=-log(v)/2-(y-v1)^2/(2*v)
  result[is.na(result)|abs(result)==Inf]=0
  return(result)
}

lik_gamma=function(y,v1,Sigma,gamma){
  v=Sigma^2
  z=(y-v1)/Sigma*sqrt(gamma/(1-gamma))
  result=-log(v)/2-(y-v1)^2/(2*v)+log(1-pnorm(z))
  result[is.na(result)|abs(result)==Inf]=0
  return(result)
}


## Score function
Q1.1=function(y, v1, v3){
  v=v3^2
  result=(y-v1)/v
  return(result)
}

Q1.3=function(y, v1, v3){
  result=-1/v3+(y-v1)^2/v3^3
  return(result)
}


#### updata theta

#### calculate bspline function ########################
cal_bs=function(Xbeta){
  a=max(abs(Xbeta))
  CC=matrix(0,df-1,df)   
  Boundary.knots=range(-a, a)
  knots <- seq.int(from =(-a), to=a, length.out = nIknots + 2L)[-c(1L, nIknots + 2L)]  
  BS1=bs(Xbeta,degree=m,knots=knots,Boundary.knots=Boundary.knots, intercept =T)  
  tBS1=bs(Xbeta,degree=m-1,knots=knots,Boundary.knots=Boundary.knots, intercept=T) 
  temp1=rep(-a,m+1);temp2=rep(a,m+1)
  Aknots <- c(temp1, knots,temp2)  
  for(t in 1:(df-1)){
    CC[t,(t:(t+1))]=c(-m/(Aknots[t+m+1]-Aknots[t+1]),m/(Aknots[t+m+1]-Aknots[t+1]))
  }
  BS=BS1    
  tBS=(tBS1%*%CC)
  return(list(BS=BS,tBS=tBS))
}



uptheta=function(data1,beta.old){
  beta=beta.old
  beta1=beta$beta1
  beta3=beta$beta3
  
  Xmat=data1$Xmat; Y=data1$Y
  Xbeta1=c(Xmat%*%beta1)
  Xbeta3=c(Xmat%*%beta3)
  bs.beta1=cal_bs(Xbeta1)
  bs.beta3=cal_bs(Xbeta3)

  theta1.old=lm(Y~bs.beta1$BS-1)$coef
  eta11=c(bs.beta1$BS%*%theta1.old)
  eta12=c(bs.beta1$tBS%*%theta1.old)
  theta3.old=lm(abs(Y-eta11)~bs.beta3$BS-1)$coef
  eta31=c(bs.beta3$BS%*%theta3.old)
  theta.old=c(theta1.old,theta3.old)
  
  lik_theta=function(theta){
    theta1=theta[1:df];theta3=theta[(df+1):(2*df)]
    eta11=c(bs.beta1$BS%*%theta1)
    eta31=c(bs.beta3$BS%*%theta3)
    Sigma=abs(eta31)
    gamma=0/Sigma^2
    result=mean(lik_gamma(Y,eta11,Sigma,gamma))
    return(-result)
  }
  
  diff=1;iter=0
  while(diff>epsilon&iter<iter.max){
    theta.new=optim(theta.old,f=lik_theta,method="BFGS")$par
    eta11=c(bs.beta1$BS%*%theta.new[1:df])
    eta31=c(bs.beta3$BS%*%theta.new[(df+1):(2*df)])
    Sigma=abs(eta31)
    gamma.new=0/Sigma^2
    diff=L2_norm(theta.new-theta.old)
    theta.old=theta.new
    gamma=gamma.new
    iter=iter+1
  }

  eta11.new=eta11
  theta1.new=lm(eta11.new~bs.beta1$BS-1)$coef
  eta12.new=c(bs.beta1$tBS%*%theta1.new)
  eta31.new=Sigma
  eta32.new=c(bs.beta3$tBS%*%theta.new[(df+1):(2*df)])
  return(list(eta11=eta11.new,eta12=eta12.new,eta31=eta31.new,eta32=eta32.new))
}

upbeta1=function(data1,beta.old,eta.old){
  beta1=beta.old$beta1
  beta3=beta.old$beta3
  eta11=eta.old$eta11
  eta12=eta.old$eta12
  eta31=eta.old$eta31
  Xmat=data1$Xmat; Y=data1$Y
  Xbeta1=c(Xmat%*%beta1)
  Xbeta3=c(Xmat%*%beta3)
  q1=Q1.1(Y,eta11,eta31)*eta12*Xmat
  u1=-apply(q1,2,mean)
  H1=diag(var(q1))
  XX1=sqrt((H1));XX=diag(XX1)
  YY=XX1*beta1-u1/(XX1+rep(1e-5,p))
  fit=grpreg(XX,YY,penalty="cMCP",family="gaussian")
  fit.beta=fit$beta[-1,]
  re.beta=fit.beta/matrix(apply(fit.beta,2,L2_norm),nr=p,nc=100,byrow=T)
  re.beta[is.na(re.beta)|is.nan(re.beta)]=0
  
  reest_beta=function(beta1){
    Xbeta1=c(as.matrix(Xmat[,id.beta1])%*%beta1);
    Xbeta1=pmax(pmin(Xbeta1,bound),-bound)
    id1=round((Xbeta1+bound)/(2*bound)*num)+1;id1[id1>num]=num;
    eta11=eta.newk$eta11[id1]
    result=mean(lik(Y,eta11,eta31))
    return(-result)
  }
  
  len.lam=ncol(re.beta);fit.loss=BIC=rep(1e+15,len.lam)
  for(j in 1:len.lam){
    re.betaj=re.beta[,j]
    id.beta1=(1:p)[re.betaj!=0]
    if(!is.na(sum(id.beta1))&sum(id.beta1)!=0){
      fit.df=sum(re.betaj!=0)
      sgn.betaj=sign(re.betaj[re.betaj!=0])
      re.betaj=re.betaj*sgn.betaj[1]
      fit.loss[j]=reest_beta(re.betaj[id.beta1])
      BIC[j]=fit.loss[j]/(1-fit.df/n)^2
    }
  }
  BIC.min1=(1:100)[BIC==min(BIC)]
  BIC.min=BIC.min1[1]
  
  re.beta1=rep(0,p)
  id.beta1=(1:p)[re.beta[,BIC.min]!=0]
  len.l=length(id.beta1)
  if(len.l>0) re.beta1=re.beta[,BIC.min]
  sgn.beta1=sign(re.beta1[re.beta1!=0])
  re.beta1=re.beta1*sgn.beta1[1]
  return(re.beta1)
}

upbeta3=function(data1,beta.old,eta.old){
  beta=beta.old
  beta3=beta$beta3
  eta11=eta.old$eta11;eta12=eta.old$eta12
  eta31=eta.old$eta31 ;eta32=eta.old$eta32
  Xmat=data1$Xmat; Y=data1$Y
  Xbeta3=c(Xmat%*%beta3)
  q3=Q1.3(Y,eta11,eta31)*eta32*Xmat
  u3=-apply(q3,2,mean)
  H3=diag(var(q3))  
  XX3=sqrt(H3);XX=diag(XX3)
  YY=XX3*beta3-u3/(XX3+rep(1e-5,p))
  fit=grpreg(XX,YY,penalty="cMCP",family="gaussian")
  fit.beta=fit$beta[-1,]
  re.beta=fit.beta/matrix(apply(fit.beta,2,L2_norm),nr=p,nc=100,byrow=T)*ss
  re.beta[is.na(re.beta)|is.nan(re.beta)]=0
  
  reest_beta=function(beta3){
    Xbeta3=c(as.matrix(Xmat[,id.beta3])%*%beta3);
    Xbeta3=pmax(pmin(Xbeta3,bound),-bound)
    id3=round((Xbeta3+bound)/(2*bound)*num)+1;id3[id3>num]=num;
    eta31=eta.newk$eta31[id3]
    result=mean(lik(Y,eta11,eta31))
    return(-result)
  }
  
  len.lam=ncol(re.beta);fit.loss=BIC=rep(1e+15,len.lam)
  for(j in 1:len.lam){
    re.betaj=re.beta[,j]
    id.beta3=(1:p)[re.betaj!=0]
    if(!is.na(sum(id.beta3))&sum(id.beta3)!=0){
      fit.df=sum(re.betaj!=0)
      sgn.betaj=sign(re.betaj[re.betaj!=0])
      re.betaj=-re.betaj*sgn.betaj[1]
      fit.loss[j]=reest_beta(re.betaj[id.beta3])
      BIC[j]=2*fit.loss[j]/(1-fit.df/n)^2
    }
  }
  BIC.min3=(1:100)[BIC==min(BIC)]
  BIC.min=BIC.min3[1]
  
  re.beta3=rep(0,p)
  id.beta3=(1:p)[re.beta[,BIC.min]!=0]
  len.l=length(id.beta3)
  if(len.l>0) re.beta3=re.beta[,BIC.min]
  sgn.beta3=sign(re.beta3[re.beta3!=0])
  re.beta3=re.beta3*sgn.beta3[1]
  return(-re.beta3)
}

upbeta=function(data1,beta.old,eta.old){
  beta=beta.old
  beta1.old=beta$beta1
  beta3.old=beta$beta3
  eta11=eta.old$eta11;eta12=eta.old$eta12
  eta31=eta.old$eta31 ;eta32=eta.old$eta32
  Xmat=data1$Xmat; Y=data1$Y

  beta.old=c(beta1.old,beta3.old)
  
  reest_beta=function(beta){
    beta1=beta[1:p];beta3=beta[(p+1):(2*p)]
    Xbeta1=c(Xmat%*%beta1)
    Xbeta3=c(Xmat%*%beta3)
 
    bound=3;
    Xbeta1=pmax(pmin(Xbeta1,bound),-bound)
    id1=round((Xbeta1+bound)/(2*bound)*num)+1;id1[id1>num]=num
    eta11=eta11[id1]
    bound=3;Xbeta3=pmax(pmin(Xbeta3,bound),-bound)
    id3=round((Xbeta3+bound)/(2*bound)*num)+1;id3[id3>num]=num
    eta31=eta31[id3]
    
    result=mean(lik(Y,eta11,eta31))
    return(-result)
  }
  
  beta.new=optim(beta.old,f=reest_beta,method="BFGS")$par

  re.beta1=beta.new[1:p]
  re.beta1=re.beta1/L2_norm(re.beta1)
  sgn.beta1=sign(re.beta1[re.beta1!=0])
  re.beta1=re.beta1*sgn.beta1[1]
  
  re.beta3=beta.new[(p+1):(2*p)]
  re.beta3=re.beta3/L2_norm(re.beta3)
  sgn.beta3=sign(re.beta3[re.beta3!=0])
  re.beta3=-re.beta3*sgn.beta3[1]
  return(list(beta1=re.beta1,beta3=re.beta3))
}

uptheta.ker=function(data1,beta.old,eta.old){
  beta=beta.old
  beta1=beta$beta1
  beta3=beta$beta3
  eta11=eta.old$eta11;eta12=eta.old$eta12
  eta31=eta.old$eta31
  Xmat=data1$Xmat; Y=data1$Y
  Xbeta1=c(Xmat%*%beta1)
  Xbeta3=c(Xmat%*%beta3)
  Xbetaij.mat1=matrix(Xbeta1,n,n)-matrix(Xbeta1,n,n,byrow=T)
  Xbetaij.mat3=matrix(Xbeta3,n,n)-matrix(Xbeta3,n,n,byrow=T)
  bs.beta1=cal_bs(Xbeta1)
  bs.beta3=cal_bs(Xbeta3)
  Kij.mat1=K(Xbetaij.mat1,h)
  Kij.mat3=K(Xbetaij.mat3,h)
  
  lik_theta=function(theta){
    theta1=theta[1:df];theta3=theta[(df+1):(2*df)]
    eta11=c(bs.beta1$BS%*%theta1)
    eta12.new=c(bs.beta1$tBS%*%theta1)
    eta31=c(bs.beta3$BS%*%theta3)
    eta32.new=c(bs.beta3$tBS%*%theta3)
    v1j=t((eta11+eta12.new*t(Xbetaij.mat1))*t(Kij.mat1))++eta11*(1-Kij.mat1)
    v3j=t((eta31+eta32.new*t(Xbetaij.mat3))*t(Kij.mat3))+eta31*(1-Kij.mat3)
    result=mean(lik(Y,v1j,abs(v3j)))
    return(-result)
  }
  theta1.init=lm(eta11~bs.beta1$BS-1)$coef
  theta3.init=lm(abs(eta31)~bs.beta3$BS-1)$coef
  theta.init=c(theta1.init,theta3.init)
  
  theta.new=optim(theta.init,f=lik_theta,method="BFGS")$par
  eta11.new=c(bs.beta1$BS%*%theta.new[1:df])
  eta12.new=c(bs.beta1$tBS%*%theta.new[1:df])
  eta31.new=c(bs.beta3$BS%*%theta.new[(df+1):(2*df)])
  eta32.new=c(bs.beta3$tBS%*%theta.new[(df+1):(2*df)])
  return(list(eta11=eta11.new,eta12=eta12.new,eta31=eta31.new,eta32=eta32.new))
}

sortxy=function(x1,y1){
  result=numeric(lenx)
  sortx = sort(x1)
  sorty = y1[order(x1)]
  for(i in 1:lenx) {
    if(pointx[i]< sortx[1])  result[i] = sorty[1] else
      if(pointx[i] >= sortx[n])  result[i] = sorty[n] else
        for(j in 2:n) 
          if(pointx[i] >= sortx[j-1] & pointx[i] < sortx[j])
            result[i] = sorty[j-1]
  }
  return(result)
}

library("mgcv")
library("grpreg")
library("splines")
set.seed(88650)
NSIM=10;epsilon=0.001;iter.max=100;eps=1e-6
rho=0.5
lam.seq=c(0.05,(1:10)/10);len.lam=length(lam.seq)
n=600; 
p=8;q1=4
beta1=c(0.8,0.4,-0.4,0.2,rep(0,p-q1));beta1=beta1/L2_norm(beta1)
beta3=c(-0.45,rep(0,6),0.9);beta3=beta3/L2_norm(beta3)
beta=list(beta1=beta1,beta3=beta3)
epi=1e-3
m=3;qn=4
df=m+qn
order=m+1  
nIknots=df-order 
bound=3
pointx=seq(-bound,bound,length.out=1000);lenx=length(pointx)
h=n^{-1/5}
num=1000

beta1.est2=beta2.est2=beta3.est2=matrix(0,NSIM,p)
eta11.est2=eta31.est2=matrix(0,NSIM,lenx)
for(isim in 1:NSIM){
  data1=gdata(n,rho,p)
  Xmat=data1$Xmat;Y=data1$Y
  beta1.old=beta1+runif(p,0,0.1)
  beta3.old=beta3+runif(p,0,0.1) 
  beta.old=list(beta1=beta1.old,beta3=beta3.old)

  for(j in 1:30){
    Xbeta1.old=c(Xmat%*%beta.old$beta1)
    Xbeta3.old=c(Xmat%*%beta.old$beta3)
    eta.new=uptheta(data1,beta.old)
    eta11=sortxy(Xbeta1.old,eta.new$eta11)
    eta31=sortxy(Xbeta3.old,eta.new$eta31)
    eta.newk=list(eta11=eta11,eta31=eta31)
    beta.new=upbeta(data1,beta.old,eta.newk)
    diff=max(abs(beta.new$beta1-beta.old$beta1))
    if(diff<1e-3) break 
    beta.old=beta.new
  }
  eta.new2=uptheta.ker(data1,beta.new,eta.new)
  beta1.est2[isim,]=upbeta1(data1,beta.new,eta.new2)
  beta3.est2[isim,]=upbeta3(data1,beta.new,eta.new2)
  eta11.est2[isim,]=sortxy(Xbeta1.old,eta.new2$eta11)
  eta31.est2[isim,]=sortxy(Xbeta3.old,eta.new2$eta31)
}
