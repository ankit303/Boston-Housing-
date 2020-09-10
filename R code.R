#-------------------------------------------------------(a)-------------------------------------------------------------------------------
library("matlib")
library("Matrix")
library("ggplot2")
my<-read.delim("abhi.txt",sep=",")
mdata<-data.frame(my)
D<-as.matrix(mdata)
X<-D[,c(-which(colnames(mdata)=="zn"),-which(colnames(mdata)=="chas"),-which(colnames(mdata)=="rad"),-which(colnames(mdata)=="medv"))]

for(i in 1:ncol(X)){                                    ## Standarized regressor 
  X[,i]<-(X[ ,i]-mean(X[,i]))/sd(X[,i])
  }
Y<-D[,which(colnames(mdata)=="medv")]
Y<-as.matrix(Y)
set.seed(1230)
I<-sample(1:506,100)
Xtest<-X[I,]
Ytest<-Y[I,]
##View(Ytest)
Xtrain<-X[-I,]
Ytrain<-Y[-I,]
x_1<-matrix(data=c(rep(1,406)),ncol=1)
Xtrain<-cbind(x_1,Xtrain)
colnames(Xtrain)<-c("intersept","crim","indus","nox","rm","age","dis","tax","ptratio","b","lstat")
x_1<-matrix(data=c(rep(1,100)),ncol=1)
Xtest<-cbind(x_1,Xtest)
colnames(Xtest)<-c("intersept","crim","indus","nox","rm","age","dis","tax","ptratio","b","lstat")
b_ols<-solve(t(Xtrain)%*%Xtrain)%*%t(Xtrain)%*%Ytrain                  ##ols estimate of beta 
ytest_ols<-Xtest%*%b_ols                                             ## ols estimate of y 
##--------------------------------calculating RMSE---------------------------------------
Rmse<-sqrt((t(Ytest-ytest_ols)%*%(Ytest-ytest_ols)))/10                ## RMSE 
Rmse
##-----------------------------------(A)----------------------------------------------------
## calculating hii
H<- Xtrain%*%solve(t(Xtrain)%*%Xtrain)%*%t(Xtrain)
H<-as.matrix(H)
hii<-cbind(rep(0,nrow(H)))
i=1
while(i <=406){
  hii[i,]<-H[i,i]
  i<-i+1
}
p<-sum(hii)
q<-2*p/nrow(hii)
P<-cbind(rep(1,nrow(hii)))*(q+.02)
S<-which(hii>P)
s<-c(S)
A<-data.frame(hii,P,t=c(1:406))
ggplot(A,aes(t,hii))+geom_point(col="red",size=2)+scale_y_continuous(limits=c(0,0.5))+
  ggtitle("Hii vs Observations")+
xlab("Observation")+ylab("Hii")+theme(axis.title=element_text(colour="Green",size=20),axis.text=element_text(size=20),plot.title=element_text(colour="Black",size=20))

###-----------------------------------------------1(b)------------------------------
## cook's distance
Ytrain<-as.matrix(Ytrain)
Xtrain<-as.matrix(Xtrain)
Ytrain<-as.matrix(Ytrain[-s,])
Xtrain<-Xtrain[-s,]
i=1
beta_i<-matrix(c(rep(0,11*nrow(Xtrain))),nrow=11,ncol=nrow(Xtrain))
while(i<=nrow(Xtrain)){
  beta_i[ ,i]<-solve(t(Xtrain[-i,])%*%Xtrain[-i,])%*%t(Xtrain[-i,])%*%(Ytrain[-i,])
  i<-i+1
}
## here beta_i denotes estimated value of beta after deleting the ith observation  
H1<- Xtrain%*%solve(t(Xtrain)%*%Xtrain)%*%t(Xtrain)
H1<-as.matrix(H1)
i=1
t<-nrow(Xtrain)
t1<-rankMatrix(H1)
In<-diag(x=1,ncol=t,nrow=t)
Msres<-(t(Ytrain)%*%(In-H1)%*%(Ytrain))/(t-t1)
## here MSres represent ms residual
Di<-matrix(c(rep(0,t)),nrow=t,ncol=1)
for(i in 1:t){
Di[i,]<-(t(beta_i[,i]-b_ols)%*%t(Xtrain)%*%Xtrain%*%(beta_i[,i]-b_ols))/(t1*Msres)
}
A3<- matrix(c(rep(0,t)),nrow=t,ncol=1)
A3[which(Di>1)]<-Di[which(Di>1) ]
influential_cook<-c(which(Di>1))

B<-data.frame(Di,ti=c(1:t))
ggplot(B,aes(ti,Di))+geom_point()+ ggtitle("Cook's vs Observations")+
  xlab("Observation")+ylab("Cook's")+theme(axis.title=element_text(colour="Green",size=20),axis.text=element_text(size=20),plot.title=element_text(colour="Red",size=20))

col_name<-colnames(Xtrain)
#--------------------------Calculating Covratio------------------------------------
b_ols_D<-solve(t(Xtrain)%*%Xtrain)%*%t(Xtrain)%*%Ytrain
Ytrain_hat<-Xtrain%*%(as.matrix(b_ols_D))
i<-1
ei2<-(Ytrain-Ytrain_hat)*(Ytrain-Ytrain_hat)
Si2<-matrix(c(rep(0,t)),nrow=t,ncol=1)
for(i in 1:t){
  Si2[i]<-((t-t1)*Msres-(ei2[i]/(1-hii[i])))/(t-t1-1)
}
covr<-matrix(c(rep(0,t)),nrow=t,ncol=1)
i<-1
for(i in 1:t){
  covr[i]<-(Si2[i]^t1)/((Msres^t1)*(1-hii[i]))
}
##  calculating inluencial point 
A1<-matrix(c(rep(0,t)),nrow=t,ncol=1)

A1[which(covr>1+(3*t1)/t |covr<1-(3*t1)/t )]<-covr[which(covr>1+(3*t1)/t |covr<1-(3*t1)/t )]
influential_cov<-c(which(covr>1+(3*t1)/t |covr<1-(3*t1)/t ))

#-------------------calculating Dffit------------------------
student<-matrix(c(rep(0,t)),nrow=t,ncol=1)
i<-1
for(i in 1:t){
  student[i]<-sqrt(ei2[i]/Si2[i]*(1-hii[i]))
}
i<-1

Dfit_i<-matrix(c(rep(0,t)),nrow=t,ncol=1)

for(i in 1:t){
  Dfit_i[i]<-sqrt(hii[i]/(1-hii[i]))*student[i]
}
A2<-matrix(c(rep(0,t)),nrow=t,ncol=1)
A2[which(Dfit_i > 2*sqrt(t1/t) |Dfit_i< -2*sqrt(t1/t))]<-Dfit_i[which(Dfit_i > 2*sqrt(t1/t) |Dfit_i< -2*sqrt(t1/t)) ]
inflential_dffit<-c(which(Dfit_i > 2*sqrt(t1/t) |Dfit_i< -2*sqrt(t1/t)))

##------------------------------Calculating DFBETA-----------------------------------------------------------
DFBETA<-matrix(data=NA,ncol=t1,nrow=t)
c_matrix<-solve(t(Xtrain)%*%(Xtrain))
for(j in 1:t1){
  for(i in 1:t ){
    DFBETA[i,j]<-(b_ols_D[j]-beta_i[j,i])/sqrt(Si2[i]*c_matrix[j,j])
  }
  
}
colnames(DFBETA)<-c("intersept","crim","indus","nox","rm","age","dis","tax","ptratio","b","lstat")
##calculating influential point
influential_df<-c(rep(0,t))
for(j in 1:t){
  
  influential_df[j]<-length(as.numeric(c(which(abs(DFBETA[j,])>2/sqrt(t)))))
}
influencial_dfbeta<-c(which(influential_df>=8))

#------------------------------1(C)--------------------------------------------------------------
ti=c(1:t)
plot(ti,covr,col=ifelse(covr==A1,"Red","black"),pch=ifelse(covr==A1,19,1),cex=ifelse(covr==A1,2,1),main="Covratio vs Observation",xlab="Observation",ylab="Covratio",col.lab="Green") 
plot(ti,Dfit_i,col=ifelse(Dfit_i==A2,"Red","black"),pch=ifelse(Dfit_i==A2,19,1),cex=ifelse(Dfit_i==A2,2,1),main="Dffit vs Observation",xlab="Observation",ylab="Dffit",col.lab="Green")
plot(ti,Di,col=ifelse(Di==A3,"Red","black"),pch=ifelse(Di==A3,19,1),cex=ifelse(Di==A3,2,1),main="Cooks vs Observations",xlab="Observation",ylab="Cooks",col.lab="Green")

influential<-(union(union(influential_cook,influencial_dfbeta),intersect(inflential_dffit,influential_cov)))
Xtrain<-as.matrix(Xtrain)
Ytrain<-as.matrix(Ytrain)
Xtrain<-Xtrain[-c(influential),]
Ytrain<-Ytrain[-c(influential),]

#(d)----------------------------------------------------------------
b_ols1<-solve(t(Xtrain)%*%Xtrain)%*%t(Xtrain)%*%Ytrain                  ##ols estimate of beta 

dim(b_ols1)
dim(as.matrix(Ytrain))
ytest_ols1<-Xtest%*%b_ols1                                             ## ols estimate of y 
##--------------------------------calculating RMSE---------------------------------------
Rmse1<-sqrt((t(Ytest-ytest_ols1)%*%(Ytest-ytest_ols1)))/(rankMatrix(Xtrain))                ## RMSE 
Rmse1

 #Question 2---------------------------------------------------------------------------------------------------------------------
#-------------------------------(a)----------------------------------------------------------------------------------------------

ytrain_hat1<-(Xtrain)%*%solve(t(Xtrain)%*%Xtrain)%*%t(Xtrain)%*%Ytrain
e_1<-Ytrain-ytrain_hat1
col_name<-colnames(Xtrain)
par(mfrow=c(2,5))
for(i in 2:11){
  plot(x=Xtrain[,i],y=e_1,xlab=col_name[i],ylab="Regressor",col="Red",cex=0.8,col.lab="Blue")
}

#-----------------------------2(b)----------------------------
Beta_D<-solve(t(Xtrain)%*%Xtrain)%*%t(Xtrain)%*%Ytrain  # calculating beta cap asfter deleting some influencial point
i<-1
par(mfrow=c(2,2))
for(i in 2:3){
  Xtrain_ap<-cbind(Xtrain,(as.vector(Xtrain[,i])^2))
  ytrain_hat2<-(Xtrain_ap)%*%solve(t(Xtrain_ap)%*%Xtrain_ap)%*%t(Xtrain_ap)%*%Ytrain
  gama<-solve(t(Xtrain_ap)%*%Xtrain_ap)%*%t(Xtrain_ap)%*%Ytrain
  e_2<-Ytrain-ytrain_hat2
  rs<-as.vector(e_1)+as.vector(Beta_D[i]%*%Xtrain[,i])     # for calculating Cpr
  rs1<-as.vector(e_2)+as.vector(gama[i]%*%Xtrain[,i])+as.vector(gama[12]*as.vector((Xtrain[,i])^2))   #for calculating Apr
  plot(x=Xtrain[,i],y=rs,xlab=col_name[i],ylab="CPR",col="Red",cex=0.8,col.lab="Blue")
  plot(x=Xtrain[,i],y=rs1,xlab=col_name[i],ylab="APR",col="Red",cex=0.8,col.lab="Blue")
}

par(mfrow=c(2,2))
for(i in 4:5){
  Xtrain_ap<-cbind(Xtrain,(as.vector(Xtrain[,i])^2))
  ytrain_hat2<-(Xtrain_ap)%*%solve(t(Xtrain_ap)%*%Xtrain_ap)%*%t(Xtrain_ap)%*%Ytrain
  gama<-solve(t(Xtrain_ap)%*%Xtrain_ap)%*%t(Xtrain_ap)%*%Ytrain
  e_2<-Ytrain-ytrain_hat2
  rs<-as.vector(e_1)+as.vector(Beta_D[i]%*%Xtrain[,i])     # for calculating Cpr
  rs1<-as.vector(e_2)+as.vector(gama[i]%*%Xtrain[,i])+as.vector(gama[12]*as.vector((Xtrain[,i])^2))   #for calculating Apr
  plot(x=Xtrain[,i],y=rs,xlab=col_name[i],ylab="CPR",col="Red",cex=0.8,col.lab="Blue")
  plot(x=Xtrain[,i],y=rs1,xlab=col_name[i],ylab="APR",col="Red",cex=0.8,col.lab="Blue")
}

par(mfrow=c(2,2))
for(i in 6:7){
  Xtrain_ap<-cbind(Xtrain,(as.vector(Xtrain[,i])^2))
  ytrain_hat2<-(Xtrain_ap)%*%solve(t(Xtrain_ap)%*%Xtrain_ap)%*%t(Xtrain_ap)%*%Ytrain
  gama<-solve(t(Xtrain_ap)%*%Xtrain_ap)%*%t(Xtrain_ap)%*%Ytrain
  e_2<-Ytrain-ytrain_hat2
  rs<-as.vector(e_1)+as.vector(Beta_D[i]%*%Xtrain[,i])     # for calculating Cpr
  rs1<-as.vector(e_2)+as.vector(gama[i]%*%Xtrain[,i])+as.vector(gama[12]*as.vector((Xtrain[,i])^2))   #for calculating Apr
  plot(x=Xtrain[,i],y=rs,xlab=col_name[i],ylab="CPR",col="Red",cex=0.8,col.lab="Blue")
  plot(x=Xtrain[,i],y=rs1,xlab=col_name[i],ylab="APR",col="Red",cex=0.8,col.lab="Blue")
}
par(mfrow=c(2,2))
for(i in 8:9){
  Xtrain_ap<-cbind(Xtrain,(as.vector(Xtrain[,i])^2))
  ytrain_hat2<-(Xtrain_ap)%*%solve(t(Xtrain_ap)%*%Xtrain_ap)%*%t(Xtrain_ap)%*%Ytrain
  gama<-solve(t(Xtrain_ap)%*%Xtrain_ap)%*%t(Xtrain_ap)%*%Ytrain
  e_2<-Ytrain-ytrain_hat2
  rs<-as.vector(e_1)+as.vector(Beta_D[i]%*%Xtrain[,i])     # for calculating Cpr
  rs1<-as.vector(e_2)+as.vector(gama[i]%*%Xtrain[,i])+as.vector(gama[12]*as.vector((Xtrain[,i])^2))   #for calculating Apr
  plot(x=Xtrain[,i],y=rs,xlab=col_name[i],ylab="CPR",col="Red",cex=0.8,col.lab="Blue")
  plot(x=Xtrain[,i],y=rs1,xlab=col_name[i],ylab="APR",col="Red",cex=0.8,col.lab="Blue")
}
par(mfrow=c(2,2))
for(i in 10:11){
  Xtrain_ap<-cbind(Xtrain,(as.vector(Xtrain[,i])^2))
  ytrain_hat2<-(Xtrain_ap)%*%solve(t(Xtrain_ap)%*%Xtrain_ap)%*%t(Xtrain_ap)%*%Ytrain
  gama<-solve(t(Xtrain_ap)%*%Xtrain_ap)%*%t(Xtrain_ap)%*%Ytrain
  e_2<-Ytrain-ytrain_hat2
  rs<-as.vector(e_1)+as.vector(Beta_D[i]%*%Xtrain[,i])     # for calculating Cpr
  rs1<-as.vector(e_2)+as.vector(gama[i]%*%Xtrain[,i])+as.vector(gama[12]*as.vector((Xtrain[,i])^2))   #for calculating Apr
  plot(x=Xtrain[,i],y=rs,xlab=col_name[i],ylab="CPR",col="Red",cex=0.8,col.lab="Blue")
  plot(x=Xtrain[,i],y=rs1,xlab=col_name[i],ylab="APR",col="Red",cex=0.8,col.lab="Blue")
}
#----------------------2(c)---------------------------------------------------------------------------------

###ploting residual plot of lstat and rm
I_1<-diag(1,nrow=nrow(Xtrain),ncol=nrow(Xtrain))

par(mfrow=c(1,2))
for(i in c(which(colnames(Xtrain)=="rm"),which(colnames(Xtrain)=="dis"))){
hat_matrix<-Xtrain[,i]%*%solve(t(Xtrain[,i])%*%Xtrain[,i])%*%t(Xtrain[,i])
plot(y=(I_1-hat_matrix)%*%Ytrain,x=(I_1-hat_matrix)%*%Xtrain[,i],xlab=col_name[i],ylab="Partial y",col.lab="Blue",col="Red",cex=0.8)
}

###---------------------------------------2(d)--------------------------------------------------
#transform rm variable
cor(Xtrain[,which(colnames(Xtrain)=="rm")],Ytrain)
plot(x=Xtrain[,which(colnames(Xtrain)=="rm")],Ytrain,xlab=col_name[5],col="Red")
y_trans<-log(Ytrain)
X_trans<-cbind(c(rep(1,nrow(Xtrain))),Xtrain[,which(colnames(Xtrain)=="rm")])
bita<-solve(t(X_trans)%*%(X_trans))%*%t(X_trans)%*%y_trans
bita[1]
bita[2]
Xtrain_transe<-exp(bita[1])*exp(bita[2]*Xtrain[,which(colnames(Xtrain)=="rm")])
cor(Xtrain_transe,Ytrain)
Xtrain_tran<-Xtrain
Xtrain_tran[,which(colnames(Xtrain)=="rm")]<-c(Xtrain_transe)

plot(Xtrain_tran[,which(colnames(Xtrain)=="rm")],Ytrain,xlab=col_name[5],col="Red")
ytest_trans<-Xtest%*%solve(t(Xtrain_tran)%*%Xtrain_tran)%*%t(Xtrain_tran)%*%Ytrain

#----------------------------------2(e)----------------------------------------------
Rmse_trans<-sqrt((t(Ytest-ytest_trans)%*%(Ytest-ytest_trans)))/(rankMatrix(Xtrain_tran))               ## RMSE 
Rmse_trans
# here RMSE IS no decreasing---------------------------------------------------

#--------------------------------(3)------------------------------------------------------
plot(x=ytrain_hat1,y=e_1,xlab="fit",ylab="residual")
#----a u shaped pattern is observed from the graph------------------------
#again we know from the residual vs rgressor in question 2 that regressor rm has a u shaped pattern
#------There is a pattern in the error vs rm graph--------------
d=rep(0,length(e_1))
for(i in 1:length(e_1))
{d[i]=length(e_1)*(e_1[i]^2)/sum(e_1^2)}
#the graph of rm is actually flatter than u shaped , so we will use approximately err=c1+c2.x+c3.x^2+eps
X=cbind(rep(1,length(e_1)),Xtrain[,5],(Xtrain[,5])^2)
c_cap=solve(t(X)%*%X)%*%t(X)%*%e_1
# so transfomation is x----> -1.06-0.83*x+1.29*x^2
x_tr=c_cap[1]+c_cap[2]*Xtrain[,5]+c_cap[3]*(Xtrain[,5])^2
plot(x_tr,e_1,pch=20)
#graph shows approximately no relation
#again for breusch pagan test we choos our z1=x_tr
#d=alpha_0+alpha_1.z1+eps1
X1=cbind(rep(1,385),x_tr)
alpha_cap=solve(t(X1)%*%X1)%*%t(X1)%*%d
d_cap=X1%*%alpha_cap
R2=((cov(d,d_cap))^2)/(var(d)*var(d_cap))
R2
Q=385*R2
Q
#here we test at level 0.05
Q>qchisq(0.95,1) #true #as here k=1, only one regressor is changed
#there is heteroschedastticity

#sigmai^2 are not same
#so var(eps) in the original matrix is not In, let  be E,diagonal
#let sigmai^2=h(z,alpha,beta)=alpha_0+alpha_1*z1, as we have k=1
#clearly sigmai^2=d_cap[i]


E=diag(as.vector(d_cap))
betahat_updated1=solve(t(Xtrain)%*%solve(E)%*%Xtrain)%*%t(Xtrain)%*%solve(E)%*%Ytrain
betahat_updated0=b_ols1
y_updated1=Xtrain%*%betahat_updated1
e_updated1=Ytrain-y_updated1
d1=e_updated1^2/(sum((e_updated1)^2)/385)
alpha_cap0=alpha_cap
alpha_cap1=solve(t(X1)%*%X1)%*%t(X1)%*%d1
d_cap1=X1%*%alpha_cap1  #1st iteration sigmai^2=d_cap1
while (max(sum((alpha_cap0-alpha_cap1)^2)>0.0001, sum((betahat_updated0-betahat_updated1)^2))>0.0001)  {
  betahat_updated0=betahat_updated1
  alpha_cap0=alpha_cap1
  E=diag(as.vector(d_cap1))
  betahat_updated1=solve(t(Xtrain)%*%solve(E)%*%Xtrain)%*%t(Xtrain)%*%solve(E)%*%Ytrain
  y_updated1=Xtrain%*%betahat_updated1
  e_updated1=Ytrain-y_updated1
  d1=e_updated1^2/(sum((e_updated1)^2)/385)
  alpha_cap1=solve(t(X1)%*%X1)%*%t(X1)%*%d1
  d_cap1=X1%*%alpha_cap1  
  
}
beta_new=betahat_updated0
ytest_hat_new=Xtest%*%beta_new
rmse=sqrt(sum((Ytest-ytest_hat_new)^2))/10
rmse
#-----------------------------------------(4)---------------------------------------------------------------
n_updated<-nrow(Xtrain)
qp.archak=function(Z,y)    # Z=Design Matrix including 1-vector column, y=response vector
{
  n=NROW(Z)
  p=NCOL(Z)-1
  
  P=Z%*%solve(t(Z)%*%Z)%*%t(Z) 
  
  y_h=P%*%y   # pridicted value of y
  e=y-y_h   # residuals vector
  
  RSS=sum(e*e)  
  
  RSS_i=rep(0,n)  
  for(i in 1:n)
  {
    RSS_i[i]=RSS-((e[i]^2)/(1-P[i,i]))
  }
  
  r_i=rep(0,n)      #R-Student residuals
  for(i in 1:n)
  {
    r_i[i]=e[i]/(sqrt(RSS_i[i]*(1-P[i,i])/(n-p-1)))
  }
  
  r_i.ord=sort(r_i)  #ordered R-Student residuals
  p_i=rep(0,n)       #vector for population quantiles of t(n-p-1) distn which is followed by the R-Student residuals
  for(i in 1:n)
  {
    p_i[i]=qt(i/n,n-p-1)
  }
  plot(p_i,r_i.ord,main="Q-Q PLOT",xlab="Population",ylab="Sample")
  abline(a=0,b=1)  #straight line with intercept 0, slope 1 drawn to make inferences
}
qp.archak(Xtrain,Ytrain)

P=Xtrain%*%solve(t(Xtrain)%*%Xtrain)%*%t(Xtrain)
u_lambda=function(y,lambda)
{
  if(lambda!=0)
  {return ((y^lambda-1)/lambda)}
  else
  {return (log(y))}
} 
u=rep(0,n_updated)
RSS_lambda=function(lambda)
{
  for(i in 1:n_updated)
    u[i]=u_lambda(Ytrain[i],lambda)
  u_h=P%*%u
  e_u=u-u_h
  return (sum(e_u^2))
}
L=function(lambda)
{
  const=(-n_updated/2)*log(2*pi)-n_updated/2
  var=(-n_updated/2)*log(RSS_lambda(lambda)/n_updated)+(lambda-1)*sum(log(Ytrain))
  return (const+var)
}
lambda=seq(-3,3,0.08)
L_val=rep(0,length(lambda))
for(i in 1:length(lambda))
{
  L_val[i]=L(lambda[i])
}

plot(lambda,L_val,xlab="lambda",ylab="Likelihood")
lambda[L_val==max(L_val)]
y_transformed=u_lambda(Ytrain,0.36)
qp.archak(Xtrain,Ytrain)


















