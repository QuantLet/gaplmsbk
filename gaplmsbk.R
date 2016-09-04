# Spline-backfitted kernel estimation for component function in Generalized Addtive Partially Linear Model
# x: n*p matrix
# y: n*1 vector
# alpha: the component function to be estimated
# xfixed: the x values for the estimation of component function
# initial: initial value for the component function to be estimated
# c1,c2: adjust the number of knots for spline
# ch: adjust the bandwidth of kernel
# family: provide a convenient way to specify the details of the models used, same as as glm in R.

gplmsbk <- function(x,t,y,alpha,xfixed,initial,c1,c2,ch,family,thetaalphaoracle) {

library("splines")


d  = ncol(x)
n  = nrow(x)
dt = ncol(t)
B  = array(1,c(n,1))
N  = floor(c1*(n^(1/4))*log(n))+c2    #number of knots
N  = min(N,floor((n/4-1)/d-1))

i  = 1
while (i<=d) {
  Bi = bs(x[,i],knots = min(x[,i])+(max(x[,i])-min(x[,i]))*c(1:(N-1))/N,degree = 1)
  B  = cbind(B,Bi)	
	i = i+1
}

Nd     = ncol(B)
result = glm(y~cbind(B[,2:Nd],t),family)
lamda  = data.matrix(result$coefficients)
beta   = lamda[(Nd+1):(Nd+dt)]
lamda  = lamda[1:Nd]
mhat   = data.matrix(lamda[1]*B[,1])
chat   = lamda[1]

i = 1
while (i<=d) {
	mhati = B[,((i-1)*N+2):(i*N+1)]%*%lamda[((i-1)*N+2):(i*N+1)]
	chat  = chat + data.matrix(colMeans(mhati))
	mhat  = cbind(mhat,mhati-colMeans(mhati))
	i     = i+1
}

mhat[,1] = chat
xfixed   = t(xfixed)
dfixed   = ncol(xfixed)
initial  = t(initial)
mhatsbk  = initial

xalpha   = data.matrix(x[,alpha])
hoptalpha= sqrt(var(xalpha))*(n^(-0.2))
halpha   = ch*hoptalpha
ualpha   = xalpha%*%array(1,c(1,dfixed))-array(1,c(n,1))%*%xfixed
kalpha   = 15/16*((1-(ualpha/halpha[1])^2+abs(1-(ualpha/halpha[1])^2))/2)^2/halpha[1]

theta      = rowSums(mhat)+t%*%beta
mhatalpha  = mhat[,alpha+1]
thetaalpha = theta-mhatalpha
difference = 1
delta      = 0.0000001

i = 1
while (i<=20 & difference>delta) {
	thetanew   = thetaalpha%*%array(1,c(1,dfixed))+array(1,c(n,1))%*%mhatsbk
	lderiv1    = colSums((y%*%array(1,c(1,dfixed))-1/(1+exp(-thetanew)))*kalpha)
	lderiv2    = colSums(-((exp(-thetanew)/((1+exp(-thetanew))^2)))*kalpha)
	mhatsbki   = lderiv1/lderiv2
	mhatsbk    = mhatsbk-mhatsbki
	difference = max(abs(t(mhatsbki)))
	i          = i+1
}

mhatsbkoracle    = initial
differenceoracle = 1
deltaoracle      = 0.0000001

i = 1
while (i<=20 & differenceoracle>deltaoracle) {
	thetaneworacle = thetaalphaoracle%*%array(1,c(1,dfixed))+array(1,c(n,1))%*%mhatsbkoracle
	lderiv1oracle  = colSums((y%*%array(1,c(1,dfixed))-1/(1+exp(-thetaneworacle)))*kalpha)
	lderiv2oracle  = colSums(-((exp(-thetaneworacle)/((1+exp(-thetaneworacle))^2)))*kalpha)
	mhatsbkioracle = lderiv1oracle/lderiv2oracle
	mhatsbkoracle = mhatsbkoracle-mhatsbkioracle
	differenceoracle = max(abs(t(mhatsbkioracle)))
	i = i+1
}
mhatsbkalpha=cbind(t(mhatsbk),t(mhatsbkoracle))
return(mhatsbkalpha)

}
