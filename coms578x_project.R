
#############################################################
# Two norm square
# y^T y 
# @param column vector y
# @return Scalar value
##############################################################
norm <- function(y) {
  r = t(y)%*%y
  r = as.numeric(r)
  return (r)
}
#############################################################
# Backtracking step
# @param x0: current position
# @param func: Function subroutine
# @param derv: derivative function subroutine
# @param alpha: between 0 and 0.5
# @param beta: between 0 and 1
# @param iteration: count
##############################################################
backtrack<-function(x0, fn, derv, s, alpha, beta, count) {
  while (fn(x0 - s*derv(x0)) > fn(x0) -  alpha * s * norm(derv(x0)) ) {
    s = s*beta
    count = count+1
  }
  return (s)
}

#############################################################
# Gradient Descent with fixed step size of 1/sqrt(k)
# where k is the iteration number
# @param func: Function subroutine
# @param derv: derivative function subroutine
# @param start: start position
# @param s : step size
##############################################################
gradient_descent<-function(func, derv, start, s, tol=0.001, maxiters=5000){
  x1 = start
  k=1
  ###########
  ftrace <- c()
  xtrace <- list()
  ftrace <- rbind(ftrace, func(x1))
  xtrace[[k]] <- x1
  k=k+1
  s=0.029
  ###########
  while (norm(derv(x1)) > tol) {
    #s = 1/sqrt(k)
    df <- derv(x1) #cost of finding gradient is O(n)
    x1 <- x1 - s * df
    ##########
    ftrace <- rbind(ftrace, func(x1))
    xtrace[[k]] <- x1
    ############
    k=k+1
    if (k==maxiters) {
      break;
    }
  }
  
  data.frame(
    "f_x" = ftrace,
    "x" = I(xtrace),
    "k" = c(1:length(ftrace))
  )
}


#############################################################
# Stochastic Gradient Descent
# @param func: Function subroutine
# @param derv: derivative function subroutine
# @param start: start position
# @param s : step size
##############################################################
stochastic_gradient_descent<-function(func, derv, derv_i, start, n, m=1, tol=0.001, maxiters=1000){
  x1 = start
  k=1
  ###########
  ftrace <- c()
  vtrace <- c()
  xtrace <- list()
  ftrace <- rbind(ftrace, func(x1))
  vtrace<-rbind(vtrace, var(derv(x1)))
  xtrace[[k]] <- x1
  k=k+1
  ###########
  #s=0.03
  while (norm(derv(x1)) > tol) {
    s = 1/sqrt(k)
    
    ik = sample(1:n,m,replace=F)
    #for (i in index) { #cost of finding gradient is O(m) where m<n
    x1 <- x1 - s * derv(x1,ik) #cost of finding gradient is O(1)
    #}
    ##########
    ftrace <- rbind(ftrace, func(x1))
    vtrace<-rbind(vtrace, var(derv(x1)))
    #print(func(x1))
    xtrace[[k]] <- x1
    ############
    k=k+1
    if (k==maxiters) {
      break;
    }
  }
  
  data.frame(
    "f_x" = ftrace,
    "v_x" = ftrace,
    "x" = I(xtrace),
    "k" = c(1:length(ftrace))
  )
}

#############################################################
# Stochastic Gradient Descent for SAG start value tuning
# @param func: Function subroutine
# @param derv: derivative function subroutine
# @param start: start position
# @param s : step size
##############################################################
stochastic_gradient_descent_sag<-function(func, derv, derv_i, start, n, m=1, tol=0.001, maxiters=1000){
  x1 = start
  k=1
  ###########
  ftrace <- c()
  xtrace <- list()
  ftrace <- rbind(ftrace, func(x1))
  xtrace[[k]] <- x1
  k=k+1
  ###########
  
  while (k<maxiters) {
    s = 1/sqrt(k)
    
    ik = sample(1:n,m,replace=F)
    #for (i in index) { #cost of finding gradient is O(m) where m<n
    x1 <- x1 - s * derv(x1,ik) #cost of finding gradient is O(1)
    #}
    ##########
    ftrace <- rbind(ftrace, func(x1))
    #print(func(x1))
    xtrace[[k]] <- x1
    ############
    k=k+1
  }
  
  data.frame(
    "f_x" = ftrace,
    "x" = I(xtrace),
    "k" = c(1:length(ftrace))
  )
}

#############################################################
# Stochastic Average Gradient: SAG
# @param func: Function subroutine
# @param derv: derivative function subroutine
# @param start: start position
# @param s : step size
##############################################################
sgd_sag<-function(func, derv, derv_i, start, n, m=1, s = 0.01, tol=0.01, maxiters=2000){
  x1 = start
  k=1
  g = list()
  for  (i in 1:n) {
    g[[i]] <- x1
  }
  ###########
  ftrace <- c()
  vtrace <- c()
  xtrace <- list()
  ftrace <- rbind(ftrace, func(x1))
  vtrace<-rbind(vtrace, var(derv(x1)))
  xtrace[[k]] <- x1
  k=k+1
  ###########
  
  while (norm(derv(x1)) > tol) {
    ik = sample(1:n,m,replace=T)
    
    sumg0<- matrix(rep(0, length.out=length(x1)), nrow=length(x1),ncol=1)
    for  (i in 1:n) {
      sumg0 = sumg0+g[[i]]
    }
    g0 = g[[ik]]
    g1 = derv_i(ik, x1)
    g[[ik]]=g1
    
    #for (i in index) { #cost of finding gradient is O(m) where m<n
    x1 <- x1 - s * (g1 - g0 + sumg0) #cost of finding gradient is O(1)
    #}
    ##########
    ftrace <- rbind(ftrace, func(x1))
    vtrace<-rbind(vtrace, var(derv(x1)))
    #print(abs(func(x1)-fstar))
    #print(func(x1))
    #print(var(derv(x1)))
    xtrace[[k]] <- x1
    ############
    k=k+1
    if (k==maxiters) {
      break;
    }
  }
  
  data.frame(
    "f_x" = ftrace,
    "v_x" = ftrace,
    "x" = I(xtrace),
    "k" = c(1:length(ftrace))
  )
}

#############################################################
# Stochastic Gradient Descent: SAGA
# @param func: Function subroutine
# @param derv: derivative function subroutine
# @param start: start position
# @param s : step size
##############################################################
sgd_saga<-function(func, derv, derv_i, start, n, m=1, s=0.01, tol=0.001, maxiters=1000){
  x1 = start
  k=1
  g = list()
  for  (i in 1:n) {
    g[[i]] <- x1
  }
  
  ###########
  ftrace <- c()
  vtrace <- c()
  xtrace <- list()
  ftrace <- rbind(ftrace, func(x1))
  vtrace<-rbind(vtrace, var(derv(x1)))
  xtrace[[k]] <- x1
  k=k+1
  ###########
  
  while (norm(derv(x1)) > tol) {
    #s = 0.01
    ik = sample(1:n,m,replace=F)
    
    sumg0<- matrix(rep(0, length.out=length(x1)), nrow=length(x1),ncol=1)
    for  (i in 1:n) {
      sumg0 = sumg0+g[[i]]
    }
    g0 = g[[ik]]
    g1 = derv_i(ik, x1)
    g[[ik]]=g1
    
    #for (i in index) { #cost of finding gradient is O(m) where m<n
    x1 <- x1 - s * (g1*n - g0*n + sumg0) #cost of finding gradient is O(1)
    #}
    ##########
    ftrace <- rbind(ftrace, func(x1))
    vtrace<-rbind(vtrace, var(derv(x1)))
    #print(abs(func(x1)-fstar))
    #print(func(x1))
    #print(var(derv(x1)))
    xtrace[[k]] <- x1
    ############
    k=k+1
    if (k==maxiters) {
      break;
    }
  }
  
  data.frame(
    "f_x" = ftrace,
    "v_x" = ftrace,
    "x" = I(xtrace),
    "k" = c(1:length(ftrace))
  )
}

#############################################################
# Stochastic Gradient Descent: SVRG
# @param func: Function subroutine
# @param derv: derivative function subroutine
# @param start: start position
# @param s : step size
##############################################################
sgd_svrg<-function(func, derv, derv_i, start, n, m=1, s=0.01, tol=0.001, maxiters=1000){
  x1 = start
  k=1
  ###########
  ftrace <- c()
  vtrace <- c()
  xtrace <- list()
  ftrace <- rbind(ftrace, func(x1))
  vtrace<-rbind(vtrace, var(derv(x1)))
  xtrace[[k]] <- x1
  ###########
  k=k+1
  #s=0.029
  while (norm(derv(x1)) > tol) {
    #s = 1/sqrt(k)
    g0 = derv(x1)
    ibatch = sample(1:n, m, replace = F)
    x0 = x1
    for (i in ibatch) {
      x1 <- x1 - s*(n*derv_i(i, x1) - n*derv_i(i, x0) + g0)
    }
    ##########
    ftrace <- rbind(ftrace, func(x1))
    #print(abs(func(x1)-fstar))
    vtrace<-rbind(vtrace, var(derv(x1)))
    xtrace[[k]] <- x1
    ############
    k=k+1
    if (k==maxiters) {
      break;
    }
  }
  
  data.frame(
    "f_x" = ftrace,
    "v_x" = vtrace,
    "x" = I(xtrace),
    "k" = c(1:length(ftrace))
  )
}

###########################################################################
#########################################
# Function: squared error cost function
#########################################
cost <- function(theta) {
  c <- 0
  for (i in 1:length(y)) {
    c<- c+cost_i(i,theta)
  }
  c
}

cost_i <- function(i, theta) {
  #print(c(i,X[i,],theta))
  ( (X[i,] %*% theta - y[i])^2 ) / (2*length(y))
}

#######################################################
# Gradient: Derivative of squared error cost function
#######################################################
delta <- function(theta, B=c()) {
  d <- matrix(rep(0,length.out=length(theta)),nrow=length(theta),ncol=1)
  if (length(B)==0) {
    B = 1:length(y)
  }
  for (i in B) {
    d <- d + delta_i(i,theta, B)
  }
  d
}

delta_i <- function(i, theta, B=c()) {
  if (length(B)==0) {
    N = length(y)
  }
  else {
    N = length(B)
  }
  error<-as.numeric(X[i,] %*% theta - y[i])
  t(error*t(X[i,])/N)
}
#########################################################################

#Simulate data
x <- runif(100, -5, 5) #try 100 works well
y <- x + rnorm(100) + 3

# fit a linear model
res <- lm( y ~ x )
print(res$coefficients)
# plot the data and the model
plot(x,y, col="grey")
abline(res, col='blue')

# learning rate and iteration limit
s <- 0.01

# initialize coefficients
beta0 <- matrix(c(0,0), nrow=2)

# add a column of 1's for the intercept coefficient
X <- cbind(1, matrix(x))

thetastar = matrix(c(res$coefficients),ncol=1,nrow=2)
fstar = cost(theta = thetastar)
#norm(delta(thetastar))
#####################################################
df = data.frame(
  method=character(),
  time = double(),
  steps = integer(),
  est_beta0 = double(),
  est_beta0 = double()
)


df = rbind(df, data.frame(method="LR", time=0, steps=0, beta0=round(res$coefficients[[1]],4), beta1=round(res$coefficients[[2]],4) ) )

a=Sys.time()
gdf = gradient_descent(cost,delta,beta0+c(2.5,0.5), s)
b=Sys.time()-a
max(gdf$k)
beta_est <- c(gdf$x[length(gdf$k)][[1]])
beta_est
df = rbind(df, data.frame(method="GD", time=round(b,4), steps=max(gdf$k), beta0=round(beta_est[1],4), beta1=round(beta_est[2],4) ) )
y1<-gdf

a=Sys.time()
gdf = stochastic_gradient_descent(cost,delta, delta_i, beta0+c(2.5,0.5), length(y))
b=Sys.time()-a
max(gdf$k)
beta_est <- c(gdf$x[length(gdf$k)][[1]])
beta_est
df = rbind(df, data.frame(method="SGD-1", time=round(b,4), steps=max(gdf$k), beta0=round(beta_est[1],4), beta1=round(beta_est[2],4) ) )
y2<-gdf

a=Sys.time()
gdf = stochastic_gradient_descent(cost,delta, delta_i, beta0, length(y), 10)
b=Sys.time()-a
max(gdf$k)
beta_est <- c(gdf$x[length(gdf$k)][[1]])
beta_est
df = rbind(df, data.frame(method="SGD-10", time=round(b,4), steps=max(gdf$k), beta0=round(beta_est[1],4), beta1=round(beta_est[2],4) ) )
y3<-gdf

a=Sys.time()
gdf = stochastic_gradient_descent(cost,delta, delta_i, beta0, length(y), 100)
b=Sys.time()-a
b
max(gdf$k)
beta_est <- c(gdf$x[length(gdf$k)][[1]])
beta_est
df = rbind(df, data.frame(method="SGD-500", time=round(b,4), steps=max(gdf$k), beta0=round(beta_est[1],4), beta1=round(beta_est[2],4) ) )
y4<-gdf
df

#sgdburn <- stochastic_gradient_descent_sag(cost,delta, delta_i, beta0, length(y))
#beta0_sag <- sgdburn$x[[length(y)-1]]

a=Sys.time()
gdf = sgd_sag(cost,delta, delta_i, beta0, length(y), m=1, s=0.029, tol=0.001) #step size as high as possible before diverging
b=Sys.time()-a
beta_est <- c(gdf$x[length(gdf$k)][[1]])
beta_est
max(gdf$k)
df = rbind(df, data.frame(method="SAG", time=round(b,4), steps=max(gdf$k), beta0=round(beta_est[1],4), beta1=round(beta_est[2],4) ) )
y5<-gdf$f_x

a=Sys.time()
gdf = sgd_saga(cost,delta, delta_i, beta0, length(y), m=1, s=0.03)
b=Sys.time()-a
beta_est <- c(gdf$x[length(gdf$k)][[1]])
beta_est
max(gdf$k)
df = rbind(df, data.frame(method="SAGA", time=round(b,4), steps=max(gdf$k), beta0=round(beta_est[1],4), beta1=round(beta_est[2],4) ) )
y6<-gdf$f_x
df
a=Sys.time()
gdf = sgd_svrg(cost,delta, delta_i, beta0, length(y), m=5, s=0.02, tol=0.001)
b=Sys.time()-a
beta_est <- c(gdf$x[length(gdf$k)][[1]])
beta_est
max(gdf$k)
df = rbind(df, data.frame(method="SVRG", time=round(b,4), steps=max(gdf$k), beta0=round(beta_est[1],4), beta1=round(beta_est[2],4) ) )
y6<-gdf$f_x
df

#Plot variance updates for SVRG with changes to update frequency
y4 = sgd_svrg(cost,delta, delta_i, beta0, length(y), m=20, s=0.02, tol=0.001)
y3 = sgd_svrg(cost,delta, delta_i, beta0, length(y), m=15, s=0.02, tol=0.001)
y2 = sgd_svrg(cost,delta, delta_i, beta0, length(y), m=10, s=0.02, tol=0.001)
y1 = sgd_svrg(cost,delta, delta_i, beta0, length(y), m=5, s=0.02, tol=0.001)
plot(y4$v_x, type="l", lty=1, col="black", xlab = "iteration", ylab="variance of gradient")
lines(y3$v_x, lty=1, col="red")
lines(y2$v_x, lty=1, col="green")
lines(y1$v_x, lty=1, col="blue")
legend("topright", c("m=20","m=15","m=10","m=5"), col=c("black","red", "green","blue"), lty=c(1,1,1,1))

#Plot variance updates for SGD and compare with mini-batch setting
y3 = stochastic_gradient_descent(cost,delta, delta_i, beta0, length(y), 100)
y2 = stochastic_gradient_descent(cost,delta, delta_i, beta0, length(y), 10)
y1 = stochastic_gradient_descent(cost,delta, delta_i, beta0, length(y))
plot(y3$v_x, type="l", lty=1, col="red", xlab = "iteration", ylab="variance of gradient")
lines(y2$v_x, lty=1, col="green")
lines(y2$v_x, lty=1, col="blue")
legend("topright", c("SGD-100","SGD-10","SGD"), col=c("red","green", "blue"), lty=c(1,1,1))

#Plot variance updates for SGD and compare with mini-batch setting
y4 = sgd_svrg(cost,delta, delta_i, start = beta0, n=length(y), m=5, s=0.02, tol=0.001)
y3 = sgd_sag(cost,delta, delta_i, beta0, length(y), m=1, s=0.029, tol=0.001) 
y2 = sgd_saga(cost,delta, delta_i, beta0, length(y), m=1, s=0.03, tol=0.001)
y1 = stochastic_gradient_descent(cost,delta, delta_i, beta0, length(y), m=1)
plot(y4$v_x, type="l", lty=1, col="black", xlab = "iteration", ylab="variance of gradient")
lines(y3$v_x, lty=2, col="red")
lines(y2$v_x, lty=3, col="green")
lines(y1$v_x, lty=4, col="blue")
legend("topright", c("SVRG-5","SAG","SAGA","SGD"), col=c("black","red", "green","blue"), lty=c(1,2,3,4))

y3 = sgd_sag(cost,delta, delta_i, beta0, length(y), m=1, s=0.01, tol=0.001) 
y2 = sgd_saga(cost,delta, delta_i, beta0, length(y), m=1, s=0.03, tol=0.001)
y1 = stochastic_gradient_descent(cost,delta, delta_i, beta0, length(y), m=1)
plot(y3$v_x, type="l", lty=1, col="black", xlab = "iteration", ylab="variance of gradient")
lines(y2$v_x, lty=2, col="red")
lines(y1$v_x, lty=4, col="blue")
legend("topright", c("SAG","SAGA","SGD"), col=c("black","red", "blue"), lty=c(1,2,4))
