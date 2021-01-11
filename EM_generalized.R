# Handmade EM4MoG ---------------------------------------------------------

?matrix
?dnorm
?rowSums

require(ramify)

d.init <- function(y, p, mu, sigma) {
  k <- length(p)
  len.y <- length(y)
  
  d <- matrix(NA, k, len.y)
  
  for(i in 1:k) {
    d[i,] <- p[i]*dnorm(y, mu[i], sigma[i])
  }
  
  return(d)
}


likelihood <- function(y, p, mu, sigma) {
  kk <- length(p)
  
  like <- 0
  for(i in 1:kk) {
    like <- like + p[i]*dnorm(y, mu[i], sigma[i])
  }
  
  return(like)
}

assign.color <- function(cols, d) {
  argm <- argmax(d, rows = FALSE)
  return(cols[argm])
}

# p is a vector of weights (pi)
handmade.em <- function(y, p, mu, sigma, n.iter, plot.flag = T)
{
  k <- length(p)
  len.y <- length(y)
  
  # set k random colors
  # cols <- randomColor(count = k, 
  #                     hue = c(" ", "random", "red", "orange", "yellow", "green", "blue", "purple", "pink", "monochrome"), 
  #                     luminosity = c(" ", "random", "light", "bright", "dark"))
  cols <- c(rgb(32/255, 74/255, 135/255, 0.7),
            rgb(204/255, 0, 0, 0.7),
            rgb(200/255, 141/255, 0, 0.7),
            rgb(78/255, 153/255, 6/255, 0.7),
            rgb(32/255, 74/255, 135/255, 0.3),
            rgb(204/255, 0, 0, 0.3),
            rgb(200/255, 141/255, 0, 0.3),
            rgb(78/255, 153/255, 6/255, 0.3))
  
  # compute the likelihood for k Gaussians
  like <- likelihood(y, p, mu, sigma)
  
  deviance <- -2*sum(log(like)) 
  res      <- matrix(NA,n.iter + 1, 3*k+2) # number of columns
  res[1,]  <- c(0, p, mu, sigma, deviance)
  
  
  d <- d.init(y, p, mu, sigma) # E step - part 1
  for (iter in 1:n.iter) {
    
    # E step - part 2
    r <- matrix(NA, k, len.y)
    for(i in 1:k) {
      r[i,] <- d[i,] / colSums(d)
    }
    
    # M step
    for(i in 1:k) {
      p[i] <- mean(r[i,])
      mu[i]    <- sum(r[i,]*y)/sum(r[i,])
      sigma[i] <- sqrt( sum(r[i,]*(y^2))/sum(r[i,]) - (mu[i])^2 )
    }
    
    # -2 x log-likelihood (a.k.a. deviance)
    like <- likelihood(y, p, mu, sigma)
    
    deviance <- -2*sum( log(like) )
    
    # Save
    res[iter+1,] <- c(iter, p, mu, sigma, deviance)
    
    # E step - part 1.1
    d <- d.init(y, p, mu, sigma) 
    
    # Plot
    if (plot.flag){
      hist(y, prob = T, breaks = 50, col = gray(.8), border = NA, 
           main = "", xlab = paste("EM Iteration: ", iter, "/", n.iter, sep = ""))
      set.seed(123)
      points(jitter(y), rep(0,length(y)), 
             pch = 19, cex = .6, 
             col = assign.color(cols, d))
      curve(likelihood(x, p, mu, sigma),
            lwd = 4, col = rgb(0,0,0,.5), add = TRUE)
      Sys.sleep(1.5)
    }
  }
  res <- data.frame(res)
  
  names(res) <- c("iteration", paste("p", 1:k, sep=""), paste("mu", 1:k, sep=""), paste("sigma", 1:k, sep=""), "deviance")
  out <- list(parameters = c(p = p, mu = mu, sigma = sigma), deviance = deviance, res = res)
  return(out)
}


# Show training MoG with k = 2 and k = 3 ----------------------------------

# require(randomcoloR)

data("faithful")
y <- faithful$eruptions

p <- c(.5,.5)
mu <- c(1.10,5.20)
sigma <- c(1.5,1.5)

handmade.em(y, p, mu, sigma, n.iter = 20, plot.flag = T)

# try with another famous dataset called "iris"
data("iris")
y2 <- iris$Sepal.Length

p2 <- c(.33, .33, .34)
mu2<- c(1.10, 2.50, 7.20)
sigma2 <- c(1.5, 1.5, 2.5)

handmade.em(y2, p2, mu2, sigma2, n.iter = 20, plot.flag = T)

# Bart samples

suppressMessages(require(mixtools, quietly = T))

n <- 250 # Sample size
XX <- rnormmix(n,
               lambda = c(0.5, rep(0.1,5)),
               mu = c(0, ((0:4)/2)-1),
               sigma = c(1, rep(0.1,5)) )

parameters <- initialize.parameters(XX, 6)
p3 <- parameters$p # c(0.40, 0.10, 0.20, 0.15, 0.05, 0.10)
mu3 <- parameters$mu # c(0, ((0:4)/2)-1)
sigma3 <- parameters$sigma # c(1, rep(0.1,5))

handmade.em(XX, p3, mu3, sigma3, n.iter = 600, plot.flag = T)


# Justify n1 < n2 .. ------------------------------------------------------

n1 <- 30 # Sample size that follows a non-asymptotic way
XX <- rnormmix(n1,
               lambda = c(0.5, rep(0.1,5)),
               mu = c(0, ((0:4)/2)-1),
               sigma = c(1, rep(0.1,5)) )

p3 <- c(0.40, 0.10, 0.20, 0.15, 0.05, 0.10)
mu3 <- c(0, ((0:4)/2)-1)
sigma3 <- c(1, rep(0.1,5))

handmade.em(XX, p3, mu3, sigma3, n.iter = 20, plot.flag = T)

n2 <- 200 # Sample size that follows a reasonably asymptotic way
XX <- rnormmix(n2,
               lambda = c(0.5, rep(0.1,5)),
               mu = c(0, ((0:4)/2)-1),
               sigma = c(1, rep(0.1,5)) )

p3 <- c(0.40, 0.10, 0.20, 0.15, 0.05, 0.10)
mu3 <- c(0, ((0:4)/2)-1)
sigma3 <- c(1, rep(0.1,5))

handmade.em(XX, p3, mu3, sigma3, n.iter = 20, plot.flag = T)


# Model selections ... ---------------------------------------------

initialize.parameters <- function(y, k) {
  p <- rep(1/k, k)
  mu    <- runif(k, min = 0.1, max= 0.9)
  sigma <- runif(k, min = 0.1, max = 0.9)
  
  return(list(p = p, mu = mu, sigma = sigma))
}

# AIC

AIC <- function(y, kmax, n.iter) {
  aic.j <- c()
  
  for (k in 1:kmax){
    parameters <- initialize.parameters(y, k)
    p <- parameters$p # runif(k,0.1,0.9)
    # p <- p/sum(p)
    
    mu <- parameters$mu # runif(k, min = -1*sd(y), max = 1*sd(y))
    sigma <- parameters$sigma # runif(k, min = 2*sd(y), max = 4*sd(y))
    
    out.opt <- handmade.em(y, p, mu, sigma, n.iter = n.iter, plot.flag = F)$parameters
    
    p <- out.opt[1:k]
    mu <- out.opt[(k+1):(2*k)]
    sigma <- out.opt[(2*k+1):(3*k)]
      
    like <- likelihood(y, p, mu, sigma)

    aic.j <- c(aic.j, 2 * sum(log(like)) - 2 * (3*k-1))
  }
  
  best.k <- which.max(aic.j) 
  return(best.k)
}

# BIC

BIC <- function(y, kmax, n.iter) {
  bic.j <- c()
  len.y <- length(y)
  
  for (k in 1:kmax){
    parameters <- initialize.parameters(y, k)
    p <- parameters$p 
    mu <- parameters$mu 
    sigma <- parameters$sigma
      
    out.opt <- handmade.em(y, p, mu, sigma, n.iter = n.iter, plot.flag = F)$parameters
    
    p <- out.opt[1:k]
    mu <- out.opt[(k+1):(2*k)]
    sigma <- out.opt[(2*k+1):(3*k)]
    
    like <- likelihood(y, p, mu, sigma)
    
    bic.j <- c(bic.j, sum(log(like)) - log(len.y)/2 * (3*k-1))
  }
  
  best.k <- which.max(bic.j) 
  return(best.k)
}

# Sample Splitting with 50% in train and 50% in test

split.data <- function(y, perc) {
  size.s <- floor(perc * length(y))
  trainIndex <- sample(seq_len(length(y)), size = size.s)
  return (trainIndex)
}

# Train model with a sample splitted at x% of training data

training.model <- function(y.train, y.test, kmax, n.iter) {
  like.j <- c()
  
  for (k in 1:kmax){
    parameters <- initialize.parameters(y.train, k)
    p <- parameters$p 
    mu <- parameters$mu 
    sigma <- parameters$sigma
    
    out.opt <- handmade.em(y.train, p, mu, sigma, n.iter = n.iter, plot.flag = F)$parameters
    
    p <- out.opt[1:k]
    mu <- out.opt[(k+1):(2*k)]
    sigma <- out.opt[(2*k+1):(3*k)]
    
    like <- likelihood(y.test, p, mu, sigma)
    
    like.j <- c(like.j, sum(log(like)))
  }
  
  best.k <- which.max(like.j) 
  return(best.k)
}

# K-fold Cross-Validation

library(caret)

k.fold.cross.validation <- function(y, k.folds, kmax, n.iter) {
  best.cross.v <- c() 
  expectations <- c()
  like.j <- c()
  
  lst.k.folds <- createFolds(y, k = k.folds, list = TRUE, returnTrain = FALSE)
  
  fold.names <- names(lst.k.folds)
  for (k in 1:kmax) {
    for (x in 1:k.folds) {
      n.k.fold <- fold.names[x]
      
      test.set <- y[unlist(lst.k.folds[n.k.fold], use.names=FALSE)]            
      train.set <- y[-unlist(lst.k.folds[n.k.fold], use.names = FALSE)]          
      
      parameters <- initialize.parameters(train.set, k)
      p <- parameters$p 
      mu <- parameters$mu 
      sigma <- parameters$sigma
      
      out.opt <- handmade.em(train.set, p, mu, sigma, n.iter = n.iter, plot.flag = F)$parameters
      
      p <- out.opt[1:k]
      mu <- out.opt[(k+1):(2*k)]
      sigma <- out.opt[(2*k+1):(3*k)]
      
      like <- likelihood(test.set, p, mu, sigma)
      
      like.j <- c(like.j, mean(log(like)))
    }
    
    expectations <- c(expectations, mean(like.j))
  }

  best.cross.v <- which.max(expectations)
  return(best.cross.v)
}

# Wasserstein method score

library(KScorrect)

wass.score <- function(y.train, y.test, kmax, n.iter) {
  f <- function(z, p, mu, sigma, y.test) {
    return(abs(qmixnorm(p = z, mean = mu, sd = sigma, pro = p) - quantile(y.test, probs = z)))
  }
  
  wass.sc <- c();
  
  for (k in 1:kmax){
    parameters <- initialize.parameters(y.train, k)
    p <- parameters$p 
    mu <- parameters$mu 
    sigma <- parameters$sigma
    
    out.opt <- handmade.em(y.train, p, mu, sigma, n.iter = n.iter, plot.flag = F)$parameters
    
    p <- out.opt[1:k]
    mu <- out.opt[(k+1):(2*k)]
    sigma <- out.opt[(2*k+1):(3*k)]
    
    wass.sc <- c(wass.sc, integrate(f, lower = 0, upper = 1,
                                    p = p, mu = mu, sigma = sigma, y.test = y.test, 
                                    rel.tol=.Machine$double.eps^.05)$value)
  }
  
  best.wass <- which.min(wass.sc)
  return(best.wass)
}

# return the max indexes of k components

max.k <- function(res) {
  k <- 0
  max <- -1
  for(i in 1:length(res)) {
    if(max <= res[i] && k < i) {
      k <- i
      max <- res[i]
    }
  }
  return(k)
}


# Define the Bart's sample

suppressMessages(require(mixtools, quietly = T))

n <- 500 # define the sample size

M <- 10 
n.iter <- 20
kmax <- 6
best.score50and50 <- c(); best.score70and30 <- c(); best.score30and70 <- c();

best.aic.k <- c(); best.bic.k <- c();

best.cross.validation5 <- c(); best.cross.validation10 <- c();

wasses <- c();
for (i in 1:M) {
  XX <- rnormmix(n,
                 lambda = c(0.5, rep(0.1,5)),
                 mu = c(0, ((0:4)/2)-1),
                 sigma = c(1, rep(0.1,5)) )
  
  # find the best k components with AIC
  best.aic.k <- c(best.aic.k, AIC(XX, kmax=kmax, n.iter))
  
  best.bic.k <- c(best.bic.k, BIC(XX, kmax=kmax, n.iter))
  
  # 50% Training-Test
  indexes50.50 <- split.data(XX, perc=0.5)
  XX.Train50 <- XX[indexes50.50]
  XX.Test50 <- XX[-indexes50.50]
  
  best.score50and50 <- c(best.score50and50, training.model(XX.Train50, XX.Test50, kmax=kmax, n.iter))
  
  # 70%-30% Training-Test
  indexes70.30 <- split.data(XX, perc=0.7)
  XX.Train70 <- XX[indexes70.30]
  XX.Test30 <- XX[-indexes70.30]
  
  best.score70and30 <- c(best.score70and30, training.model(XX.Train70, XX.Test30, kmax=kmax, n.iter))
  
  # 30%-70% Training-Test
  indexes30.70 <- split.data(XX, perc=0.3)
  XX.Train30 <- XX[indexes30.70]
  XX.Test70 <- XX[-indexes30.70]
  
  best.score30and70 <- c(best.score30and70, training.model(XX.Train30, XX.Test70, kmax=kmax, n.iter))
  
  # K-Fold Cross-Validation
  best.cross.validation5 <- c(best.cross.validation5, k.fold.cross.validation(XX, k.folds=5, kmax=kmax, n.iter))
  best.cross.validation10 <- c(best.cross.validation10, k.fold.cross.validation(XX, k.folds=5, kmax=kmax, n.iter))
  
  #Wass score
  wasses <- c(wasses, wass.score(XX.Train50, XX.Test50, kmax=kmax, n.iter))
}

paste("The best k component with the AIC method is: ", max.k(tabulate(best.aic.k)))

paste("The best k component with the BIC method is: ", max.k(tabulate(best.bic.k)))

paste("The best k component with the sample splitting at 50% and 50% is: ", max.k(tabulate(best.score50and50)))
paste("The best k component with the sample splitting at 70% and 30% method is: ", max.k(tabulate(best.score70and30)))
paste("The best k component with the sample splitting at 30% and 70% method is: ", max.k(tabulate(best.score30and70)))

paste("The best k component with the k-fold 5 cross-validation method is: ", max.k(tabulate(best.cross.validation5)))
paste("The best k component with the k-fold 10 cross-validation method is: ", max.k(tabulate(best.cross.validation10)))

paste("The best k component with the wass score method is: ", max.k(tabulate(wasses)))
