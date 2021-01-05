# Handmade EM4MoG ---------------------------------------------------------

?matrix
?dnorm
?rowSums

require(ramify)

d_init <- function(y, p, mu, sigma) {
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

assign_color <- function(cols, d) {
  argm <- argmax(d, rows = FALSE)
  return(cols[argm])
}

# p is a vector of weights (pi)
handmade.em <- function(y, p, mu, sigma, n_iter, plot_flag = T)
{
  k <- length(p)
  len.y <- length(y)
  
  # set k random colors
  cols <- randomColor(count = k, 
                      hue = c(" ", "random", "red", "orange", "yellow", "green", "blue", "purple", "pink", "monochrome"), 
                      luminosity = c(" ", "random", "light", "bright", "dark"))
  
  # compute the likelihood for k Gaussians
  like <- likelihood(y, p, mu, sigma)
  
  deviance <- -2*sum(log(like)) 
  res      <- matrix(NA,n_iter + 1, 3*k+2) # number of columns
  res[1,]  <- c(0, p, mu, sigma, deviance)
  
  
  d <- d_init(y, p, mu, sigma) # E step - part 1
  for (iter in 1:n_iter) {
    
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
    d <- d_init(y, p, mu, sigma) 
    
    # Plot
    if (plot_flag){
      hist(y, prob = T, breaks = 50, col = gray(.8), border = NA, 
           main = "", xlab = paste("EM Iteration: ", iter, "/", n_iter, sep = ""))
      set.seed(123)
      points(jitter(y), rep(0,length(y)), 
             pch = 19, cex = .6, 
             col = assign_color(cols, d))
      curve(likelihood(x, p, mu, sigma),
            lwd = 4, col = rgb(0,0,0,.5), add = TRUE)
      Sys.sleep(1)
    }
  }
  res <- data.frame(res)
  
  names(res) <- c("iteration", paste("p", 1:k, sep=""), paste("mu", 1:k, sep=""), paste("sigma", 1:k, sep=""), "deviance")
  out <- list(parameters = c(p = p, mu = mu, sigma = sigma), deviance = deviance, res = res)
  return(out)
}


# Show training MoG with k = 2 and k = 3 ----------------------------------

require(randomcoloR)

data("faithful")
y <- faithful$eruptions

p <- c(.5,.5)
mu <- c(1.10,5.20)
sigma <- c(1.5,1.5)

handmade.em(y, p, mu, sigma, n_iter = 20, plot_flag = T)

# try with another famous dataset called "iris"
data("iris")
y2 <- iris$Sepal.Length

p2 <- c(.33, .33, .34)
mu2<- c(1.10, 2.50, 7.20)
sigma2 <- c(1.5, 1.5, 2.5)

handmade.em(y2, p2, mu2, sigma2, n_iter = 20, plot_flag = T)

# Bart samples

suppressMessages(require(mixtools, quietly = T))

n <- 250 # Sample size
XX <- rnormmix(n,
               lambda = c(0.5, rep(0.1,5)),
               mu = c(0, ((0:4)/2)-1),
               sigma = c(1, rep(0.1,5)) )

p3 <- c(0.40, 0.10, 0.20, 0.15, 0.05, 0.10)
mu3 <- c(0, ((0:4)/2)-1)
sigma3 <- c(1, rep(0.1,5))

handmade.em(XX, p3, mu3, sigma3, n_iter = 20, plot_flag = T)


# Justify n1 < n2 .. ------------------------------------------------------

n1 <- 30 # Sample size that follows a non-asymptotic way
XX <- rnormmix(n1,
               lambda = c(0.5, rep(0.1,5)),
               mu = c(0, ((0:4)/2)-1),
               sigma = c(1, rep(0.1,5)) )

p3 <- c(0.40, 0.10, 0.20, 0.15, 0.05, 0.10)
mu3 <- c(0, ((0:4)/2)-1)
sigma3 <- c(1, rep(0.1,5))

handmade.em(XX, p3, mu3, sigma3, n_iter = 20, plot_flag = T)

n2 <- 200 # Sample size that follows a reasonably asymptotic way
XX <- rnormmix(n2,
               lambda = c(0.5, rep(0.1,5)),
               mu = c(0, ((0:4)/2)-1),
               sigma = c(1, rep(0.1,5)) )

p3 <- c(0.40, 0.10, 0.20, 0.15, 0.05, 0.10)
mu3 <- c(0, ((0:4)/2)-1)
sigma3 <- c(1, rep(0.1,5))

handmade.em(XX, p3, mu3, sigma3, n_iter = 20, plot_flag = T)


# Model selections ... ---------------------------------------------

# AIC

AIC <- function(y, kmax) {
  aic_j <- c()
  
  for (k in 1:kmax){
    p <- runif(k)
    p <- p/sum(p)
    
    mu <- runif(k, min = -1*sd(y), max = 1*sd(y))
    sigma <- runif(k, min = 2*sd(y), max = 4*sd(y))
    
    out_opt <- handmade.em(y, p, mu, sigma, n_iter = k*100, plot_flag = F)$parameters
    
    p <- out_opt[1:k]
    mu <- out_opt[(k+1):(2*k)]
    sigma <- out_opt[(2*k+1):(3*k)]
      
    like <- likelihood(y, p, mu, sigma)

    aic_j <- c(aic_j, 2 * sum(log(like)) - 2 * k)
  }
  
  best_k <- which.max(aic_j) 
  return(best_k)
}

# BIC

# Train model with a sample splitted at x% of training data

training_model <- function(y_train, y_test, kmax) {
  aic_j <- c()
  
  for (k in 1:kmax){
    p <- runif(k)
    p <- p/sum(p)
    
    mu <- runif(k, min = -1*sd(y_train), max = 1*sd(y_train))
    sigma <- runif(k, min = 2*sd(y_train), max = 4*sd(y_train))
    
    out_opt <- handmade.em(y_train, p, mu, sigma, n_iter = k*100, plot_flag = F)$parameters
  }
  
  p <- out_opt[1:k]
  mu <- out_opt[(k+1):(2*k)]
  sigma <- out_opt[(2*k+1):(3*k)]
  
  like <- likelihood(y_test, p, mu, sigma)
  
  aic_j <- c(aic_j, 2 * sum(log(like)) - 2 * k)
  
  best_k <- which.max(aic_j) 
  return(best_k)
}

# Sample Splitting with 50% in train and 50% in test

split_50_50 <- function(y) {
  size_s <- floor(0.50 * length(y))
  trainIndex <- sample(seq_len(length(y)), size = size_s)
  return (trainIndex)
}

# Sample Splitting with 70% in train and 30% in test

split_70_30 <- function(y) {
  size_s <- floor(0.70 * length(y))
  trainIndex <- sample(seq_len(length(y)), size = size_s)
  return (trainIndex)
}

# Sample Splitting with 30% in train and 70% in test

split_30_70 <- function(y) {
  size_s <- floor(0.30 * length(y))
  trainIndex <- sample(seq_len(length(y)), size = size_s)
  return (trainIndex)
}

# Define the Bart's sample

suppressMessages(require(mixtools, quietly = T))

n <- 250 # define the sample size

M <- 10 
best_score50and50 <- c(); best_score70and30 <- c(); best_score30and70 <- c(); best_aic_k <- c();
for (i in 1:M) {
  XX <- rnormmix(n,
                 lambda = c(0.5, rep(0.1,5)),
                 mu = c(0, ((0:4)/2)-1),
                 sigma = c(1, rep(0.1,5)) )
  
  # 50% Training-Test
  indexes50_50 <- split_50_50(XX)
  XX_Train50 <- XX[indexes50_50]
  XX_Test50 <- XX[-indexes50_50]
  
  best_score50and50 <- c(best_score50and50, training_model(XX_Train50, XX_Test50, 6))
  
  # 70%-30% Training-Test
  indexes70_30 <- split_70_30(XX)
  XX_Train70 <- XX[indexes70_30]
  XX_Test30 <- XX[-indexes70_30]
  
  best_score70and30 <- c(best_score70and30, training_model(XX_Train70, XX_Test30, 6))
  
  # 30%-70% Training-Test
  indexes30_70 <- split_30_70(XX)
  XX_Train30 <- XX[indexes30_70]
  XX_Test70 <- XX[-indexes30_70]
  
  best_score30and70 <- c(best_score30and70, training_model(XX_Train30, XX_Test70, 6))
  
  # find the best k components with AIC
  best_aic_k <- c(best_aic_k, AIC(XX, 6))
}

paste("The best k component with the AIC method is: ", which.max(tabulate(best_aic_k)))

paste("The best k component with the sample splitting at 50% and 50% is: ", which.max(tabulate(best_score50and50)))
paste("The best k component with the sample splitting at 70% and 30% method is: ", which.max(tabulate(best_score70and30)))
paste("The best k component with the sample splitting at 30% and 70% method is: ", which.max(tabulate(best_score30and70)))
