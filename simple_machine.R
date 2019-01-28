# Note!! The code has to be run in chunks and some chunks has to be ran
# multiple times to run successfully. Sometimes you get unlucky and get numerical problems,
# but not often.

# packages
library(ggplot2)
#library(MASS)  # multivariate normal distributions
library(mvtnorm) # multivariate normal distributions
library(nimble)  # inverse-gamma distributions
library(truncdist)  # truncated distritubtions
library(fields)  # calculating euclidean distance
library(gridExtra)
library(grid)
library(plyr)

# Importing dataset
data_simple_machine <- read.csv("~/Documents/R/Prosjektoppgaven/data_simple_machine.txt", header=FALSE)


# Assigning dataset
X11 <- as.double(data_simple_machine[1,!is.na(data_simple_machine[1,])])
X31 <- as.double(data_simple_machine[2,!is.na(data_simple_machine[2,])])
X61 <- as.double(data_simple_machine[3,!is.na(data_simple_machine[3,])])
true11 <- as.double(data_simple_machine[4,!is.na(data_simple_machine[4,])])
true31 <- as.double(data_simple_machine[5,!is.na(data_simple_machine[5,])])
true61 <- as.double(data_simple_machine[6,!is.na(data_simple_machine[6,])])
interpolation_x <- as.double(data_simple_machine[7,!is.na(data_simple_machine[7,])])
interpolation_zeta <- as.double(data_simple_machine[8,!is.na(data_simple_machine[8,])])
extrapolation_x <- 6
extrapolation_zeta <- 3
# Errors genereated by zero mean normal dustribution with variance 0.01^2
error11 <- c(5.836698e-05, -1.294391e-04, 8.636004e-05, 2.545971e-05, -1.899284e-04, 1.973852e-05, 1.626383e-04,-4.632656e-05, 1.047955e-05,-1.613296e-05, 6.475140e-05)
error31 <- c(-4.186443e-05,-9.484696e-05,-2.630421e-04, 1.979543e-05,-1.340001e-04, 7.022427e-05,-1.829711e-05,-1.725726e-05, 8.821428e-05,-6.535563e-05, 3.569571e-05, 1.084952e-04,
             1.666466e-04, 1.690780e-04,-1.743096e-04, 7.644912e-05, 1.561225e-05, 3.874042e-05,-4.357500e-05,-4.972623e-06, 1.487953e-04,-2.183937e-05, 7.796819e-05, 1.090586e-05,
             -2.449725e-06,-1.031338e-04, 1.342935e-04, 2.182190e-05,-7.859128e-05, 3.167328e-05,-6.178424e-05)
error61 <- c(1.135047e-05,-1.353582e-04, 1.998339e-05, 6.537765e-05, 8.548965e-05,-9.822860e-05, 1.050306e-04,-1.852564e-05,-9.424297e-05, 4.670077e-05,-1.302927e-05,-8.294267e-05,
             1.150562e-05, 1.276892e-05, 8.515714e-05, 1.058882e-04,-1.877829e-05,-9.967227e-06 ,1.247205e-04, 5.449387e-05,-3.027834e-06, 2.076499e-04,-2.841888e-05,-1.452947e-04,
             -4.399744e-05, 1.504915e-04, 1.197973e-04,-1.021390e-04, 7.924050e-05, 7.449373e-05,-2.522495e-05, 7.893844e-05, 4.634988e-05, 5.132954e-05,-5.177323e-05,-1.175639e-04,
             4.903257e-06,-1.635492e-04, 7.889011e-05,-2.409988e-05,-1.266366e-08, 1.485976e-04,-8.939738e-05, 1.727700e-04, 1.501395e-04,-1.241436e-04, 1.544317e-04, 1.094756e-04,
             -3.972290e-05, 1.851934e-04,-1.706003e-04,-9.863664e-05,-2.091363e-06,-9.621987e-05, 1.030752e-04, 1.413960e-04, 2.927164e-06,-2.497293e-04,-1.343687e-04, 7.243032e-05,
             -3.356590e-05)
Z11 <- true11 + error11
Z31 <- true31 + error31
Z61 <- true61 + error61


# Bayesian linear regression
blr <- function(X, Z, interpolation, extrapolation){
  n = length(X)
  theta_hat <- solve(t(X) %*% X) %*% t(X) %*% Z
  S2 <- (1/(n-1))%*%t(Z-X%*%theta_hat)%*%(Z-X%*%theta_hat)
  
  mean_theta <- (sqrt(S2)/sqrt(t(X) %*% X)) %*% qt(0.5,n-1) + theta_hat
  lower_bound_theta <- (sqrt(S2)/sqrt(t(X) %*% X)) %*% qt(0.05,n-1) + theta_hat
  upper_bound_theta <- (sqrt(S2)/sqrt(t(X) %*% X)) %*% qt(0.95,n-1) + theta_hat
  interpolation <- c(interpolation,1.5)
  mean_ip <- as.numeric(mean_theta)*interpolation
  lower_bound_ip <- as.numeric(lower_bound_theta)*interpolation
  upper_bound_ip <- as.numeric(upper_bound_theta)*interpolation
  mean_ep <- as.numeric(mean_theta)*extrapolation
  lower_bound_ep <- as.numeric(lower_bound_theta)*extrapolation
  upper_bound_ep <- as.numeric(upper_bound_theta)*extrapolation
  
  
  
  # calculating values for plotting the posterior density of theta, interpolation and extrapolation
  points <- 100
  x_tdistribution <- seq(qt(0.001,n-1),qt(0.999,n-1),length = points)
  plot_theta_x <- ( sqrt(S2) / c(sqrt(t(X) %*% X)) ) %*% (x_tdistribution) + c(theta_hat)
  plot_theta_y <- dt(x_tdistribution,n-1) * (qt(0.001,n-1)-qt(0.999,n-1))/((plot_theta_x[1] - plot_theta_x[points]) )
  plot_ipx <- matrix(double(length(interpolation) * points), nrow = length(interpolation))
  plot_ipy <- matrix(double(length(interpolation) * points), nrow = length(interpolation))
  plot_epx <- matrix(double(length(extrapolation) * points), nrow = length(extrapolation))
  plot_epy <- matrix(double(length(extrapolation) * points), nrow = length(extrapolation))
  for (i in 1:length(interpolation)){
    plot_ipx[i,] <- plot_theta_x * interpolation[i]
    plot_ipy[i,] <- dt(x_tdistribution,n-1)*(qt(0.001,n-1)-qt(0.999,n-1))/((plot_ipx[i,1] - plot_ipx[i,points]))
  }
  for (i in 1:length(extrapolation)){
    plot_epx[i,] <- plot_theta_x * extrapolation[i]
    plot_epy[i,] <- dt(x_tdistribution,n-1)*(qt(0.001,n-1)-qt(0.999,n-1))/((plot_epx[i,1] - plot_epx[i,points]))
  }
  
  
  return(list("mean_theta" = mean_theta, "lower_bound_theta" = lower_bound_theta, "upper_bound_theta" = upper_bound_theta, 
              "mean_ip" = mean_ip, "lower_bound_ip" = lower_bound_ip, "upper_bound_ip" = upper_bound_ip, 
              "mean_ep" = mean_ep, "lower_bound_ep" = lower_bound_ep, "upper_bound_ep" = upper_bound_ep, 
              "theta" = plot_theta_x, "probability_theta"  = plot_theta_y, "ipx" = plot_ipx, "ipy" = plot_ipy, 
              "epx" = plot_epx, "epy" = plot_epy))
}


# model with model discrepancy
md <- function(X,Z, interpolation, extrapolation){
  Nsim <- 10000
  burn_in <- 1000
  #mcmc <-function(){
  
  
  # parameters for prior
  a_e <- 181/19
  b_e <- 81/95000
  a <- 13/5
  b <- 18/125
  a_psi <- 5
  b_psi <- 5
  
  
  # sampling from prior
  sigma2_e <- sigma2 <- psi <- double(Nsim)
  sigma2_e[1] <- rinvgamma(1, shape = a_e, scale = b_e)
  sigma2[1] <- rinvgamma(1, shape = a, scale = b)
  psi[1] <- rtrunc(1, spec = "gamma", a = 0, b = 4, shape = a_psi, rate = b_psi)
  
  
  # creating lambda and H function
  lambda <- function(x1, x2, psi)exp(- (rdist(x1,x2)/psi)^2)
  H <- function(x_o, x_r, psi)rbind(diag(length(x_o)), lambda(x_r, x_o, psi)%*%solve(lambda(x_o,x_o,psi)))
  
  
  # creating x_o, x_r, x_p and x_op
  x_o <- c(0.2, 0.96, 1.72, 2.48, 3.24, 4.00)
  if (length(X)==11){
    x_r <- X[-c(1,3,5,7,9,11)] }
  else if (length(X)==31){
    x_r <- X[-c(1,7,13,19,25,31)] }
  else{
    x_r <- X[-c(1,13,25,37,49,61)] }
  x_p <- c(interpolation, extrapolation)
  x_op <- c(x_o, x_p)
  n_o <- length(x_o)
  n_r <- length(x_r)
  n_p <- length(x_p)
  n_op <- length(x_op)
  n <- length(X)
  
  # Sampling for theta, delta_o and delta_p
  theta_delta <- matrix(double(9*Nsim), nrow = 9)
  A <- cbind(H(x_o, x_r, psi[1]), 0,0)
  m <- matrix(double(9), ncol = 1)
  m[1] <- (1/sigma2_e[1])*t(X)%*%Z
  m[2:7] <- (1/sigma2_e[1])*t(H(x_o, x_r, psi[1]))%*%Z
  Cov <- matrix(double(9*9), nrow = 9)
  Cov[1,1] <- (1/sigma2_e[1])*t(X)%*%X
  Cov[1,2:9] <- (1/sigma2_e[1])*t(X)%*%A
  Cov[2:9,1] <- (1/sigma2_e[1])*t(A)%*%X
  Cov[2:9,2:9] <- (1/sigma2_e[1])*t(A)%*%A + (1/sigma2[1])*solve(lambda(x_op, x_op, psi[1]))
  theta_delta[,1] <- rmvnorm(1, mean = solve(Cov) %*% m, sigma = solve(Cov))  # sd or Cov!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  tau <- 0.1
  # iterations
  for (i in 2:Nsim){
    sigma2[i] <- rinvgamma(1, shape = a + n_o/2, scale = b + (1/2)*t(theta_delta[2:7, i-1])%*% solve(lambda(x_o,x_o,psi[i-1])) %*%theta_delta[2:7, i-1])
    sigma2_e[i] <- rinvgamma(1, shape = a_e + n/2, scale = b_e + (1/2) * t(Z - X*theta_delta[1, i-1] - H(x_o, x_r, psi[i-1]) %*% theta_delta[2:7, i-1]) %*% (Z - X*theta_delta[1, i-1] - H(x_o, x_r, psi[i-1]) %*% theta_delta[2:7, i-1]))
    omega_old <- log(psi[i-1])
    omega_new <- rtrunc(1, spec = "norm", a = -Inf, b = log(4), mean = omega_old, sd = tau)
    c_yx <- dnorm(omega_new, mean = omega_old, sd = tau)
    c_xy <- dnorm(omega_old, mean = omega_new, sd = tau)
    log_prob <- -1/2*log(det(lambda(x_o,x_o,exp(omega_new)))) + 1/2*log(det(lambda(x_o,x_o,exp(omega_old)))) + 
      (omega_new - omega_old)*(a_psi-1) - b_psi*(exp(omega_new)-exp(omega_old)) - 
      (1/(2*sigma2[i])) * t(theta_delta[2:7,i-1]) %*% (solve(lambda(x_o,x_o,exp(omega_new))) - solve(lambda(x_o,x_o,exp(omega_old)))) %*% theta_delta[2:7,i-1] -
      (1/(2*sigma2_e[i])) * t(Z - X*theta_delta[1,i-1] - H(x_o,x_r,exp(omega_new))%*%theta_delta[2:7,i-1]) %*% (Z - X*theta_delta[1,i-1] - H(x_o,x_r,exp(omega_new))%*%theta_delta[2:7,i-1]) +
      (1/(2*sigma2_e[i])) * t(Z - X*theta_delta[1,i-1] - H(x_o,x_r,exp(omega_old))%*%theta_delta[2:7,i-1]) %*% (Z - X*theta_delta[1,i-1] - H(x_o,x_r,exp(omega_old))%*%theta_delta[2:7,i-1]) +
      log(c_xy/c_yx)
    if (runif(1) < exp(log_prob)){
      psi[i] <- exp(omega_new)
    }else{
      psi[i] <- psi[i-1]
    }
    A <- cbind(H(x_o, x_r, psi[i]), 0, 0)
    m[1] <- (1/sigma2_e[i]) * t(X)%*%Z
    m[2:7] <- (1/sigma2_e[i]) * t(H(x_o,x_r,psi[i])) %*% Z
    Cov[1,1] <- (1/sigma2_e[i]) * t(X) %*% X
    Cov[1,2:9] <- (1/sigma2_e[i]) * t(X) %*% A
    Cov[2:9,1] <- (1/sigma2_e[i]) * t(A) %*% X
    Cov[2:9,2:9] <- (1/sigma2_e[i]) * t(A) %*% A + (1/sigma2[i]) * solve(lambda(x_op,x_op,psi[i]))
    theta_delta[,i] <- rmvnorm(1, mean = solve(Cov) %*% m, sigma = solve(Cov))
  }
  
  
  # interpolation and extrapolation
  ip <- ep <- double(Nsim-burn_in)
  for (i in 1:(Nsim-burn_in)){
    ip[i] <- rnorm(1, mean = interpolation*theta_delta[1,i] + theta_delta[8,i], sd = sqrt(sigma2_e[i]))
    ep[i] <- rnorm(1, mean = extrapolation*theta_delta[1,i] + theta_delta[9,i], sd = sqrt(sigma2_e[i]))
  }
  ip <- sort(ip)
  ep <- sort(ep)
  theta<-sort(theta_delta[1,(burn_in+1):Nsim])
  d0 <- sort(theta_delta[2:7,(burn_in+1):Nsim])
  dp <- sort(theta_delta[8:9,(burn_in+1):Nsim])
  return(list("theta" = theta, "d0" = d0, "dp" = dp, "ip" = ip, "ep" = ep))
}

# Applying function
blr11 <- blr(X11, Z11, interpolation_x, extrapolation_x)
blr31 <- blr(X31, Z31, interpolation_x, extrapolation_x)
blr61 <- blr(X61, Z61, interpolation_x, extrapolation_x)
md11 <- md(X11, Z11, 1.5, extrapolation_x)
md31 <- md(X31, Z31, 1.5, extrapolation_x)
md61 <- md(X61, Z61, 1.5, extrapolation_x)


# confidence intervals and means
data.frame(blr11$lower_bound_theta, blr11$mean_theta, blr11$upper_bound_theta, blr11$lower_bound_ip[21], blr11$mean_ip[21], blr11$upper_bound_ip[21], blr11$lower_bound_ep, blr11$mean_ep, blr11$upper_bound_ep)
data.frame(blr31$lower_bound_theta, blr31$mean_theta, blr31$upper_bound_theta, blr31$lower_bound_ip[21], blr31$mean_ip[21], blr31$upper_bound_ip[21], blr31$lower_bound_ep, blr31$mean_ep, blr31$upper_bound_ep)
data.frame(blr61$lower_bound_theta, blr61$mean_theta, blr61$upper_bound_theta, blr61$lower_bound_ip[21], blr61$mean_ip[21], blr61$upper_bound_ip[21], blr61$lower_bound_ep, blr61$mean_ep, blr61$upper_bound_ep)
data.frame(md11$theta[round(0.05*9000)], md11$theta[round(0.5*9000)], md11$theta[round(0.95*9000)], md11$ip[round(0.05*9000)], md11$ip[round(0.5*9000)], md11$ip[round(0.95*9000)], md11$ep[round(0.05*9000)], md11$ep[round(0.5*9000)], md11$ep[round(0.95*9000)])
data.frame(md31$theta[round(0.05*9000)], md31$theta[round(0.5*9000)], md31$theta[round(0.95*9000)], md31$ip[round(0.05*9000)], md31$ip[round(0.5*9000)], md31$ip[round(0.95*9000)], md31$ep[round(0.05*9000)], md31$ep[round(0.5*9000)], md31$ep[round(0.95*9000)])
data.frame(md61$theta[round(0.05*9000)], md61$theta[round(0.5*9000)], md61$theta[round(0.95*9000)], md61$ip[round(0.05*9000)], md61$ip[round(0.5*9000)], md61$ip[round(0.95*9000)], md61$ep[round(0.05*9000)], md61$ep[round(0.5*9000)], md61$ep[round(0.95*9000)])


## plots ##
###########


# Plot output (work) for simple machine, observation and true process:
X_new <- c(X11,6)
true_new <- c(true11, 3)
simple_machine_output <- 0.65*X_new
predx <- c(1.5, extrapolation_x)
predy <- c(0.65*1.5/(1+1.5/20), extrapolation_zeta)

work_effort = as.data.frame(cbind(X_new, simple_machine_output, true_new))
work_effort2 = as.data.frame(cbind(X11, Z11))
work_effort3 = as.data.frame(cbind(predx, predy))

plot1 <- ggplot() +
  theme_bw() +
  theme(legend.title=element_blank(),
        panel.grid.major = element_line(size = 0.5),
        panel.grid.minor = element_line(size = 0.5),
        legend.background = element_rect(linetype = 1, size = 0.2, colour = 1),
        plot.title = element_text(hjust = 0.5, vjust = 3, size = 18, margin = margin(t = 10, b = 10)),
        axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black", size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.text=element_text(size=14),
        legend.position = c(0.2, 0.8)) + 
  geom_point(data = work_effort2, aes(x = X11,y = Z11, color = "Observations"), size = 3) +
  geom_point(data = work_effort3, aes(x = predx, y = predy, color = "Prediction points"), size = 3 ) +
  geom_line(data = work_effort, aes(x = X_new, y = simple_machine_output , color = "Simulator"), size = 0.8) +
  geom_line(data = work_effort, aes(x = X_new,y = true_new, color = "True process"), size = 0.8) +
  scale_colour_manual(values=c("black", "red", "red", "blue"),
                      guide = guide_legend(override.aes = list(
                        linetype = c("blank", "blank", "solid", "solid"),
                        shape = c(16, 16, NA, NA)))) +
  labs(title="Work delivered by machine",x="x (effort)", y = "Work")
plot1


# plot estimation of theta:
theta1 <- as.data.frame(cbind(as.numeric(blr11$theta), as.numeric(blr31$theta), as.numeric(blr61$theta), 
                              as.numeric(blr11$probability_theta), as.numeric(blr31$probability_theta), as.numeric(blr61$probability_theta)))
theta2 <-as.data.frame(cbind(md11$theta, md31$theta, md61$theta))

plot2 <- ggplot() +
  theme_bw() +
  theme(legend.title=element_blank(),
        panel.grid.major = element_line(size = 0.5),
        panel.grid.minor = element_line(size = 0.5),
        legend.background = element_rect(linetype = 1, size = 0.2, colour = 1),
        plot.title = element_text(hjust = 0.5, vjust = 3, size = 18, margin = margin(t = 10, b = 10)),
        axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black", size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.text=element_text(size=14),
        legend.position = c(0.7, 0.8),
        legend.key.width = unit(2,"cm")) + 
  geom_density(data = theta2, aes(x = V1, color = "11 observations"), size = 0.8, show.legend = FALSE) +
  stat_density(data = theta2, aes(x = V1, color= "11 observations"), geom="line", size = 0.8, position="identity") +
  geom_density(data = theta2, aes(x = V2, color = "31 observations"), size = 0.8, show.legend = FALSE) +
  stat_density(data = theta2, aes(x = V2, color= "31 observations"), geom="line", size = 0.8, position="identity") +
  geom_density(data = theta2, aes(x = V3, color = "61 observations"), size = 0.8, show.legend = FALSE) + 
  stat_density(data = theta2, aes(x = V3, color= "61 observations"), geom="line", size = 0.8, position="identity") +
  geom_line(data = theta1, aes(x = V1, y = V4, color = "11 observations, m1"), size = 0.8, linetype = "dashed") +
  geom_line(data = theta1, aes(x = V2, y = V5, color = "31 observations, m1"), size = 0.8, linetype = "dashed") +
  geom_line(data = theta1, aes(x = V3, y = V6, color = "61 observations, m1"), size = 0.8, linetype = "dashed") + 
  geom_vline(aes(xintercept=0.65, color = "True value"), size = 0.8, show.legend=FALSE) +
  geom_hline(aes(yintercept=0, color = "True value"), size = 0.8, show.legend=FALSE) +
  scale_colour_manual(values=c("dodgerblue", "dodgerblue", "red", "red", "green", "green", "black"),
                      guide = guide_legend(override.aes = list(
                        linetype = c("dashed", "solid", "dashed", "solid", "dashed", "solid", "solid"))),
                      labels = c("11 obs. model 1","11 obs. model 2", "31 obs. model 1",
                                 "31 obs. model 2","61 obs. model 1", "61 obs. model 2",
                                 expression(paste("True value, ", theta, " = ", 0.65)))) +
  labs(title=expression(paste("Estimation of ",theta," for simple machine")), x=expression(theta), y = "Posterior density")
plot2


# plot interpolation x = 1.5
ip1_1.5 <- as.data.frame(cbind(blr11$ipx[length(interpolation_x) + 1,], blr31$ipx[length(interpolation_x) + 1,], blr61$ipx[length(interpolation_x) + 1,], 
                              blr11$ipy[length(interpolation_x) + 1,], blr31$ipy[length(interpolation_x) + 1,], blr61$ipy[length(interpolation_x) + 1,]))
ip2_1.5 <- as.data.frame(cbind(md11$ip, md31$ip, md61$ip))

plot3 <- ggplot() +
  theme_bw() +
  theme(legend.title=element_blank(),
        panel.grid.major = element_line(size = 0.5),
        panel.grid.minor = element_line(size = 0.5),
        legend.background = element_rect(linetype = 1, size = 0.2, colour = 1),
        plot.title = element_text(hjust = 0.5, vjust = 3, size = 18, margin = margin(t = 10, b = 10)),
        axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black", size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.text=element_text(size=14),
        legend.position = c(0.8, 0.8),
        legend.key.width = unit(2,"cm")) + 
  geom_density(data = ip2_1.5, aes(x = V1, color = "11 observations"), size = 0.8, show.legend = FALSE) +
  stat_density(data = ip2_1.5, aes(x = V1, color= "11 observations"), geom="line", size = 0.8, position="identity") +
  geom_density(data = ip2_1.5, aes(x = V2, color = "31 observations"), size = 0.8, show.legend = FALSE) +
  stat_density(data = ip2_1.5, aes(x = V2, color= "31 observations"), geom="line", size = 0.8, position="identity") +
  geom_density(data = ip2_1.5, aes(x = V3, color = "61 observations"), size = 0.8, show.legend = FALSE) + 
  stat_density(data = ip2_1.5, aes(x = V3, color= "61 observations"), geom="line", size = 0.8, position="identity") +
  geom_vline(data = ip2_1.5, aes(xintercept=0.65 * 1.5 /(1 + 1.5/20), color = "True value"), size = 0.8, show.legend=FALSE) +
  geom_hline(data = ip2_1.5, aes(yintercept=0, color = "True value"), size = 0.8, show.legend=FALSE) +
  geom_line(data = ip1_1.5, aes(x = V1, y = V4, color = "11 observations, m1"), size = 0.8, linetype = "dashed") +
  geom_line(data = ip1_1.5, aes(x = V2, y = V5, color = "31 observations, m1"), size = 0.8, linetype = "dashed") +
  geom_line(data = ip1_1.5, aes(x = V3, y = V6, color = "61 observations, m1"), size = 0.8, linetype = "dashed") + 
  xlim(0.6, 1.2) +
  scale_colour_manual(values=c("dodgerblue", "dodgerblue", "red", "red", "green", "green","black"),
                      guide = guide_legend(override.aes = list(
                        linetype = c("dashed", "solid", "dashed", "solid", "dashed", "solid", "solid"))),
                      labels = c("11 obs. model 1","11 obs. model 2", "31 obs. model 1",
                                "31 obs. model 2","61 obs. model 1", "61 obs. model 2",
                                expression(paste("True value, ", zeta(1.5))))) +
  labs(title="Interpolation at point x = 1.5 for simple machine", x="Work", y = "Posterior density")
plot3


# plot extrapolation x = 6
ep1_6 <- as.data.frame(cbind(as.numeric(blr11$epx), as.numeric(blr31$epx), as.numeric(blr61$epx), 
                            as.numeric(blr11$epy), as.numeric(blr31$epy), as.numeric(blr61$epy)))
ep2_6 <- as.data.frame(cbind(md11$ep, md31$ep, md61$ep))



plot4 <- ggplot() +
  theme_bw() +
  theme(legend.title=element_blank(),
        panel.grid.major = element_line(size = 0.5),
        panel.grid.minor = element_line(size = 0.5),
        legend.background = element_rect(linetype = 1, size = 0.2, colour = 1),
        plot.title = element_text(hjust = 0.5, vjust = 3, size = 18, margin = margin(t = 10, b = 10)),
        axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black", size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.text=element_text(size=14),
        legend.position = c(0.8, 0.8),
        legend.key.width = unit(2,"cm")) + 
  geom_density(data = ep2_6, aes(x = V1, color = "11 observations"), size = 0.8, show.legend = FALSE) +
  stat_density(data = ep2_6, aes(x = V1, color= "11 observations"), geom="line", size = 0.8, position="identity") +
  geom_density(data = ep2_6, aes(x = V2, color = "31 observations"), size = 0.8, show.legend = FALSE) +
  stat_density(data = ep2_6, aes(x = V2, color= "31 observations"), geom="line", size = 0.8, position="identity") +
  geom_density(data = ep2_6, aes(x = V3, color = "61 observations"), size = 0.8, show.legend = FALSE) + 
  stat_density(data = ep2_6, aes(x = V3, color= "61 observations"), geom="line", size = 0.8, position="identity") +
  geom_vline(data = ep2_6, aes(xintercept=3, color = "True value"), size = 0.8, show.legend=FALSE) +
  geom_hline(data = ep2_6, aes(yintercept=0, color = "True value"), size = 0.8, show.legend=FALSE) +
  geom_line(data = ep1_6, aes(x = V1, y = V4, color = "11 observations, m1"), size = 0.8, linetype = "dashed") +
  geom_line(data = ep1_6, aes(x = V2, y = V5, color = "31 observations, m1"), size = 0.8, linetype = "dashed") +
  geom_line(data = ep1_6, aes(x = V3, y = V6, color = "61 observations, m1"), size = 0.8, linetype = "dashed") + 
  scale_colour_manual(values=c("dodgerblue", "dodgerblue", "red", "red", "green", "green","black"),
                      guide = guide_legend(override.aes = list(
                        linetype = c("dashed", "solid", "dashed", "solid", "dashed", "solid", "solid"))),
                      labels = c("11 obs. model 1","11 obs. model 2", "31 obs. model 1",
                                 "31 obs. model 2","61 obs. model 1", "61 obs. model 2",
                                 expression(paste("True value, ", zeta(6))))) +
  labs(title="Extrapolation at point x = 6 for simple machine", x="Work", y = "Posterior density")
plot4


# Model 1, Interpolation at all interpolation points (takes 8-9 min to run):
p<-list()
intercept <- 0.65 * interpolation_x /(1 + interpolation_x/20)
for (i in 1:(length(interpolation_x))){
  ip <- as.data.frame(cbind(blr11$ipx[i,], blr31$ipx[i,], blr61$ipx[i,], 
                            blr11$ipy[i,], blr31$ipy[i,], blr61$ipy[i,], intercept[i]))
  p[[i]] <- ggplot(ip) +
    theme_bw() +
    theme(legend.title=element_blank(),
          panel.grid.major = element_line(size = 0.5),
          panel.grid.minor = element_line(size = 0.5),
          legend.background = element_rect(linetype = 1, size = 0.2, colour = 1),
          plot.title = element_text(hjust = 0.5, vjust = 3, size = 18, margin = margin(t = 10, b = 10)),
          axis.text.x = element_text(color = "black", size = 6),
          axis.text.y = element_text(color = "black", size = 6),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.text=element_text(size=12),
          legend.position = c(0.6, 0.8)) + 
    geom_line(aes(x = V1, y = V4), color = "dodgerblue", size = 0.8) +
    geom_line(aes(x = V2, y = V5), color = "red", size = 0.8) +
    geom_line(aes(x = V3, y = V6), color = "green", size = 0.8) + 
    geom_vline(aes(xintercept=V7), color = "black", size = 0.8, show.legend=FALSE)
}
do.call("grid.arrange", c(p, ncol = 4, top = "Distribution of interpolation points for simple machine, model 1"))


# Model 2, Interpolation at all interpolation points:
md_ip61 <- list()
intercept <- 0.65 * interpolation_x /(1 + interpolation_x/20)

start_time <- Sys.time()
for (i in 1:(length(interpolation_x))){
  md_ip61[[i]] <- md(X61, Z61, interpolation_x[i], extrapolation_x)
  print(i)
}
end_time <- Sys.time()
time <- end_time - start_time
print(time)
p<-list()
for (i in 1:(length(interpolation_x))){
  ip <- as.data.frame(cbind(md_ip61[[i]]$ip,intercept[i]))
  p[[i]] <- ggplot(ip) +
    theme_bw() +
    theme(legend.title=element_blank(),
          panel.grid.major = element_line(size = 0.5),
          panel.grid.minor = element_line(size = 0.5),
          legend.background = element_rect(linetype = 1, size = 0.2, colour = 1),
          plot.title = element_text(hjust = 0.5, vjust = 3, size = 18, margin = margin(t = 10, b = 10)),
          axis.text.x = element_text(color = "black", size = 6),
          axis.text.y = element_text(color = "black", size = 6),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.text=element_text(size=12),
          legend.position = c(0.6, 0.8)) + 
    geom_density(aes(x = V1), color = "green", size = 0.8) +
    geom_vline(aes(xintercept=V2), color = "black", size = 0.8, show.legend=FALSE) +
    geom_hline(aes(yintercept=0), color = "black", size = 0.8, show.legend=FALSE)
}
do.call("grid.arrange", c(p, ncol = 4, top = "Distribution of interpolation points for simple machine, model 2"))


#comparing interpolation for model 1 and 2
p<-list()
ip1 <- list()
ip2 <- list()
for (i in 1:(length(interpolation_x))){
  ip1[[i]] <- as.data.frame(cbind(md_ip61[[i]]$ip,intercept[i]))
  ip2[[i]] <- as.data.frame(cbind(as.numeric(blr61$ipx[i,]), as.numeric(blr61$ipy[i,])))
  p[[i]] <- ggplot() +
    theme_bw() +
    theme(legend.title=element_blank(),
          panel.grid.major = element_line(size = 0.5),
          panel.grid.minor = element_line(size = 0.5),
          legend.background = element_rect(linetype = 1, size = 0.2, colour = 1),
          plot.title = element_text(hjust = 0.5, vjust = 3, size = 18, margin = margin(t = 10, b = 10)),
          axis.text.x = element_text(color = "black", size = 6),
          axis.text.y = element_text(color = "black", size = 6),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.text=element_text(size=12),
          legend.position = c(0.6, 0.8)) + 
    geom_density(data = ip1[[i]], aes(x = V1), color = "green", size = 0.8) + 
    geom_line(data = ip2[[i]], aes(x = V1, y = V2), color = "dodgerblue", size = 0.8) +
    geom_vline(data = ip1[[i]], aes(xintercept=V2), color = "black", size = 0.8, show.legend=FALSE) +
    geom_hline(aes(yintercept=0), color = "black", size = 0.8, show.legend=FALSE)
}
do.call("grid.arrange", c(p, ncol = 4, top = "Post. distribution of interpolation points for simple machine"))


# plot with 90% confidence interval for all interpolation points.
lower_ip<-c()
upper_ip<-c()
for (i in 1:20){
  lower_ip[i]<-md_ip61[[i]]$ip[round(0.05*9000)]
  upper_ip[i]<-md_ip61[[i]]$ip[round(0.95*9000)]
}
interpolation1_90 <- list(as.data.frame(cbind(X61, Z61, true61)), 
                        as.data.frame(cbind( as.numeric(interpolation_x),
                                             as.numeric(blr61$lower_bound_theta%*%interpolation_x), 
                                             as.numeric(blr61$upper_bound_theta%*%interpolation_x), 
                                             as.numeric(interpolation_zeta))))
interpolation1_90 <- do.call(rbind.fill, interpolation1_90)
interpolation2_90 <- list(as.data.frame(cbind(X61, Z61, true61)), 
                        as.data.frame(cbind( as.numeric(interpolation_x),as.numeric(lower_ip), 
                                             as.numeric(upper_ip), as.numeric(interpolation_zeta))))
interpolation2_90 <- do.call(rbind.fill, interpolation2_90)

plot5 <- ggplot() +
  theme_bw() +
  theme(legend.title=element_blank(),
        panel.grid.major = element_line(size = 0.5),
        panel.grid.minor = element_line(size = 0.5),
        legend.background = element_rect(linetype = 1, size = 0.2, colour = 1),
        plot.title = element_text(hjust = 0.5, vjust = 3, size = 18, margin = margin(t = 10, b = 10)),
        axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black", size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.text=element_text(size=14),
        legend.position = c(0.2, 0.8)) + 
  geom_line(data = interpolation2_90, aes(x = X61, y = true61, color = "True process"), size = 0.5) + 
  geom_line(data = interpolation2_90, aes(x = V1, y = V2, color = "90% Cred. interval, model 2"), size = 0.5) + 
  geom_line(data = interpolation2_90, aes(x = V1, y = V3, color = "90% Cred. interval, model 2"), size = 0.5, show.legend = FALSE) + 
  geom_point(data = interpolation2_90, aes(x = V1, y = V4, color = "Interpolation points"), size = 1) +
  geom_line(data = interpolation1_90, aes(x = V1, y = V2, color = "90% Cred. interval, model 1"), size = 0.5) + 
  geom_line(data = interpolation1_90, aes(x = V1, y = V3, color = "90% Cred. interval, model 1"), size = 0.5, show.legend = FALSE) + 
  scale_colour_manual(values=c("dodgerblue", "black", "black", "red"),
                      guide = guide_legend(override.aes = list(
                        linetype = c("solid", "solid", "blank", "solid"),
                        shape = c(NA, NA, 16, NA)))) +
  labs(title="Work delivered by simple machine with 90% CI interpolation",x="x (effort)", y = "Work")
plot5