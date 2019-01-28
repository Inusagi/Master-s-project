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

# loading datasets for golf model
data_golfmodel_standard <- read.csv("~/Documents/R/Prosjektoppgaven/data_golfmodel_standard.txt", header=FALSE)


# velocity^2 and landingposition for golf model for 11, 31 and 61 observations + interpolation and extrpolation points
X11 <- as.double(data_golfmodel_standard[1,!is.na(data_golfmodel_standard[1,])])
X31 <- as.double(data_golfmodel_standard[2,!is.na(data_golfmodel_standard[2,])])
X61 <- as.double(data_golfmodel_standard[3,!is.na(data_golfmodel_standard[3,])])
true11 <- as.double(data_golfmodel_standard[4,!is.na(data_golfmodel_standard[4,])])
true31 <- as.double(data_golfmodel_standard[5,!is.na(data_golfmodel_standard[5,])])
true61 <- as.double(data_golfmodel_standard[6,!is.na(data_golfmodel_standard[6,])])
interpolation_v2 <- as.double(data_golfmodel_standard[7,!is.na(data_golfmodel_standard[7,])])
interpolation_x <- as.double(data_golfmodel_standard[8,!is.na(data_golfmodel_standard[8,])])
extrapolation_v2 <- 6400
extrapolation_x <- 146.82771316504866
variance <- 0.01^2*((123.05649-0.351517)/(2.1666667-0.1287129))
set.seed(5)
error11 <- rnorm(11,0,variance)
error31 <- rnorm(31,0,variance)
error61 <- rnorm(61,0,variance)
Z11 <- true11 + error11
Z31 <- true31 + error31
Z61 <- true61 + error61

# function using bayesian linear regressino to calculate 90% confidence interval for g, landing position for 
# interpolating and extrapolating points, points for plotting posterior distribution of g and 
# landing position for interpolating and extrapolating points. 
blr <- function(Z, X, interpolation, extrapolation, angle, g = NA){
  
  
  n = length(X)
  theta_hat <- solve(t(X) %*% X) %*% t(X) %*% Z
  S2 <- as.numeric((1/(n-1))%*%t(Z-X%*%theta_hat)%*%(Z-X%*%theta_hat))
  
  
  if (is.na(g)){  # if g is "unknown"
    
    C <- 2*tan(angle) * (cos(angle))^2  # constant
    
    
    # calulating 90% confidence interval for posterior and mean of g, x(4200) and x(6400)
    mean_g <- C/((sqrt(S2)/sqrt(t(X) %*% X)) %*% qt(0.50,n-1) + theta_hat)
    upper_bound_g <- C/((sqrt(S2)/sqrt(t(X) %*% X)) %*% qt(0.05,n-1) + theta_hat)
    lower_bound_g <- C/((sqrt(S2)/sqrt(t(X) %*% X)) %*% qt(0.95,n-1) + theta_hat)
    interpolation <- c(interpolation,1700)
    mean_ip <- as.numeric(C/mean_g) * interpolation  # interpolation
    lower_bound_ip <- as.numeric(C/upper_bound_g) * interpolation
    upper_bound_ip <- as.numeric(C/lower_bound_g) * interpolation
    mean_ep <- as.numeric(C/mean_g) * extrapolation  # extrapolation
    lower_bound_ep <- as.numeric(C/upper_bound_g) * extrapolation
    upper_bound_ep <- as.numeric(C/lower_bound_g) * extrapolation
    
    
    
    # calculating values for plotting the posterior density of g, x(4200) and x(6400)
    points <- 100
    x_tdistribution <- seq(qt(0.001,n-1),qt(0.999,n-1),length = points)
    plot_gx <- C/((sqrt(S2)/c(sqrt(t(X) %*% X))) * (x_tdistribution) + c(theta_hat))
    plot_gy <- dt(x_tdistribution,n-1)*(qt(0.999,n-1)-qt(0.001,n-1))/((plot_gx[1] - plot_gx[points]))
    plot_ipx <- matrix(double(length(interpolation) * points), nrow = length(interpolation))
    plot_ipy <- matrix(double(length(interpolation) * points), nrow = length(interpolation))
    plot_epx <- matrix(double(length(extrapolation) * points), nrow = length(extrapolation))
    plot_epy <- matrix(double(length(extrapolation) * points), nrow = length(extrapolation))
    for (i in 1:length(interpolation)){
      plot_ipx[i,] <- C/rev(plot_gx) * interpolation[i]
      plot_ipy[i,] <- dt(x_tdistribution,n-1)*(qt(0.999,n-1)-qt(0.001,n-1))/((plot_ipx[i,1] - plot_ipx[i,points]))
    }
    for (i in 1:length(extrapolation)){
      plot_epx[i,] <- C/rev(plot_gx) * extrapolation[i]
      plot_epy[i,] <- dt(x_tdistribution,n-1)*(qt(0.999,n-1)-qt(0.001,n-1))/((plot_epx[i,1] - plot_epx[i,points]))
    }
    
    # returns a list with mean, lower_bound, upper_bound, g and density values for g
    return(list("mean_g" = mean_g, "lower_bound_g" = lower_bound_g, "upper_bound_g" = upper_bound_g, 
                "mean_ip" = mean_ip, "lower_bound_ip" = lower_bound_ip, "upper_bound_ip" = upper_bound_ip, 
                "mean_ep" = mean_ep, "lower_bound_ep" = lower_bound_ep, "upper_bound_ep" = upper_bound_ep, 
                "g" = plot_gx, "probability_g"  = plot_gy, "ipx" = plot_ipx, "ipy" = plot_ipy, 
                "epx" = plot_epx, "epy" = plot_epy))
  }
  
}


## Gibbs sampling ##
####################

md <- function(X,Z, interpolation, extrapolation){
  Nsim <- 10000
  burn_in <- 1000
  #mcmc <-function(){
  
  
  # parameters for prior
  a_e <- 181/19
  b_e <- 81/95000
  a <- 13/5
  b <- 18/125
  constant <- (4900-10)/(4*(4-0.2))  # simple machine to golf model
  a_psi <- 5
  b_psi <- 5/constant
  trunc_psi <- 4*constant
  
  
  # sampling from prior
  sigma2_e <- sigma2 <- psi <- double(Nsim)
  sigma2_e[1] <- rinvgamma(1, shape = a_e, scale = b_e)
  sigma2[1] <- rinvgamma(1, shape = a, scale = b)
  psi[1] <- rtrunc(1, spec = "gamma", a = 0, b = trunc_psi, shape = a_psi, rate = b_psi)
  
  
  # creating lambda and H function
  lambda <- function(x1, x2, psi)exp(- (rdist(x1,x2)/psi)^2)
  H <- function(x_o, x_r, psi)rbind(diag(length(x_o)), lambda(x_r, x_o, psi)%*%solve(lambda(x_o,x_o,psi)))
  
  
  # creating x_o, x_r, x_p and x_op
  x_o <- c(10, 988, 1966, 2944, 3922, 4900)
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
    omega_new <- rtrunc(1, spec = "norm", a = -Inf, b = log(trunc_psi), mean = omega_old, sd = tau)
    c_yx <- integrate(dnorm, mean = omega_old, sd = tau, -Inf, log(trunc_psi))$value
    c_xy <- integrate(dnorm, mean = omega_new, sd = tau, -Inf, log(trunc_psi))$value
    log_prob <- -1/2*log(det(lambda(x_o,x_o,exp(omega_new)))) + 1/2*log(det(lambda(x_o,x_o,exp(omega_old)))) + 
      (omega_new - omega_old)*(a_psi-1) - b_psi*(exp(omega_new)-exp(omega_old)) - 
      (1/(2*sigma2[i])) * t(theta_delta[2:7,i-1]) %*% (solve(lambda(x_o,x_o,exp(omega_new))) - solve(lambda(x_o,x_o,exp(omega_old)))) %*% theta_delta[2:7,i-1] -
      (1/(2*sigma2_e[i])) * t(Z - X*theta_delta[1,i-1] - H(x_o,x_r,exp(omega_new))%*%theta_delta[2:7,i-1]) %*% (Z - X*theta_delta[1,i-1] - H(x_o,x_r,exp(omega_new))%*%theta_delta[2:7,i-1]) +
      (1/(2*sigma2_e[i])) * t(Z - X*theta_delta[1,i-1] - H(x_o,x_r,exp(omega_old))%*%theta_delta[2:7,i-1]) %*% (Z - X*theta_delta[1,i-1] - H(x_o,x_r,exp(omega_old))%*%theta_delta[2:7,i-1]) +
      log(c_yx/c_xy)
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
  
  return(list("theta" = theta_delta[1,(burn_in+1):Nsim], "d0" = theta_delta[2:7,(burn_in+1):Nsim], "dp" = theta_delta[8:9,(burn_in+1):Nsim], "ip" = ip, "ep" = ep))
}


# applying functions
launch_angle_deg <- 11 * pi/180  # launch angle degree
blr11 <- blr(Z11, X11, interpolation_v2, extrapolation_v2, angle = launch_angle_deg)
blr31 <- blr(Z31, X31, interpolation_v2, extrapolation_v2, angle = launch_angle_deg)
blr61 <- blr(Z61, X61, interpolation_v2, extrapolation_v2, angle = launch_angle_deg)
md11 <- md(X11, Z11, 1700, extrapolation_v2)
md31 <- md(X31, Z31, 1700, extrapolation_v2)
md61 <- md(X61, Z61, 1700, extrapolation_v2)


#Printing confidence intervals and mean:
C <- 2*tan(launch_angle_deg) * (cos(launch_angle_deg))^2
data.frame(blr11$lower_bound_g, blr11$mean_g, blr11$upper_bound_g, blr11$lower_bound_ip[21], blr11$mean_ip[21], blr11$upper_bound_ip[21], blr11$lower_bound_ep, blr11$mean_ep, blr11$upper_bound_ep)
data.frame(blr31$lower_bound_g, blr31$mean_g, blr31$upper_bound_g, blr31$lower_bound_ip[21], blr31$mean_ip[21], blr31$upper_bound_ip[21], blr31$lower_bound_ep, blr31$mean_ep, blr31$upper_bound_ep)
data.frame(blr61$lower_bound_g, blr61$mean_g, blr61$upper_bound_g, blr61$lower_bound_ip[21], blr61$mean_ip[21], blr61$upper_bound_ip[21], blr61$lower_bound_ep, blr61$mean_ep, blr61$upper_bound_ep)
data.frame(C/sort(md11$theta)[round(0.05*9000)], C/sort(md11$theta)[round(0.5*9000)], C/sort(md11$theta)[round(0.95*9000)], sort(md11$ip)[round(0.05*9000)], sort(md11$ip)[round(0.5*9000)], sort(md11$ip)[round(0.95*9000)], sort(md11$ep)[round(0.05*9000)], sort(md11$ep)[round(0.5*9000)], sort(md11$ep)[round(0.95*9000)])
data.frame(C/sort(md31$theta)[round(0.05*9000)], C/sort(md31$theta)[round(0.5*9000)], C/sort(md31$theta)[round(0.95*9000)], sort(md31$ip)[round(0.05*9000)], sort(md31$ip)[round(0.5*9000)], sort(md31$ip)[round(0.95*9000)], sort(md31$ep)[round(0.05*9000)], sort(md31$ep)[round(0.5*9000)], sort(md31$ep)[round(0.95*9000)])
data.frame(C/sort(md61$theta)[round(0.05*9000)], C/sort(md61$theta)[round(0.5*9000)], C/sort(md61$theta)[round(0.95*9000)], sort(md61$ip)[round(0.05*9000)], sort(md61$ip)[round(0.5*9000)], sort(md61$ip)[round(0.95*9000)], sort(md61$ep)[round(0.05*9000)], sort(md61$ep)[round(0.5*9000)], sort(md61$ep)[round(0.95*9000)])


# plot of difference between true process and golf ball model without air resistance
C <- 2*tan(launch_angle_deg) * (cos(launch_angle_deg))^2
X_new <- c(X11,6400)
true_new <- c(true11, extrapolation_x)
simple <- X_new*C/9.81
predx <- c(1700, 6400)
predy <- c(54.427808530863516, extrapolation_x)

landing_velocitysq <- as.data.frame(cbind(X_new, simple, true_new))
landing_velocitysq2 <- as.data.frame(cbind(X11, Z11))
landing_velocitysq3 <- as.data.frame(cbind(predx, predy))


ggplot() +
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
  geom_point(data = landing_velocitysq2, aes(x = X11,y = Z11, color = "Observations"), size = 3) +
  geom_point(data = landing_velocitysq3, aes(x = predx, y = predy, color = "Prediction points"), size = 3 ) +
  geom_line(data = landing_velocitysq, aes(x = X_new, y = simple, color = "Simulator"), size = 0.8) +
  geom_line(data = landing_velocitysq, aes(x = X_new,y = true_new, color = "True process"), size = 0.8) +
  scale_colour_manual(values=c("black", "red", "red", "blue"),
                      guide = guide_legend(override.aes = list(
                        linetype = c("blank", "blank", "solid", "solid"),
                        shape = c(16, 16, NA, NA)))) +
  labs(title="Landing position of golf ball", x = expression(paste("Initial velocity squared ", v^2," ", (m^2/s^2))), y = "Landing position (m)")


# plot estimation of g:
g1 <- as.data.frame(cbind(as.numeric(blr11$g), as.numeric(blr31$g), as.numeric(blr61$g), 
                               as.numeric(blr11$probability_g), as.numeric(blr31$probability_g), as.numeric(blr61$probability_g)))
g2 <-as.data.frame(cbind(C/md11$theta, C/md31$theta, C/md61$theta))

ggplot() +
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
        legend.position = c(0.3, 0.8),
        legend.key.width = unit(2,"cm")) + 
  geom_density(data = g2, aes(x = V1, color = "11 observations"), size = 0.8, show.legend = FALSE) +
  stat_density(data = g2, aes(x = V1, color= "11 observations"), geom="line", size = 0.8, position="identity") +
  geom_density(data = g2, aes(x = V2, color = "31 observations"), size = 0.8, show.legend = FALSE) +
  stat_density(data = g2, aes(x = V2, color= "31 observations"), geom="line", size = 0.8, position="identity") +
  geom_density(data = g2, aes(x = V3, color = "61 observations"), size = 0.8, show.legend = FALSE) + 
  stat_density(data = g2, aes(x = V3, color= "61 observations"), geom="line", size = 0.8, position="identity") +
  geom_line(data = g1, aes(x = V1, y = V4, color = "11 observations, m1"), size = 0.8, linetype = "dashed") +
  geom_line(data = g1, aes(x = V2, y = V5, color = "31 observations, m1"), size = 0.8, linetype = "dashed") +
  geom_line(data = g1, aes(x = V3, y = V6, color = "61 observations, m1"), size = 0.8, linetype = "dashed") + 
  geom_vline(aes(xintercept=9.81, color = "True value"), size = 0.8, show.legend=FALSE) +
  geom_hline(aes(yintercept=0, color = "True value"), size = 0.8, show.legend=FALSE) +
  scale_colour_manual(values=c("dodgerblue", "dodgerblue", "red", "red", "green", "green", "black"),
                      guide = guide_legend(override.aes = list(
                        linetype = c("dashed", "solid", "dashed", "solid", "dashed", "solid", "solid"))),
                      labels = c("11 obs. model 1","11 obs. model 2", "31 obs. model 1",
                                 "31 obs. model 2","61 obs. model 1", "61 obs. model 2",
                                 expression(paste("True value, ", g, " = ", 9.81)))) +
  labs(title=expression(paste("Estimation of ",g," for golf shots")), x=expression(paste("g ",(m/s^2))), y = "Posterior density")


# plot interpolation v2 = 1700
ip1_1700 <- as.data.frame(cbind(blr11$ipx[length(interpolation_v2) + 1,], blr31$ipx[length(interpolation_v2) + 1,], blr61$ipx[length(interpolation_v2) + 1,], 
                               blr11$ipy[length(interpolation_v2) + 1,], blr31$ipy[length(interpolation_v2) + 1,], blr61$ipy[length(interpolation_v2) + 1,]))

ip2_1700 <- as.data.frame(cbind(md11$ip, md31$ip, md61$ip))

ggplot() +
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
        legend.position = c(0.2, 0.8),
        legend.key.width = unit(2,"cm")) + 
  geom_density(data = ip2_1700, aes(x = V1, color = "11 observations"), size = 0.8, show.legend = FALSE) +
  stat_density(data = ip2_1700, aes(x = V1, color= "11 observations"), geom="line", size = 0.8, position="identity") +
  geom_density(data = ip2_1700, aes(x = V2, color = "31 observations"), size = 0.8, show.legend = FALSE) +
  stat_density(data = ip2_1700, aes(x = V2, color= "31 observations"), geom="line", size = 0.8, position="identity") +
  geom_density(data = ip2_1700, aes(x = V3, color = "61 observations"), size = 0.8, show.legend = FALSE) + 
  stat_density(data = ip2_1700, aes(x = V3, color= "61 observations"), geom="line", size = 0.8, position="identity") +
  geom_vline(data = ip2_1700, aes(xintercept=54.427808530863516, color = "True value"), size = 0.8, show.legend=FALSE) +
  geom_hline(data = ip2_1700, aes(yintercept=0, color = "True value"), size = 0.8, show.legend=FALSE) +
  geom_line(data = ip1_1700, aes(x = V1, y = V4, color = "11 observations, m1"), size = 0.8, linetype = "dashed") +
  geom_line(data = ip1_1700, aes(x = V2, y = V5, color = "31 observations, m1"), size = 0.8, linetype = "dashed") +
  geom_line(data = ip1_1700, aes(x = V3, y = V6, color = "61 observations, m1"), size = 0.8, linetype = "dashed") + 
  scale_colour_manual(values=c("dodgerblue", "dodgerblue", "red", "red", "green", "green","black"),
                      guide = guide_legend(override.aes = list(
                        linetype = c("dashed", "solid", "dashed", "solid", "dashed", "solid", "solid"))),
                      labels = c("11 obs. model 1","11 obs. model 2", "31 obs. model 1",
                                 "31 obs. model 2","61 obs. model 1", "61 obs. model 2",
                                 expression(paste("True value, ", zeta(1700))))) +
  labs(title=expression(paste("Interpolation at point ", v^2 ," = 1700 for golf shots")), x="Landing position (m)", y = "Posterior density")


# plot extrapolation v2 = 6400
ep1_6400 <- as.data.frame(cbind(as.numeric(blr11$epx), as.numeric(blr31$epx), as.numeric(blr61$epx), 
                               as.numeric(blr11$epy), as.numeric(blr31$epy), as.numeric(blr61$epy)))
ep2_6400 <- as.data.frame(cbind(md11$ep, md31$ep, md61$ep))



ggplot() +
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
        legend.position = c(0.82, 0.8),
        legend.key.width = unit(2,"cm")) + 
  geom_density(data = ep2_6400, aes(x = V1, color = "11 observations"), size = 0.8, show.legend = FALSE) +
  stat_density(data = ep2_6400, aes(x = V1, color= "11 observations"), geom="line", size = 0.8, position="identity") +
  geom_density(data = ep2_6400, aes(x = V2, color = "31 observations"), size = 0.8, show.legend = FALSE) +
  stat_density(data = ep2_6400, aes(x = V2, color= "31 observations"), geom="line", size = 0.8, position="identity") +
  geom_density(data = ep2_6400, aes(x = V3, color = "61 observations"), size = 0.8, show.legend = FALSE) + 
  stat_density(data = ep2_6400, aes(x = V3, color= "61 observations"), geom="line", size = 0.8, position="identity") +
  geom_vline(data = ep2_6400, aes(xintercept=146.82771316504866, color = "True value"), size = 0.8, show.legend=FALSE) +
  geom_hline(data = ep2_6400, aes(yintercept=0, color = "True value"), size = 0.8, show.legend=FALSE) +
  geom_line(data = ep1_6400, aes(x = V1, y = V4, color = "11 observations, m1"), size = 0.8, linetype = "dashed") +
  geom_line(data = ep1_6400, aes(x = V2, y = V5, color = "31 observations, m1"), size = 0.8, linetype = "dashed") +
  geom_line(data = ep1_6400, aes(x = V3, y = V6, color = "61 observations, m1"), size = 0.8, linetype = "dashed") + 
  scale_colour_manual(values=c("dodgerblue", "dodgerblue", "red", "red", "green", "green","black"),
                      guide = guide_legend(override.aes = list(
                        linetype = c("dashed", "solid", "dashed", "solid", "dashed", "solid", "solid"))),
                      labels = c("11 obs. model 1","11 obs. model 2", "31 obs. model 1",
                                 "31 obs. model 2","61 obs. model 1", "61 obs. model 2",
                                 expression(paste("True value, ", zeta(6400))))) +
  labs(title=expression(paste("Extrapolation at point ", v^2 ," = 6400 for golf shots")), x="Landing position (m)", y = "Posterior density")


# Model 1, Interpolation at all interpolation points:
p<-list()
for (i in 1:(length(interpolation_x))){
  ip <- as.data.frame(cbind(blr11$ipx[i,], blr31$ipx[i,], blr61$ipx[i,], 
                            blr11$ipy[i,], blr31$ipy[i,], blr61$ipy[i,], interpolation_x[i]))
  p[[i]] <- ggplot(ip) +
    theme_bw() +
    theme(legend.title=element_blank(),
          panel.grid.major = element_line(size = 0.5),
          panel.grid.minor = element_line(size = 0.5),
          legend.background = element_rect(linetype = 1, size = 0.2, colour = 1),
          plot.title = element_text(hjust = 0.5, vjust = 3, size = 18, margin = margin(t = 10, b = 10)),
          axis.text.x = element_text(color = "black", size = 7),
          axis.text.y = element_text(color = "black", size = 7),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.text=element_text(size=12),
          legend.position = c(0.6, 0.8)) + 
    geom_line(aes(x = V1, y = V4), color = "dodgerblue", size = 0.8) +
    geom_line(aes(x = V2, y = V5), color = "red", size = 0.8) +
    geom_line(aes(x = V3, y = V6), color = "green", size = 0.8) + 
    geom_vline(aes(xintercept=V7), color = "black", size = 0.8, show.legend=FALSE)
}
do.call("grid.arrange", c(p, ncol = 4, top = "Distribution of interpolation points for golf shots, model 1"))


# Model 2, Interpolation at all interpolation points, (takes 8-9 min to run):
md_ip61 <- list()

start_time <- Sys.time()
for (i in 1:(length(interpolation_x))){
  md_ip61[[i]] <- md(X61, Z61, interpolation_v2[i], extrapolation_x)
  print(i)
}
end_time <- Sys.time()
time <- end_time - start_time
print(time)
p<-list()
for (i in 1:(length(interpolation_x))){
  ip <- as.data.frame(cbind(md_ip61[[i]]$ip,interpolation_x[i]))
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
do.call("grid.arrange", c(p, ncol = 4, top = "Distribution of interpolation points for golf shots, model 2"))


#comparing interpolation for model 1 and 2
p<-list()
ip1 <- list()
ip2 <- list()
for (i in 1:(length(interpolation_x))){
  ip1[[i]] <- as.data.frame(cbind(md_ip61[[i]]$ip,interpolation_x[i]))
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
do.call("grid.arrange", c(p, ncol = 4, top = "Post. distribution of interpolation points for golf shots"))


# plot with 90% confidence interval for all interpolation points.
lower_ip<-c()
upper_ip<-c()
for (i in 1:20){
  lower_ip[i]<-sort(md_ip61[[i]]$ip)[round(0.05*9000)]
  upper_ip[i]<-sort(md_ip61[[i]]$ip)[round(0.95*9000)]
}
interpolation1_90 <- list(as.data.frame(cbind(X61, Z61, true61)), 
                          as.data.frame(cbind( as.numeric(interpolation_v2),
                                               as.numeric((C/blr61$lower_bound_g)%*%interpolation_v2), 
                                               as.numeric((C/blr61$upper_bound_g)%*%interpolation_v2), 
                                               as.numeric(interpolation_x))))
interpolation1_90 <- do.call(rbind.fill, interpolation1_90)
interpolation2_90 <- list(as.data.frame(cbind(X61, Z61, true61)), 
                          as.data.frame(cbind( as.numeric(interpolation_v2),as.numeric(lower_ip), 
                                               as.numeric(upper_ip), as.numeric(interpolation_x))))
interpolation2_90 <- do.call(rbind.fill, interpolation2_90)

ggplot() +
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
  labs(title="Landing position of golf ball with 90% CI interpolation",x=expression(paste("Initial velocity squared ", v^2," ", (m^2/s^2))),y = "Landing position (m)")
