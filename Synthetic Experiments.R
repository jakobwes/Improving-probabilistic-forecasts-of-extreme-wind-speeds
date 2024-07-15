library("latex2exp")

CRPS <- function(f, X, Y_X){
  f_X <- f(X)
  mean(abs(f_X - Y_X)) - 0.5*mean(dist(f_X, method = "manhattan"))
}

twCRPS <- function(f, X, Y_X, t){
  f_X <- f(X)
  
  f_X[f_X < t] <- t
  if(Y_X < t) Y_X <- t
  
  mean(abs(f_X - Y_X)) - 0.5*mean(dist(f_X, method = "manhattan"))
}


# Experiment 1 ------------------------------------------------------------
# GPD normal mixture: identify xi

xi <- 0.5
mu <- 2

n <- 1000

target <- 0.4 *rnorm(n, mean = mu) + 0.6 * (runif(n)^(-xi) - 1)/xi

candid_vals <- seq(0.01, 0.5, length.out = 100)
candid_vals <- seq(0.25, 0.75, length.out = 100)

get_crps <- function(xi_val, n_generations = 250){
  f <- function(xi_val, n_generations = 250) {0.4 *rnorm(n_generations, mean = mu) + 0.6 * (runif(n_generations)^(-xi_val) - 1)/xi_val}
  mean(sapply(1:length(target), function(i){CRPS(f, xi_val, target[i])}))
}

get_twcrps <- function(xi_val, n_generations = 250, q = 0.9){
  f <- function(xi_val, n_generations = 250) {0.4 *rnorm(n_generations, mean = mu) + 0.6 * (runif(n_generations)^(-xi_val) - 1)/xi_val}
  mean(sapply(1:length(target), function(i){twCRPS(f, xi_val, target[i], quantile(target, q))}))
}

crps_vals <- sapply(candid_vals, get_crps) 
twcrps_vals <- sapply(candid_vals, get_twcrps) 
twcrps_vals_095 <- sapply(candid_vals, get_twcrps, q = 0.95)

png(filename="plots/synthetic_experiments/ex_1_crps.png", width = 500, height = 300)
plot(candid_vals, crps_vals, xlab = TeX("Candidate values for $\\xi$"), ylab = "CRPS"); abline(v=xi, col = "blue")
dev.off()

png(filename="plots/synthetic_experiments/ex_1_twcrps_1.png", width = 500, height = 300)
plot(candid_vals, twcrps_vals, xlab = TeX("Candidate values for $\\xi$"), ylab = "twCRPS"); abline(v=xi, col = "blue")
dev.off()

png(filename="plots/synthetic_experiments/ex_1_twcrps_2.png", width = 500, height = 300)
plot(candid_vals, twcrps_vals_095, xlab = TeX("Candidate values for $\\xi$"), ylab = "twCRPS"); abline(v=xi, col = "blue")
dev.off()


# Experiment 2 ------------------------------------------------------------
# GPD-normal mixture: identify xi and mu

xi <- 0.5
mu <- 2

n <- 1000
m <- 250

target <- 0.4 *rnorm(n, mean = mu) + 0.6 * (runif(n)^(-xi) - 1)/xi

candid_vals_xi <- seq(0.25, 0.75, length.out = 20)
candid_vals_mu <- seq(1, 3, length.out = 20)


CRPS <- function(f, Y_X, ...){
  f_X <- f(...)
  mean(abs(f_X - Y_X)) - 0.5*mean(dist(f_X, method = "manhattan"))
}

twCRPS <- function(f, Y_X, t, ...){
  f_X <- f(...)
  
  f_X[f_X < t] <- 0
  if(Y_X < t) Y_X <- 0
  
  mean(abs(f_X - Y_X)) - 0.5*mean(dist(f_X, method = "manhattan"))
}

get_crps <- function(mu_val, xi_val, n_generations = 250){
  f <- function(mu_val, xi_val, n_generations = 250) {0.4 *rnorm(n_generations, mean = mu_val) + 0.6 * (runif(n_generations)^(-xi_val) - 1)/xi_val}
  mean(sapply(1:length(target), function(i){CRPS(f, target[i], mu_val, xi_val)}))
}

get_twcrps <- function(mu_val, xi_val, n_generations = 250, q = 0.9){
  f <- function(mu_val, xi_val, n_generations = 250) {0.4 *rnorm(n_generations, mean = mu_val) + 0.6 * (runif(n_generations)^(-xi_val) - 1)/xi_val}
  mean(sapply(1:length(target), function(i){twCRPS(f, target[i], quantile(target, q),  mu_val, xi_val)}))
}

grid <- expand.grid(candid_vals_mu, candid_vals_xi)

crps_vals <- do.call(mapply, c(function(a, b) get_crps(mu_val = a, xi_val = b), unname(grid)))
twcrps_vals <- do.call(mapply, c(function(a, b) get_twcrps(mu_val = a, xi_val = b), unname(grid)))

grid$crps <- crps_vals
grid$twcrps <- twcrps_vals


library("lattice")
levelplot(crps ~ Var1*Var2, data = grid, xlab = "mu", ylab = "xi")
levelplot(twcrps ~ Var1*Var2, data = grid, xlab = "mu", ylab = "xi")

png(filename="plots/synthetic_experiments/ex_2_crps_1.png", width = 400, height = 300)
plot(grid$Var1, crps_vals, xlab = TeX("Candidate values for $\\mu$"), ylab = "CRPS"); abline(v=mu, col = "blue")
dev.off()

png(filename="plots/synthetic_experiments/ex_2_crps_2.png", width = 400, height = 300)
plot(grid$Var2, crps_vals, xlab = TeX("Candidate values for $\\xi$"), ylab = "CRPS"); abline(v=xi, col = "blue")
dev.off()

plot(grid$Var1, (crps_vals-mean(crps_vals))/sd(crps_vals), xlab = "mu", ylab = "standardized CRPS"); abline(v=xi, col = "blue")
plot(grid$Var2, (crps_vals-mean(crps_vals))/sd(crps_vals), xlab = "xi", ylab = "standardized CRPS"); abline(v=xi, col = "blue")


png(filename="plots/synthetic_experiments/ex_2_twcrps_1.png", width = 400, height = 300)
plot(grid$Var1, twcrps_vals, xlab = TeX("Candidate values for $\\mu$"), ylab = "twCRPS"); abline(v=mu, col = "blue")
dev.off()

png(filename="plots/synthetic_experiments/ex_2_twcrps_2.png", width = 400, height = 300)
plot(grid$Var2, twcrps_vals, xlab = TeX("Candidate values for $\\xi$"), ylab = "twCRPS"); abline(v=xi, col = "blue")
dev.off()

plot(grid$Var1, (twcrps_vals-mean(twcrps_vals))/sd(twcrps_vals), xlab = "mu", ylab = "standardized twCRPS"); abline(v=xi, col = "blue")
plot(grid$Var2, (twcrps_vals-mean(twcrps_vals))/sd(twcrps_vals), xlab = "xi", ylab = "standardized twCRPS"); abline(v=xi, col = "blue")

