# ==============================================================================
# Figure 2: Representative population trajectories
#
# Purpose:
#   Reproduce the four panels of Fig. 2 for the manuscript
#   "Demographic stochasticity can drive population variability even in large
#   populations".
#
# Description:
#   Simulates female age-structured population dynamics under four combinations
#   of reproductive overdispersion (phi) and environmental variability
#   (sigma_E), and saves each panel as a separate PDF.
#
# Output:
#   figure/Fig2a_phi100_sigmaE0.pdf
#   figure/Fig2b_phi100_sigmaE06.pdf
#   figure/Fig2c_phi1e-5_sigmaE0.pdf
#   figure/Fig2d_phi1e-5_sigmaE06.pdf
#
# Dependencies:
#   Base R only.
# ==============================================================================

## Poisson-lognormal draw with specified arithmetic mean and variance
rPLN <- function(mu, var) {
  if (mu <= 0) return(0L)
  if (var <= mu) return(rpois(1, mu))
  sigma2 <- log(1 + (var - mu) / mu^2)
  lambda <- rlnorm(1, log(mu) - 0.5 * sigma2, sqrt(sigma2))
  rpois(1, lambda)
}

## Run one age-structured population simulation
simulate_population <- function(phi, sigmaE, seed = 1, years = 1000) {
  set.seed(seed)

  age <- 0:9
  maturity <- as.integer(age >= 3)
  survival <- rep(0.5, length(age))
  weight <- (1 - exp(-0.25 * (age + 1.5)))^3
  alpha <- 2000
  beta <- 0.001

  A <- matrix(0L, nrow = years, ncol = length(age))
  A[1, ] <- c(10000L, rep(0L, length(age) - 1))

  N <- numeric(years)
  Nb <- numeric(years)
  recruitment <- numeric(years)

  for (t in seq_len(years)) {

    ## Binomial survival at each age
    A[t, ] <- rbinom(length(age), A[t, ], survival)

    ## Mature females contributing to reproduction
    mature <- rbinom(length(age), A[t, ], maturity)
    N[t] <- sum(mature)

    if (N[t] >= 2) {
      freq <- mature / N[t]
      m1 <- sum(weight * freq)
      m2 <- sum(weight^2 * freq)

      ## Effective number of female breeders
      cfac <- (1 + 1 / phi) * m2 / m1^2
      Nb[t] <- min((N[t] - 1) / cfac + 1, N[t])

      ## Mean recruitment from the Beverton-Holt function
      SSB <- sum(mature * weight)
      Rmean <- alpha * SSB / (1 + beta * SSB)

      ## Recruitment variance
      if (Nb[t] > 1 && Rmean > 0) {
        nonPoisson <- 1 + (1 - Nb[t] / N[t]) / (Nb[t] - 1)
        varR <- Rmean + Rmean^2 * (exp(sigmaE^2) * nonPoisson - 1)
        recruitment[t] <- rPLN(Rmean, max(varR, Rmean))
      }
    }

    if (t < years) {
      A[t + 1, ] <- c(as.integer(recruitment[t]), A[t, -length(age)])
    }
  }

  list(N = N, Nb = Nb)
}

## Four parameter combinations corresponding to Fig. 2a-d
scenarios <- data.frame(
  panel  = c("a", "b", "c", "d"),
  phi    = c(1e2, 1e2, 1e-5, 1e-5),
  sigmaE = c(0.0, 0.6, 0.0, 0.6),
  seed   = c(1, 2, 3, 4)
)

## Plot years after the 500-year burn-in
plot_years <- 501:900

## Create output directory in the repository root
dir.create("figure", showWarnings = FALSE)

for (i in seq_len(nrow(scenarios))) {

  sim <- simulate_population(
    phi = scenarios$phi[i],
    sigmaE = scenarios$sigmaE[i],
    seed = scenarios$seed[i]
  )

  phi_label <- if (scenarios$phi[i] == 1e2) "phi100" else "phi1e-5"
  sigma_label <- if (scenarios$sigmaE[i] == 0) "sigmaE0" else "sigmaE06"

  outfile <- file.path(
    "figure",
    paste0("Fig2", scenarios$panel[i], "_", phi_label, "_", sigma_label, ".pdf")
  )

  pdf(outfile, width = 5, height = 5, pointsize = 12, bg = "white")
  plot(
    seq_along(plot_years),
    sim$N[plot_years],
    type = "l",
    lwd = 2,
    ylim = c(0, 1e6),
    xlab = "",
    ylab = ""
  )
  dev.off()
}
