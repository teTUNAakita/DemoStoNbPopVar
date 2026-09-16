# ==============================================================================
# Figure 3 and Figure S2: Population and recruitment variability vs. Nb
#
# Purpose:
#   Reproduce Fig. 3 and Fig. S2 for the manuscript
#   "Demographic stochasticity can drive population variability even in large
#   populations".
#
# Description:
#   Runs stochastic female age-structured simulations across reproductive
#   overdispersion (phi) and environmental variability (sigma_E), then compares
#   simulated CV[N] and CV[R] with the analytical approximation.
#
# Output:
#   figure/Fig3_Nb_vs_CVN.pdf
#   figure/FigS2_Nb_vs_CVR.pdf
#
# Dependencies:
#   Base R only.
# ==============================================================================

set.seed(777)

## Simulation settings
phi_values    <- c(1e-5, 2e-5, 5e-5, 1e-4, 2e-4, 5e-4,
                   1e-3, 2e-3, 5e-3, 1e-2, 2e-2, 5e-2, 100)
sigmaE_values <- c(0, 0.3, 0.6)
n_rep <- 20
T <- 10000
burn <- 500

## Life-history and recruitment parameters
age <- 0:9
maturity <- as.integer(age >= 3)
survival <- rep(0.5, length(age))
lambda <- (1 - exp(-0.25 * (age + 1.5)))^3
alpha <- 2000
beta <- 0.001

h <- cumprod(survival) * maturity
gamma <- sqrt(sum(h^2)) / sum(h)

CV <- function(x) sd(x) / mean(x)

## Poisson-lognormal draw with specified arithmetic mean and variance
rPLN <- function(mu, var) {
  if (mu <= 0) return(0L)
  if (var <= mu) return(rpois(1, mu))
  sigma2 <- log(1 + (var - mu) / mu^2)
  lambdaPLN <- rlnorm(1, log(mu) - sigma2 / 2, sqrt(sigma2))
  rpois(1, lambdaPLN)
}

## Run one stochastic age-structured simulation
simulate_population <- function(phi, sigmaE) {

  A <- matrix(0L, nrow = T, ncol = length(age))
  A[1, 1] <- 10000L

  R <- integer(T)
  N <- integer(T)
  Nb <- rep(NA_real_, T)

  for (t in seq_len(T)) {

    ## Binomial survival and mature female abundance
    A[t, ] <- rbinom(length(age), A[t, ], survival)
    adults <- A[t, ] * maturity
    N[t] <- sum(adults)

    if (N[t] >= 2) {
      f <- adults / N[t]
      m1 <- sum(lambda * f)
      m2 <- sum(lambda^2 * f)

      ## Effective number of female breeders
      Nb[t] <- (N[t] - 1) / ((1 + 1 / phi) * m2 / m1^2) + 1

      ## Beverton-Holt mean recruitment
      total_output <- sum(adults * lambda)
      Rbar <- alpha * total_output / (1 + beta * total_output)

      ## Recruitment variance
      if (Nb[t] > 1 && Rbar > 0) {
        D <- 1 + (1 - Nb[t] / N[t]) / (Nb[t] - 1)
        varR <- Rbar + Rbar^2 * (exp(sigmaE^2) * D - 1)
        R[t] <- rPLN(Rbar, max(varR, Rbar))
      }
    }

    if (t < T) {
      A[t + 1, ] <- c(as.integer(R[t]), A[t, 1:(length(age) - 1)])
      if (sum(A[t + 1, ]) == 0) return(c(Nb = NA, CVR = NA, CVN = NA))
    }
  }

  keep <- (burn + 1):T
  c(
    Nb  = mean(Nb[keep], na.rm = TRUE),
    CVR = CV(R[keep]),
    CVN = CV(N[keep])
  )
}

## Run all parameter combinations
design <- expand.grid(
  phi = phi_values,
  sigmaE = sigmaE_values,
  rep = seq_len(n_rep)
)

sim <- t(mapply(
  simulate_population,
  phi = design$phi,
  sigmaE = design$sigmaE
))

results <- cbind(design, as.data.frame(sim))

## Analytical approximation
Nb_grid <- exp(seq(log(2), log(25000), length.out = 200))

CVR_theory <- function(Nb, sigmaE) {
  sqrt(exp(sigmaE^2) * (1 + 1 / (Nb - 1)) - 1)
}

CVN_theory <- function(Nb, sigmaE) {
  gamma * CVR_theory(Nb, sigmaE)
}

cols <- c("black", grey(0.6), grey(0.8))

## Reusable plotting function
plot_variability <- function(yvar, theory_fun, ylim, outfile) {

  pdf(outfile, width = 5, height = 6, pointsize = 12, bg = "white")

  plot(
    NA,
    log = "xy",
    xlim = c(2, 25000),
    ylim = ylim,
    xlab = "",
    ylab = ""
  )

  for (j in seq_along(sigmaE_values)) {
    d <- results[results$sigmaE == sigmaE_values[j], ]

    points(
      d$Nb,
      d[[yvar]],
      pch = 16,
      cex = 0.5,
      col = cols[j]
    )

    lines(
      Nb_grid,
      theory_fun(Nb_grid, sigmaE_values[j]),
      lwd = 2,
      col = cols[j]
    )
  }

  legend(
    "bottomleft",
    legend = c(
      expression(sigma[E] == 0),
      expression(sigma[E] == 0.3),
      expression(sigma[E] == 0.6)
    ),
    pch = 16,
    col = cols,
    pt.cex = 0.8,
    bty = 1
  )

  dev.off()
}

## Save figures in the repository-level figure directory
dir.create("figure", showWarnings = FALSE)

plot_variability(
  yvar = "CVN",
  theory_fun = CVN_theory,
  ylim = c(0.005, 1.1),
  outfile = "figure/Fig3_Nb_vs_CVN.pdf"
)

plot_variability(
  yvar = "CVR",
  theory_fun = CVR_theory,
  ylim = c(0.01, 1.5),
  outfile = "figure/FigS2_Nb_vs_CVR.pdf"
)
