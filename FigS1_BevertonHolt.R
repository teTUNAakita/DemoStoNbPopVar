# ==============================================================================
# Figure S1: Beverton-Holt recruitment function and replacement line
#
# Purpose:
#   Reproduce Fig. S1 for the manuscript
#   "Demographic stochasticity can drive population variability even in large
#   populations".
#
# Description:
#   This script plots the Beverton-Holt recruitment function used in the
#   age-structured simulations together with the replacement line. Their
#   intersection gives the deterministic equilibrium.
#
# Model:
#   Mean recruitment:
#       R = alpha * x / (1 + beta * x)
#
#   where x is the total expected reproductive output. In the fish-based
#   simulation used in the manuscript, x corresponds to spawning stock biomass
#   because individual expected reproductive output is represented by body
#   weight.
#
#   At demographic equilibrium:
#       x = SPR * R
#
#   where SPR is spawning biomass per recruit under the specified survival,
#   maturation, and body-weight schedules.
#
# Output:
#   figures/FigS1_BevertonHolt.pdf
#
# Dependencies:
#   Base R only.
#
# Notes:
#   Parameter values are identical to those used in the simulation analysis.
# ==============================================================================


## -----------------------------------------------------------------------------
## 1. Model parameters
## -----------------------------------------------------------------------------

alpha <- 2000
beta  <- 0.001

age      <- 0:9
maturity <- as.integer(age >= 3)
survival <- rep(0.5, length(age))
weight   <- (1 - exp(-0.25 * (age + 1.5)))^3


## -----------------------------------------------------------------------------
## 2. Spawning biomass per recruit and deterministic equilibrium
## -----------------------------------------------------------------------------

# Spawning biomass per recruit (SPR) under the specified life-history schedule.
SPR <- sum(weight * maturity * cumprod(survival))

# Deterministic equilibrium obtained from the intersection of the
# Beverton-Holt recruitment curve and the replacement line.
x_star <- (alpha * SPR - 1) / beta
R_star <- x_star / SPR


## -----------------------------------------------------------------------------
## 3. Recruitment and replacement functions
## -----------------------------------------------------------------------------

# x denotes total expected reproductive output and corresponds to spawning
# stock biomass in this implementation.
x <- seq(0, x_star * 2, length.out = 500)

# Beverton-Holt recruitment function.
BH <- alpha * x / (1 + beta * x)

# Replacement line: recruitment required to maintain spawning biomass x.
replacement <- x / SPR


## -----------------------------------------------------------------------------
## 4. Plot Figure S1
## -----------------------------------------------------------------------------

# Create an output directory if it does not already exist.
dir.create("figures", showWarnings = FALSE)

pdf(
  "figures/FigS1_BevertonHolt.pdf",
  width = 5,
  height = 6,
  pointsize = 12,
  bg = "white"
)

plot(
  x, BH,
  type = "l",
  lwd = 2,
  xlab = "",
  ylab = ""
)

lines(
  x, replacement,
  lty = 2,
  lwd = 2
)

points(
  x_star, R_star,
  pch = 16
)

legend(
  "bottomright",
  legend = c(
    "Beverton-Holt",
    "Replacement line",
    "Equilibrium"
  ),
  lty = c(1, 2, NA),
  pch = c(NA, NA, 16),
  lwd = c(2, 2, NA),
  bty = "n"
)

dev.off()
