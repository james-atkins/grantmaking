TOLERANCE <- sqrt(.Machine$double.eps)

DISTRIBUTIONS <- list(
  NORMAL = list(density = dnorm, distribution = pnorm, quantile = qnorm, support_lower = -Inf, support_upper = Inf),
  UNIFORM = list(density = dunif, distribution = punif, quantile = qunif, support_lower = 0, support_upper = 1),
  EXPONENTIAL = list(density = dexp, distribution = pexp, quantile = qexp, support_lower = 0, support_upper = Inf)
)

checked_integrate <- function(f, lower, upper) {
  integral <- stats::integrate(f, lower, upper, rel.tol = TOLERANCE, subdivisions = 1000L)
  stopifnot(integral$message == "OK")
  integral$value
}

# A field is parameterised by
# - type distribution
# - noise distribution
# - gamma (cost-benefit ratio)
# - sigma
new_field <- function(type_dist, noise_dist, gamma, sigma) {
  list(
    type_dist = type_dist,
    noise_dist = noise_dist,
    gamma = gamma,
    sigma = sigma
  )
}


# Convenience function that "unpacks" the field into the environment of the calling function.
# In particular, the type and noise distribution functions are assigned to variables using the
# same notation as in the paper, for example the density function of the noise distribution is
# unpacked as "f".
unpack_field <- function(field) {
  env <- parent.frame()

  env$f <- field$noise_dist$density
  env$F <- field$noise_dist$distribution
  env$F_inv <- field$noise_dist$quantile
  env$F_support <- list(lower = field$noise_dist$support_lower, upper = field$noise_dist$support_upper)

  env$g <- field$type_dist$density
  env$G <- field$type_dist$distribution
  env$G_inv <- field$type_dist$quantile
  env$G_support <- list(lower = field$type_dist$support_lower, upper = field$type_dist$support_upper)

  env$gamma <- field$gamma
  env$sigma <- field$sigma
}


# Application demand
a_D <- function(x, field) {
  unpack_field(field)
  stopifnot(sigma > 0)
  1 - G(x - sigma * F_inv(1 - gamma))
}

# Inverse demand
x_D <- function(a, field) {
  unpack_field(field)
  stopifnot(sigma > 0)
  G_inv(1-a) + sigma * F_inv(1 - gamma)
}


.supply_integral <- function(x, a, field) {
  unpack_field(field)
  integrand <- function(v) 1 - F((x - G_inv(1-v)) / sigma)
  checked_integrate(integrand, 0, a)
}


x_P <- function(a, p, field) {
  unpack_field(field)
  stopifnot(0 <= a && a <= 1)
  stopifnot(sigma > 0)
  stopifnot(0 < p && p < 1)

  if (a == 0) {
    return(sigma * F_inv(1-p) + G_inv(1))
  }

  func <- function(x) .supply_integral(x, a, field) - (p * a)

  # TODO: Be better at finding lower and upper bounds for the root finder.
  # How about Newton-Raphson instead?
  uniroot(func, c(-1000, 1000), extendInt = "downX", check.conv = TRUE)$root
}

x_P <- Vectorize(x_P, c("a", "p"))


# Equilibria

.system_equations_p <- function(fields, payline, as) {
  ps <- payline(as)

  purrr::pmap_dbl(list(fields, as, ps), function(field, a, p) {
    unpack_field(field)

    x <- x_D(a, field)

    # Handle non-interior equilibria where supply intersects with the vertical
    # part of demand
    if (abs(a - 0) <= TOLERANCE && x < x_P(0, p, field)) {
      return(0)
    } else if (abs(a - 1) <= TOLERANCE && x > x_P(1, p, field)) {
      return(0)
    }

    val <- 1 - F((x - G_inv(1 - a)) / sigma)

    if (a > 0) {
      integrand <- function(theta) (1 - G(theta)) * f((x - theta) / sigma)
      val <- val + (checked_integrate(integrand, G_inv(1 - a), G_support$upper) / (a * sigma))
    }

    val - p
  })
}


compute_equilibrium <- function(fields, payline) {
  objective <- function(as) {
    sum(.system_equations_p(fields, payline, as) ^ 2)
  }

  # Constrained optimisation: 0 <= a <= 1
  ui <- as.matrix(Matrix::bdiag(rep(list(c(1, -1)), length(fields))))
  ci <- rep(c(0, -1), length(fields))

  result <- constrOptim(
    theta = rep(0.5, length(fields)),
    f = objective,
    ui = ui,
    ci = ci,
    method = "Nelder-Mead",
    control = list(warn.1d.NelderMead = FALSE)
  )

  stopifnot(result$convergence == 0)

  as <- result$par
  as[abs(as - 0) <= TOLERANCE] <- 0
  as[abs(as - 1) <= TOLERANCE] <- 1

  ps <- payline(as)

  list(
    as = as,
    xs = purrr::pmap_dbl(list(as, ps, fields), x_P)
  )
}

constant_payline <- function(ps) {
  function(as) {
    stopifnot(length(as) == length(ps))
    ps
  }
}

qpa <- function(rhos, total_budget) {
  function(as) {
    stopifnot(length(as) == length(rhos))
    as[abs(as - 0) <= TOLERANCE] <- 0
    as[abs(as - 1) <= TOLERANCE] <- 1
    total_budget * as ^ (rhos - 1) / sum(as ^ rhos)
  }
}
