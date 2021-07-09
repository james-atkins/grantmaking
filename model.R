
DISTRIBUTIONS <- list(
  NORMAL = list(density = dnorm, distribution = pnorm, quantile = qnorm, support_lower = -Inf, support_upper = Inf),
  UNIFORM = list(density = dunif, distribution = punif, quantile = qunif, support_lower = 0, support_upper = 1),
  EXPONENTIAL = list(density = dexp, distribution = pexp, quantile = qexp, support_lower = 0, support_upper = Inf)
)

checked_integrate <- function(f, lower, upper) {
  integral <- stats::integrate(f, lower, upper, rel.tol = .Machine$double.eps^0.5, subdivisions = 1000L)
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
  
  func <- function(x) .supply_integral(x, a, field) - (p * a)
  
  # TODO: Be better at finding lower and upper bounds for the root finder.
  # How about Newton-Raphson instead?
  uniroot(func, c(-1000, 1000), extendInt = "downX", check.conv = TRUE)$root
}

x_P <- Vectorize(x_P, "a")


# Equilibria

.field_equations <- function(field, x, a) {
  unpack_field(field)
  
  # Remember to handle the non-interior parts of the demand curve
  x_0 <- x_D(0, field)
  x_1 <- x_D(1, field)

  if (a == 0 && x >= x_0) {
    demand <- 0
  } else if (a == 1 && x <= x_1) {
    demand <- 0
  } else {
    demand <- 1 - F((x - G_inv(1-a)) / sigma) - gamma
  }

  supply <- .supply_integral(x, a, field) / a
  
  c(demand, supply)
}

.system_equations <- function(fields, payline, xs, as) {
  # Compute vector of demand and supply for each field
  equations <- purrr::flatten_dbl(purrr::pmap(list(fields, xs, as), .field_equations))
  
  # Subtract p from the supply equations
  ps <- payline$func(as)
  for (i in 1:length(ps)) {
    equations[i*2] <- equations[i*2] - ps[i]
  }
  
  equations
}

# Compute the Jacobian for each field
.field_jacobian <- function(field, x, a) {
  unpack_field(field)
  
  D_x <- -(1 / sigma) * f((x - G_inv(1-a)) / sigma)
  D_a <- -(1 / sigma) * f((x - G_inv(1-a)) / sigma) / g(G_inv(1-a))
  
  # TODO: Integration by substitution here?
  P_x <- -(1 / (a * sigma)) * checked_integrate(
    function(theta) f((x-theta) / sigma) * g(theta),
    G_inv(1-a),
    G_support$upper
  )
  
  P_a <- -(1/a^2 * sigma) * checked_integrate(
    function(theta) (1 - G(theta)) * f((x-theta) / sigma),
    G_inv(1-a),
    G_support$upper
  )
  
  matrix(c(D_x, D_a, P_x, P_a), nrow = 2, ncol = 2, byrow = TRUE)
}


.system_jacobian <- function(fields, payline, xs, as) {
  # First build a block diagonal matrix of each field's demand/supply system
  jacobian <- as.matrix(Matrix::bdiag(purrr::pmap(list(fields, xs, as), .field_jacobian)))
  
  # Now subtract the payline Jacobian
  payline_jacobian <- payline$jacobian(as)
  for (i in 1:nrow(payline_jacobian)) {
    for (j in 1:ncol(payline_jacobian)) {
      jacobian[i*2, j*2] <- jacobian[i*2, j*2] - payline_jacobian[i, j]
    }
  }

  jacobian
}

compute_equilibrium <- function(fields, payline) {
  get_xs <- function(vars) vars[1:(length(vars)/2)]
  get_as <- function(vars) vars[(length(vars)/2+1):length(vars)]
  
  func <- function(vars) {
    .system_equations(fields, payline, get_xs(vars), get_as(vars))
  }
  
  jacobian <- function(vars) {
    .system_jacobian(fields, payline, get_xs(vars), get_as(vars))
  }
  
  # Start on the demand curve at a=1/2
  start_as <- rep(0.5, length(fields))
  start_xs <- purrr::map2_dbl(start_as, fields, x_D)
  
  solution <- rootSolve::multiroot(
    f = func,
    start = c(start_xs, start_as),
    jactype = "fullusr",
    jacfunc = jacobian
  )
  
  # TODO: Check that we actually found an equilibrium
  
  list(
    xs = get_xs(solution$root),
    as = get_as(solution$root)
  )
}

constant_payline <- function(ps) {
  list(
    func = function(as) {
      stopifnot(length(as) == length(ps))
      ps
    },
    
    jacobian = function(as) {
      stopifnot(length(as) == length(ps))
      matrix(0, nrow = length(ps), ncol = length(ps))
    }
  )
}

qpa <- function(rhos, total_budget) {
  list(
    func = function(as) {
      stopifnot(length(as) == length(rhos))
      total_budget * as ^ (rhos - 1) / sum(as ^ rhos)
    },
    
    jacobian = function(as) {
      stopifnot(length(as) == length(rhos))
      (diag((rhos-1) * as^(rhos-2) * sum(as^rhos)) - (outer(as^(rhos-1), as^(rhos-1)) %*% diag(rhos))) / sum(as^rhos)^2
    }
  )
}

