# [ GLM checks ]

hinkley_test <- function (m, test) {
  md <- if (class(m$data) == "data.frame") within(m$data, eta <- predict(m))
    else data.frame(eta = predict(m))
  if (missing(test)) {
    mf <- m$family$family
    test <- ifelse(mf == "binomial" || mf == "poisson", "Chisq", "F")
  }
  anova(m, update(m, ~ . + I(eta ^ 2), data = md), test = test)
}

residual_plot <- function (m, type = "studentized", nl = 3L) {
  mu <- fitted(m)
  r <- if (type == "studentized") rstudent(m) else rstandard(m)
  ii <- sort.list(abs(r), decreasing = TRUE)[1:nl]
  plot(mu, r, xlab = "fitted values", ylab = paste(type, "residuals"))
  abline(h = 0, lty = 3)
  lines(lowess(mu, r), col = "red")
  text(mu[ii], r[ii], labels = ii, cex = 0.8, offset = 0.25,
       pos = ifelse(mu[ii] < mean(mu), 4, 2)) # left (2) or right (4)?
  invisible(r)
}

cook_plot <- function (m) {
  plot(m, 5)
  abline(v = 1 - df.residual(m) / length(m$y), lty = 3, col = "gray")
}

partial_residual_plot <- function (m, i) {
  tl <- attr(m$terms, "term.labels")
  if (i < 1 || i > length(tl)) stop("invalid term index: ", i)
  x <- model.matrix(m)
  ii <- which(attr(x, "assign") == i)
  z <- residuals(m, type = "working") + x[, ii, drop = FALSE] %*% coef(m)[ii]
  xi <- m$model[[tl[i]]]
  mr <- glm(z ~ xi)
  if (class(xi) == "factor") {
    boxplot(z ~ xi, xlab = tl[i], ylab = "partial residuals")
    points(xi, fitted(mr), pch = 19)
  } else {
    plot(xi, z, xlab = tl[i], ylab = "partial residuals")
    abline(mr, lty = 3)
    lines(lowess(xi, z), col = "red")
  }
  invisible(mr)
}

# partial regression plots (added variable plots):
# if variable shouldn't be in the model, slope ~ 0
chol.solve <- function (C, y) backsolve(C, backsolve(C, y, transpose = TRUE))
partial_regression_plot <- function (m, i, alpha = .05) {
  tl <- attr(m$terms, "term.labels")
  if (i < 1 || i > length(tl)) stop("invalid term index: ", i)
  x <- sqrt(m$weights) * model.matrix(m)
  ii <- which(attr(x, "assign") == i)
  if (length(ii) > 1) stop("factor terms are not allowed")
  z <- residuals(m, type = "working") + predict(m) # adjusted dep variable
  z <- sqrt(m$weights) * z # weight-scaled
  xi <- x[, -ii]
  # adjust response and covariate
  C <- chol(crossprod(xi))
  zr <- z - xi %*% chol.solve(C, crossprod(xi, z))
  ur <- x[, ii] - xi %*% chol.solve(C, crossprod(xi, x[, ii]))
  plot(ur, zr, xlab = paste0(tl[i], " | ", paste(tl[-i], collapse = ", ")),
       ylab = paste0("adj resp | ", paste(tl[-i], collapse = ", ")))
  abline(h = 0, lty = 3); abline(v = 0, lty = 3)
  # plot marginal fit
  mr <- glm(zr ~ ur - 1)
  abline(mr, lty = 2)
  mf <- m$family$family
  df <- ifelse(mf == "binomial" || mf == "poisson", Inf, df.residual(m))
  crit_value <- qt(1 - alpha / 2, df)
  sr <- sqrt(diag(vcov(m))[ii])
  abline(0, coef(mr) - crit_value * sr, lty = 2, col = "gray") # lower conf
  abline(0, coef(mr) + crit_value * sr, lty = 2, col = "gray") # upper conf
  lines(lowess(ur, zr), col = "red")
  invisible(mr)
}

# Deviance curves for parametric shape
# `df.dev` = Inf is Chi-squared test, i.e., assuming dispersion = 1
deviance_curve <- function (param, dev, df.dev = Inf, alpha = .95) {
  n <- length(param)
  if (length(dev) != n) stop("inputs have different lenghts")
  plot(param, dev, type = "l", ylab = "deviance", xlab = "parameter")
  imin <- which.min(dev)
  crit_value <- qf(alpha, 1, df.dev)
  phi <- ifelse(is.infinite(df.dev), 1, dev[imin] / df.dev)
  dev_conf <- dev[imin] + phi * crit_value
  abline(h = dev_conf, lty = 2)
  lines(rep(param[imin], 2), c(0, dev[imin]), lty = 2)
  param_lower <- param[1:imin][which.min(abs(dev[1:imin] - dev_conf))]
  param_upper <- param[imin:n][which.min(abs(dev[imin:n] - dev_conf))]
  lines(rep(param_lower, 2), c(0, dev_conf), lty = 3)
  lines(rep(param_upper, 2), c(0, dev_conf), lty = 3)
  invisible(list(lower = param_lower, upper = param_upper,
                 param = param[imin], dev = dev[imin], dev_conf =dev_conf))
}

# Power link and variance
power_link <- function (lambda = 1) {
  if (lambda == 0) return(make.link("log"))
  structure(list(
    name = paste0("power(", lambda, ")"),
    linkfun = function (mu) (mu ^ lambda - 1) / lambda,
    linkinv = function (eta)
      pmax((lambda * eta + 1) ^ (1 / lambda), .Machine$double.eps),
    mu.eta = function (eta)
      pmax((lambda * eta + 1) ^ (1 / lambda - 1), .Machine$double.eps),
    valideta = function (eta)
      all(is.finite(eta)) && all(eta > 0) # -1 / lambda
  ), class = "link-glm")
}

int_p <- function (p) { # int_y^mu dt / t^p
  if (p == 1) # corner case?
    function (y, mu) log(mu / y)
  else
    function (y, mu) (mu ^ (1 - p) - y ^ (1 - p)) / (1 - p)
}
power_var <- function (p) {
  list(
    name = paste0("power(", p, ")"),
    varfun = function(mu) mu ^ p,
    validmu = function(mu) all(mu > 0),
    dev.resids = function(y, mu, wt)
      -2 * wt * (y * int_p(p)(y, mu) - int_p(p - 1)(y, mu)),
    initialize = expression({mustart <- y + .1 * (y == 0)})
  )
}

