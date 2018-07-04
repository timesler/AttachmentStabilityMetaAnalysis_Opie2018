trimfill.default <- function(yi, vi, weights = NULL, ni = NULL, maxiter = 100, estimator = "L0", verbose = F)
{
  idix <- sort(yi, index.return = TRUE)$ix
  yi <- yi[idix]
  vi <- vi[idix]
  weights <- weights[idix]
  ni <- ni[idix]
  k <- length(yi)
  k0.sav <- -1
  k0 <- 0
  side = "left"
  iter <- 0
  while (abs(k0 - k0.sav) > 0) {
    k0.sav <- k0
    iter <- iter + 1
    if (iter > maxiter) 
      stop("Trim and fill algorithm did not converge.")
    yi.t <- yi[1:(k - k0)]
    vi.t <- vi[1:(k - k0)]
    weights.t <- weights[1:(k - k0)]
    res <- suppressWarnings(rma.uni(yi.t, vi.t, weights = weights.t))
    b <- c(res$b)
    yi.c <- yi - b
    yi.c.r <- rank(abs(yi.c), ties.method = "first")
    yi.c.r.s <- sign(yi.c) * yi.c.r
    if (estimator == "R0") {
      k0 <- (k - max(-1 * yi.c.r.s[yi.c.r.s < 0])) - 1
      se.k0 <- sqrt(2 * max(0, k0) + 2)
    }
    if (estimator == "L0") {
      Sr <- sum(yi.c.r.s[yi.c.r.s > 0])
      k0 <- (4 * Sr - k * (k + 1))/(2 * k - 1)
      varSr <- 1/24 * (k * (k + 1) * (2 * k + 1) + 10 * 
                         k0^3 + 27 * k0^2 + 17 * k0 - 18 * k * k0^2 - 
                         18 * k * k0 + 6 * k^2 * k0)
      se.k0 <- 4 * sqrt(varSr)/(2 * k - 1)
    }
    if (estimator == "Q0") {
      Sr <- sum(yi.c.r.s[yi.c.r.s > 0])
      k0 <- k - 1/2 - sqrt(2 * k^2 - 4 * Sr + 1/4)
      varSr <- 1/24 * (k * (k + 1) * (2 * k + 1) + 10 * 
                         k0^3 + 27 * k0^2 + 17 * k0 - 18 * k * k0^2 - 
                         18 * k * k0 + 6 * k^2 * k0)
      se.k0 <- 2 * sqrt(varSr)/sqrt((k - 1/2)^2 - k0 * 
                                      (2 * k - k0 - 1))
    }
    k0 <- max(0, round(k0))
    se.k0 <- max(0, se.k0)
    if (verbose) 
      cat("Iteration:", iter, "\tmissing =", k0, "\t  b =", 
          ifelse(side == "right", -1 * b, b), "\n")
  }
  if (k0 > 0) {
    if (side == "right") {
      yi.c <- -1 * (yi.c - b)
    }
    else {
      yi.c <- yi.c - b
    }
    yi.fill <- c(yi, -1 * yi.c[(k - k0 + 1):k])
    vi.fill <- c(vi, vi[(k - k0 + 1):k])
    weights.fill <- c(weights, weights[(k - k0 + 1):k])
    ni.fill <- c(ni, ni[(k - k0 + 1):k])
    res <- suppressWarnings(rma.uni(yi.fill, vi.fill, weights = weights.fill, 
                                    ni = ni.fill))
    res$fill <- c(rep(FALSE, k), rep(TRUE, k0))
    res$slab.null <- FALSE
  }
  else {
    res <- list()
    res$fill <- rep(FALSE, k)
  }
  res$k0 <- k0
  res$se.k0 <- se.k0
  res$side <- side
  res$k0.est <- estimator
  if (estimator == "R0") {
    m <- -1:(k0 - 1)
    res$p.k0 <- 1 - sum(choose(0 + m + 1, m + 1) * 0.5^(0 + 
                                                          m + 2))
  }
  else {
    res$p.k0 <- NA
  }
  class(res) <- c("default.trimfill", class(res))
  return(res)
}