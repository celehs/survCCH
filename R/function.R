#' sCCH resampling
#' 
#' @param J ...
#' @param D ...
#' @param R ...
#' @param V ...
#' @param B ...
#' 
#' @export
ipw.scch <- function(J, D, R, V, B) {
  # browser()
  J <- factor(J)
  N <- length(J)
  M <- t(matrix(rep(levels(J), N), ncol = N))
  I <- 1 * (J == M)
  tmp <- colSums(B * R * I) / colSums(B * I)
  p.hat <- D + (1 - D) * (I %*% tmp)
  V * B / p.hat
}


#' Calculating the baseline cumulative hazard 
#' 
#' @param X ...
#' @param Z ...
#' @param beta ...
#' @param t ...
#' @param order ...
#' @param W ...
#' 
#' @export
PI.fun <- function(X, Z, beta, t = Inf, order = c(0, 1, 2), W = NULL) {
  # browser()
  if (is.null(W)) W <- rep(1, length(X))
  Z <- as.matrix(Z)
  theta <- c(exp(Z %*% beta))
  M <- outer(X, t, ">=")
  PI0 <- PI1 <- PI2 <- NULL
  if (0 %in% order) PI0 <- c(t(W * theta) %*% M)
  if (1 %in% order) PI1 <- t(t(W * theta * Z) %*% M)
  if (2 %in% order) {
    p <- NCOL(Z)
    Z2 <- Z[, rep(1:p, each = p)] * Z[, rep(1:p, p)]
    PI2 <- t(t(W * theta * Z2) %*% M)
  }
  list(PI0 = PI0, PI1 = PI1, PI2 = PI2)  
}


#' Main function for sCCH analysis
#' 
#' @param N ...
#' @param X ...
#' @param D ...
#' @param Z ...
#' @param R ...
#' @param J ...
#' @param n.ptb ...
#' @param Z.new ... 
#' 
#' @export
survcox.scch <- function(N, X, D, Z, R, J, n.ptb = 1000, Z.new = NULL) {
  # browser()
  n <- table(J)
  size <- cbind(N, n)
  M1 <- cbind(J = as.integer(J), D, R, V = 1)
  M0 <- cbind(J = as.integer(rep(names(n), N - n)), D = 0, R = 0, V = 0)
  M.all <- rbind(M1, M0)
  W <- ipw.scch(M.all[, "J"], M.all[, "D"], M.all[, "R"], M.all[, "V"], rep(1, NROW(M.all)))[1:NROW(M1)]
  fit <- coxph(Surv(X, D) ~ Z, weights = W, ties = "breslow")
  if (is.null(Z.new)) Z.new <- fit$mean
  coef <- fit$coef
  t <- seq(0.6, 2.8, length = 200)
  PI <- PI.fun(X, Z, coef, X, order = c(0, 1), W = W)
  Z.bar <- (1/PI$PI0) * PI$PI1
  M <- outer(t, X, ">=")
  # browser()
  log.cumhaz <- log(M %*% (W * D / PI$PI0)) + sum(Z.new * coef)
  U <- M %*% (W * D * (Z - Z.bar))
  se <- sqrt(diag(fit$var))
  U.std <- t(se * t(U))  
  Q <- apply(abs(U.std), 2, max)
  Q.max <- max(Q)
  p <- length(coef)
  m <- length(t)
  log.cumhaz.ptb <- matrix(NA, n.ptb, m)
  U.std.ptb <- array(NA, dim = c(m, p, n.ptb))
  Q.ptb <- matrix(NA, n.ptb, p)
  Q.max.ptb <- rep(NA, n.ptb)
  coef.ptb <- matrix(NA, n.ptb, p)  
  # perturbation resampling
  for (i in 1:n.ptb) {
    W.ptb <- ipw.scch(M.all[, "J"], M.all[, "D"], M.all[, "R"], M.all[, "V"], rexp(NROW(M.all)))[1:NROW(M1)]
    fit.ptb <- coxph(Surv(X, D) ~ Z, weights = W.ptb, ties = "breslow", init = coef)
    coef.ptb[i, ] <- fit.ptb$coef
    PI.ptb <- PI.fun(X, Z, coef.ptb[i, ], X, order = c(0, 1), W = W.ptb)
    Z.bar.ptb <- (1/PI.ptb$PI0) * PI.ptb$PI1
    log.cumhaz.ptb[i, ] <- log(M %*% (W.ptb * D / PI.ptb$PI0)) + sum(Z.new * coef.ptb[i, ])
    U.ptb <- M %*% (W.ptb * D * (Z - Z.bar.ptb))
    U.std.ptb[, , i] <- t(se * t(U.ptb))
    Q.ptb[i, ] <- apply(abs(U.std.ptb[, , i] - U.std), 2, max)
    Q.max.ptb[i] <- max(Q.ptb[i, ])
  }  
  p.value.PH <- colMeans((Q.ptb > outer(rep(1, n.ptb), Q)))
  pval.omni <- mean(Q.max.ptb > Q.max)
  coef.se <- apply(coef.ptb, 2, sd)
  summary <- data.frame(coef = coef, 
                        coef.se = coef.se, 
                        HR = exp(coef),
                        HR.lower = exp(coef - 1.96 * coef.se),
                        HR.upper = exp(coef + 1.96 * coef.se),
                        p.value.HR = 2 * pnorm(-abs(coef / coef.se)),
                        p.value.PH = p.value.PH)  
  list(size = size, coef = coef, coef.ptb = coef.ptb, summary = round(summary, 3), 
       t = t, log.cumhaz = log.cumhaz, log.cumhaz.ptb = log.cumhaz.ptb, 
       U.std = U.std, U.std.ptb = U.std.ptb, pval.omni = pval.omni)
}

