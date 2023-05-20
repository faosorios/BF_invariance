## Id: invar_LW.R, last updated May 5, 2023
## Author: originally coded by Federico Crudu, with contributions of Felipe Osorio

simul.LW <- function(Nsize = 5000, nobs = 500, coef, sigma = 1, k = 2, alpha = 0.05, trace = TRUE, msg = NULL)
{ ## function to perform the simulation experiment (Section 4 of the manuscript)
  constr.A <- function(coef) {
    coef[2] - 1.0
  }
  constr.B <- function(coef) {
    coef[2]^k - 1.0
  }
  deriv.A <- function(coef) {
    c(0.0, 1.0)
  }
  deriv.B <- function(coef) {
    c(0.0, k * coef[2]^(k - 1))
  }
  objective.fnc <- function(coef) {
    res <- y - coef[1] * x1 - coef[2] * x2
    sum(c(res)^2)
  }
  Gplus <- function(coef) {
    g <- c(0.0, (1 / k) * coef[2]^(1 - k))
    matrix(g, ncol = 1)
  }
  Kmat <- function(coef) {
    mat <- matrix(0, nrow = 2, ncol = 2)
    mat[1,1] <- 1.0
    mat[2,2] <- k * coef[2]^(k - 1)
    mat
  }

  if (is.null(msg))
    msg <- "Progress"

  # results containers
  ok <- matrix(FALSE, nrow = Nsize, ncol = 7)
  now <- proc.time()

  if (trace) {
    cat(" ", paste(msg, ":", sep = ""), "\n")
    pb <- txtProgressBar(min = 0, max = Nsize, style = 3)
  }

  # setting 'x' matrix
  set.seed(29)
  x1 <- runif(nobs)
  x2 <- runif(nobs)
  x  <- cbind(x1, x2)
  xx <- crossprod(x)
  V <- solve(xx)

  # Monte Carlo iterations
  for (i in 1:Nsize) {
    set.seed(123 + i)

    # model building
    coef.true <- coef
    eps <- rnorm(nobs, 0, sigma)
    y <- coef.true[1] *x1 + coef.true[2] * x2 + eps

    # fitting unconstrained model
    fm0 <- lsfit(x, y, intercept = FALSE)
    cf0 <- fm0$coef 
    res0 <- fm0$resid

    # fitting constrained models
    fm1 <- auglag(par = coef.true, fn = objective.fnc, heq = constr.A, control.outer = list(trace = FALSE))
    cf1 <- fm1$par
    res1 <- y - cf1[1] * x1 - cf1[2] * x2
    fm2 <- auglag(par = coef.true, fn = objective.fnc, heq = constr.B, control.outer = list(trace = FALSE))
    cf2 <- fm2$par
    res2 <- y - cf2[1] * x1 - cf2[2] * x2

    # distance metric statistic
    DM <- objective.fnc(cf1) - objective.fnc(cf0)

    # Wald statistic (hypothesis A)
    a1 <- constr.A(cf0)
    a1.dot <- matrix(deriv.A(cf0), nrow = 1)
    Wald.A <- c(a1^2 / (a1.dot %*% V %*% t(a1.dot)))

    # LM statistic (hypothesis A)
    U1 <- crossprod(x, res1)
    LM <- c(crossprod(U1, V %*% U1))

    # BF statistic (hypothesis A)
    a1.plus <- t(a1.dot) %*% solve(a1.dot %*% t(a1.dot))
    prod <- c(crossprod(U1, a1.plus))
    BF.A <- prod * a1

    # Wald statistic (hypothesis B)
    a2 <- constr.B(cf0)
    a2.dot <- matrix(deriv.B(cf0), nrow = 1)
    Wald.B <- c(a2^2 / (a2.dot %*% V %*% t(a2.dot)))

    # BF statistic (hypothesis B)
    a2.plus <- t(a2.dot) %*% solve(a2.dot %*% t(a2.dot))
    prod <- c(crossprod(U2, a2.plus))
    BF.B <- prod * a2

    # BF statistic (corrected)
    K2 <- Kmat(cf2)
    U2 <- crossprod(K2, U2)
    a.plus <- Gplus(cf2)
    prod <- c(crossprod(U2, a.plus))
    BF <- prod * a1

    # saving results
    cutoff <- qchisq(1 - alpha, 1)

    ok[i,1] <- Wald.A > cutoff
    ok[i,2] <- Wald.B > cutoff
    ok[i,3] <- BF.A > cutoff
    ok[i,4] <- BF.B > cutoff
    ok[i,5] <- BF > cutoff
    ok[i,6] <- LM > cutoff
    ok[i,7] <- DM > cutoff

    # update progress bar
    if (trace)
      setTxtProgressBar(pb, i)
  }
  if (trace)
    close(pb)

  pnames <- c("Wald.A","Wald.B","BF.A","BF.B","BF","LM","D")
  percentage <- rep(0, 7)
  percentage[1] <- sum(ok[,1], na.rm = TRUE) / Nsize
  percentage[2] <- sum(ok[,2], na.rm = TRUE) / Nsize
  percentage[3] <- sum(ok[,3], na.rm = TRUE) / Nsize
  percentage[4] <- sum(ok[,4], na.rm = TRUE) / Nsize
  percentage[5] <- sum(ok[,5], na.rm = TRUE) / Nsize
  percentage[6] <- sum(ok[,6], na.rm = TRUE) / Nsize
  percentage[7] <- sum(ok[,7], na.rm = TRUE) / Nsize
  #percentage <- apply(ok, 2, sum) / Nsize
  bad <- rep(0, 7)
  if (any(is.na(ok[,1]))) bad[1] <- sum(is.na(ok[,1]))
  if (any(is.na(ok[,2]))) bad[2] <- sum(is.na(ok[,2]))
  if (any(is.na(ok[,3]))) bad[3] <- sum(is.na(ok[,3]))
  if (any(is.na(ok[,4]))) bad[4] <- sum(is.na(ok[,4]))
  if (any(is.na(ok[,5]))) bad[5] <- sum(is.na(ok[,5]))
  if (any(is.na(ok[,6]))) bad[6] <- sum(is.na(ok[,6]))
  if (any(is.na(ok[,7]))) bad[7] <- sum(is.na(ok[,7]))
  names(percentage) <- pnames
  names(bad) <- pnames
  speed <- proc.time() - now
  out <- list(percentage = 100 * percentage, bad = bad, cutoff = cutoff, speed = speed)
  out
}

summary.LW <- function(Nsize = 10000, coef, sigma, k = 2, alpha = 0.05, trace = TRUE) {
  out <- matrix(0, nrow = 4, ncol = 7)
  bad <- matrix(0, nrow = 4, ncol = 7)
  now <- proc.time()

  o <- simul.LW(Nsize, nobs = 25,  coef, sigma, k, alpha, trace, msg = "Stage 1 of 4")
  out[1,] <- o$percentage
  bad[1,] <- o$bad
  o <- simul.LW(Nsize, nobs = 50,  coef, sigma, k, alpha, trace, msg = "Stage 2 of 4")
  out[2,] <- o$percentage
  bad[2,] <- o$bad
  o <- simul.LW(Nsize, nobs = 100, coef, sigma, k, alpha, trace, msg = "Stage 3 of 4")
  out[3,] <- o$percentage
  bad[3,] <- o$bad
  o <- simul.LW(Nsize, nobs = 500, coef, sigma, k, alpha, trace, msg = "Stage 4 of 4")
  out[4,] <- o$percentage
  bad[4,] <- o$bad
  cutoff <- o$cutoff

  rownames(out) <- c("25","50","100","500")
  colnames(out) <- c("Wald.A","Wald.B","BF.A","BF.B","BF","LM","D")
  rownames(bad) <- c("25","50","100","500")
  colnames(bad) <- c("Wald.A","Wald.B","BF.A","BF.B","BF","LM","D")

  list(output = out, bad = bad, cutoff = cutoff, speed = proc.time() - now)
}
