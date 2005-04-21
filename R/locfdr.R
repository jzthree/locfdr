locfdr <- function(zz, bre = 120, df = 7, pct = 1/1000, pct0 = 1/3, nulltype = 1,
                   type = 0, plot = 1) {
  ## use  help(locfdr) for definitions and suggestions
  require("splines")
  if(length(bre) > 1) {
    lo <- min(bre)
    up <- max(bre)
    bre <- length(bre)
  }
  else {
    if(length(pct) > 1) {
      lo <- pct[1]
      up <- pct[2]
    }
    else {
      if(pct == 0) {
        lo <- min(zz)
        up <- max(zz)
      }
      else {
        v <- quantile(zz, c(pct, 1 - pct))
        lo <- v[1]
        up <- v[2]
      }
    }
  }
  zzz <- pmax(pmin(zz, up), lo)
  breaks <- seq(lo, up, length = bre)
  zh <- hist(zzz, breaks = breaks, plot = F)
  x <- (breaks[-1] + breaks[ - length(breaks)])/2
  yall <- y <- zh$counts
  K <- length(y)
  N <- length(zz)
  if(pct > 0) {
    y[1.] <- min(y[1.], 1.)
    y[K] <- min(y[K], 1.)
  }
  if(type == 0) {
    f <- glm(y ~ ns(x, df = df), poisson)$fit
  }
  else {
    f <- glm(y ~ poly(x, df = df), poisson)$fit
  }
  l <- log(f)
  Fl <- cumsum(f)
  Fr <- cumsum(rev(f))
  D <- (y - f)/(f + 1)^0.5
  D <- sum(D[2:(K - 1)]^2)/(K - 2 - df)
  if(D > 1.5) {
    print(paste("CHECK FIT, INCREASE DF?  MISFIT=", round(
                                                          D, 1)))
  }
  ## ..............begin f0 calcs................................
  imax <- seq(l)[l == max(l)][1]
  xmax <- x[imax]
  if(length(pct0) == 1) {
    pctup <- 1 - pct0
    pctlo <- pct0
  }
  else {
    pctlo <- pct0[1]
    pctup <- pct0[2]
  }
  lo0 <- quantile(zz[zz < xmax], pctlo)
  hi0 <- quantile(zz[zz > xmax], pctup)
  nx <- length(x)
  i0 <- (1.:nx)[x > lo0 & x < hi0]
  x0 <- x[i0]
  y0 <- l[i0]
  if(nulltype == 2) {
    X00 <- cbind((x0 - xmax)^2, pmax(x0 - xmax, 0)^2)
  }
  else {
    X00 <- cbind(x0 - xmax, (x0 - xmax)^2)
  }
  lr <- lm(y0 ~ X00)
  co <- lr$coef
  if(nulltype == 2) {
    X0 <- cbind(1, (x - xmax)^2, pmax(x - xmax, 0)^2)
    sigs <- 1/sqrt(-2 * (c(co[2], co[2] + co[3])))
    fp0 <- c(xmax, sigs)
  }
  else {
    X0 <- cbind(1, x - xmax, (x - xmax)^2)
    sighat <- 1./sqrt(-2. * co[3.])
    xmaxx <-  - co[2.]/(2. * co[3.]) + xmax
    fp0 <- c(xmaxx, sighat)
  }
  l0 <- as.vector(X0 %*% co)
  f0 <- exp(l0)
  fdr <- pmin(f0/f, 1)
  p0 <- sum(f0)/sum(f)
  f0 <- f0/p0
  fp0 <- c(fp0, p0)
  if(nulltype == 2)
    names(fp0) <- c("zmax", "sigleft", "sigright", "p0")
  else names(fp0) <- c("zmax", "sig", "p0")
  f00 <- exp( - x^2/2)
  f00 <- (f00 * sum(f[i0]))/sum(f00[i0])
  fdr0 <- pmin(f00/f, 1)
  p0theo <- sum(f00)/sum(f)
  f00 <- f00/p0theo
  fp0 <- c(fp0, p0theo)
  names(fp0)[length(fp0)] <- "p0theo"
  f0p <- p0 * f0
  if(nulltype == 0)
    f0p <- p0theo * f00
  F0l <- cumsum(f0p)
  F0r <- cumsum(rev(f0p))
  Fdrl <- F0l/Fl
  Fdrr <- rev(F0r/Fr)
  Int <- (1 - fdr) * f * (fdr < 0.90000000000000002)
  ## raise fdr to 1 near xmax, do Efdr calcs...............................
  xlo <- min(x[x <= xmax & fdr > 0.5])
  xhi <- max(x[x >= xmax & fdr > 0.5])
  xxlo <- min(x[x <= xmax & fdr == 1])
  xxhi <- max(x[x >= xmax & fdr == 1])
  fdr[x >= xxlo & x <= xxhi] <- 1
  xlo <- min(x[x <= xmax & fdr0 > 0.5])
  xhi <- max(x[x >= xmax & fdr0 > 0.5])
  xxlo <- min(x[x <= xmax & fdr0 == 1])
  xxhi <- max(x[x >= xmax & fdr0 == 1])
  fdr0[x >= xxlo & x <= xxhi] <- 1
  p1 <- sum((1 - fdr) * f)/N
  p1theo <- sum((1 - fdr0) * f)/N
  fall <- f + (yall - y)
  Efdr <- sum((1 - fdr) * fdr * fall)/sum((1 - fdr) * fall)
  Efdrtheo <- sum((1 - fdr0) * fdr0 * fall)/sum((1 - fdr0) * fall)
  iup <- (1:K)[x >= xmax]
  ido <- (1:K)[x <= xmax]
  Eleft <- sum((1 - fdr[ido]) * fdr[ido] * fall[ido])/
    sum((1 - fdr[ido]) * fall[ido])
  Eleft0 <- sum((1 - fdr0[ido]) * fdr0[ido] * fall[ido])/
    sum((1 - fdr0[ido]) * fall[ido])
  Eright <- sum((1 - fdr[iup]) * fdr[iup] * fall[iup])/
    sum((1 - fdr[iup]) * fall[iup])
  Eright0 <- sum((1 - fdr0[iup]) * fdr0[iup] * fall[iup])/
    sum((1 - fdr0[iup]) * fall[iup])
  Efdr <- c(Efdr, Eleft, Eright, Efdrtheo, Eleft0, Eright0)
  names(Efdr) <- c("Efdr", "Eleft", "Eright", "Efdrtheo", "Eleft0",
                   "Eright0")
  if(nulltype == 0)
    f1 <- (1 - fdr0) * fall
  else f1 <- (1 - fdr) * fall
  ##................. Accuracy Calcs .................................
  if(type == 0) X <- cbind(1, ns(x, df = df)) else X <- cbind(1,
                                      poly(x, df = df))
  if(nulltype == 0)
    X0 <- matrix(1, length(x), 1)
  Xtil <- X[i0,  ]
  X0til <- X0[i0,  ]
  G <- t(X) %*% (f * X)
  M <- solve(G) %*% t(X)
  A <- X %*% M
  G0 <- t(X0til) %*% X0til
  B0 <- X0 %*% (solve(G0) %*% t(X0til)) %*% Xtil
  C <- B0 - X
  Cov <- C %*% solve(G) %*% t(C)
  lfdrse <- diag(Cov)^0.5
  ## find cdf1, the cdf of fdr according to f1 density..................
  p1 <- seq(0.01, 0.98999999999999999, 0.01)
  cdf1 <- rep(0, 99)
  fd <- fdr
  if(nulltype == 0)
    fd <- fdr0
  for(i in 1:99)
    cdf1[i] <- sum(f1[fd <= p1[i]])
  cdf1 <- cbind(p1, cdf1/cdf1[99])
  mat <- cbind(x, fdr, Fdrl, Fdrr, f, f0, f00, fdr0, yall, lfdrse,
               f1)
  namat <- c("x", "fdr", "Fdrleft", "Fdrright", "f", "f0", "f0theo",
             "fdrtheo", "counts", "lfdrse", "f1")
  if(nulltype == 0)
    namat[c(3, 4, 10)] <- c("Fdrltheo", "Fdrrtheo", 
                            "lfdrsetheo")
  dimnames(mat) <- list(NULL, namat)
  if(plot > 0) {
    if (plot > 1) {
      par(mfrow= c(1,2))
    }

    hist(zzz, breaks = breaks, xlab = " ")
    yt <- yall * (1 - fd)
    for(k in 1:K)
      lines(c(x[k], x[k]), c(0, yt[k]), lwd = 2, col = 6
            )
    if(nulltype == 2)
      title(xlab = paste("zmax=", round(xmax, 3), 
              "sigleft=", round(sigs[1], 3), 
              " sigright=", round(sigs[2], 3), "p0=",
              round(fp0[4], 3)))
    if(nulltype == 1)
      title(xlab = paste("delhat=", round(xmaxx, 3),
              " sighat=", round(sighat, 3), "p0=", round(
                                                         fp0[3], 3)))
    lines(x, f, lwd = 3, col = 3)
    lines(x, f0, lwd = 2, lty = 2, col = 4)
    z2hi <- approx(fd[x > 0], x[x > 0], 0.2)$y
    z2lo <- approx(fd[x < 0], x[x < 0], 0.2)$y
    if(!is.na(z2hi))
      points(z2hi, 0, pch = 15)
    if(!is.na(z2lo))
      points(z2lo, 0, pch = 15)
    if(nulltype == 1)
      Ef <- Efdr[1]
    else Ef <- Efdr[4]

    if(plot == 2) {
      matplot(x, cbind(fdr, Fdrl, Fdrr), type = "l",
              lwd = 3, xlab = " ", ylim = c(0, 
                                     1.1000000000000001), main = 
              "fdr (solid); Fdr's (dashed)")
      title(xlab = paste("Efdr= ", round(Ef, 3)))
      abline(0, 0, lty = 3, col = 2)
      lines(c(0, 0), c(0, 1), lty = 3, col = 2)
    }
    if(plot == 3) {
      plot(cdf1[, 1], cdf1[, 2], type = "l", lwd = 3,
           xlab = "fdr level", ylim = c(0, 1), ylab
           = "f1 proportion < fdr level", main = 
           "f1 cdf of estimated fdr")
      title(sub = paste("Efdr= ", round(Ef, 3)))
      lines(c(0.20000000000000001, 0.20000000000000001),
            c(0, cdf1[20, 2]), col = 4, lty = 2)
      lines(c(0, 0.20000000000000001), rep(cdf1[20,
                                                2], 2), col = 4, lty = 2)
      text(0.050000000000000003, cdf1[20, 2], round(
                                                    cdf1[20, 2], 2))
      abline(0, 0, col = 2)
      lines(c(0, 0), c(0, 1), col = 2)
    }
  }
  if(nulltype == 0) {
    ffdr <- approx(x, fdr0, zz, rule = 2)$y
  }
  else ffdr <- approx(x, fdr, zz, rule = 2)$y
  list(fdr = ffdr, fp0 = fp0, Efdr = Efdr, cdf1 = cdf1, mat = mat)
}
