"locfdr" <-
function(zz,
                     bre = 120,
                     df = 7,
                     pct = 1/1000,
                     pct0 = 2/3,
                     nulltype = 1,
                     type = 0) {
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
  zh <- hist(zzz, breaks = breaks, plot = FALSE)
  x <- (breaks[-1] + breaks[ - length(breaks)])/2
  y <- zh$counts
  N <- length(y)
  if(pct > 0) {
    y[1.] <- min(y[1.], 1.)
    y[N] <- min(y[N], 1.)
  }
  if(type == 0) {
    f <- glm(y ~ ns(x, df = df), poisson)$fit
  } else {
    f <- glm(y ~ poly(x, df = df), poisson)$fit
  }
  l <- log(f)
  Fl <- cumsum(f)
  Fr <- cumsum(rev(f))
  ## ..............begin f0 calcs................................
  imax <- seq(l)[l == max(l)][1]
  xmax <- x[imax]
  lo0 <- quantile(zz[zz < xmax], 1 - pct0)
  hi0 <- quantile(zz[zz > xmax], pct0)
  nx <- length(x)
  i0 <- (1.:nx)[x > lo0 & x < hi0]
  x0 <- x[i0]
  y0 <- l[i0]
  if(nulltype == 2) {
    X0 <- cbind((x0 - xmax)^2, pmax(x0 - xmax, 0)^2)
  } else {
    X0 <- cbind(x0 - xmax, (x0 - xmax)^2)
  }
  lr <- lm(y0 ~ X0)
  co <- lr$coef
  if(nulltype == 2) {
    X00 <- cbind(1, (x - xmax)^2, pmax(x - xmax, 0)^2)
    sigs <- 1/sqrt(-2 * (c(co[2], co[2] + co[3])))
    f0.p0 <- c(xmax, sigs)
  } else {
    X00 <- cbind(1, x - xmax, (x - xmax)^2)
    sighat <- 1./sqrt(-2. * co[3.])
    xmaxx <-  - co[2.]/(2. * co[3.]) + xmax
    f0.p0 <- c(xmaxx, sighat)
  }
  l0 <- as.vector(X00 %*% co)
  f0 <- exp(l0)
  fdr <- pmin(f0/f, 1)
  F0l <- cumsum(f0)
  F0r <- cumsum(rev(f0))
  Fdrl <- F0l/Fl
  Fdrr <- rev(F0r/Fr)
  p0 <- sum(f0)/sum(f)
  f0 <- f0/p0
  f0.p0 <- c(f0.p0, p0)
  if(nulltype == 2) {
    names(f0.p0) <- c("zmax", "sigleft", "sigright", "p0")
  } else {
    names(f0.p0) <- c("zmax", "sig", "p0")
  }
  namat <- c("z.", "fdr", "fdrtheo", "Fdrleft", "Fdrright", "f",
             "f0", "f0theo", "counts")
  if(nulltype == 0) {
    f00 <- exp( - x^2/2)
    f00 <- (f00 * sum(f[i0]))/sum(f00[i0])
    fdr0 <- pmin(f00/f, 1)
    p00 <- sum(f00)/sum(f)
    f00 <- f00/p00
    mat <- cbind(x, fdr, fdr0, Fdrl, Fdrr, f, f0, f00, y)
    f0.p0 <- c(0, 1, p00)
    dimnames(mat) <- list(NULL, namat)
  } else {
    mat <- cbind(x, fdr, Fdrl, Fdrr, f, f0, y)
    dimnames(mat) <- list(NULL, namat[ - c(3, 8)])
  }
  if(nulltype == 0) {
    ffdr <- approx(x, fdr0, zz, rule = 2)$y
  } else {
    ffdr <- approx(x, fdr, zz, rule = 2)$y
  }
  list(fdr = ffdr, f0.p0 = f0.p0, mat = mat)
}
