#adopted from 

t.test.partial <- function(paired, unpaired.x, unpaired.y,
                           alternative=c("two.sided", "less", "greater"), conf.level=0.95) {
  
  # Also known as "Optimal pooled t-test"
  # Beibei Guo and Ying Yuan (2015/2017)
  # DOI: 10.1177/0962280215577111
  
  if (any(is.na(paired)) & NCOL(paired) > 1) {
    compc <- complete.cases(paired)
    unpaired.x <- paired[!compc, 1]
    unpaired.x <- unpaired.x[!is.na(unpaired.x)]
    unpaired.y <- paired[!compc, 2]
    unpaired.y <- unpaired.y[!is.na(unpaired.y)]
    paired <- paired[compc, 1:2]
  } else {
    if (missing(unpaired.x) & missing(unpaired.y)) {
      stop("missing unpaired samples")
    }
  }
  
  # paired samples
  if (NCOL(paired) > 1) {
    xpd <- paired[, 1] - paired[, 2]
  } else xpd <- paired
  
  np <- length(xpd)
  
  # unpaired samples
  if (NCOL(unpaired.x) > 1) {
    unpaired.y <- unpaired.x[, 2]
    unpaired.x <- unpaired.x[, 1]
  }    
  
  if (length(unpaired.x) < 2 | length(unpaired.y) < 2) {
    stop("Too few unpaired samples")
  }
  
  xu1 <- unpaired.x
  xu2 <- unpaired.y
  
  xu1 <- xu1[!is.na(xu1)]
  xu2 <- xu2[!is.na(xu2)]
  
  nu1 <- length(xu1)
  nu2 <- length(xu2)
  
  vu1 <- var(xu1)
  vu2 <- var(xu2)
  
  # estimate variances
  vp <- var(xpd)/np
  vu <- vu1/nu1 + vu2/nu2
  
  # weights
  wt <- (1/vp) + (1/vu)
  wp <- 1 / (wt * vp)
  wu <- 1 / (wt * vu)
  
  # estimates
  ep <- mean(xpd)
  eu <- mean(xu1) - mean(xu2)
  et <- (ep*wp) + (eu*wu)
  
  # degrees of freedom
  dfp <- np - 1
  dfu <- vu^2 / ( ((vu1/nu1)^2 / (nu1 - 1)) + 
                    ((vu2/nu2)^2 / (nu2 - 1)) )
  dfd <- 1 / ((wp^2)/dfp + (wu^2)/dfu) # df of t-dist
  
  # total estimate variance (standard error)
  vt <- (1 + (4*wp*(1 - wp)/dfp) + 
           (4*wu*(1 - wu)/dfu))/wt
  se <- sqrt(vt)
  
  # t statistic
  Ts <- et/se
  
  # confidence interval
  alpha <- 1 - conf.level
  alt <- match.arg(alternative)
  
  cint <- switch(alt,
                 "less" = c(-Inf, Ts + qt(conf.level, dfd)),
                 "greater" = c(Ts - qt(conf.level, dfd), Inf),
                 "two.sided" = {
                   ci <- qt(1 - alpha/2, dfd)
                   Ts + c(-ci, ci)
                 }
  )
  cint <- cint * se
  
  # p-value    
  pval <- switch(alt,
                 "less" = pt(Ts, dfd),
                 "greater" = pt(Ts, dfd, lower.tail=FALSE),
                 "two.sided" = 2 * pt(-abs(Ts), dfd)
  )
  
  names(Ts) <- "t"
  names(dfd) <- "df"
  names(et) <- "mean of the differences"
  attr(cint, "conf.level") <- conf.level
  rval <- list(statistic=Ts, parameter=dfd, p.value=pval, 
               conf.int=cint, estimate=et, alternative=alt,
               method="Partially paired t-test")
  class(rval) <- "htest"
  rval
}

list(paired = structure(list(test = c(52.87, 52.37, 47.39,
                                      52.56, 54.43, 50.79), 
                             cont = c(52.51, 57.87, 51.41,
                                      57.37, 56.84, 52.47)), 
                        class = "data.frame", 
                        row.names = c(NA, -6L)), 
     unpaired = structure(list(test = c(52.32, 53.26,53.16, 
                                        50.23, 52.34, 
                                        NA, NA, NA), 
                               cont = c(47.21, 46.38, 54.34,                                                                                                                     60.15, 55.08, 55.4, 58.05, 55.88)), class = "data.frame",
                              row.names = c(NA, -8L))) -> m2

set.seed(1)

n <- 40
p <- 0.3
r <- 0.5

m <- MASS::mvrnorm(n, c(0, 0), matrix(c(1, r, r, 1), 2))
m <- t(t(m) * c(1, 2))
m <- t(t(m) + c(51, 50))
m <- round(m, 1)
m[sample(n*2, p*n*2)] <- NA
colnames(m) <- c("test", "cont")
