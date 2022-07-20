#defined values for tidyness
sv.colnames<-c("UCLFilename","Channel","L(C)","L(C).Units",
  "Min(C)","Min(C).Units","Min(C).Time","Min(C).Time.Units",
  "Max(C)","Max(C).Units","Max(C).Time","Max(C).Time.Units",
  "L5(C)","L5(C).Units","L10(C)","L10(C).Units",
  "L50(C)","L50(C).Units","L90(C)","L90(C).Units",
  "L95(C)","L95(C).Units","L(A)","L(A).Units",
  "Min(A)","Min(A).Units","Min(A).Time","Min(A).Time.Units",
  "Max(A)","Max(A).Units","Max(A).Time","Max(A).Time.Units",
  "L5(A)","L5(A).Units","L10(A)","L10(A).Units",
  "L50(A)","L50(A).Units","L90(A)","L90(A).Units",
  "L95(A)","L95(A).Units","N5","N5.Units",
  "Min(N)","Min(N).Units","Min(N).Time","Min(N).Time.Units",
  "Max(N)","Max(N).Units","Max(N).Time","Max(N).Time.Units",
  "N10","N10.Units","N50","N50.Units",
  "N90","N90.Units","N95","N95.Units","S","S.Units",
  "Min(S)","Min(S).Units","Min(S).Time","Min(S).Time.Units",
  "Max(S)","Max(S).Units","Max(S).Time","Max(S).Time.Units",
  "S5","S5.Units","S10","S10.Units",
  "S50","S50.Units","S90","S90.Units",
  "S95","S95.Units","R","R.Units",
  "Min(R)","Min(R).Units","Min(R).Time","Min(R).Time.Units",
  "Max(R)","Max(R).Units","Max(R).Time","Max(R).Time.Units",
  "R5","R5.Units","R10","R10.Units",
  "R50","R50.Units","R90","R90.Units",
  "R95","R95.Units", "Tone","Tone.Units",
  "Min(Tone)","Min(Tone).Units","Min(Tone).Time","Min(Tone).Time.Units",
  "Max(Tone)","Max(Tone).Units","Max(Tone).Time","Max(Tone).Time.Units",
  "Tone5","Tone5.Units","Tone10","Tone10.Units",
  "Tone50","Tone50.Units","Tone90","Tone90.Units",
  "Tone95","Tone95.Units","Min(Tonef)","Min(Tonef).Units",
  "Min(Tonef).Time","Min(Tonef).Time.Units",
  "Max(Tonef)","Max(Tonef).Units","Max(Tonef).Time","Max(Tonef).Time.Units",
  "Tonef5","Tonef5.Units","Tonef10","Tonef10.Units",
  "Tonef50","Tonef50.Units","Tonef90","Tonef90.Units",
  "Tonef95","Tonef95.Units","Fluc","Fluc.Units",
  "Min(Fluc)","Min(Fluc).Units","Min(Fluc).Time","Min(Fluc).Time.Units",
  "Max(Fluc)","Max(Fluc).Units","Max(Fluc).Time","Max(Fluc).Time.Units",
  "Fluc5","Fluc5.Units","Fluc10","Fluc10.Units",
  "Fluc50","Fluc50.Units","Fluc90","Fluc90.Units",
  "Fluc95","Fluc95.Units")

params.udrTest<-
  c("L(C)","L5(C)","L10(C)","L50(C)","L90(C)","L95(C)",
"L(A)","L5(A)","L10(A)","L50(A)","L90(A)","L95(A)",
"N5","N10","N50","N90","N95",
"S","S5","S10","S50","S90","S95",
"R","R5","R10","R50","R90","R95",
"Tone5","Tone10","Tone50","Tone90","Tone95",
"Fluc5","Fluc10","Fluc50","Fluc90","Fluc95")

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

#function for BA plots

baplotOpts<-function(data,X,Y,color, #ggplot aes
                     fg.XY, #facet grid formula
                     np,#no. of params
                     #hline for mean diff, upper LoA, lower LoA
                     hl.mean, hl.upper, hl.lower, 
                     #lower/upper 95% limits of lwr/upr LoA & mean
                     lwr.ymin,lwr.ymax,upr.ymin,upr.ymax,m.ymin,m.ymax,
                     geom.text.size,theme.size,geom.point.size,colorp,
                     xlim.low,xlim.upp,ylim.low,ylim.upp,ylabel,xlabel){
  ggplot(data = data,
         aes(x={{X}},y={{Y}},color={{color}})) +
    facet_grid({{fg.XY}}) +
    geom_point(size=geom.point.size,alpha=0.5) + #add differences
    #average difference line (mean)
    geom_hline(aes(yintercept = {{hl.mean}}), color=set1clr[5]) +
    geom_text(aes(xlim.upp,{{hl.mean}},
                  label = paste("M =",as.character(round({{hl.mean}},2))), 
                  vjust = -0.5,hjust="right"), 
              color=set1clr[5],size=geom.text.size, 
              check_overlap = T) +
    #95% cf. int. mean
    geom_rect(data = data[1:(np*3),], #prevent drawing rect every repeated row
              aes(xmin = -Inf, xmax = Inf, 
                  ymin = {{m.ymin}}, ymax = {{m.ymax}}),
              fill = set1clr[5], alpha = 0.2,color = NA) +
    #lower bound line
    geom_hline(aes(yintercept = {{hl.lower}}), 
               color = set1clr[4], linetype="dashed") +
    geom_text(aes(xlim.upp,{{hl.lower}},
                  label = paste("M - 1.96SD:",as.character(round({{hl.lower}},2))), 
                  vjust = 2.5,hjust="right"), 
              color=set1clr[4],size=geom.text.size,
              check_overlap = T, parse = F) +
    geom_rect(data = data[1:(np*3),],aes(xmin = -Inf, xmax = Inf, 
                                         ymin = {{lwr.ymin}}, ymax = {{lwr.ymax}}),
              fill = set1clr[4], alpha = 0.2,color = NA) +
    #upper bound line
    geom_hline(aes(yintercept = {{hl.upper}}), 
               color = set1clr[4], linetype="dashed")  +
    geom_text(aes(xlim.upp,{{hl.upper}},
                  label = paste("M + 1.96SD:",as.character(round({{hl.upper}},2))), 
                  vjust = -1.5 ,hjust="right"), size=geom.text.size,
              color=set1clr[4],check_overlap = T, parse = F) +
    geom_rect(data = data[1:(np*3),],aes(xmin = -Inf, xmax = Inf, 
                                         ymin = {{upr.ymin}}, ymax = {{upr.ymax}}),
              fill = set1clr[4], alpha = 0.2,color = NA) +
    xlim(xlim.low, xlim.upp) + ylim(ylim.low, ylim.upp) +
    scale_color_brewer(palette = "Set1") +
    ylab(label = ylabel) + xlab(label = xlabel) +
    theme(text = element_text(size=theme.size),legend.position="none",
          axis.text.x = element_text(angle = 30, 
                                     vjust = 0.75, hjust=0.5))
}

#rounding function
round_any = function(x, accuracy, f=round){f(x/ accuracy) * accuracy}