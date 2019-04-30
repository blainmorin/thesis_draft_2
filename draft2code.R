library(metafor)

########## Copas functions for Benefit Outcomes
b.lik = function(theta, tau2) {
  
  QH = sum(log(pnorm((1.96*sH - theta)/sqrt(sH^2 + tau2)) - pnorm((-1.96*sH - theta)/sqrt(sH^2 + tau2))))
  
  -1 * sum((y - theta)^2/(s^2 + tau2) + log(s^2 + tau2))/2 + QH
  
}


b.CI = function(){
  
  lik1 = function(par, theta) { 
    -1 * b.lik(theta, par[1]) 
  }
  
  lik2 = function(theta, lik1){ 
    
    parhat = c(0)
    
    nlminb(parhat, lik1, lower = c(0), upper = c(Inf), theta = theta)$objective 
    
  }
  
  temp1 = nlminb(sum(y/s^2)/sum(1/s^2), lik2, lik1 = lik1)
  
  thetahat = temp1$par
  
  likmax = -1 * temp1$objective
  
  lik3 = function(theta, likmax, lik1, lik2){
    (likmax + lik2(theta, lik1) - 1.96)^2 
  }
  
  thetahatL = nlminb(thetahat - 0.01, lik3, upper = thetahat, likmax = likmax, lik1 = lik1, lik2 = lik2)$par
  
  thetahatU = nlminb(thetahat + 0.01, lik3, lower = thetahat, likmax = likmax, lik1 = lik1, lik2 = lik2)$par
  
  c(thetahatL, thetahat, thetahatU)
  
}

CIRE = function(){
  
  likRE = function(tau2, theta){ 
    sum((y - theta)^2/(s^2 + tau2) + log(s^2 + tau2))/2 
  }
  
  likRE1 = function(theta, likRE){
    nlminb(0, likRE, theta = theta, lower = 0)$objective 
  }
  
  temp = nlminb(sum(y/s^2)/sum(1/s^2), likRE1, likRE = likRE)
  
  thetaREhat = temp$par
  
  likREmax = -1 * temp$objective
  
  likRE2 = function(theta, likREmax, likRE, likRE1){ 
    (likREmax + likRE1(theta, likRE) - 1.96)^2 
  }
  
  thetaREhatL = nlminb(thetaREhat - 0.01, likRE2, upper = thetaREhat, likREmax = likREmax, likRE = likRE, likRE1 = likRE1)$par
  
  thetaREhatU = nlminb(thetaREhat + 0.01, likRE2, lower = thetaREhat, likREmax = likREmax, likRE = likRE, likRE1 = likRE1)$par
  
  c(thetaREhatL, thetaREhat, thetaREhatU)
  
}

###### Copas functions for harm outcomes

h.lik = function(theta, tau2) {
  
  QH = sum(log(pnorm(theta/sqrt(sH^2 + tau2))))
  
  -1 * sum((y - theta)^2/(s^2 + tau2) + log(s^2 + tau2))/2 + QH
  
}


h.CI = function(){
  
  lik1 = function(par, theta) { 
    -1 * h.lik(theta, par[1]) 
  }
  
  lik2 = function(theta, lik1){ 
    
    parhat = c(0)
    
    nlminb(parhat, lik1, lower = c(0), theta = theta, upper = c(Inf))$objective
    
  }
  
  temp1 = nlminb(sum(y/s^2)/sum(1/s^2), lik2, lik1 = lik1)
  
  thetahat = temp1$par
  
  likmax = -1 * temp1$objective
  
  lik3 = function(theta, likmax, lik1, lik2){
    (likmax + lik2(theta, lik1) - 1.96)^2 
  }
  
  thetahatL = nlminb(thetahat - 0.01, lik3, upper = thetahat, likmax = likmax, lik1 = lik1, lik2 = lik2)$par
  
  thetahatU = nlminb(thetahat + 0.01, lik3, lower = thetahat, likmax = likmax, lik1 = lik1, lik2 = lik2)$par
  
  c(thetahatL, thetahat, thetahatU)
  
}

#############################################################################
### Review Benefits THR
#############################################################################

### Make author year pmid column
cites  = c(
  "Eriksson2010, 20088935",
  "Fuji2014A, 25047458",
  "Fuji2014D, 22952213",
  "Fuji2015, 26269694",
  "Lassen2010A, 21175312",
  "Raskob2010, 20589317",
  "Eriksson2014, 24136153",
  "Yokote2011, 21282767",
  "Zhang2013, EMBASE 2014592535",
  "Turpie2002, 12049860",
  "Lassen2002, 12049858",
  "Kim2016, 26790579",
  "Kakkar2008, 18582928",
  "Turpie2001, 11228275",
  "Eriksson2005, 16409461",
  "Eriksson2006A,	17292948",
  "Eriksson2008, 18579811",
  "Eriksson2006B,	17116766"
)

##############
### Total VTE
##############

### LMWH
as = c(24, 3, 2, 17, 72, 63, 48, 5, 66, 85, 24, 76, 16, 18, 18, 58, 27)
n1 = c(127, 74, 82, 248, 1917, 144, 314, 83, 797, 919, 351, 869, 171, 106, 107, 1558, 107)

### FXaI
cs = c(123, 5, 11, 6, 23, 120, 149, 6, 48, 37, 26, 15, 26, 63, 64, 17, 55)
n2 = c(581, 150, 270, 255, 1949, 630, 1132, 84, 787, 908, 350, 864, 422, 442, 359, 1595, 511)

### ni for each High risk of bias
nhigh = c(106)


### Calculate ORs and variances
ors = escalc(measure = "OR", ai = as, n1i = n1, ci = cs, n2i = n2)

y = ors$yi
s = sqrt(ors$vi)
nreported = n1 + n2

### Estimate unreported variances
kbar = sum(s^-2)/sum(nreported)

sH = sqrt(1/(kbar*nhigh))

unadjusted = CIRE()
adjusted = b.CI()
n.high = length(nhigh)

##### Make forest plot
cites  = c(
  "Eriksson2010, 20088935",
  "Fuji2014A, 25047458",
  "Fuji2014D, 22952213",
  "Fuji2015, 26269694",
  "Lassen2010A, 21175312",
  "Raskob2010, 20589317",
  "Eriksson2014, 24136153",
  "Yokote2011, 21282767",
  "Turpie2002, 12049860",
  "Lassen2002, 12049858",
  "Kim2016, 26790579",
  "Kakkar2008, 18582928",
  "Turpie2001, 11228275",
  "Eriksson2005, 16409461",
  "Eriksson2006A,	17292948",
  "Eriksson2008, 18579811",
  "Eriksson2006B,	17116766"
)

temp = data.frame(as = as, n1 = n1, cs = cs, n2 = n2, row.names = NULL)
test = rma.uni(yi = ors$yi, vi = ors$vi)
forest(test, slab = cites,
       ilab = temp,
       xlim = c(-14, 8),
       ylim = c(-3, 22),
       rows = 18,
       ilab.xpos = c(-7,-6,-5,-4),
       cex = .8,
       addfit = FALSE,
       atransf = exp)

addpoly.default(x = unadjusted[2], ci.lb = unadjusted[1], ci.ub = unadjusted[3], 
                rows = -1, 
                atransf = exp,
                cex = .8,
                mlab = "Maximum Likelihood RE Estimate",
                col = "Black")

addpoly.default(x = adjusted[2], ci.lb = adjusted[1], ci.ub = adjusted[3], 
                rows = -2, 
                atransf = exp,
                cex = .8,
                mlab = "ORB Adjusted Estimate",
                col = "Red")

abline(h = 0)
text(c(-6.5, -4.5), 22, c("LMWH", "FXaI"))
text(c(-7,-6,-5,-4), 21, c("n", "N", "n", "N"), cex = .8)
text(c(-12), 21, c("Author Year PMID"), cex = .8)
text(c(-1.2, 1.2), 21, c("Favors LMWH", "Favors FXaI"), cex = .8)
text(c(6), 21, c("Odds Ratio [95% CI]"), cex = .8)
title("Updated Forest Plot: THR, Total VTE, LMWH vs FXaI")

####################################################################

########################
### Symptomatic VTE
####################

### LMWH
as = c(0, 0, 10, 0, 1, 3, 8, 15, 0, 11, 1)
n1 = c(127, 82, 2699, 83, 1128, 1123, 251, 1207, 107, 2206, 107)

### FXaI
cs = c(4, 1, 4, 1, 10, 5, 8, 3, 4, 6, 0)
n2 = c(581, 270, 2708, 84, 1126, 1129, 350, 1212, 359, 2193, 511)

### ni for each High risk of bias
nhigh = c(sum(144, 630), sum(314, 1132), 106, sum(171, 422))

### Calculate ORs and variances
ors = escalc(measure = "OR", ai = as, n1i = n1, ci = cs, n2i = n2)

y = ors$yi
s = sqrt(ors$vi)
nreported = n1 + n2

### Estimate unreported variances
kbar = sum(s^-2)/sum(nreported)

sH = sqrt(1/(kbar*nhigh))

unadjusted = CIRE()
adjusted = b.CI()
n.high = length(nhigh)

##### Make forest plot
cites  = c(
  "Eriksson2010, 20088935",
  "Fuji2014D, 22952213",
  "Lassen2010A, 21175312",
  "Yokote2011, 21282767",
  "Turpie2002, 12049860",
  "Lassen2002, 12049858",
  "Kim2016, 26790579",
  "Kakkar2008, 18582928",
  "Eriksson2006A,	17292948",
  "Eriksson2008, 18579811",
  "Eriksson2006B,	17116766"
)

temp = data.frame(as = as, n1 = n1, cs = cs, n2 = n2, row.names = NULL)
test = rma.uni(yi = ors$yi, vi = ors$vi)
forest(test, slab = cites,
       ilab = temp,
       xlim = c(-14, 8),
       ylim = c(-3, 16),
       rows = 12,
       ilab.xpos = c(-7,-6,-5,-4),
       clim = c(-2, 3),
       cex = .8,
       addfit = FALSE,
       atransf = exp)

addpoly.default(x = unadjusted[2], ci.lb = unadjusted[1], ci.ub = unadjusted[3], 
                rows = -1, 
                atransf = exp,
                cex = .8,
                mlab = "Maximum Likelihood RE Estimate",
                col = "Black")

addpoly.default(x = adjusted[2], ci.lb = adjusted[1], ci.ub = adjusted[3], 
                rows = -2, 
                atransf = exp,
                cex = .8,
                mlab = "ORB Adjusted Estimate",
                col = "Red")
abline(h = 0)
text(c(-6.5, -4.5), 16, c("LMWH", "FXaI"))
text(c(-7,-6,-5,-4), 15, c("n", "N", "n", "N"), cex = .8)
text(c(-12), 15, c("Author Year PMID"), cex = .8)
text(c(-1.2, 1.2), 15, c("Favors LMWH", "Favors FXaI"), cex = .8)
text(c(6), 15, c("Odds Ratio [95% CI]"), cex = .8)
title("Updated Forest Plot: THR, Symptomatic VTE, LMWH vs FXaI")

####################################################################

########################
### Total PE
####################

### LMWH
as = c(5, 1, 2, 1, 5, 0, 0, 1)
n1 = c(2699, 1128, 1123, 251, 869, 171, 107, 1558)

### FXaI
cs = c(3, 5, 2, 2, 1, 2, 3, 4)
n2 = c(2708, 1126, 1129, 350, 864, 422, 359, 1595)

### ni for each High risk of bias
nhigh = c(sum(144, 630), sum(314,1132), 106)

### Calculate ORs and variances
ors = escalc(measure = "OR", ai = as, n1i = n1, ci = cs, n2i = n2)

y = ors$yi
s = sqrt(ors$vi)
nreported = n1 + n2

### Estimate unreported variances
kbar = sum(s^-2)/sum(nreported)

sH = sqrt(1/(kbar*nhigh))

unadjusted = CIRE()
adjusted = b.CI()
n.high = length(nhigh)

##### Make forest plot
cites  = c(
  "Lassen2010A, 21175312",
  "Turpie2002, 12049860",
  "Lassen2002, 12049858",
  "Kim2016, 26790579",
  "Kakkar2008, 18582928",
  "Turpie2001, 11228275",
  "Eriksson2006A,	17292948",
  "Eriksson2008, 18579811"
)

temp = data.frame(as = as, n1 = n1, cs = cs, n2 = n2, row.names = NULL)
test = rma.uni(yi = ors$yi, vi = ors$vi)
forest(test, slab = cites,
       ilab = temp,
       xlim = c(-14, 8),
       ylim = c(-3, 11),
       rows = 8,
       ilab.xpos = c(-7,-6,-5,-4),
       clim = c(-2, 3),
       addfit = FALSE,
       cex = .8,
       atransf = exp)

addpoly.default(x = unadjusted[2], ci.lb = unadjusted[1], ci.ub = unadjusted[3], 
                rows = -1, 
                atransf = exp,
                cex = .8,
                mlab = "Maximum Likelihood RE Estimate",
                col = "Black")

addpoly.default(x = adjusted[2], ci.lb = adjusted[1], ci.ub = adjusted[3], 
                rows = -2, 
                atransf = exp,
                cex = .8,
                mlab = "ORB Adjusted Estimate",
                col = "Red")
abline(h = 0)
text(c(-6.5, -4.5), 11, c("LMWH", "FXaI"))
text(c(-7,-6,-5,-4), 10, c("n", "N", "n", "N"), cex = .8)
text(c(-12), 10, c("Author Year PMID"), cex = .8)
text(c(-1.2, 1.2), 10, c("Favors LMWH", "Favors FXaI"), cex = .8)
text(c(6), 10, c("Odds Ratio [95% CI]"), cex = .8)
title("Updated Forest Plot: THR, Total PE, LMWH vs FXaI")

####################################################################


########################
### Fatal PE
####################

### LMWH
as = c(0, 1, 1, 0, 0)
n1 = c(2699, 1128, 869, 171, 107)

### FXaI
cs = c(1, 0, 0, 1, 1)
n2 = c(2708, 1126, 864, 422, 359)

### ni for each High risk of bias
nhigh = c(106)

### Calculate ORs and variances
ors = escalc(measure = "OR", ai = as, n1i = n1, ci = cs, n2i = n2)

y = ors$yi
s = sqrt(ors$vi)
nreported = n1 + n2

### Estimate unreported variances
kbar = sum(s^-2)/sum(nreported)

sH = sqrt(1/(kbar*nhigh))

unadjusted = CIRE()
adjusted = b.CI()
n.high = length(nhigh)

##### Make forest plot

cites  = c(
  "Lassen2010A, 21175312",
  "Turpie2002, 12049860",
  "Kakkar2008, 18582928",
  "Turpie2001, 11228275",
  "Eriksson2006A,	17292948"
)


temp = data.frame(as = as, n1 = n1, cs = cs, n2 = n2, row.names = NULL)
test = rma.uni(yi = ors$yi, vi = ors$vi)
forest(test, slab = cites,
       ilab = temp,
       xlim = c(-14, 8),
       ylim = c(-3, 8),
       rows = 5,
       ilab.xpos = c(-7,-6,-5,-4),
       clim = c(-2, 3),
       cex = .8,
       addfit = FALSE,
       atransf = exp)

addpoly.default(x = unadjusted[2], ci.lb = unadjusted[1], ci.ub = unadjusted[3], 
                rows = -1, 
                atransf = exp,
                cex = .8,
                mlab = "Maximum Likelihood RE Estimate",
                col = "Black")

addpoly.default(x = adjusted[2], ci.lb = adjusted[1], ci.ub = adjusted[3], 
                rows = -2, 
                atransf = exp,
                cex = .8,
                mlab = "ORB Adjusted Estimate",
                col = "Red")

abline(h = 0)
text(c(-6.5, -4.5), 8, c("LMWH", "FXaI"))
text(c(-7,-6,-5,-4), 7, c("n", "N", "n", "N"), cex = .8)
text(c(-12), 7, c("Author Year PMID"), cex = .8)
text(c(-1.2, 1.2), 7, c("Favors LMWH", "Favors FXaI"), cex = .8)
text(c(6), 7, c("Odds Ratio [95% CI]"), cex = .8)
title("Updated Forest Plot: THR, Fatal PE, LMWH vs FXaI")

####################################################################


########################
### Symptomatic PE
####################

### LMWH
as = c(5, 1, 2, 1, 5, 0, 0, 1)
n1 = c(2699, 1128, 1123, 251, 869, 171, 107, 1558)

### FXaI
cs = c(3, 5, 2, 2, 1, 2, 3, 4)
n2 = c(2708, 1126, 1129, 350, 864, 422, 359, 1595)

### ni for each High risk of bias
nhigh = c(sum(144, 630), sum(314,1132), 106)

### ni for each High risk of bias
nhigh = c(sum(144, 630), sum(314,1132), 106)

### Calculate ORs and variances
ors = escalc(measure = "OR", ai = as, n1i = n1, ci = cs, n2i = n2)

y = ors$yi
s = sqrt(ors$vi)
nreported = n1 + n2

### Estimate unreported variances
kbar = sum(s^-2)/sum(nreported)

sH = sqrt(1/(kbar*nhigh))

unadjusted = CIRE()
adjusted = b.CI()
n.high = length(nhigh)

##### Make forest plot

cites  = c(
  "Lassen2010A, 21175312",
  "Turpie2002, 12049860",
  "Lassen2002, 12049858",
  "Kim2016, 26790579",
  "Kakkar2008, 18582928",
  "Turpie2001, 11228275",
  "Eriksson2006A,	17292948",
  "Eriksson2008, 18579811"
)

temp = data.frame(as = as, n1 = n1, cs = cs, n2 = n2, row.names = NULL)
test = rma.uni(yi = ors$yi, vi = ors$vi)
forest(test, slab = cites,
       ilab = temp,
       xlim = c(-14, 8),
       ylim = c(-3, 11),
       rows = 8,
       ilab.xpos = c(-7,-6,-5,-4),
       clim = c(-2, 3),
       cex = .8,
       addfit = FALSE,
       atransf = exp)

addpoly.default(x = unadjusted[2], ci.lb = unadjusted[1], ci.ub = unadjusted[3], 
                rows = -1, 
                atransf = exp,
                cex = .8,
                mlab = "Maximum Likelihood RE Estimate",
                col = "Black")

addpoly.default(x = adjusted[2], ci.lb = adjusted[1], ci.ub = adjusted[3], 
                rows = -2, 
                atransf = exp,
                cex = .8,
                mlab = "ORB Adjusted Estimate",
                col = "Red")
abline(h = 0)
text(c(-6.5, -4.5), 11, c("LMWH", "FXaI"))
text(c(-7,-6,-5,-4), 10, c("n", "N", "n", "N"), cex = .8)
text(c(-12), 10, c("Author Year PMID"), cex = .8)
text(c(-1.2, 1.2), 10, c("Favors LMWH", "Favors FXaI"), cex = .8)
text(c(6), 10, c("Odds Ratio [95% CI]"), cex = .8)
title("Updated Forest Plot: THR, Symptomatic PE, LMWH vs FXaI")

####################################################################


##############
### Total DVT
##############

### LMWH
as = c(24, 3, 2, 17, 69, 5, 7, 65, 83, 23, 71, 23, 18, 18, 53, 27)
n1 = c(127, 74, 82, 248, 1911, 83, 53, 796, 918, 351, 869, 351, 106, 107, 1558, 107)

### FXaI
cs = c(123, 5, 10, 6, 22, 6, 0, 44, 36, 24, 14, 24, 63, 62, 12, 55)
n2 = c(581, 150, 270, 255, 1944, 84, 53, 784, 908, 350, 864, 422, 442, 359, 1595, 511)

### ni for each High risk of bias
nhigh = c(sum(144, 630), sum(314, 1132))


### Calculate ORs and variances
ors = escalc(measure = "OR", ai = as, n1i = n1, ci = cs, n2i = n2)

y = ors$yi
s = sqrt(ors$vi)
nreported = n1 + n2

### Estimate unreported variances
kbar = sum(s^-2)/sum(nreported)

sH = sqrt(1/(kbar*nhigh))

unadjusted = CIRE()
adjusted = b.CI()
n.high = length(nhigh)

##### Make forest plot
cites  = c(
  "Eriksson2010, 20088935",
  "Fuji2014A, 25047458",
  "Fuji2014D, 22952213",
  "Fuji2015, 26269694",
  "Lassen2010A, 21175312",
  "Yokote2011, 21282767",
  "Zhang2013, EMBASE 2014592535",
  "Turpie2002, 12049860",
  "Lassen2002, 12049858",
  "Kim2016, 26790579",
  "Kakkar2008, 18582928",
  "Turpie2001, 11228275",
  "Eriksson2005, 16409461",
  "Eriksson2006A,	17292948",
  "Eriksson2008, 18579811",
  "Eriksson2006B,	17116766"
)

temp = data.frame(as = as, n1 = n1, cs = cs, n2 = n2, row.names = NULL)
test = rma.uni(yi = ors$yi, vi = ors$vi)
forest(test, slab = cites,
       ilab = temp,
       xlim = c(-14, 8),
       ylim = c(-3, 20),
       rows = 16,
       ilab.xpos = c(-7,-6,-5,-4),
       cex = .8,
       addfit = FALSE,
       clim = c(-2, 4),
       atransf = exp)

addpoly.default(x = unadjusted[2], ci.lb = unadjusted[1], ci.ub = unadjusted[3], 
                rows = -1, 
                atransf = exp,
                cex = .8,
                mlab = "Maximum Likelihood RE Estimate",
                col = "Black")

addpoly.default(x = adjusted[2], ci.lb = adjusted[1], ci.ub = adjusted[3], 
                rows = -2, 
                atransf = exp,
                cex = .8,
                mlab = "ORB Adjusted Estimate",
                col = "Red")
abline(h = 0)
text(c(-6.5, -4.5), 20, c("LMWH", "FXaI"))
text(c(-7,-6,-5,-4), 19, c("n", "N", "n", "N"), cex = .8)
text(c(-12), 19, c("Author Year PMID"), cex = .8)
text(c(-1.2, 1.2), 19, c("Favors LMWH", "Favors FXaI"), cex = .8)
text(c(6), 19, c("Odds Ratio [95% CI]"), cex = .8)
title("Updated Forest Plot: THR, Total DVT, LMWH vs FXaI")

####################################################################


###################
### Symptomatic DVT
##############

### LMWH
as = c(0, 5, 0, 0, 1, 7, 10, 0, 10, 1)
n1 = c(127, 2699, 83, 1128, 1123, 351, 1207, 107, 2206, 107)

### FXaI
cs = c(4, 1, 1, 5, 3, 6, 2, 1, 2, 0)
n2 = c(581, 2708, 84, 1126, 1129, 350, 1212, 359, 2193, 511)

### ni for each High risk of bias
nhigh = c(sum(144, 630), sum(314, 1132), 106, sum(171, 422))


### Calculate ORs and variances
ors = escalc(measure = "OR", ai = as, n1i = n1, ci = cs, n2i = n2)

y = ors$yi
s = sqrt(ors$vi)
nreported = n1 + n2

### Estimate unreported variances
kbar = sum(s^-2)/sum(nreported)

sH = sqrt(1/(kbar*nhigh))

unadjusted = CIRE()
adjusted = b.CI()
n.high = length(nhigh)

##### Make forest plot
cites  = c(
  "Eriksson2010, 20088935",
  "Lassen2010A, 21175312",
  "Yokote2011, 21282767",
  "Turpie2002, 12049860",
  "Lassen2002, 12049858",
  "Kim2016, 26790579",
  "Kakkar2008, 18582928",
  "Eriksson2006A,	17292948",
  "Eriksson2008, 18579811",
  "Eriksson2006B,	17116766"
)

temp = data.frame(as = as, n1 = n1, cs = cs, n2 = n2, row.names = NULL)
test = rma.uni(yi = ors$yi, vi = ors$vi)
forest(test, slab = cites,
       ilab = temp,
       xlim = c(-14, 8),
       ylim = c(-3, 14),
       rows = 10,
       ilab.xpos = c(-7,-6,-5,-4),
       cex = .8,
       addfit = FALSE,
       clim = c(-2, 4),
       atransf = exp)

addpoly.default(x = unadjusted[2], ci.lb = unadjusted[1], ci.ub = unadjusted[3], 
                rows = -1, 
                atransf = exp,
                cex = .8,
                mlab = "Maximum Likelihood RE Estimate",
                col = "Black")

addpoly.default(x = adjusted[2], ci.lb = adjusted[1], ci.ub = adjusted[3], 
                rows = -2, 
                atransf = exp,
                cex = .8,
                mlab = "ORB Adjusted Estimate",
                col = "Red")
abline(h = 0)
text(c(-6.5, -4.5), 14, c("LMWH", "FXaI"))
text(c(-7,-6,-5,-4), 13, c("n", "N", "n", "N"), cex = .8)
text(c(-12), 13, c("Author Year PMID"), cex = .8)
text(c(-1.2, 1.2), 13, c("Favors LMWH", "Favors FXaI"), cex = .8)
text(c(6), 13, c("Odds Ratio [95% CI]"), cex = .8)
title("Updated Forest Plot: THR, Symptomatic DVT, LMWH vs FXaI")

####################################################################

