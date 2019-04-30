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

cites = c(
  'Fuji2014d	22952213',
  'Zou2014	24695091',
  'Bonneux2006	16387501',
  'Hu2015	na',
  'Lassen2010b	20206776',
  'Bauer2001	11794149',
  'Fuji2014c	25294589',
  'Weitz2010	20886185',
  'Cohen2013	23782955',
  'Lassen2009	19657123',
  'Turpie2009	19411100',
  'Lassen2008	18579812',
  'Turpie2005	16241946'
)

##############
### Total VTE
##############

### LMWH
as = c(15, 14, 146, 101, 41, 24, 34, 99, 94, 164, 31)
n1 = c(66, 102, 976, 363, 295, 109, 188, 1596, 1508, 1217, 70)

### FXaI
cs = c(33, 3, 243, 45, 22, 160, 119, 105, 66, 79, 92)
n2 = c(152, 112, 997, 361, 299, 621, 561, 1599, 1526, 1201, 296)

### ni for each High risk of bias
nhigh = c()


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
cites = c(
  'Fuji2014d	22952213',
  'Zou2014	24695091',
  'Lassen2010b	20206776',
  'Bauer2001	11794149',
  'Fuji2014c	25294589',
  'Weitz2010	20886185',
  'Cohen2013	23782955',
  'Lassen2009	19657123',
  'Turpie2009	19411100',
  'Lassen2008	18579812',
  'Turpie2005	16241946'
)

temp = data.frame(as = as, n1 = n1, cs = cs, n2 = n2, row.names = NULL)
test = rma.uni(yi = ors$yi, vi = ors$vi)
forest(test, slab = cites,
       ilab = temp,
       xlim = c(-14, 8),
       ylim = c(-3, 14),
       rows = 11,
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
text(c(-6.5, -4.5), 14, c("LMWH", "FXaI"))
text(c(-7,-6,-5,-4), 13, c("n", "N", "n", "N"), cex = .8)
text(c(-12), 13, c("Author Year PMID"), cex = .8)
text(c(-1.2, 1.2), 13, c("Favors LMWH", "Favors FXaI"), cex = .8)
text(c(6), 13, c("Odds Ratio [95% CI]"), cex = .8)
title("Updated Forest Plot: TKR, Total VTE, LMWH vs FXaI")

####################################################################

########################
### Symptomatic VTE
##############

### LMWH
as = c(1, 1, 7, 7, 1, 4, 4, 13, 18, 24, 2)
n1 = c(90, 102, 1529, 517, 295, 109, 188, 1596, 1508, 1217, 70)

### FXaI
cs = c(1, 0, 7, 3, 4, 14, 3, 19, 11, 8, 4)
n2 = c(162, 112, 1528, 517, 299, 621, 561, 1599, 1526, 1201, 296)

### ni for each High risk of bias
nhigh = c()


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
cites = c(
  'Fuji2014d	22952213',
  'Zou2014	24695091',
  'Lassen2010b	20206776',
  'Bauer2001	11794149',
  'Fuji2014c	25294589',
  'Weitz2010	20886185',
  'Cohen2013	23782955',
  'Lassen2009	19657123',
  'Turpie2009	19411100',
  'Lassen2008	18579812',
  'Turpie2005	16241946'
)

temp = data.frame(as = as, n1 = n1, cs = cs, n2 = n2, row.names = NULL)
test = rma.uni(yi = ors$yi, vi = ors$vi)
forest(test, slab = cites,
       ilab = temp,
       xlim = c(-14, 8),
       ylim = c(-3, 14),
       rows = 11,
       ilab.xpos = c(-7,-6,-5,-4),
       cex = .8,
       clim = c(-2, 3),
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
text(c(-6.5, -4.5), 14, c("LMWH", "FXaI"))
text(c(-7,-6,-5,-4), 13, c("n", "N", "n", "N"), cex = .8)
text(c(-12), 13, c("Author Year PMID"), cex = .8)
text(c(-1.2, 1.2), 13, c("Favors LMWH", "Favors FXaI"), cex = .8)
text(c(6), 13, c("Odds Ratio [95% CI]"), cex = .8)
title("Updated Forest Plot: TKR, Symptomatic VTE, LMWH vs FXaI")

####################################################################


########################
### Total PE
##############

### LMWH
as = c(1, 0, 4, 2, 3, 7, 8, 4, 0)
n1 = c(90, 1529, 517, 109, 188, 1596, 1508, 1217, 70)

### FXaI
cs = c(0, 4, 1, 5, 2, 16, 5, 0, 2)
n2 = c(162, 1528, 517, 621, 561, 1599, 1526, 1201, 296)

### ni for each High risk of bias
nhigh = c()


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
cites = c(
  'Fuji2014d	22952213',
  'Lassen2010b	20206776',
  'Bauer2001	11794149',
  'Weitz2010	20886185',
  'Cohen2013	23782955',
  'Lassen2009	19657123',
  'Turpie2009	19411100',
  'Lassen2008	18579812',
  'Turpie2005	16241946'
)

temp = data.frame(as = as, n1 = n1, cs = cs, n2 = n2, row.names = NULL)
test = rma.uni(yi = ors$yi, vi = ors$vi)
forest(test, slab = cites,
       ilab = temp,
       xlim = c(-14, 8),
       ylim = c(-3, 12),
       rows = 9,
       ilab.xpos = c(-7,-6,-5,-4),
       cex = .8,
       clim = c(-2, 3),
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
text(c(-6.5, -4.5), 12, c("LMWH", "FXaI"))
text(c(-7,-6,-5,-4), 11, c("n", "N", "n", "N"), cex = .8)
text(c(-12), 11, c("Author Year PMID"), cex = .8)
text(c(-1.2, 1.2), 11, c("Favors LMWH", "Favors FXaI"), cex = .8)
text(c(6), 11, c("Odds Ratio [95% CI]"), cex = .8)
title("Updated Forest Plot: TKR, Total PE, LMWH vs FXaI")

####################################################################

########################
### Fatal PE
##############

### LMWH
as = c(0, 0, 0)
n1 = c(1529, 1596, 1508)

### FXaI
cs = c(1, 2, 1)
n2 = c(1528, 1599, 1526)

### ni for each High risk of bias
nhigh = c()


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
cites = c(
  'Lassen2010b	20206776',
  'Lassen2009	19657123',
  'Turpie2009	19411100'
)

temp = data.frame(as = as, n1 = n1, cs = cs, n2 = n2, row.names = NULL)
test = rma.uni(yi = ors$yi, vi = ors$vi)
forest(test, slab = cites,
       ilab = temp,
       xlim = c(-14, 8),
       ylim = c(-3, 6),
       rows = 3,
       ilab.xpos = c(-7,-6,-5,-4),
       cex = .8,
       clim = c(-2, 3),
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
text(c(-6.5, -4.5), 6, c("LMWH", "FXaI"))
text(c(-7,-6,-5,-4), 5, c("n", "N", "n", "N"), cex = .8)
text(c(-12), 5, c("Author Year PMID"), cex = .8)
text(c(-1.2, 1.2), 5, c("Favors LMWH", "Favors FXaI"), cex = .8)
text(c(6), 5, c("Odds Ratio [95% CI]"), cex = .8)
title("Updated Forest Plot: TKR, Fatal PE, LMWH vs FXaI")

####################################################################


########################
### Symptomatic PE
##############

### LMWH
as = c(1, 0, 4, 2, 3, 7, 8, 4, 0)
n1 = c(90, 1529, 517, 109, 188, 1596, 1508, 1217, 70)

### FXaI
cs = c(0, 4, 1, 5, 2, 16, 5, 0, 2)
n2 = c(162, 1528, 517, 621, 561, 1599, 1526, 1201, 296)

### ni for each High risk of bias
nhigh = c()


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
cites = c(
  'Fuji2014d	22952213',
  'Lassen2010b	20206776',
  'Bauer2001	11794149',
  'Weitz2010	20886185',
  'Cohen2013	23782955',
  'Lassen2009	19657123',
  'Turpie2009	19411100',
  'Lassen2008	18579812',
  'Turpie2005	16241946'
)

temp = data.frame(as = as, n1 = n1, cs = cs, n2 = n2, row.names = NULL)
test = rma.uni(yi = ors$yi, vi = ors$vi)
forest(test, slab = cites,
       ilab = temp,
       xlim = c(-14, 8),
       ylim = c(-3, 12),
       rows = 9,
       ilab.xpos = c(-7,-6,-5,-4),
       cex = .8,
       clim = c(-2, 3),
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
text(c(-6.5, -4.5), 12, c("LMWH", "FXaI"))
text(c(-7,-6,-5,-4), 11, c("n", "N", "n", "N"), cex = .8)
text(c(-12), 11, c("Author Year PMID"), cex = .8)
text(c(-1.2, 1.2), 11, c("Favors LMWH", "Favors FXaI"), cex = .8)
text(c(6), 11, c("Odds Ratio [95% CI]"), cex = .8)
title("Updated Forest Plot: TKR, Symptomatic PE, LMWH vs FXaI")

####################################################################


########################
### Total DVT
##############

### LMWH
as = c(14, 14, 243, 98, 41, 22, 31, 92, 86, 160, 31)
n1 = c(66, 102, 997, 361, 295, 109, 188, 1122, 1508, 1217, 70)

### FXaI
cs = c(33, 3, 142, 45, 22, 156, 117, 89, 61, 79, 90)
n2 = c(152, 112, 971, 361, 299, 621, 561, 1142, 1526, 1201, 296)

### ni for each High risk of bias
nhigh = c()


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
cites = c(
  'Fuji2014d	22952213',
  'Zou2014	24695091',
  'Lassen2010b	20206776',
  'Bauer2001	11794149',
  'Fuji2014c	25294589',
  'Weitz2010	20886185',
  'Cohen2013	23782955',
  'Lassen2009	19657123',
  'Turpie2009	19411100',
  'Lassen2008	18579812',
  'Turpie2005	16241946'
)


temp = data.frame(as = as, n1 = n1, cs = cs, n2 = n2, row.names = NULL)
test = rma.uni(yi = ors$yi, vi = ors$vi)
forest(test, slab = cites,
       ilab = temp,
       xlim = c(-14, 8),
       ylim = c(-3, 14),
       rows = 11,
       ilab.xpos = c(-7,-6,-5,-4),
       cex = .8,
       clim = c(-2, 3),
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
text(c(-6.5, -4.5), 14, c("LMWH", "FXaI"))
text(c(-7,-6,-5,-4), 13, c("n", "N", "n", "N"), cex = .8)
text(c(-12), 13, c("Author Year PMID"), cex = .8)
text(c(-1.2, 1.2), 13, c("Favors LMWH", "Favors FXaI"), cex = .8)
text(c(6), 13, c("Odds Ratio [95% CI]"), cex = .8)
title("Updated Forest Plot: TKR, Total DVT, LMWH vs FXaI")

####################################################################




########################
### Symptomatic DVT
##############

### LMWH
as = c(1, 1, 1, 7, 4, 1, 2, 1, 7, 10, 20, 2)
n1 = c(90,102,54, 1529, 517, 295, 109, 188, 1596, 1508, 1217, 70)

### FXaI
cs = c(1, 0, 2, 3, 3, 4, 9, 1, 3, 6, 8, 2)
n2 = c(162, 112, 55, 1528, 517, 299, 621, 561, 1599, 1526, 1201, 296)

### ni for each High risk of bias
nhigh = c()


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
cites = c(
  'Fuji2014d	22952213',
  'Zou2014	24695091',
  'Bonneux2006	16387501',
  'Lassen2010b	20206776',
  'Bauer2001	11794149',
  'Fuji2014c	25294589',
  'Weitz2010	20886185',
  'Cohen2013	23782955',
  'Lassen2009	19657123',
  'Turpie2009	19411100',
  'Lassen2008	18579812',
  'Turpie2005	16241946'
)

temp = data.frame(as = as, n1 = n1, cs = cs, n2 = n2, row.names = NULL)
test = rma.uni(yi = ors$yi, vi = ors$vi)
forest(test, slab = cites,
       ilab = temp,
       xlim = c(-14, 8),
       ylim = c(-3, 15),
       rows = 12,
       ilab.xpos = c(-7,-6,-5,-4),
       cex = .8,
       clim = c(-2, 3),
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
text(c(-6.5, -4.5), 15, c("LMWH", "FXaI"))
text(c(-7,-6,-5,-4), 14, c("n", "N", "n", "N"), cex = .8)
text(c(-12), 14, c("Author Year PMID"), cex = .8)
text(c(-1.2, 1.2), 14, c("Favors LMWH", "Favors FXaI"), cex = .8)
text(c(6), 14, c("Odds Ratio [95% CI]"), cex = .8)
title("Updated Forest Plot: TKR, Symptomatic DVT, LMWH vs FXaI")

####################################################################


########################
### Proximal DVT
##############

### LMWH
as = c(1, 26, 20, 1, 3, 1, 11, 13, 20, 3)
n1 = c(74, 1199, 372, 295, 109, 188, 1207, 1508, 1217, 70)

### FXaI
cs = c(6, 9, 9, 0, 18, 7, 9, 3, 9, 9)
n2 = c(157, 1192, 368, 299, 621, 561, 1254, 1526, 1201, 296)

### ni for each High risk of bias
nhigh = c(sum(102, 112))


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
cites = c(
  'Fuji2014d	22952213',
  'Lassen2010b	20206776',
  'Bauer2001	11794149',
  'Fuji2014c	25294589',
  'Weitz2010	20886185',
  'Cohen2013	23782955',
  'Lassen2009	19657123',
  'Turpie2009	19411100',
  'Lassen2008	18579812',
  'Turpie2005	16241946'
)

temp = data.frame(as = as, n1 = n1, cs = cs, n2 = n2, row.names = NULL)
test = rma.uni(yi = ors$yi, vi = ors$vi)
forest(test, slab = cites,
       ilab = temp,
       xlim = c(-14, 8),
       ylim = c(-3, 13),
       rows = 10,
       ilab.xpos = c(-7,-6,-5,-4),
       cex = .8,
       clim = c(-2, 3),
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
text(c(-6.5, -4.5), 13, c("LMWH", "FXaI"))
text(c(-7,-6,-5,-4), 12, c("n", "N", "n", "N"), cex = .8)
text(c(-12), 12, c("Author Year PMID"), cex = .8)
text(c(-1.2, 1.2), 12, c("Favors LMWH", "Favors FXaI"), cex = .8)
text(c(6), 12, c("Odds Ratio [95% CI]"), cex = .8)
title("Updated Forest Plot: TKR, Proximal DVT, LMWH vs FXaI")

####################################################################


