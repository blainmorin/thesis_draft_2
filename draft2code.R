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

unadjusted = exp(CIRE())
adjusted = exp(b.CI())
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
       ilab.xpos = c(-7,-6,-5,-4), 
       atransf = exp)


