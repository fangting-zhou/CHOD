library(splines2); library(igraph); library(mvtnorm); library(matrixcalc); library(ggplot2); library(umap)

Rcpp::sourceCpp("AcyclicGraphWithConfounders.cpp")

# three-node graph
r0 = matrix(0, 3, 3); r0[1, 3] = r0[2, 3] = r0[1, 2] = 1
index = 1 : 2; lr = r0[index, index]

# four-node graph
r0 = matrix(0, 4, 4); r0[1, 2] = r0[1, 4] = r0[2, 4] = r0[3, 4] = 1
index = 1 : 3; lr = r0[index, index]

# misspecified model: case 1
generateData = function(r0, n0, s0) {
  p0 = nrow(r0); z0 = seq(-1, 1, length = n0)
  
  x0 = matrix(0, n0, p0)
  for(i in 1 : n0) {
    b0 = diag(p0)
    for(j in 1 : p0) {
      for(l in 1 : p0) {
        if(r0[j, l] == 1) b0[j, l] = b0[j, l] - (1.0 * z0[i]^2 + 0.00)
        # b0[j, l] = b0[j, l] - (0.8 * z0[i]^2 + 0.08)
        # b0[j, l] = b0[j, l] - (0.5 * z0[i]^2 + 0.18)
        # b0[j, l] = b0[j, l] - (0.2 * z0[i]^2 + 0.28)
        # b0[j, l] = b0[j, l] - (0.0 * z0[i]^2 + 0.34)
      }
    }
    
    x0[i, ] = solve(b0) %*% rnorm(p0, s0)
  }
  
  return(list(z0 = z0, x0 = x0))
}

# misspecified model: case 2
generateData = function(r0, n0, s0 = 0.5) {
  p0 = nrow(r0); z0 = seq(-1, 1, length = n0)
  
  x0 = matrix(0, n0, p0)
  for(i in 1 : n0) {
    b0 = diag(p0)
    for(j in 1 : p0) {
      for(l in 1 : p0) {
        if(r0[j, l] == 1) b0[j, l] = b0[j, l] - ceiling(10 * abs(z0[i])) / 10
      }
    }
    
    x0[i, ] = solve(b0) %*% (rep(z0[i], p0) + rnorm(p0, s0))
  }
  
  return(list(z0 = z0, x0 = x0))
}

rR = NULL; num = 1

while(num <= 50) {
  # use in case 1
  data = generateData(r0, 250, c(0.5, 0.3, 0.2))
  # s0 = c(0.5, 0.3, 0.5, 0.2) in four-node graph
  z0 = data$z0; x0 = data$x0; xl = x0[, index]
  
  # use in case 2
  # data = generateData(r0, 250)
  # z0 = data$z0; x0 = data$x0; xl = x0[, index]
  # zhat = umap(xl, n_components = 1)$layout; zhat = (2* zhat - min(zhat) - max(zhat)) / (max(zhat) - min(zhat))
  # for(j in 1 : ncol(xl)) {
  # xl[, j] = lm(xl[, j] ~ zhat)$residual
  # }
  
  n0 = nrow(xl); p0 = ncol(xl); bknot = c(-1, 1); degree = 3; len = 3; iknot = seq(-1, 1, length = len)[2 : (len - 1)]
  bs = bSpline(z0, knots = iknot, degree = degree, Boundary.knots = bknot, intercept = TRUE)
  # len = 10; iknot = quantile(zhat, seq(0, 1, length = len))[2 : (len - 1)]
  # bs = bSpline(zhat, knots = iknot, degree = degree, Boundary.knots = bknot, intercept = TRUE)
  r = matrix(0, p0, p0); b = array(0, c(p0, p0, ncol(bs))); v = 1; u = 0.5
  
  xfactor = calculatex(b, xl, bs)
  factorb = xfactor$factorb; factorx = xfactor$factorx
  
  iter = 1; maxit = 2000; burnin = 1000; interval = 5; recordr = NULL
  
  while(iter <= maxit) {
    s = updates(factorx)
    sfactor = calculates(factorx, s)
    svec = sfactor$svec; sval = sfactor$sval; factors = sfactor$factors
    R = updater(b, factorb, r, xl, factorx, factors, bs, svec, sval, v, u)
    r = R$r; b = R$b; factorb = R$factorb; factorx = R$factorx; factors = R$factors
    u = updateu(r)
    v = updatev(b)
  
    recordr = c(recordr, list(r))
    
    iter = iter + 1
  }
  
  r = matrix(0, p0, p0)
  for(k in seq(burnin + 1, maxit, interval)) r = r + recordr[[k]]
  r = r / length(seq(burnin + 1, maxit, interval))
  
  rR = c(rR, list(r))
  
  num = num + 1
}

cutoff = seq(0, 1.01, length = 200); tpr = fpr = rep(0, 200)

for(i in 1 : 200) {
  for(j in 1 : 50) {
    r = matrix(rR[[j]] >= cutoff[i], length(index), length(index)); diag(r) = 0
    
    truePositive = sum(r[lr != 0])
    trueNegative = sum(1 - r[lr == 0]) - nrow(r)
    falsePositive = sum(r[lr == 0])
    falseNegative = sum(1 - r[lr != 0])
    
    xx = truePositive / (truePositive + falseNegative)
    yy = falsePositive / (falsePositive + trueNegative)
    
    tpr[i] = tpr[i] + xx
    fpr[i] = fpr[i] + yy
  }
  
  tpr[i] = tpr[i] / 50
  fpr[i] = fpr[i] / 50
}

AUC(fpr, tpr)
