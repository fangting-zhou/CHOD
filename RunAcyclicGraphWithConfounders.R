library(splines2); library(igraph); library(mvtnorm); library(matrixcalc); library(ggplot2)

Rcpp::sourceCpp("AcyclicGraphWithConfounders.cpp")

generateGraph = function(p0) {
  r0 = matrix(0, p0, p0); r0[which(lower.tri(r0))] = runif(sum(lower.tri(r0))) < 2 / p0
  while(sum((rowSums(r0) + colSums(r0)) == 0) != 0) {
    r0 = matrix(0, p0, p0); r0[which(lower.tri(r0))] = runif(sum(lower.tri(r0))) < 2 / p0
  }
  
  f0 = r0 * sample(3, p0 * p0, replace = TRUE)
  
  s0 = matrix(0, p0, p0)
  s0[lower.tri(s0)] = (runif(sum(lower.tri(s0))) < 0.01)
  s0[lower.tri(s0)] = s0[lower.tri(s0)] * runif(sum(lower.tri(s0)), -1, 1)
  s0 = s0 + t(s0); diag(s0) = 1
  
  while(!is.positive.definite(s0)) {
    s0 = matrix(0, p0, p0)
    s0[lower.tri(s0)] = (runif(sum(lower.tri(s0))) < 0.01)
    s0[lower.tri(s0)] = s0[lower.tri(s0)] * runif(sum(lower.tri(s0)), -1, 1)
    s0 = s0 + t(s0); diag(s0) = 1
  }
  
  return(list(f0 = f0, s0 = s0))
}

G0 = generateGraph(10); f0 = G0$f0; s0 = G0$s0

b1 = function(z) z
b2 = function(z) cos(3 * z)
b3 = function(z) tanh(3 * z)

generateData = function(f0, s0, n0) {
  p0 = nrow(f0); z0 = seq(-1, 1, length = n0)
  
  x0 = matrix(0, n0, p0)
  for(i in 1 : n0) {
    b0 = diag(p0)
    for(j in 1 : p0) {
      for(l in 1 : p0) {
        if(f0[j, l] == 1) b0[j, l] = b0[j, l] - b1(z0[i])
        if(f0[j, l] == 2) b0[j, l] = b0[j, l] - b2(z0[i])
        if(f0[j, l] == 3) b0[j, l] = b0[j, l] - b3(z0[i])
      }
    }
    
    x0[i, ] = solve(b0) %*% t(rmvnorm(1, sigma = s0))
  }
  
  return(list(z0 = z0, x0 = x0))
}

res = matrix(0, 3, 50); num = 1

while(num <= 50) {
  data = generateData(f0, s0, 500)
  z0 = data$z0; x0 = data$x0
  
  n0 = nrow(x0); p0 = ncol(x0); bknot = c(-1, 1); degree = 3
  len = 9; iknot = seq(-1, 1, length = len)[2 : (len - 1)]
  bs = bSpline(z0, knots = iknot, degree = degree, Boundary.knots = bknot)
  r = matrix(0, p0, p0); b = array(0, c(p0, p0, ncol(bs))); v = 1; u = 0.1
  
  xfactor = calculatex(b, x0, bs)
  factorb = xfactor$factorb; factorx = xfactor$factorx
  
  iter = 1; maxit = 2000; burnin = 1000; interval = 5; recordb = recordr = NULL
  
  while(iter <= maxit) {
    s = updates(factorx)
    sfactor = calculates(factorx, s)
    svec = sfactor$svec; sval = sfactor$sval; factors = sfactor$factors
    R = updater(b, factorb, r, x0, factorx, factors, bs, svec, sval, v, u)
    r = R$r; b = R$b; factorb = R$factorb; factorx = R$factorx; factors = R$factors
    u = updateu(r)
    v = updatev(b)
    
    recordb = c(recordb, list(b))
    recordr = c(recordr, list(r))
    
    iter = iter + 1
  }

  
  r = matrix(0, p0, p0); b = array(0, dim = c(p0, p0, ncol(bs)))
  for(k in seq(burnin + 1, maxit, interval)) {
    r = r + recordr[[k]]; b = b + recordb[[k]]
  }
  r = r / length(seq(burnin + 1, maxit, interval))
  b = b / length(seq(burnin + 1, maxit, interval))
  
  R = matrix(as.numeric(r > 0.5), nrow(r), ncol(r))
  
  truePositive = sum(R[f0 != 0])
  trueNegative = sum(1 - R[f0 == 0])
  falsePositive = sum(R[f0 == 0])
  falseNegative = sum(1 - R[f0 != 0])
  
  res[1, num] = truePositive / (truePositive + falseNegative)
  res[2, num] = falsePositive / (falsePositive + truePositive)
  res[3, num] = (truePositive * trueNegative - falsePositive * falseNegative) / 
    sqrt((truePositive + falsePositive) * (truePositive + falseNegative) * (trueNegative + falsePositive) * (trueNegative + falseNegative))
  
  num = num + 1
}

