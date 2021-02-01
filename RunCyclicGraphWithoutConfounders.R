library(splines2); library(igraph); library(ggplot2)

Rcpp::sourceCpp("CyclicGraphWithoutConfounders.cpp")

generateGraph = function(p0) {
  r0 = matrix(0, p0, p0); r0[runif(r0) < 1 / p0] = 1
  while(sum(diag(r0)) != 0) {
    r0 = matrix(0, p0, p0); r0[runif(r0) < 1 / p0] = 1
  }
  
  f0 = r0 * sample(3, p0 * p0, replace = TRUE)
  
  return(list(r0 = r0, f0 = f0))
}

G0 = generateGraph(10); r0 = G0$r0; f0 = G0$f0

b1 = function(z) 0.8 * z
b2 = function(z) 0.8 * cos(3 * z)
b3 = function(z) 0.8 * tanh(2 * z)

generateData = function(f0, n0, s = 0.5) {
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
    
    x0[i, ] = solve(b0) %*% rnorm(p0, sd = s)
  }
  
  return(x0)
}

res = matrix(0, 3, 50)

for(num in 1 : 50) {
  x0 = generateData(f0, 500)
  
  n0 = nrow(x0); p0 = ncol(x0); z0 = seq(-1, 1, length = n0); bknot = c(-1, 1)
  degree = 3; len = 8; iknot = seq(-1, 1, length = len)[2 : (len - 1)]
  bs = bSpline(z0, knots = iknot, degree = degree, Boundary.knots = bknot, intercept = TRUE)
  r = matrix(0, p0, p0); b = array(0, c(p0, p0, ncol(bs))); u = 0.5; w = 1
  
  factors = calculatey(b, x0, bs)
  factorb = factors$factorb; factory = factors$factory; factord = factors$factord
  
  iter = 1; maxit = 2000; burn = 1000; int = 5; recordr = NULL
  
  while(iter <= maxit) {
    s = updates(factory)
    R = updater(b, factorb, r, x0, factory, bs, s, factord, w, u)
    r = R$r; b = R$b; factorb = R$factorb; factord = R$factord; factory = R$factory
    w = updatew(b)
    u = updateu(r)
    
    recordr = c(recordr, list(r))
    iter = iter + 1
  }
  
  r = matrix(0, p0, p0)
  for(k in seq(burn + 1, maxit, int)) r = r + recordr[[k]]
  r = r / length(seq(burn + 1, maxit, int))
  r = matrix(as.numeric(r > 0.5), nrow(r), ncol(r))
  
  truePositive = sum(r[r0 != 0])
  trueNegative = sum(1 - r[r0 == 0])
  falsePositive = sum(r[r0 == 0])
  falseNegative = sum(1 - r[r0 != 0])
  
  TPR = truePositive / (truePositive + falseNegative)
  FDR = falsePositive / (falsePositive + truePositive)
  MCC = (truePositive * trueNegative - falsePositive * falseNegative) / 
    sqrt((truePositive + falsePositive) * (truePositive + falseNegative) * (trueNegative + falsePositive) * (trueNegative + falseNegative))
  
  res[1, num] = TPR
  res[2, num] = FDR
  res[3, num] = MCC
}

