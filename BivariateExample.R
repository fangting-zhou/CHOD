library(mvtnorm); library(splines2)

Rcpp::sourceCpp("CyclicGraphWithConfounders.cpp")

n = 1000; s = matrix(c(1, 0.5, 0.5, 1), 2, 2); z = seq(-1, 1, length = n)

x = matrix(0, n, 2)
for(i in 1 : n) {
  b0 = diag(2); 
  # x1 -> x2 & x2 <- x1 comment out if necessary
  b0[2, 1] = b0[2, 1] - 0.5 * sin(pi * z[i]); b0[1, 2] = b0[1, 2] - 0.5 * sin(pi * z[i])
  x[i, ] = solve(b0) %*% t(rmvnorm(1, sigma = s))
}

x0 = x; z0 = z; n0 = nrow(x0); p0 = ncol(x0)
bknot = c(-1, 1); degree = 3; len = 4; iknot = seq(-1, 1, length = len)[2 : (len - 1)]
bs = bSpline(z0, knots = iknot, degree = degree, Boundary.knots = bknot)
r = matrix(0, p0, p0); b = array(0, c(p0, p0, ncol(bs))); w = 1; u = 0.1

factors = calculatey(b, x0, bs)
factorb = factors$factorb; factory = factors$factory; factord = factors$factord

iter = 1; maxit = 1000; burnin = 500; interval = 5; recordb = recordr = NULL

while(iter <= maxit) {
  s = updates(factory)
  factors = calculatex(factory, s)
  svec = factors$svec; sval = factors$sval; factorx = factors$factorx
  R = updater(b, factorb, r, x0, factory, factorx, bs, svec, sval, factord, w, u)
  r = R$r; b = R$b; factorb = R$factorb; factory = R$factory; factord = R$factord; factorx = R$factorx
  B = updateb(b, factorb, r, x0, factory, factorx, bs, svec, sval, factord, w)
  b = B$b; factorb = B$factorb; factory = B$factory; factord = B$factord; factorx = B$factorx
  u = updateu(r)
  w = updatew(b)
  
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

xx = x * x

var1 = var2 = rep(0, n)
for(i in 1 : n) {
  var1[i] = sum(exp(- 10 * (z[i] - z)^2) * xx[,1]) / sum(exp(- 10 * (z[i] - z)^2))
  var2[i] = sum(exp(- 10 * (z[i] - z)^2) * xx[,2]) / sum(exp(- 10 * (z[i] - z)^2))
}

par(mfrow = c(1, 3))

plot(x[, 1], x[, 2], pch = 20, cex = 0.8, xlab = "", ylab = "", cex.axis = 1.5, cex.lab = 2)
plot(z, var1, xlab = "", ylab = "", type = "l", ylim = c(0, 2.6), cex.axis = 1.5, cex.lab = 2)
plot(z, var2, xlab = "", ylab = "", type = "l", ylim = c(0, 2.6), cex.axis = 1.5, cex.lab = 2)
