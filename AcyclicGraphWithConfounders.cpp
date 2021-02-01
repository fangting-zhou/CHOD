// [[Rcpp::depends(RcppArmadillo)]]
# include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
bool checkDAG(mat G) {
  mat cG = G; bool isDAG = TRUE;
  
  uvec index = find(sum(cG, 0) != 0);
  while(index.n_elem < cG.n_rows) {
    cG = cG(index, index);
    if(accu(cG) == 0) break; else index = find(sum(cG, 0) != 0);
  }
  
  if(accu(cG) != 0) isDAG = FALSE;
  
  return(isDAG);
}

// [[Rcpp::export]]
List calculatex(cube b, mat x, mat bs) {
  mat factorx = zeros<mat>(x.n_rows, x.n_cols); 
  cube factorb = zeros<cube>(x.n_rows, x.n_cols, x.n_cols);
  
  for(int i = 0; i < x.n_rows; i ++) {
    for(int j = 0; j < x.n_cols; j ++) {
      for(int l = 0; l < x.n_cols; l ++) {
        if(j == l) factorb(i, j, l) = 1; else {
          factorb(i, j, l) = 0 - accu(bs.row(i) * vec(b.tube(j, l)));
        }
      }
    }
    
    factorx.row(i) = x.row(i) * trans(mat(factorb.row(i)));
  }
  
  return List::create(Named("factorb") = factorb, Named("factorx") = factorx);
}

// [[Rcpp::export]]
List calculates(mat factorx, mat s) {
  mat svec; vec sval; eig_sym(sval, svec, s);
  
  mat factors = zeros<mat>(factorx.n_rows, factorx.n_cols);
  for(int i = 0; i < factorx.n_rows; i ++) factors.row(i) = factorx.row(i) * trans(svec);
  
  return List::create(Named("svec") = svec, Named("sval") = sval, Named("factors") = factors);
}

// [[Rcpp::export]]
mat updates(mat factorx) {
  mat anew = trans(factorx) * factorx + diagmat(ones<vec>(factorx.n_cols)); 
  double bnew = factorx.n_rows + factorx.n_cols;
  
  mat s = iwishrnd(anew, bnew);
  
  return(s);
}

// [[Rcpp::export]]
List updater(cube b, cube factorb, mat r, mat x, mat factorx, mat factors, mat bs, mat svec, vec sval, double v, double u, double s0 = 1) {
  mat rnew = r, xden = factorx, xnum = factorx, sden = factors, snum = factors;
  cube bnew = b, bden = factorb, bnum = factorb;
  
  for(int j = 0; j < r.n_rows; j ++) {
    for(int l = 0; l < r.n_cols; l ++) {
      double flip = randu();
      if((flip < 0.5) & (r(j, l) == 1)) {
        for(int k = 0; k < b.n_slices; k ++) bnew(j, l, k) = bnew(j, l, k) + s0 * randn();
        
        double den = 0 - accu(pow(b.tube(j, l), 2)) / (2 * v);
        double num = 0 - accu(pow(bnew.tube(j, l), 2)) / (2 * v);
        
        for(int i = 0; i < x.n_rows; i ++) {
          bnum(i, j, l) = 0 - accu(bs.row(i) * vec(bnew.tube(j, l)));
          xnum(i, j) =  xnum(i, j) - x(i, l) * (bden(i, j, l) - bnum(i, j, l));
          
          snum.row(i) = snum.row(i) - (xden(i, j) - xnum(i, j)) * trans(svec.col(j)); 
          
          den = den - accu(pow(trans(sden.row(i)), 2) / sval) / 2;
          num = num - accu(pow(trans(snum.row(i)), 2) / sval) / 2;
        }
        
        if((num - den) >= log(randu())) {
          b = bnew; bden = bnum; xden = xnum; sden = snum;
        } else {
          bnew = b; bnum = bden; xnum = xden; snum = sden;
        }
      }
      
      if(flip >= 0.5) {
        rnew(j, l) = 1 - r(j, l);
        if(checkDAG(rnew)) {
          if(rnew(j, l) == 0) bnew.tube(j, l) = zeros<vec>(b.n_slices); else bnew.tube(j, l) = randn(b.n_slices) * sqrt(v);
          
          double den = log(u) * r(j, l) + log(1 - u) * (1 - r(j, l));
          double num = log(u) * rnew(j, l) + log(1 - u) * (1 - rnew(j, l));
          
          for(int i = 0; i < x.n_rows; i ++) {
            bnum(i, j, l) = 0 - accu(bs.row(i) * vec(bnew.tube(j, l)));
            xnum(i, j) =  xnum(i, j) - x(i, l) * (bden(i, j, l) - bnum(i, j, l));
            
            snum.row(i) = snum.row(i) - (xden(i, j) - xnum(i, j)) * trans(svec.col(j)); 
            
            den = den - accu(pow(trans(sden.row(i)), 2) / sval) / 2;
            num = num - accu(pow(trans(snum.row(i)), 2) / sval) / 2;
          }
          
          if((num - den) >= log(randu())) {
            r = rnew; b = bnew; bden = bnum; xden = xnum; sden = snum;
          } else {
            rnew = r; bnew = b; bnum = bden; xnum = xden; snum = sden;
          }
        } else rnew(j, l) = 0;
      }
    }
  }
  
  return List::create(Named("r") = r, Named("b") = b, Named("factorx") = xden, Named("factorb") = bden, Named("factors") = sden);
}

// [[Rcpp::export]]
double updateu(mat r, double a0 = 0.5, double b0 = 0.5) {
  double anew = a0 + accu(r);
  double bnew = b0 + accu(1 - r); double c0 = 1; 
    
  double x = randg(distr_param(anew, c0));
  double y = randg(distr_param(bnew, c0));
      
  double u = x / (x + y);
  
  return(u);
}

// [[Rcpp::export]]
double updatev(cube b, double a0 = 0.1, double b0 = 1) {
  double anew = a0 + accu(b != 0) / 2;
  double bnew = b0 + accu(pow(b, 2)) / 2;
  
  double v = 1 / randg(distr_param(anew, 1 / bnew));
  
  return v;
}