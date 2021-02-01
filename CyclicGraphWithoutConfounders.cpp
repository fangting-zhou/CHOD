// [[Rcpp::depends(RcppArmadillo)]]
# include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
List calculatey(cube b, mat x, mat bs) {
  vec factord = zeros<vec>(x.n_rows);
  mat factory = zeros<mat>(x.n_rows, x.n_cols); 
  cube factorb = zeros<cube>(x.n_rows, x.n_cols, x.n_cols);
  
  for(int i = 0; i < x.n_rows; i ++) {
    for(int j = 0; j < x.n_cols; j ++) {
      for(int l = 0; l < x.n_cols; l ++) {
        if(j == l) factorb(i, j, l) = 1; else {
          factorb(i, j, l) = 0 - accu(bs.row(i) * vec(b.tube(j, l)));
        }
      }
    }
    
    factord(i) = log(abs(det(mat(factorb.row(i)))));
    factory.row(i) = x.row(i) * trans(mat(factorb.row(i)));
  }
  
  return List::create(Named("factord") = factord, Named("factorb") = factorb, Named("factory") = factory);
}

// [[Rcpp::export]]
vec updates(mat factory, double a0 = 0.001, double b0 = 0.001) {
  vec s = zeros<vec>(factory.n_cols);
  
  double anew = a0 + factory.n_rows / 2;
  for(int j = 0; j < factory.n_cols; j ++) {
    double bnew = b0 + accu(pow(factory.col(j), 2)) / 2;
    
    s(j) = 1 / randg(distr_param(anew, 1 / bnew));
  }
  
  return s;
}

// [[Rcpp::export]]
List updater(cube b, cube factorb, mat r, mat x, mat factory, mat bs, vec s, vec factord, double w, double u) {
  vec dden = factord, dnum = factord;
  mat rnew = r, yden = factory, ynum = factory;
  cube bnew = b, bden = factorb, bnum = factorb;
  
  for(int j = 0; j < r.n_rows; j ++) {
    for(int l = 0; l < r.n_cols; l ++) {
      if(j != l) {
        rnew(j, l) = 1 - r(j, l);
        if(rnew(j, l) == 0) bnew.tube(j, l) = zeros<vec>(b.n_slices); else bnew.tube(j, l) = randn(b.n_slices) * sqrt(w);
        
        double den = log(u) * r(j, l) + log(1 - u) * (1 - r(j, l));
        double num = log(u) * rnew(j, l) + log(1 - u) * (1 - rnew(j, l));
        
        for(int i = 0; i < x.n_rows; i ++) {
          bnum(i, j, l) = 0 - accu(bs.row(i) * vec(bnew.tube(j, l)));
          ynum(i, j) =  ynum(i, j) - x(i, l) * (bden(i, j, l) - bnum(i, j, l));
          dnum(i) = log(abs(det(mat(bnum.row(i)))));
          
          den = den - accu(pow(trans(yden.row(i)), 2) / s) / 2 + dden(i);
          num = num - accu(pow(trans(ynum.row(i)), 2) / s) / 2 + dnum(i);
        }
        
        if((num - den) >= log(randu())) {
          r = rnew; b = bnew; bden = bnum; yden = ynum; dden = dnum;
        } else {
          rnew = r; bnew = b; bnum = bden; ynum = yden; dnum = dden;
        }
      }
    }
  }
  
  return List::create(Named("r") = r, Named("b") = b, Named("factory") = yden, Named("factorb") = bden, Named("factord") = dden);
}

// [[Rcpp::export]]
List updateb(cube b, cube factorb, mat r, mat x, mat factory, mat bs, vec s, vec factord, double w, double s0 = 0.1) {
  vec dden = factord, dnum = factord;
  mat rnew = r, yden = factory, ynum = factory;
  cube bnew = b, bden = factorb, bnum = factorb;
  
  for(int j = 0; j < r.n_rows; j ++) {
    for(int l = 0; l < r.n_cols; l ++) {
      if(r(j, l) == 1) {
        for(int k = 0; k < b.n_slices; k ++) bnew(j, l, k) = bnew(j, l, k) + s0 * randn();
        
        double den = 0 - accu(pow(b.tube(j, l), 2)) / (2 * w);
        double num = 0 - accu(pow(bnew.tube(j, l), 2)) / (2 * w);
        
        for(int i = 0; i < x.n_rows; i ++) {
          bnum(i, j, l) = 0 - accu(bs.row(i) * vec(bnew.tube(j, l)));
          ynum(i, j) =  ynum(i, j) - x(i, l) * (bden(i, j, l) - bnum(i, j, l));
          dnum(i) = log(abs(det(mat(bnum.row(i)))));
          
          den = den - accu(pow(trans(yden.row(i)), 2) / s) / 2 + dden(i);
          num = num - accu(pow(trans(ynum.row(i)), 2) / s) / 2 + dnum(i);
        }
        
        if((num - den) >= log(randu())) {
          b = bnew; bden = bnum; yden = ynum; dden = dnum;
        } else {
          bnew = b; bnum = bden; ynum = yden; dnum = dden;
        }
      }
    }
  }
  
  return List::create(Named("b") = b, Named("factory") = yden, Named("factorb") = bden, Named("factord") = dden);
}

// [[Rcpp::export]]
double updatew(cube b, double a0 = 0.1, double b0 = 1) {
  double anew = a0 + accu(b != 0) / 2;
  double bnew = b0 + accu(pow(b, 2)) / 2;
  
  double w = 1 / randg(distr_param(anew, 1 / bnew));
  
  return w;
}

// [[Rcpp::export]]
double updateu(mat r, double a0 = 0.5, double b0 = 0.5) {
  double anew = a0 + accu(r); double c0 = 1;
  double bnew = b0 + accu(1 - r);
  
  double x = randg(distr_param(anew, c0));
  double y = randg(distr_param(bnew, c0));
  
  double u = x / (x + y);
  
  return(u);
}

