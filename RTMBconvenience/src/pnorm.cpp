// [[Rcpp::depends(TMB)]]
#include <RTMB.h>

#include "fromR.hpp"


namespace adaptive {

  template<class Float>
  Float pnorm1_1x(Float x, Float lower_tail, Float log_p){
    return fromR::pnorm5_raw(x,Float(0.0),Float(1.0),lower_tail,log_p);    
  };


  template<class Float>
  Float log_ipnorm1_1x(Float x, Float y){
    if(x > y) Rf_error("y should be larger than x!");
    if(x < 0.5){ // Use lower tail
      Float v1 = fromR::pnorm5_raw(x,Float(0.0),Float(1.0),Float(1),Float(1));
      Float v2 = fromR::pnorm5_raw(y,Float(0.0),Float(1.0),Float(1),Float(1));
      return v2 + atomic::robust_utils::R_Log1_Exp(v1 - v2);
    }else{ // Use upper tail
      Float v2 = fromR::pnorm5_raw(x,Float(0.0),Float(1.0),Float(0),Float(1));
      Float v1 = fromR::pnorm5_raw(y,Float(0.0),Float(1.0),Float(0),Float(1));
      return v2 + atomic::robust_utils::R_Log1_Exp(v1 - v2);
    }    
  };

  
  template<class Float>
  Float pt_raw(Float x, Float n, Float lower_tail, Float log_p){
    return fromR::pt(x,n,(int)trunc(lower_tail),(int)trunc(log_p));    
  };

  template<class Float>
  Float ps_raw (Float Pf) {
    Float af = 0.0;
    if(Pf > 10.0){
      af = (Pf - 1.0) / Pf;
    }else if(Pf > 1e-6){
      af = (Pf-1.0 + pow(1.0-Pf/10.0,5)) / Pf;
    }else if(Pf > -1e-6){
      af = 0.5;
    }else if(Pf > -10.0){
      af = (pow(1.0+Pf/10.0,5) - 1.0) / Pf;
    }else{
      af = -1.0 / Pf;
    }
    return af;
  }


  template<class Float>
  Float logspace_add2_raw (Float logx, Float logy) {
    // Was:
    //  fmax2 (logx, logy) + log1p (exp (-fabs (logx - logy)));
    if(logx == R_NegInf && logy == R_NegInf)
      return(R_NegInf);
    if(logx == R_NegInf)
      return(logy);
    if(logy == R_NegInf)
      return(logx);
    return ( logx < logy ?
	     logy + log1p (exp (logx - logy)) :
	     logx + log1p (exp (logy - logx)) );
  }

  template<class Float>
  Float logspace_sub2_raw (Float logx, Float logy) {
    if(logx == logy)
      return R_NegInf;
    if(logy == R_NegInf)
      return(logx);
    if(logx < logy)
      Rf_error("logx < logy in logspace_sub2");
    return logx + atomic::robust_utils::R_Log1_Exp(logy - logx);
  }

  template<class Float>
  Float quantreg_raw (Float x, Float tau) {
    if(x < 0)
      return x * (tau - 1.0);
    return x * tau;
  }


  template<class Float>
  Float login_log_besselI_raw(Float logx, Float nu){
    if(logx == R_NegInf){
      if(fabs(nu) < 1e-16) return(0);
      return(R_NegInf);
    }
    if(logx == R_PosInf) return(R_PosInf);
    if(nu < 0){
      if(fabs(nu - floor(nu)) > 1e-16)
	Rf_error("Not implemented for negative fractions");
      nu = -nu;
    }
    Float r = R_NegInf;
    Float dr = R_PosInf;
    Float m = 0;
    while(dr > -30){
      dr = (m * 2.0 + nu) * (logx - log(2.0));
      m += 1.0;
      //dr -= atomic::Rmath::D_lgamma(m,0.0) + atomic::Rmath::D_lgamma(m+nu,0.0);
      dr -= fromR::lgammafn(m) + fromR::lgammafn(m+nu);
      r = logspace_add2_raw(r,dr);
    }
    return r;
  }

  template<class Float>
  Float log_besselI_raw (Float x, Float nu) {
    Float v = atomic::bessel_utils::bessel_i(x, nu, 2.0);
    return log(v) + x;
  }

  template<class Float>
  Float log_MarcumQ_raw(Float a, Float b, Float nu){
    if(fabs(b) < 1e-16 || a == R_PosInf || nu == R_PosInf) return 0.0;
    if(b == R_NegInf) return 0.0;

    Float r = R_NegInf;
    Float dr = R_PosInf;
    Float k = 1 - nu;
    Float ab = a * b;
    while(dr > -30){
      dr = k * (log(a) - log(b)) + log_besselI_raw(ab,(Float)-k);
      k += 1.0;
      r = logspace_add2_raw(r,dr);
    }
    return -(a*a + b*b) / 2.0 + r;
  }

  template<class Float>
  Float login_log_MarcumQ_raw(Float loga, Float logb, Float nu){
    if(logb == R_NegInf || loga == R_PosInf || nu == R_PosInf) return 0.0;

    Float r = R_NegInf;
    Float dr = R_PosInf;
    Float k = 1 - nu;
    Float logab = loga + logb;
    while(dr > -30){
      dr = k * (loga - logb) + login_log_besselI_raw(logab,(Float)-k);
      k += 1.0;
      r = logspace_add2_raw(r,dr);
    }
    return -exp(logspace_add2_raw((Float)(2.0 * loga), (Float)(2.0 * logb))) / 2.0 + r;
  }

  template<class Float>
  Float log_Marcum1mQ_raw(Float a, Float b, Float nu){
    if(fabs(b) < 1e-16 || a == R_PosInf || nu == R_PosInf) return R_NegInf;
    if(b == R_NegInf) return R_PosInf;

    Float r = R_NegInf;
    Float dr = R_PosInf;
    Float alpha = nu;
    Float ab = a * b;
    while(dr > -30){
      dr = (Float)(alpha) * (log(b) - log(a)) + log_besselI_raw(ab,alpha);
      alpha += 1.0;
      r = logspace_add2_raw(r,dr);
    }
    return -(a*a + b*b) / 2.0 + r;
  }

  template<class Float>
  Float login_log_Marcum1mQ_raw(Float loga, Float logb, Float nu){
    if(logb == R_NegInf || loga == R_PosInf || nu == R_PosInf) return R_NegInf;

    Float r = R_NegInf;
    Float dr = R_PosInf;
    Float alpha = nu;
    Float logab = loga + logb;
    while(dr > -30){
      dr = (Float)(alpha) * (logb - loga) + login_log_besselI_raw(logab,alpha);
      alpha += 1.0;
      r = logspace_add2_raw(r,dr);
    }
    return -exp(logspace_add2_raw((Float)(2.0 * loga), (Float)(2.0 * logb))) / 2.0 + r;
  }

  template<class Float>
  Float logdrice_raw(Float logx, Float lognu, Float logsigma){
    return logx - 2.0 * logsigma - exp(logspace_add2_raw((Float)(2.0 * logx), (Float)(2.0 * lognu)) - log(2.0) - 2.0 * logsigma) + login_log_besselI_raw((Float)(logx+lognu-2.0 * logsigma), (Float)0);
  }

  template<class Float>
  Float logprice_raw(Float logx, Float lognu, Float logsigma, Float lower_tail){
    Float v1 = (lognu - logsigma);
    Float v2 = (logx - logsigma);
    Float ut = login_log_MarcumQ_raw(v1,v2,(Float)1);
    if((int)trunc(lower_tail))
      return logspace_sub2_raw((Float)0.0,ut);
    return ut;
  }


  
}

TMB_BIND_ATOMIC(pnorm1_2x,100,adaptive::pnorm1_1x(x[0], x[1], x[2]))  
template<class Type>
Type pnorm5_(Type x, Type mu, Type sigma, int lower_tail, int log_p){
  vector<Type> tx(4);
  tx[0] = (x-mu) / sigma;
  // tx[1] = mu;
  // tx[2] = sigma;
  tx[1] = (Type)lower_tail;
  tx[2] = (Type)log_p;
  tx[3] = 0; // extra argument for derivative order
  Type res = pnorm1_2x(CppAD::vector<Type>(tx))[0];
  return res;
}

// [[Rcpp::export]]
Rcpp::ComplexVector pnorm5_ad(Rcpp::ComplexVector x, Rcpp::ComplexVector mu, Rcpp::ComplexVector sigma, bool lower_tail, bool log_p ) {
  CHECK_INPUT(x); CHECK_INPUT(mu); CHECK_INPUT(sigma);
  size_t n1 = x.size();
  size_t n2 = mu.size();
  size_t n3 = sigma.size();
  int nmax = std::max({n1, n2, n3});
  int nmin = std::min({n1, n2, n3});
  int n = (nmin == 0 ? 0 : nmax);
  Rcpp::ComplexVector ans(n);
  const ad* X1 = adptr(x); const ad* X2 = adptr(mu); const ad* X3 = adptr(sigma);
  ad* Y = adptr(ans);
  for (int i=0; i<n; i++) Y[i] = pnorm5_(X1[i % n1], X2[i % n2], X3[i % n3], lower_tail, log_p);
  return as_advector(ans);
}

TMB_BIND_ATOMIC(log_ipnorm1_2x,11,adaptive::log_ipnorm1_1x(x[0], x[1]))  
template<class Type>
Type log_ipnorm_(Type x, Type y, Type mu, Type sigma){
  vector<Type> tx(3);
  tx[0] = (x-mu) / sigma;
  tx[1] = (y-mu) / sigma;
  tx[2] = 0; // extra argument for derivative order
  Type res = log_ipnorm1_2x(CppAD::vector<Type>(tx))[0];
  return res;
}

// [[Rcpp::export]]
Rcpp::ComplexVector log_ipnorm_ad(Rcpp::ComplexVector x, Rcpp::ComplexVector y, Rcpp::ComplexVector mu, Rcpp::ComplexVector sigma) {
  CHECK_INPUT(x); CHECK_INPUT(y); CHECK_INPUT(mu); CHECK_INPUT(sigma);
  size_t n1 = x.size();
  size_t n2 = y.size();
  size_t n3 = mu.size();
  size_t n4 = sigma.size();
  int nmax = std::max({n1, n2, n3, n4});
  int nmin = std::min({n1, n2, n3, n4});
  int n = (nmin == 0 ? 0 : nmax);
  Rcpp::ComplexVector ans(n);
  const ad* X1 = adptr(x); const ad* X2 = adptr(y); const ad* X3 = adptr(mu); const ad* X4 = adptr(sigma);
  ad* Y = adptr(ans);
  for (int i=0; i<n; i++) Y[i] = log_ipnorm_(X1[i % n1], X2[i % n2], X3[i % n3], X4[i % n4]);
  return as_advector(ans);
}


TMB_BIND_ATOMIC(pt_at,1100,adaptive::pt_raw(x[0], x[1], x[2], x[3]))  
 template<class Type>
  Type pt_(Type x, Type n, int lower_tail, int log_p){
    vector<Type> tx(5);
    tx[0] = x;
    tx[1] = n;
    // tx[2] = sigma;
    tx[2] = (Type)lower_tail;
    tx[3] = (Type)log_p;
    tx[4] = 0; // extra argument for derivative order
    Type res = pt_at(CppAD::vector<Type>(tx))[0];
    return res;
  }

// [[Rcpp::export]]
Rcpp::ComplexVector pt_ad(Rcpp::ComplexVector x, Rcpp::ComplexVector df, bool lower_tail, bool log_p ) {
  CHECK_INPUT(x); CHECK_INPUT(df);
  size_t n1 = x.size();
  size_t n2 = df.size();
  int nmax = std::max({n1, n2});
  int nmin = std::min({n1, n2});
  int n = (nmin == 0 ? 0 : nmax);
  Rcpp::ComplexVector ans(n);
  const ad* X1 = adptr(x); const ad* X2 = adptr(df);
  ad* Y = adptr(ans);
  for (int i=0; i<n; i++) Y[i] = pt_(X1[i % n1], X2[i % n2], lower_tail, log_p);
  return as_advector(ans);
}

TMB_BIND_ATOMIC(ps0,1,adaptive::ps_raw(x[0]) )
template<class Type>
Type ps_(Type x) {  
  vector<Type> tx(2);
  tx[0] = x;
  tx[1] = 0; // order
  return ps0(CppAD::vector<Type>(tx))[0];
}

// [[Rcpp::export]]
Rcpp::ComplexVector pde_scheme_ad(Rcpp::ComplexVector x) {
  CHECK_INPUT(x);
  int n = x.size();
  Rcpp::ComplexVector ans(n);
  const ad* X1 = adptr(x);
  ad* Y = adptr(ans);
  for (int i=0; i<n; i++) Y[i] = ps_(X1[i]);
  return as_advector(ans);
}


TMB_BIND_ATOMIC(logspace_add2x,
		11,
		adaptive::logspace_add2_raw(x[0], x[1]) )

template<class Type>
Type logspace_add2(Type logx, Type logy) {
  if ( !CppAD::Variable(logx) && logx == Type(R_NegInf) && !CppAD::Variable(logy) && logy == Type(R_NegInf))
    return Type(R_NegInf);
  if ( !CppAD::Variable(logx) && logx == Type(R_NegInf) )
    return logy;
  if ( !CppAD::Variable(logy) && logy == Type(R_NegInf) )
    return logx;
  vector<Type> tx(3);
  tx[0] = logx;
  tx[1] = logy;
  tx[2] = 0; // order
  return logspace_add2x(CppAD::vector<Type>(tx))[0];
}

// [[Rcpp::export]]
Rcpp::ComplexVector logspace_add_ad(Rcpp::ComplexVector x, Rcpp::ComplexVector y) {
  CHECK_INPUT(x); CHECK_INPUT(y);
  size_t n1 = x.size();
  size_t n2 = y.size();
  int nmax = std::max({n1, n2});
  int nmin = std::min({n1, n2});
  int n = (nmin == 0 ? 0 : nmax);
  Rcpp::ComplexVector ans(n);
  const ad* X1 = adptr(x); const ad* X2 = adptr(y);
  ad* Y = adptr(ans);
  for (int i=0; i<n; i++) Y[i] = logspace_add2(X1[i % n1], X2[i % n2]);
  return as_advector(ans);
}

// [[Rcpp::export]]
Rcpp::ComplexVector logspace_sum_ad(Rcpp::ComplexVector x) {
  CHECK_INPUT(x);
  int n = x.size();
  Rcpp::ComplexVector ans(1);  
  const ad* X1 = adptr(x);
  ad* Y = adptr(ans);
  Y[0] = X1[0];
  for (int i=1; i<n; i++) Y[0] = logspace_add2(X1[i], Y[0]);
  return as_advector(ans);
}



TMB_BIND_ATOMIC(logspace_sub2x,
		11,
		adaptive::logspace_sub2_raw(x[0], x[1]) )  
template<class Type>
Type logspace_sub2(Type logx, Type logy) {
  if ( !CppAD::Variable(logy) && logy == Type(R_NegInf) )
    return logx;
  vector<Type> tx(3);
  tx[0] = logx;
  tx[1] = logy;
  tx[2] = 0; // order
  return logspace_sub2x(CppAD::vector<Type>(tx))[0];
}
// [[Rcpp::export]]
Rcpp::ComplexVector logspace_sub_ad(Rcpp::ComplexVector x, Rcpp::ComplexVector y) {
  CHECK_INPUT(x); CHECK_INPUT(y);
  size_t n1 = x.size();
  size_t n2 = y.size();
  int nmax = std::max({n1, n2});
  int nmin = std::min({n1, n2});
  int n = (nmin == 0 ? 0 : nmax);
  Rcpp::ComplexVector ans(n);
  const ad* X1 = adptr(x); const ad* X2 = adptr(y);
  ad* Y = adptr(ans);
  for (int i=0; i<n; i++) Y[i] = logspace_sub2(X1[i % n1], X2[i % n2]);
  return as_advector(ans);
}


 TMB_BIND_ATOMIC(quantreg_x,11,adaptive::quantreg_raw(x[0], x[1]) )
  template<class Type>
  Type quantreg_(Type x, Type tau) {  
    vector<Type> tx(3);
    tx[0] = x;
    tx[1] = tau;
    tx[2] = 0; // order
    return quantreg_x(CppAD::vector<Type>(tx))[0];
  }
// [[Rcpp::export]]
Rcpp::ComplexVector quantreg_ad(Rcpp::ComplexVector x, Rcpp::ComplexVector tau) {
  CHECK_INPUT(x); CHECK_INPUT(tau);
  size_t n1 = x.size();
  size_t n2 = tau.size();
  int nmax = std::max({n1, n2});
  int nmin = std::min({n1, n2});
  int n = (nmin == 0 ? 0 : nmax);
  Rcpp::ComplexVector ans(n);
  const ad* X1 = adptr(x); const ad* X2 = adptr(tau);
  ad* Y = adptr(ans);
  for (int i=0; i<n; i++) Y[i] = quantreg_(X1[i % n1], X2[i % n2]);
  return as_advector(ans);
}


TMB_BIND_ATOMIC(login_log_besselIx,
		11,
		adaptive::login_log_besselI_raw(x[0], x[1]) )
template<class Type>
Type login_log_besselI_(Type logx, Type nu) {  
  vector<Type> tx(3);
  tx[0] = logx;
  tx[1] = nu;
  tx[2] = 0; // order
  return login_log_besselIx(CppAD::vector<Type>(tx))[0];
}
// [[Rcpp::export]]
Rcpp::ComplexVector login_log_besselI_ad(Rcpp::ComplexVector logx, Rcpp::ComplexVector nu) {
  CHECK_INPUT(logx); CHECK_INPUT(nu);
  size_t n1 = logx.size();
  size_t n2 = nu.size();
  int nmax = std::max({n1, n2});
  int nmin = std::min({n1, n2});
  int n = (nmin == 0 ? 0 : nmax);
  Rcpp::ComplexVector ans(n);
  const ad* X1 = adptr(logx); const ad* X2 = adptr(nu);
  ad* Y = adptr(ans);
  for (int i=0; i<n; i++) Y[i] = login_log_besselI_(X1[i % n1], X2[i % n2]);
  return as_advector(ans);
}

TMB_BIND_ATOMIC(log_besselIx,
		11,
		adaptive::log_besselI_raw(x[0], x[1]) )
template<class Type>
Type log_besselI_(Type x, Type nu) {  
  vector<Type> tx(3);
  tx[0] = x;
  tx[1] = nu;
  tx[2] = 0; // order
  return log_besselIx(CppAD::vector<Type>(tx))[0];
}
// [[Rcpp::export]]
Rcpp::ComplexVector log_besselI_ad(Rcpp::ComplexVector x, Rcpp::ComplexVector nu) {
  CHECK_INPUT(x); CHECK_INPUT(nu);
  size_t n1 = x.size();
  size_t n2 = nu.size();
  int nmax = std::max({n1, n2});
  int nmin = std::min({n1, n2});
  int n = (nmin == 0 ? 0 : nmax);
  Rcpp::ComplexVector ans(n);
  const ad* X1 = adptr(x); const ad* X2 = adptr(nu);
  ad* Y = adptr(ans);
  for (int i=0; i<n; i++) Y[i] = log_besselI_(X1[i % n1], X2[i % n2]);
  return as_advector(ans);
}


TMB_BIND_ATOMIC(log_MarcumQx,
		110,
		adaptive::log_MarcumQ_raw(x[0], x[1], x[2]) )
template<class Type>
Type log_MarcumQ_(Type a, Type b, Type nu) {  
  vector<Type> tx(4);
  tx[0] = a;
  tx[1] = b;
  tx[2] = nu;
  tx[3] = 0; // order
  return log_MarcumQx(CppAD::vector<Type>(tx))[0];
}
// [[Rcpp::export]]
Rcpp::ComplexVector log_MarcumQ_ad(Rcpp::ComplexVector a, Rcpp::ComplexVector b, Rcpp::ComplexVector nu) {
  CHECK_INPUT(a); CHECK_INPUT(b); CHECK_INPUT(nu);
  size_t n1 = a.size();
  size_t n2 = b.size();
  size_t n3 = nu.size();
  int nmax = std::max({n1, n2, n3});
  int nmin = std::min({n1, n2, n3});
  int n = (nmin == 0 ? 0 : nmax);
  Rcpp::ComplexVector ans(n);
  const ad* X1 = adptr(a); const ad* X2 = adptr(b); const ad* X3 = adptr(nu);
  ad* Y = adptr(ans);
  for (int i=0; i<n; i++) Y[i] = log_MarcumQ_(X1[i % n1], X2[i % n2], X3[i % n3]);
  return as_advector(ans);
}


TMB_BIND_ATOMIC(login_log_MarcumQx,
		110,
		adaptive::login_log_MarcumQ_raw(x[0], x[1], x[2]) )
template<class Type>
Type login_log_MarcumQ_(Type loga, Type logb, Type nu) {  
  vector<Type> tx(4);
  tx[0] = loga;
  tx[1] = logb;
  tx[2] = nu;
  tx[3] = 0; // order
  return login_log_MarcumQx(CppAD::vector<Type>(tx))[0];
}
// [[Rcpp::export]]
Rcpp::ComplexVector login_log_MarcumQ_ad(Rcpp::ComplexVector loga, Rcpp::ComplexVector logb, Rcpp::ComplexVector nu) {
  CHECK_INPUT(loga); CHECK_INPUT(logb); CHECK_INPUT(nu);
  size_t n1 = loga.size();
  size_t n2 = logb.size();
  size_t n3 = nu.size();
  int nmax = std::max({n1, n2, n3});
  int nmin = std::min({n1, n2, n3});
  int n = (nmin == 0 ? 0 : nmax);
  Rcpp::ComplexVector ans(n);
  const ad* X1 = adptr(loga); const ad* X2 = adptr(logb); const ad* X3 = adptr(nu);
  ad* Y = adptr(ans);
  for (int i=0; i<n; i++) Y[i] = login_log_MarcumQ_(X1[i % n1], X2[i % n2], X3[i % n3]);
  return as_advector(ans);
}

TMB_BIND_ATOMIC(log_Marcum1mQx,
		110,
		adaptive::log_Marcum1mQ_raw(x[0], x[1], x[2]) )
template<class Type>
Type log_Marcum1mQ_(Type a, Type b, Type nu) {  
  vector<Type> tx(4);
  tx[0] = a;
  tx[1] = b;
  tx[2] = nu;
  tx[3] = 0; // order
  return log_Marcum1mQx(CppAD::vector<Type>(tx))[0];
}
// [[Rcpp::export]]
Rcpp::ComplexVector log_Marcum1mQ_ad(Rcpp::ComplexVector a, Rcpp::ComplexVector b, Rcpp::ComplexVector nu) {
  CHECK_INPUT(a); CHECK_INPUT(b); CHECK_INPUT(nu);
  size_t n1 = a.size();
  size_t n2 = b.size();
  size_t n3 = nu.size();
  int nmax = std::max({n1, n2, n3});
  int nmin = std::min({n1, n2, n3});
  int n = (nmin == 0 ? 0 : nmax);
  Rcpp::ComplexVector ans(n);
  const ad* X1 = adptr(a); const ad* X2 = adptr(b); const ad* X3 = adptr(nu);
  ad* Y = adptr(ans);
  for (int i=0; i<n; i++) Y[i] = log_Marcum1mQ_(X1[i % n1], X2[i % n2], X3[i % n3]);
  return as_advector(ans);
}


TMB_BIND_ATOMIC(login_log_Marcum1mQx,
		110,
		adaptive::login_log_Marcum1mQ_raw(x[0], x[1], x[2]) )
template<class Type>
Type login_log_Marcum1mQ_(Type loga, Type logb, Type nu) {  
  vector<Type> tx(4);
  tx[0] = loga;
  tx[1] = logb;
  tx[2] = nu;
  tx[3] = 0; // order
  return login_log_Marcum1mQx(CppAD::vector<Type>(tx))[0];
}
// [[Rcpp::export]]
Rcpp::ComplexVector login_log_Marcum1mQ_ad(Rcpp::ComplexVector loga, Rcpp::ComplexVector logb, Rcpp::ComplexVector nu) {
  CHECK_INPUT(loga); CHECK_INPUT(logb); CHECK_INPUT(nu);
  size_t n1 = loga.size();
  size_t n2 = logb.size();
  size_t n3 = nu.size();
  int nmax = std::max({n1, n2, n3});
  int nmin = std::min({n1, n2, n3});
  int n = (nmin == 0 ? 0 : nmax);
  Rcpp::ComplexVector ans(n);
  const ad* X1 = adptr(loga); const ad* X2 = adptr(logb); const ad* X3 = adptr(nu);
  ad* Y = adptr(ans);
  for (int i=0; i<n; i++) Y[i] = login_log_Marcum1mQ_(X1[i % n1], X2[i % n2], X3[i % n3]);
  return as_advector(ans);
}


TMB_BIND_ATOMIC(lodricex,
		111,
		adaptive::logdrice_raw(x[0], x[1], x[2]) )
template<class Type>
Type logdrice_(Type logx, Type lognu, Type logsigma) {  
  vector<Type> tx(4);
  tx[0] = logx;
  tx[1] = lognu;
  tx[2] = logsigma;
  tx[3] = 0; // order
  return lodricex(CppAD::vector<Type>(tx))[0];
}
// [[Rcpp::export]]
Rcpp::ComplexVector logdrice_ad(Rcpp::ComplexVector logx, Rcpp::ComplexVector lognu, Rcpp::ComplexVector logsigma) {
  CHECK_INPUT(logx); CHECK_INPUT(lognu); CHECK_INPUT(logsigma);
  size_t n1 = logx.size();
  size_t n2 = lognu.size();
  size_t n3 = logsigma.size();
  int nmax = std::max({n1, n2, n3});
  int nmin = std::min({n1, n2, n3});
  int n = (nmin == 0 ? 0 : nmax);
  Rcpp::ComplexVector ans(n);
  const ad* X1 = adptr(logx); const ad* X2 = adptr(lognu); const ad* X3 = adptr(logsigma);
  ad* Y = adptr(ans);
  for (int i=0; i<n; i++) Y[i] = logdrice_(X1[i % n1], X2[i % n2], X3[i % n3]);
  return as_advector(ans);
}

TMB_BIND_ATOMIC(lopricex,
		1110,
		adaptive::logprice_raw(x[0], x[1], x[2], x[3]) )
template<class Type>
Type logprice_(Type logx, Type lognu, Type logsigma, int lower_tail) {  
  vector<Type> tx(5);
  tx[0] = logx;
  tx[1] = lognu;
  tx[2] = logsigma;
  tx[3] = (Type)lower_tail;
  tx[4] = 0; // order
  return lopricex(CppAD::vector<Type>(tx))[0];
}
// [[Rcpp::export]]
Rcpp::ComplexVector logprice_ad(Rcpp::ComplexVector logx, Rcpp::ComplexVector lognu, Rcpp::ComplexVector logsigma, bool lower_tail) {
  CHECK_INPUT(logx); CHECK_INPUT(lognu); CHECK_INPUT(logsigma);
  size_t n1 = logx.size();
  size_t n2 = lognu.size();
  size_t n3 = logsigma.size();
  int nmax = std::max({n1, n2, n3});
  int nmin = std::min({n1, n2, n3});
  int n = (nmin == 0 ? 0 : nmax);
  Rcpp::ComplexVector ans(n);
  const ad* X1 = adptr(logx); const ad* X2 = adptr(lognu); const ad* X3 = adptr(logsigma);
  ad* Y = adptr(ans);
  for (int i=0; i<n; i++) Y[i] = logprice_(X1[i % n1], X2[i % n2], X3[i % n3], lower_tail);
  return as_advector(ans);
}
