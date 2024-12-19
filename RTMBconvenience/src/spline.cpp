// [[Rcpp::depends(TMB)]]
#include <RTMB.h>

namespace fromR {
    template<class Float>
    Float pnorm5_raw(Float x, Float mu, Float sigma, Float lower_tail, Float log_p);
}

namespace spline_atomic {

  template<class Float>
  struct pkwnorm1_t {
    typedef Float Scalar; // Required by integrate
    Float gam;         // Parameters
    Float x0;
    // Evaluate joint density of (u, x)
    Float operator() (Float u) {
      Float a = 1.0 + 0.5 * (gam + sqrt(gam * gam));
      Float b = 1.0 - 0.5 * (gam - sqrt(gam * gam));
      Float lpv = fromR::pnorm5_raw(u,Float(0.0), Float(1.0),Float(1.0), Float(1.0)); //log(pnorm_approx(u)); //
      Float lGa = a * lpv;
      return 1.0 - exp(b * atomic::robust_utils::logspace_sub(Float(0.0),lGa));
    }
    // Integrate latent variable (u) out
    Float integrate(Float x) {
      using gauss_kronrod::integrate;
      Float ans = integrate(*this, x0, x);
      return ans;
    }
  };

   template<class Float>
   Float eval_ipkwnorm1(Float x, Float gam, Float x0) {
     pkwnorm1_t<Float> f = {gam, x0};
    return f.integrate(x);
  }

  TMB_BIND_ATOMIC(fun_ipkwnorm1, 110, eval_ipkwnorm1(x[0], x[1], x[2]))

}


template<class Type>
Type spline_ipkwnorm(Type x, Type mu, Type sd, Type gam, double x0) {
  vector<Type> args(4); // Last index reserved for derivative order
  args << (x - mu) / sd, gam, (Type)(x0), 0;
  return spline_atomic::fun_ipkwnorm1(CppAD::vector<Type>(args))[0] * sd;
};

// Implemented in pnorm.cpp
template<class Type>
Type pnorm5_(Type x, Type mu, Type sigma, int lower_tail, int log_p);

template<class Type>
Type dkwnorm_(Type x, Type mu, Type sig, Type a, Type b, bool give_log){
  Type lpv = pnorm5_(x,mu,sig,1, 1); // log(pnorm(x,mu,sig));
  Type lGa = a * lpv;
  Type log_res = log(a) + log(b) + dnorm(x,mu,sig, true) + (a-1.0) * lpv + (b-1.0) * logspace_sub(Type(0.0),(Type)lGa);
  if(give_log)
    return log_res;
  return exp(log_res);
}
// [[Rcpp::export]]
ADrep distr_dkwnorm(ADrep x, ADrep mu, ADrep sigma, ADrep a, ADrep b, bool give_log ) {
  size_t n1 = x.size();
  size_t n2 = mu.size();
  size_t n3 = sigma.size();
  size_t n4 = sigma.size();
  size_t n5 = sigma.size();
  int nmax = std::max({n1, n2, n3, n4, n5});
  int nmin = std::min({n1, n2, n3, n4, n5});
  int n = (nmin == 0 ? 0 : nmax);
  ADrep ans(n);
  const ad* X1 = adptr(x); const ad* X2 = adptr(mu); const ad* X3 = adptr(sigma);
  const ad* X4 = adptr(a); const ad* X5 = adptr(b);
  ad* Y = adptr(ans);
  for (int i=0; i<n; i++) Y[i] = dkwnorm_(X1[i % n1], X2[i % n2], X3[i % n3], X4[i % n4], X5[i % n5], give_log);
  return ans;
}

// Kumaraswamy-normal (Kw-normal) density function with special choice of a and b
template<class Type>
Type spline_dkwnorm(Type x, Type mu, Type sig, Type gam, bool give_log){
  Type a = 1.0 + 0.5 * (gam + sqrt(gam * gam));
  Type b = 1.0 - 0.5 * (gam - sqrt(gam * gam));  
  return dkwnorm_(x,mu,sig,a,b,give_log);
}

// Kumaraswamy-normal (Kw-normal) distribution function with special choice of a and b
template<class Type>
Type pkwnorm_(Type x, Type mu, Type sig, Type a, Type b, bool lower_tail, bool give_log){
  Type lpv = pnorm5_(x,mu,sig,1, 1); // log(pnorm(x,mu,sig)); //
  Type lGa = a * lpv;
  if(lower_tail){
    if(give_log)
      return logspace_sub(Type(0.0), (Type)(b * logspace_sub(Type(0.0),lGa)));
    return 1.0 - exp(b * logspace_sub(Type(0.0),lGa));
  }
  if(give_log)
    return (Type)(b * logspace_sub(Type(0.0),lGa));
  return exp(b * logspace_sub(Type(0.0),lGa));
};

// [[Rcpp::export]]
ADrep distr_pkwnorm(ADrep x, ADrep mu, ADrep sigma, ADrep a, ADrep b, bool lower_tail,bool give_log ) {
  size_t n1 = x.size();
  size_t n2 = mu.size();
  size_t n3 = sigma.size();
  size_t n4 = sigma.size();
  size_t n5 = sigma.size();
  int nmax = std::max({n1, n2, n3, n4, n5});
  int nmin = std::min({n1, n2, n3, n4, n5});
  int n = (nmin == 0 ? 0 : nmax);
  ADrep ans(n);
  const ad* X1 = adptr(x); const ad* X2 = adptr(mu); const ad* X3 = adptr(sigma);
  const ad* X4 = adptr(a); const ad* X5 = adptr(b);
  ad* Y = adptr(ans);
  for (int i=0; i<n; i++) Y[i] = pkwnorm_(X1[i % n1], X2[i % n2], X3[i % n3], X4[i % n4], X5[i % n5], lower_tail, give_log);
  return ans;
}


// Kumaraswamy-normal (Kw-normal) density function with special choice of a and b
template<class Type>
Type spline_pkwnorm(Type x, Type mu, Type sig, Type gam){
  Type a = 1.0 + 0.5 * (gam + sqrt(gam * gam));
  Type b = 1.0 - 0.5 * (gam - sqrt(gam * gam));  
  return pkwnorm_(x,mu,sig,a,b,1,0);
}

  
template<class Type>
matrix<Type> spline_getSigAndGam(vector<Type> knots){
  // if(CppAD::Variable(knots(0)))
  //   Rf_error("Knots can not be parameters");
  if(knots.size() < 3)
    Rf_error("The spline must have at least three knots.");
  matrix<Type> res(knots.size(),2); // Sigma, Gamma
  res.setZero();
  for(int i = 1; i < knots.size() - 1; ++i){
    // Sigma
    if(fabs(knots(i) - knots(i-1)) < 0.001){
      res(i,0) = 0.001 / 3.0;
    }else{
      res(i,0) = (knots(i+1) - knots(i-1)) / 3.0;
    }
    // Gamma
    if(fabs(knots(i) - knots(i-1)) < 1e-8){
      res(i,1) = res(i,0);
    }else if(fabs(knots(i) - knots(i+1)) < 1e-8){
      res(i,1) = -res(i,0);
    }else{
      res(i,1) = log(knots(i+1) - knots(i)) - log(knots(i) - knots(i-1));
    }
  }
  res(0,0) = (knots(2) - knots(0)) / 3.0;
  res(knots.size()-1,0) = (knots(knots.size()-1) - knots(knots.size()-1-2)) / 3.0;
  if(fabs(knots(1) - knots(0)) < 1e-8){
    res(0,1) = 0.0;
  }else{
    res(0,1) = res(0,0);
  }
  if(fabs(knots(knots.size()-1) - knots(knots.size()-1-1)) < 1e-8){
    res(knots.size()-1,1) = 0.0;
  }else{
    res(knots.size()-1,1) = -res(knots.size()-1,0);
  }
  return res;
};


// Spline using kw-normals as basis functions
// [[Rcpp::export]]
ADrep spline_bcspline(ADrep x, Rcpp::NumericVector knots, ADrep pars) {
  if((int)knots.size() != (int)pars.size())
    Rf_error("Pars must have the same number of elements as knots");
  size_t n1 = x.size();
  int n = n1;
  ADrep ans(n);
  const ad* X1 = adptr(x);
  ad* Y = adptr(ans);
  const ad* P = adptr(pars);
  // Do calculations here to save some time
  matrix<double> sg = spline_getSigAndGam(asVector<double>(knots));
  for (int i=0; i<n; i++){
    ad tmp = 0.0;
    for(int j = 0; j < (int)knots.size(); ++j)
      tmp += P[j] * spline_dkwnorm(X1[i], (ad)knots[j], (ad)sg(j,0), (ad)sg(j,1), false) * sg(j,0);
    Y[i] = tmp;
  }
  return ans;
}

// Integrated spline using kw-normals as basis functions
// [[Rcpp::export]]
ADrep spline_ibcspline(ADrep x, Rcpp::NumericVector knots, ADrep pars) {
  if((int)knots.size() + 1 != (int)pars.size())
     Rf_error("Pars must have one more element than knots");
  size_t n1 = x.size();
  int n = n1;
  ADrep ans(n);
  const ad* X1 = adptr(x);
  ad* Y = adptr(ans);
  const ad* P = adptr(pars);
  // Do calculations here to save some time
  matrix<double> sg = spline_getSigAndGam(asVector<double>(knots));
  for (int i=0; i<n; i++){
    ad res = 0.0;
    for(int j = 0; j < (int)knots.size(); ++j){
      double v0 = spline_pkwnorm(knots[0], knots[j], sg(j,0), sg(j,1));
      ad tmp = spline_pkwnorm(X1[i], (ad)knots[j], (ad)sg(j,0), (ad)sg(j,1));
      res += P[j] * (tmp - (ad)v0);
    }
    Y[i] = res + P[pars.size()-1];
  }
  return ans;
}


// Monotonically non-increasing spline using kw-normals as basis functions
// [[Rcpp::export]]
ADrep spline_ibcdspline(ADrep x, Rcpp::NumericVector knots, ADrep pars) {
  if((int)knots.size() + 1 != (int)pars.size())
     Rf_error("Pars must have one more element than knots");
  size_t n1 = x.size();
  int n = n1;
  ADrep ans(n);
  const ad* X1 = adptr(x);
  ad* Y = adptr(ans);
  const ad* P = adptr(pars);
  ADrep pars2(pars.size()-1);
  ad* P2 = adptr(pars2);
  for(int i = 0; i < (int)pars2.size(); ++i){
    P2[i] = -exp(P[i]);
  }
  // Do calculations here to save some time
  matrix<double> sg = spline_getSigAndGam(asVector<double>(knots));
  for (int i=0; i<n; i++){
    ad res = 0.0;
    for(int j = 0; j < (int)knots.size(); ++j){
      double v0 = spline_pkwnorm(knots[0], knots[j], sg(j,0), sg(j,1));
      ad tmp = spline_pkwnorm(X1[i], (ad)knots[j], (ad)sg(j,0), (ad)sg(j,1));
      res += P2[j] * (tmp - (ad)v0);
    }
    Y[i] = res + P[pars.size()-1];
  }
  return ans;
}


// Monotonically non-decreasing spline using kw-normals as basis functions
// [[Rcpp::export]]
ADrep spline_ibcispline(ADrep x, Rcpp::NumericVector knots, ADrep pars) {
  if((int)knots.size() + 1 != (int)pars.size())
     Rf_error("Pars must have one more element than knots");
  size_t n1 = x.size();
  int n = n1;
  ADrep ans(n);
  const ad* X1 = adptr(x);
  ad* Y = adptr(ans);
  const ad* P = adptr(pars);
  ADrep pars2(pars.size()-1);
  ad* P2 = adptr(pars2);
  for(int i = 0; i < (int)pars2.size(); ++i){
    P2[i] = exp(P[i]);
  }
  // Do calculations here to save some time
  matrix<double> sg = spline_getSigAndGam(asVector<double>(knots));
  for (int i=0; i<n; i++){
    ad res = 0.0;
    for(int j = 0; j < (int)knots.size(); ++j){
      double v0 = spline_pkwnorm(knots[0], knots[j], sg(j,0), sg(j,1));
      ad tmp = spline_pkwnorm(X1[i], (ad)knots[j], (ad)sg(j,0), (ad)sg(j,1));
      res += P2[j] * (tmp - (ad)v0);
    }
    Y[i] = res + P[pars.size()-1];
  }
  return ans;
}



// Integrated integrated spline using kw-normals as basis functions
// [[Rcpp::export]]
ADrep spline_iibcspline(ADrep x, Rcpp::NumericVector knots, ADrep pars) {
  if((int)knots.size() + 1!= (int)pars.size())
     Rf_error("Pars must have one more element than knots");
  size_t n1 = x.size();
  int n = n1;
  ADrep ans(n);
  const ad* X1 = adptr(x);
  ad* Y = adptr(ans);
  const ad* P = adptr(pars);
  // Do calculations here to save some time
  matrix<double> sg = spline_getSigAndGam(asVector<double>(knots));
  for (int i=0; i<n; i++){
    ad res = 0.0;
    for(int j = 0; j < (int)knots.size(); ++j){
      // double v0 = spline_ipkwnorm(knots[0], knots[j], sg(j,0), sg(j,1));
      ad tmp = spline_ipkwnorm(X1[i], (ad)knots[j], (ad)sg(j,0), (ad)sg(j,1), knots[0]);
      res += P[j] * tmp; //(tmp - (ad)v0 * X1[i]);
    }
    Y[i] = res + P[pars.size()-1];
  }
  return ans;
}


// Monotonically non-increasing spline using kw-normals as basis functions
// [[Rcpp::export]]
ADrep spline_iibcdspline(ADrep x, Rcpp::NumericVector knots, ADrep pars) {
  if((int)knots.size() + 1 != (int)pars.size())
     Rf_error("Pars must have one more element than knots");
  size_t n1 = x.size();
  int n = n1;
  ADrep ans(n);
  const ad* X1 = adptr(x);
  ad* Y = adptr(ans);
  const ad* P = adptr(pars);
  ADrep pars2(pars.size()-1);
  ad* P2 = adptr(pars2);
  ad p2s = 0;
  for(int i = 0; i < (int)pars2.size(); ++i){
    P2[i] = -exp(P[i]);
    p2s += P2[i];
  }
  // Do calculations here to save some time
  matrix<double> sg = spline_getSigAndGam(asVector<double>(knots));
  for (int i=0; i<n; i++){
    ad res = 0.0;
    for(int j = 0; j < (int)knots.size(); ++j){
      //double v0 = spline_pkwnorm(knots[0], knots[j], sg(j,0), sg(j,1));
      ad tmp = spline_ipkwnorm(X1[i], (ad)knots[j], (ad)sg(j,0), (ad)sg(j,1), knots[0]);
      res += P2[j] * tmp; //(tmp - (ad)v0 * X1[i]);
    }
    Y[i] = res  + P[pars.size()-1];// - p2s * X1[i];
  }
  return ans;
}


// Monotonically non-decreasing spline using kw-normals as basis functions
// [[Rcpp::export]]
ADrep spline_iibcispline(ADrep x, Rcpp::NumericVector knots, ADrep pars) {
  if((int)knots.size() + 1 != (int)pars.size())
     Rf_error("Pars must have one more element than knots");
  size_t n1 = x.size();
  int n = n1;
  ADrep ans(n);
  const ad* X1 = adptr(x);
  ad* Y = adptr(ans);
  const ad* P = adptr(pars);
  ADrep pars2(pars.size()-1);
  ad* P2 = adptr(pars2);
  ad p2s = 0;
  for(int i = 0; i < (int)pars2.size(); ++i){
    P2[i] = -exp(P[i]);
    p2s += P2[i];
  }
  // Do calculations here to save some time
  matrix<double> sg = spline_getSigAndGam(asVector<double>(knots));
  for (int i=0; i<n; i++){
    ad res = 0.0;
    for(int j = 0; j < (int)knots.size(); ++j){
      //double v0 = spline_pkwnorm(knots[0], knots[j], sg(j,0), sg(j,1));
      ad tmp = spline_ipkwnorm(X1[i], (ad)knots[j], (ad)sg(j,0), (ad)sg(j,1), knots[0]);
      res += P2[j] * tmp; //(tmp - (ad)v0 * X1[i]);
    }
    Y[i] = res  + P[pars.size()-1];// - p2s * X1[i];
  }
  return ans;
}
