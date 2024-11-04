  /*
   *  Mathlib : A C Library of Special Functions
   *  Copyright (C) 1998	    Ross Ihaka
   *  Copyright (C) 2000-2013 The R Core Team
   *  Copyright (C) 2003	    The R Foundation
   *
   *  This program is free software; you can redistribute it and/or modify
   *  it under the terms of the GNU General Public License as published by
   *  the Free Software Foundation; either version 2 of the License, or
   *  (at your option) any later version.
   *
   *  This program is distributed in the hope that it will be useful,
   *  but WITHOUT ANY WARRANTY; without even the implied warranty of
   *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   *  GNU General Public License for more details.
   *
   *  You should have received a copy of the GNU General Public License
   *  along with this program; if not, a copy is available at
   *  https://www.R-project.org/Licenses/
   *
   *  SYNOPSIS
   *
   *   #include <Rmath.h>
   *
   *   double pnorm5(double x, double mu, double sigma, int lower_tail,int log_p);
   *	   {pnorm (..) is synonymous and preferred inside R}
   *
   *   void   pnorm_both(double x, double *cum, double *ccum,
   *		       int i_tail, int log_p);
   *
   *  DESCRIPTION
   *
   *	The main computation evaluates near-minimax approximations derived
   *	from those in "Rational Chebyshev approximations for the error
   *	function" by W. J. Cody, Math. Comp., 1969, 631-637.  This
   *	transportable program uses rational functions that theoretically
   *	approximate the normal distribution function to at least 18
   *	significant decimal digits.  The accuracy achieved depends on the
   *	arithmetic system, the compiler, the intrinsic functions, and
   *	proper selection of the machine-dependent constants.
   *
   *  REFERENCE
   *
   *	Cody, W. D. (1993).
   *	ALGORITHM 715: SPECFUN - A Portable FORTRAN Package of
   *	Special Function Routines and Test Drivers".
   *	ACM Transactions on Mathematical Software. 19, 22-32.
   *
   *  EXTENSIONS
   *
   *  The "_both" , lower, upper, and log_p  variants were added by
   *  Martin Maechler, Jan.2000;
   *  as well as log1p() and similar improvements later on.
   *
   *  James M. Rath contributed bug report PR#699 and patches correcting SIXTEN
   *  and if() clauses {with a bug: "|| instead of &&" -> PR #2883) more in line
   *  with the original Cody code.
   */

  /*
   * Modified 2019 Christoffer Moesgaard Albertsen, Technical University of Denmark
   */




  namespace fromR {

#ifndef WITH_SAM_LIB

    template<class T> int R_finite(T x) { return std::isfinite(asDouble(x)); }
    template<class T> int isnan(T x) { return std::isnan(asDouble(x)); }

#undef ML_ERROR
#undef MATHLIB_ERROR
#undef MATHLIB_WARNING
#undef MATHLIB_WARNING2
#undef MATHLIB_WARNING3
#undef MATHLIB_WARNING4
#undef MATHLIB_WARNING5
#undef ML_POSINF
#undef ML_NEGINF
#undef ML_NAN
#undef M_SQRT_2dPI
#undef ISNAN
# define ML_ERROR(x, s) /* nothing */
#define ML_WARNING(x,s)
# define MATHLIB_ERROR(fmt,x) /* nothing */
# define MATHLIB_WARNING(fmt,x) /* nothing */
# define MATHLIB_WARNING2(fmt,x,x2) /* nothing */
# define MATHLIB_WARNING3(fmt,x,x2,x3) /* nothing */
# define MATHLIB_WARNING4(fmt,x,x2,x3,x4) /* nothing */
# define MATHLIB_WARNING5(fmt,x,x2,x3,x4,x5) /* nothing */
#define ML_POSINF	R_PosInf
#define ML_NEGINF	R_NegInf
#define ML_NAN		R_NaN
#define M_SQRT_2dPI	0.797884560802865355879892119869	/* sqrt(2/pi) */
#define M_LN_SQRT_2PI	0.918938533204672741780329736406	/* log(sqrt(2*pi)) */
#define M_LN_SQRT_PId2	0.225791352644727432363097614947	/* log(sqrt(pi/2))
								   == log(pi/2)/2 */
#define M_LN_2PI	1.837877066409345483560659472811	/* log(2*pi) */
#define ISNAN(x) (isnan(x)!=0)

#define ML_ERR_return_NAN return R_NaN
#define ML_WARN_return_NAN return R_NaN

#define give_log log_p
    
#define R_D__0	(log_p ? ML_NEGINF : 0.)		/* 0 */
#define R_D__1	(log_p ? 0. : 1.)			/* 1 */
#define R_DT_0	(lower_tail ? R_D__0 : R_D__1)		/* 0 */
#define R_DT_1	(lower_tail ? R_D__1 : R_D__0) /* 1 */
# define attribute_hidden __attribute__ ((visibility ("hidden")))

#define R_D_val(x)	(log_p	? log(x) : (x))		/*  x  in pF(x,..) */
#define R_D_Clog(p)	(log_p	? log1p(-(p)) : (0.5 - (p) + 0.5)) /* [log](1-p) */
#define R_DT_val(x)	(lower_tail ? R_D_val(x)  : R_D_Clog(x))
    
  
#define SIXTEN	16 /* Cutoff allowing exact "*" and "/" */

#define R_D_Cval(p) (lower_tail ? (0.5 - (p) + 0.5) : (p)) /*  1 - p */

#define M_SQRT_32	5.656854249492380195206754896838	/* sqrt(32) */
#define M_1_SQRT_2PI	0.398942280401432677939946059934	/* 1/sqrt(2pi) */

#define R_D_log(p)	(log_p	?  (p)	 : log(p))	/* log(p) */
#define R_D_LExp(x)     (log_p ? atomic::robust_utils::R_Log1_Exp(x) : log1p(-x))
#define R_D_exp(x)	(log_p	?  (x)	 : exp(x))	/* exp(x) */

#define R_P_bounds_01(x, x_min, x_max)	\
    if(x <= x_min) return R_DT_0;		\
    if(x >= x_max) return R_DT_1
    
    
#define R_Q_P01_boundaries(p, _LEFT_, _RIGHT_)		\
    if (log_p) {					\
	if(p > 0)					\
	    ML_WARN_return_NAN;				\
	if(p == 0) /* upper bound*/			\
	    return lower_tail ? _RIGHT_ : _LEFT_;	\
	if(p == ML_NEGINF)				\
	    return lower_tail ? _LEFT_ : _RIGHT_;	\
    }							\
    else { /* !log_p */					\
	if(p < 0 || p > 1)				\
	    ML_WARN_return_NAN;				\
	if(p == 0)					\
	    return lower_tail ? _LEFT_ : _RIGHT_;	\
	if(p == 1)					\
	    return lower_tail ? _RIGHT_ : _LEFT_;	\
    }

    
#define R_D_Lval(p)	(lower_tail ? (p) : (0.5 - (p) + 0.5))	/*  p  */
/*#define R_DT_qIv(p)	R_D_Lval(R_D_qIv(p))		 *  p  in qF ! */
#define R_DT_qIv(p)	(log_p ? (lower_tail ? exp(p) : - expm1(p)) \
			       : R_D_Lval(p))

/*#define R_DT_CIv(p)	R_D_Cval(R_D_qIv(p))		 *  1 - p in qF */
#define R_DT_CIv(p)	(log_p ? (lower_tail ? -expm1(p) : exp(p)) \
			       : R_D_Cval(p))




template<class Float>
static Float myfmod(Float x1, Float x2)
{
    Float q = x1 / x2;
    return x1 - floor(q) * x2;
}

template<class Float>
static Float myfmin2(Float x1, Float x2)
{
  return 0.5 * (x1 + x2) - 0.5 * fabs(x1 - x2);
}

    template<class Float>
static Float myfmax2(Float x1, Float x2)
{
  return 0.5 * (x1 + x2) + 0.5 * fabs(x1 - x2);
}

    
    template<class Float>
    attribute_hidden
    void pnorm_both_raw(Float x, Float *cum, Float *ccum, int i_tail, int log_p)
    {
      /* i_tail in {0,1,2} means: "lower", "upper", or "both" :
	 if(lower) return  *cum := P[X <= x]
	 if(upper) return *ccum := P[X >  x] = 1 - P[X <= x]
      */
      const static Float a[5] = {
	2.2352520354606839287,
	161.02823106855587881,
	1067.6894854603709582,
	18154.981253343561249,
	0.065682337918207449113
      };
      const static Float b[4] = {
	47.20258190468824187,
	976.09855173777669322,
	10260.932208618978205,
	45507.789335026729956
      };
      const static Float c[9] = {
	0.39894151208813466764,
	8.8831497943883759412,
	93.506656132177855979,
	597.27027639480026226,
	2494.5375852903726711,
	6848.1904505362823326,
	11602.651437647350124,
	9842.7148383839780218,
	1.0765576773720192317e-8
      };
      const static Float d[8] = {
	22.266688044328115691,
	235.38790178262499861,
	1519.377599407554805,
	6485.558298266760755,
	18615.571640885098091,
	34900.952721145977266,
	38912.003286093271411,
	19685.429676859990727
      };
      const static Float p[6] = {
	0.21589853405795699,
	0.1274011611602473639,
	0.022235277870649807,
	0.001421619193227893466,
	2.9112874951168792e-5,
	0.02307344176494017303
      };
      const static Float q[5] = {
	1.28426009614491121,
	0.468238212480865118,
	0.0659881378689285515,
	0.00378239633202758244,
	7.29751555083966205e-5
      };

      Float xden, xnum, temp, del, eps, xsq, y;
#ifdef NO_DENORMS
      Float min = DBL_MIN;
#endif
      int i, lower, upper;

#ifdef IEEE_754
      if(ISNAN(x)) { *cum = *ccum = x; return; } 
#endif

      /* Consider changing these : */
      eps = DBL_EPSILON * 0.5;

      /* i_tail in {0,1,2} =^= {lower, upper, both} */
      lower = i_tail != 1;
      upper = i_tail != 0;

      y = fabs(x);
      if (y <= 0.67448975) { /* qnorm(3/4) = .6744.... -- earlier had 0.66291 */
	if (y > eps) {
	  xsq = x * x;
	  xnum = a[4] * xsq;
	  xden = xsq;
	  for (i = 0; i < 3; ++i) {
	    xnum = (xnum + a[i]) * xsq;
	    xden = (xden + b[i]) * xsq;
	  }
	} else xnum = xden = 0.0;

	temp = x * (xnum + a[3]) / (xden + b[3]);
	if(lower)  *cum = 0.5 + temp;
	if(upper) *ccum = 0.5 - temp;
	if(log_p) {
	  if(lower)  *cum = log(*cum);
	  if(upper) *ccum = log(*ccum);
	}
      }
      else if (y <= M_SQRT_32) {

	/* Evaluate pnorm for 0.674.. = qnorm(3/4) < |x| <= sqrt(32) ~= 5.657 */

	xnum = c[8] * y;
	xden = y;
	for (i = 0; i < 7; ++i) {
	  xnum = (xnum + c[i]) * y;
	  xden = (xden + d[i]) * y;
	}
	temp = (xnum + c[7]) / (xden + d[7]);

#define do_del(X)							\
	xsq = trunc(X * SIXTEN) / SIXTEN;				\
	del = (X - xsq) * (X + xsq);					\
	if(log_p) {							\
	  *cum = (-xsq * xsq * 0.5) + (-del * 0.5) + log(temp);		\
	  if((lower && x > 0.) || (upper && x <= 0.))			\
	    *ccum = log1p(-exp(-xsq * xsq * 0.5) *			\
			  exp(-del * 0.5) * temp);			\
	}								\
	else {								\
	  *cum = exp(-xsq * xsq * 0.5) * exp(-del * 0.5) * temp;	\
	  *ccum = 1.0 - *cum;						\
	}

#define swap_tail						\
	if (x > 0.) {/* swap  ccum <--> cum */			\
	  temp = *cum; if(lower) *cum = *ccum; *ccum = temp;	\
	}

	do_del(y);
	swap_tail;
      }

      /* else	  |x| > sqrt(32) = 5.657 :
       * the next two case differentiations were really for lower=T, log=F
       * Particularly	 *not*	for  log_p !

       * Cody had (-37.5193 < x  &&  x < 8.2924) ; R originally had y < 50
       *
       * Note that we do want symmetry(0), lower/upper -> hence use y
       */
      else if((log_p && y < 1e170) /* avoid underflow below */
	      /*  ^^^^^ MM FIXME: can speedup for log_p and much larger |x| !
	       * Then, make use of  Abramowitz & Stegun, 26.2.13, something like

	       xsq = x*x;

	       if(xsq * DBL_EPSILON < 1.)
	       del = (1. - (1. - 5./(xsq+6.)) / (xsq+4.)) / (xsq+2.);
	       else
	       del = 0.;
	       *cum = -.5*xsq - M_LN_SQRT_2PI - log(x) + log1p(-del);
	       *ccum = log1p(-exp(*cum)); /.* ~ log(1) = 0 *./

	       swap_tail;

	       [Yes, but xsq might be infinite.]

	      */
	      || (lower && -37.5193 < x  &&  x < 8.2924)
	      || (upper && -8.2924  < x  &&  x < 37.5193)
	      ) {

	/* Evaluate pnorm for x in (-37.5, -5.657) union (5.657, 37.5) */
	xsq = 1.0 / (x * x); /* (1./x)*(1./x) might be better */
	xnum = p[5] * xsq;
	xden = xsq;
	for (i = 0; i < 4; ++i) {
	  xnum = (xnum + p[i]) * xsq;
	  xden = (xden + q[i]) * xsq;
	}
	temp = xsq * (xnum + p[4]) / (xden + q[4]);
	temp = (M_1_SQRT_2PI - temp) / y;

	do_del(x);
	swap_tail;
      } else { /* large x such that probs are 0 or 1 */
	if(x > 0) {	*cum = R_D__1; *ccum = R_D__0;	}
	else {	        *cum = R_D__0; *ccum = R_D__1;	}
      }


#ifdef NO_DENORMS
      /* do not return "denormalized" -- we do in R */
      if(log_p) {
	if(*cum > -min)	 *cum = -0.;
	if(*ccum > -min)*ccum = -0.;
      }
      else {
	if(*cum < min)	 *cum = 0.;
	if(*ccum < min)	*ccum = 0.;
      }
#endif
      return;
    }


    template<class Float>
    attribute_hidden
    Float pnorm5_raw(Float x, Float mu, Float sigma, Float lower_tail_, Float log_p_)
    {
      int lower_tail = (int)trunc(lower_tail_);
      int log_p = (int)trunc(log_p_);
      Float p, cp;

      /* Note: The structure of these checks has been carefully thought through.
       * For example, if x == mu and sigma == 0, we get the correct answer 1.
       */
#ifdef IEEE_754
      if(ISNAN(x) || ISNAN(mu) || ISNAN(sigma))
	return x + mu + sigma;
#endif
      if(!R_FINITE(x) && mu == x) return ML_NAN;/* x-mu is NaN */
      if (sigma <= 0) {
	if(sigma < 0) ML_ERR_return_NAN;
	/* sigma = 0 : */
	return (x < mu) ? R_DT_0 : R_DT_1;
      }
      p = (x - mu) / sigma;
      // if(!R_FINITE(p))
      //   return (x < mu) ? R_DT_0 : R_DT_1;
      x = p;

      pnorm_both_raw(x, &p, &cp, (lower_tail ? 0 : 1), log_p);

      return(lower_tail ? p : cp);
    }


    // template<class Float>
    // Float pnorm5_1(Float x, Float mu, Float sigma, Float lower_tail, Float log_p)
    // {
    //   return pnorm5_raw(x,mu,sigma,lower_tail,log_p);
    // }

    // TMB_BIND_ATOMIC(pnorm5_2,11100,pnorm5_1(x[0], x[1], x[2], x[3], x[4]))  

#endif

#define R_D_qIv(p)	(log_p	? exp(p) : (p))		/*  p  in qF(p,..) */


template<class Float>
Float frexp2(Float x, int* e){
  // Naïve
  *e = (x == 0) ? 0 : (int)trunc(1 + log(x) / log(2.0));
  return x * pow((Float)2.0, -(*e));
}

    template<class Float>
    Float ldexp2(Float num, Float exp){
      return num * pow((Float)2.0, exp);
    }

    
    template<class Float>
    Float qnorm5(Float p, Float mu, Float sigma, int lower_tail, int log_p)
{
    Float p_, q, r, val;

#ifdef IEEE_754
    if (ISNAN(p) || ISNAN(mu) || ISNAN(sigma))
	return p + mu + sigma;
#endif
    R_Q_P01_boundaries(p, ML_NEGINF, ML_POSINF);

    if(sigma  < 0)	ML_WARN_return_NAN;
    if(sigma == 0)	return mu;

    p_ = R_DT_qIv(p);/* real lower_tail prob. p */
    q = p_ - 0.5;

#ifdef DEBUG_qnorm
    REprintf("qnorm(p=%10.7g, m=%g, s=%g, l.t.= %d, log= %d): q = %g\n",
	     p,mu,sigma, lower_tail, log_p, q);
#endif


/*-- use AS 241 --- */
/* double ppnd16_(double *p, long *ifault)*/
/*      ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3

        Produces the normal deviate Z corresponding to a given lower
        tail area of P; Z is accurate to about 1 part in 10**16.

        (original fortran code used PARAMETER(..) for the coefficients
         and provided hash codes for checking them...)
*/
    if (fabs(q) <= .425) {/* |p~ - 0.5| <= .425  <==> 0.075 <= p~ <= 0.925 */
        r = .180625 - q * q; // = .425^2 - q^2  >= 0
	val =
            q * (((((((r * 2509.0809287301226727 +
                       33430.575583588128105) * r + 67265.770927008700853) * r +
                     45921.953931549871457) * r + 13731.693765509461125) * r +
                   1971.5909503065514427) * r + 133.14166789178437745) * r +
                 3.387132872796366608)
            / (((((((r * 5226.495278852854561 +
                     28729.085735721942674) * r + 39307.89580009271061) * r +
                   21213.794301586595867) * r + 5394.1960214247511077) * r +
                 687.1870074920579083) * r + 42.313330701600911252) * r + 1.);
    }
    else { /* closer than 0.075 from {0,1} boundary :
	    *  r := log(p~);  p~ = min(p, 1-p) < 0.075 :  */
	Float lp;
	if(log_p && ((lower_tail && q <= 0) || (!lower_tail && q > 0))) {
	    lp = p;
	} else {
	    lp = log( (q > 0) ? R_DT_CIv(p) /* 1-p */ : p_ /* = R_DT_Iv(p) ^=  p */);
	}
	// r = sqrt( - log(min(p,1-p)) )  <==>  min(p, 1-p) = exp( - r^2 ) :
        r = sqrt(-lp);
#ifdef DEBUG_qnorm
	REprintf("\t close to 0 or 1: r = %7g\n", r);
#endif
        if (r <= 5.) { /* <==> min(p,1-p) >= exp(-25) ~= 1.3888e-11 */
            r += -1.6;
            val = (((((((r * 7.7454501427834140764e-4 +
                       .0227238449892691845833) * r + .24178072517745061177) *
                     r + 1.27045825245236838258) * r +
                    3.64784832476320460504) * r + 5.7694972214606914055) *
                  r + 4.6303378461565452959) * r +
                 1.42343711074968357734)
                / (((((((r *
                         1.05075007164441684324e-9 + 5.475938084995344946e-4) *
                        r + .0151986665636164571966) * r +
                       .14810397642748007459) * r + .68976733498510000455) *
                     r + 1.6763848301838038494) * r +
                    2.05319162663775882187) * r + 1.);
        }
	else if(r <= 27) { /* p is very close to  0 or 1: r in (5, 27] :
		*  r >   5 <==> min(p,1-p)  < exp(-25) = 1.3888..e-11
		*  r <= 27 <==> min(p,1-p) >= exp(-27^2) = exp(-729) ~= 2.507972e-317
		* i.e., we are just barely in the range where min(p, 1-p) has not yet underflowed to zero.
		*/
	    // Wichura, p.478: minimax rational approx R_3(t) is for 5 <= t <= 27  (t :== r)
            r += -5.;
            val = (((((((r * 2.01033439929228813265e-7 +
                       2.71155556874348757815e-5) * r +
                      .0012426609473880784386) * r + .026532189526576123093) *
                    r + .29656057182850489123) * r +
                   1.7848265399172913358) * r + 5.4637849111641143699) *
                 r + 6.6579046435011037772)
                / (((((((r *
                         2.04426310338993978564e-15 + 1.4215117583164458887e-7)*
                        r + 1.8463183175100546818e-5) * r +
                       7.868691311456132591e-4) * r + .0148753612908506148525)
                     * r + .13692988092273580531) * r +
                    .59983220655588793769) * r + 1.);
        }
        else { // r > 27: p is *really* close to 0 or 1 .. practically only when log_p =TRUE
	    if(r >= 6.4e8) { // p is *very extremely* close to 0 or 1
		// Using the asymptotical formula ("0-th order"): qn = sqrt(2*s)
		val = r * M_SQRT2;
	    } else {
	      Float s2 = -ldexp2((Float)lp, (Float)1), // = -2*lp = 2s
		    x2 = s2 - log(M_2PI * s2); // = xs_1
		// if(r >= 36000.)  # <==> s >= 36000^2   use x2 = xs_1  above
		if(r < 36000.) {
		    x2 = s2 - log(M_2PI * x2) - 2./(2. + x2); // == xs_2
		    if(r < 840.) { // 27 < r < 840
			x2 = s2 - log(M_2PI * x2) + 2*log1p(- (1 - 1/(4 + x2))/(2. + x2)); // == xs_3
			if(r < 109.) { // 27 < r < 109
			  x2 = s2 - log(M_2PI * x2) +
			      2*log1p(- (1 - (1 - 5/(6 + x2))/(4. + x2))/(2. + x2)); // == xs_4
			  if(r < 55.) { // 27 < r < 55
			    x2 = s2 - log(M_2PI * x2) +
			      2*log1p(- (1 - (1 - (5 - 9/(8. + x2))/(6. + x2))/(4. + x2))/(2. + x2)); // == xs_5
			  }
			}
		    }
		}
                val = sqrt(x2);
	    }
	}
	if(q < 0.0)
	    val = -val;
    }
    return mu + sigma * val;
}




    

    
template<class Float>
 Float sinpi(Float x) {
#ifdef IEEE_754
    if (ISNAN(x)) return x;
#endif
    if(!R_FINITE(x)) ML_WARN_return_NAN;

    x = myfmod((Float)x, (Float)2.); // sin(pi(x + 2k)) == sin(pi x)  for all integer k
    // map (-2,2) --> (-1,1] :
    if(x <= -1) x += 2.; else if (x > 1.) x -= 2.;
    if(x == 0. || x == 1.) return 0.;
    if(x ==  0.5)	return  1.;
    if(x == -0.5)	return -1.;
    // otherwise
    return sin(M_PI * x);
}


    /*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 1998 Ross Ihaka
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 *
 *  SYNOPSIS
 *
 *    int chebyshev_init(double *dos, int nos, double eta)
 *    double chebyshev_eval(double x, double *a, int n)
 *
 *  DESCRIPTION
 *
 *    "chebyshev_init" determines the number of terms for the
 *    double precision orthogonal series "dos" needed to insure
 *    the error is no larger than "eta".  Ordinarily eta will be
 *    chosen to be one-tenth machine precision.
 *
 *    "chebyshev_eval" evaluates the n-term Chebyshev series
 *    "a" at "x".
 *
 *  NOTES
 *
 *    These routines are translations into C of Fortran routines
 *    by W. Fullerton of Los Alamos Scientific Laboratory.
 *
 *    Based on the Fortran routine dcsevl by W. Fullerton.
 *    Adapted from R. Broucke, Algorithm 446, CACM., 16, 254 (1973).
 */

/* NaNs propagated correctly */


template<class Float>
int attribute_hidden chebyshev_init(Float *dos, int nos, Float eta)
{
    int i, ii;
    Float err;

    if (nos < 1)
	return 0;

    err = 0.0;
    i = 0;			/* just to avoid compiler warnings */
    for (ii=1; ii<=nos; ii++) {
	i = nos - ii;
	err += fabs(dos[i]);
	if (err > eta) {
	    return i;
	}
    }
    return i;
}


template<class Float>
Float attribute_hidden chebyshev_eval(Float x, const Float *a, const int n)
{
    Float b0, b1, b2, twox;
    int i;

    if (n < 1 || n > 1000) ML_WARN_return_NAN;

    if (x < -1.1 || x > 1.1) ML_WARN_return_NAN;

    twox = x * 2;
    b2 = b1 = 0;
    b0 = 0;
    for (i = 1; i <= n; i++) {
	b2 = b1;
	b1 = b0;
	b0 = twox * b1 - b2 + a[n - i];
    }
    return (b0 - b2) * 0.5;
}

template<class Float>
Float attribute_hidden chebyshev_eval(Float x, std::vector<Float> a, const int n)
{
    Float b0, b1, b2, twox;
    int i;

    if (n < 1 || n > 1000) ML_WARN_return_NAN;

    if (x < -1.1 || x > 1.1) ML_WARN_return_NAN;

    twox = x * 2;
    b2 = b1 = 0;
    b0 = 0;
    for (i = 1; i <= n; i++) {
	b2 = b1;
	b1 = b0;
	b0 = twox * b1 - b2 + a[n - i];
    }
    return (b0 - b2) * 0.5;
}

template<class Float>
Float lgammafn(Float x);

/*
 *  AUTHOR
 *    Catherine Loader, catherine@research.bell-labs.com.
 *    October 23, 2000.
 *
 *  Merge in to R:
 *	Copyright (C) 2000-2024, The R Core Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 *
 *
 *  DESCRIPTION
 *
 *    Computes the log of the error term in Stirling's formula.
 *      For n > 15, uses the series 1/12n - 1/360n^3 + ...
 *      For n <=15, integers or half-integers, uses stored values.
 *      For other n < 15, uses lgamma directly (don't use this to
 *        write lgamma!)
 *
 * Merge in to R:
 * Copyright (C) 2000, The R Core Team
 * R has lgammafn, and lgamma is not part of ISO C
 */

/* stirlerr(n) = log(n!) - log( sqrt(2*pi*n)*(n/e)^n )
 *             = log Gamma(n+1) - 1/2 * [log(2*pi) + log(n)] - n*[log(n) - 1]
 *             = log Gamma(n+1) - (n + 1/2) * log(n) + n - log(2*pi)/2
 *
 * see also lgammacor() in ./lgammacor.c  which computes almost the same!
 *
 * NB: stirlerr(n/2) & stirlerr((n+1)/2) are called from dt(x,n) for all real n > 0 ;
 *     stirlerr(x) from gammafn(x) when |x| > 10,  2|x| is integer, but |x| is *not* in {11:50}
 *     stirlerr(x) from dpois_raw(x, lam) for any x > 0  which itself is called by many,
 *                                        including pgamma(), hence ppois(), ..

 *     stirlerr(n), stirlerr(x), stirlerr(n-x) from binom_raw(x, n, ..) for all possible 0 < x < n
 */


 template<class Float>
 Float lgamma1p (Float a);


template<class Float>
Float attribute_hidden stirlerr(Float n)
{

#define S0 0.083333333333333333333       /* 1/12 */
#define S1 0.00277777777777777777778     /* 1/360 */
#define S2 0.00079365079365079365079365  /* 1/1260 */
#define S3 0.000595238095238095238095238 /* 1/1680 */
#define S4 0.0008417508417508417508417508/* 1/1188 */
#define S5 0.0019175269175269175269175262 // 691/360360
#define S6 0.0064102564102564102564102561 // 1/156
#define S7 0.029550653594771241830065352  // 3617/122400
#define S8 0.17964437236883057316493850   // 43867/244188
#define S9 1.3924322169059011164274315    // 174611/125400
#define S10 13.402864044168391994478957 // 77683/5796
#define S11 156.84828462600201730636509 // 236364091/1506960
#define S12 2193.1033333333333333333333 // 657931/300
#define S13 36108.771253724989357173269 // 3392780147/93960
#define S14 691472.26885131306710839498 // 1723168255201/2492028
#define S15 15238221.539407416192283370 // 7709321041217/505920
#define S16 382900751.39141414141414141 // 151628697551/396
/* #define S17 10882266035.784391089015145 // 26315271553053477373/2418179400 */

/*
  exact values for 0, 0.5, 1.0, 1.5, ..., 14.5, 15.0.
*/
    const static Float sferr_halves[31] = {
	0.0, /* n=0 - wrong, place holder only */
	0.1534264097200273452913848,  /* 0.5 */
	0.0810614667953272582196702,  /* 1.0 */
	0.0548141210519176538961390,  /* 1.5 */
	0.0413406959554092940938221,  /* 2.0 */
	0.03316287351993628748511048, /* 2.5 */
	0.02767792568499833914878929, /* 3.0 */
	0.02374616365629749597132920, /* 3.5 */
	0.02079067210376509311152277, /* 4.0 */
	0.01848845053267318523077934, /* 4.5 */
	0.01664469118982119216319487, /* 5.0 */
	0.01513497322191737887351255, /* 5.5 */
	0.01387612882307074799874573, /* 6.0 */
	0.01281046524292022692424986, /* 6.5 */
	0.01189670994589177009505572, /* 7.0 */
	0.01110455975820691732662991, /* 7.5 */
	0.010411265261972096497478567, /* 8.0 */
	0.009799416126158803298389475, /* 8.5 */
	0.009255462182712732917728637, /* 9.0 */
	0.008768700134139385462952823, /* 9.5 */
	0.008330563433362871256469318, /* 10.0 */
	0.007934114564314020547248100, /* 10.5 */
	0.007573675487951840794972024, /* 11.0 */
	0.007244554301320383179543912, /* 11.5 */
	0.006942840107209529865664152, /* 12.0 */
	0.006665247032707682442354394, /* 12.5 */
	0.006408994188004207068439631, /* 13.0 */
	0.006171712263039457647532867, /* 13.5 */
	0.005951370112758847735624416, /* 14.0 */
	0.005746216513010115682023589, /* 14.5 */
	0.005554733551962801371038690  /* 15.0 */
    };
    Float nn;

    if (n <= 23.5) {
	nn = n + n;
	if (n <= 15. && (nn == (int)trunc(nn))) return sferr_halves[(int)trunc(nn)];
	// else:
	if (n <= 5.25) {
	    if(n >= 1.) { // "MM2"; slightly more accurate than direct form
		Float l_n = log(n);	      // ldexp(u, -1) == u/2
		return lgammafn(n) + n*(1 - l_n) + ldexp2((Float)(l_n - M_LN_2PI), (Float)-1);
	    }
	    else // n < 1
		return lgamma1p(n) - (n + 0.5)*log(n) + n - M_LN_SQRT_2PI;
	}
	// else  5.25 < n <= 23.5
	nn = n*n;
	if (n > 12.8)	return (S0-(S1-(S2-(S3-(S4-(S5 -S6/nn)/nn)/nn)/nn)/nn)/nn)/n;		// k = 7
	if (n > 12.3)	return (S0-(S1-(S2-(S3-(S4-(S5-(S6 -S7/nn)/nn)/nn)/nn)/nn)/nn)/nn)/n;	// k = 8
	if (n >  8.9)	return (S0-(S1-(S2-(S3-(S4-(S5-(S6-(S7 -S8/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/n;	// k = 9
	/* if (n >  7.9)	return (S0-(S1-(S2-(S3-(S4-(S5-(S6-(S7-(S8 -S9/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/n; skip k = 10 */
	if (n >  7.3)	return (S0-(S1-(S2-(S3-(S4-(S5-(S6-(S7-(S8-(S9-S10/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/n; // 11
	/* if (n >  6.5)	return (S0-(S1-(S2-(S3-(S4-(S5-(S6-(S7-(S8-(S9-(S10-S11/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/n; skip k=12*/
	if (n >  6.6)	return (S0-(S1-(S2-(S3-(S4-(S5-(S6-(S7-(S8-(S9-(S10-(S11-S12/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/n;
	/* if (n >  5.7)	return (S0-(S1-(S2-(S3-(S4-(S5-(S6-(S7-(S8-(S9-(S10-(S11-(S12-S13/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/n; skip k= 14 */
	if (n >  6.1)	return (S0-(S1-(S2-(S3-(S4-(S5-(S6-(S7-(S8-(S9-(S10-(S11-(S12-(S13-S14/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/n; // k = 15
	/* ....		return (S0-(S1-(S2-(S3-(S4-(S5-(S6-(S7-(S8-(S9-(S10-(S11-(S12-(S13-(S14-S15/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/n;
	 * skip order k=16 : never "good" for double prec */
	/* 6.1 >= n > 5.25 */
	return (S0-(S1-(S2-(S3-(S4-(S5-(S6-(S7-(S8-(S9-(S10-(S11-(S12-(S13-(S14-(S15-S16/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/n;
	/* return (S0-(S1-(S2-(S3-(S4-(S5-(S6-(S7-(S8-(S9-(S10-(S11-(S12-(S13-(S14-(S15-(S16-S17/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/n; */

     } else { // n > 23.5
	nn = n*n;
	if (n > 15.7e6)	return S0/n;
	if (n > 6180)	return (S0 -S1/nn)/n;
	if (n > 205)    return (S0-(S1 -S2/nn)/nn)/n;
	if (n > 86)	return (S0-(S1-(S2 -S3/nn)/nn)/nn)/n;
	if (n > 27)	return (S0-(S1-(S2-(S3 -S4/nn)/nn)/nn)/nn)/n;
	/* 23.5 < n <= 27 */
	return (S0-(S1-(S2-(S3-(S4 -S5/nn)/nn)/nn)/nn)/nn)/n;
    }
}

/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 2000-2021 The R Core Team
 *  Copyright (C) 1998 Ross Ihaka
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 *
 *  SYNOPSIS
 *
 *    #include <Rmath.h>
 *    double lgammacor(double x);
 *
 *  DESCRIPTION
 *
 *    Compute the log gamma correction factor for x >= 10 so that
 *
 *    log(gamma(x)) = .5*log(2*pi) + (x-.5)*log(x) -x + lgammacor(x)
 *
 *    [ lgammacor(x) is called	Del(x)	in other contexts (e.g. dcdflib)]
 *
 *  NOTES
 *
 *    This routine is a translation into C of a Fortran subroutine
 *    written by W. Fullerton of Los Alamos Scientific Laboratory.
 *
 *  SEE ALSO
 *
 *    Loader(1999)'s stirlerr() {in ./stirlerr.c} is *very* similar in spirit,
 *    is faster and cleaner, but is only defined "fast" for half integers.
 */

    template<class Float>
Float attribute_hidden lgammacor(Float x)
{
    const static Float algmcs[15] = {  // below, nalgm = 5 ==> only the first 5 are used!
	+.1666389480451863247205729650822e+0,
	-.1384948176067563840732986059135e-4,
	+.9810825646924729426157171547487e-8,
	-.1809129475572494194263306266719e-10,
	+.6221098041892605227126015543416e-13,
	-.3399615005417721944303330599666e-15,
	+.2683181998482698748957538846666e-17,
	-.2868042435334643284144622399999e-19,
	+.3962837061046434803679306666666e-21,
	-.6831888753985766870111999999999e-23,
	+.1429227355942498147573333333333e-24,
	-.3547598158101070547199999999999e-26,
	+.1025680058010470912000000000000e-27,
	-.3401102254316748799999999999999e-29,
	+.1276642195630062933333333333333e-30
    };

    std::vector<Float> algmcs_vec;
    for(int i = 0; i < 15; ++i)
      algmcs_vec.push_back(algmcs[i]);

/* For IEEE double precision DBL_EPSILON = 2^-52 = 2.220446049250313e-16 :
 *   xbig = 2 ^ 26.5
 *   xmax = DBL_MAX / 48 =  2^1020 / 3 */
#define nalgm 5
#define xbig  94906265.62425156
#define xmax  3.745194030963158e306

    if (x < 10) // possibly consider stirlerr()
      ML_WARN_return_NAN;
    else if (x >= xmax) {
	// ML_WARNING(ME_UNDERFLOW, "lgammacor");
	/* allow to underflow below */
    }
    else if (x < xbig) {
	Float tmp = 10 / x;
	return chebyshev_eval((Float)(tmp * tmp * 2 - 1), algmcs_vec, nalgm) / x;
    }
    // else, xbig <= x < xmax :
    return 1 / (x * 12);
}
    

    
/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 2000-2024 The R Core Team
 *  Copyright (C) 2002-2024 The R Foundation
 *  Copyright (C) 1998 Ross Ihaka
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 *
 *  SYNOPSIS
 *
 *    #include <Rmath.h>
 *    double gammafn(double x);
 *
 *  DESCRIPTION
 *
 *    This function computes the value of the gamma function.
 *
 *  NOTES
 *
 *    This function is a translation into C of a Fortran subroutine
 *    by W. Fullerton of Los Alamos Scientific Laboratory.
 *    (e.g. http://www.netlib.org/slatec/fnlib/gamma.f)
 *
 *    The accuracy of this routine compares (very) favourably
 *    with those of the Sun Microsystems portable mathematical
 *    library.
 *
 *    MM specialized the case of  n!  for n < 50 - for even better precision
 */


 template<class Float>
 Float gammafn(Float x)
{
    const static Float gamcs[42] = {
	+.8571195590989331421920062399942e-2,
	+.4415381324841006757191315771652e-2,
	+.5685043681599363378632664588789e-1,
	-.4219835396418560501012500186624e-2,
	+.1326808181212460220584006796352e-2,
	-.1893024529798880432523947023886e-3,
	+.3606925327441245256578082217225e-4,
	-.6056761904460864218485548290365e-5,
	+.1055829546302283344731823509093e-5,
	-.1811967365542384048291855891166e-6,
	+.3117724964715322277790254593169e-7,
	-.5354219639019687140874081024347e-8,
	+.9193275519859588946887786825940e-9,
	-.1577941280288339761767423273953e-9,
	+.2707980622934954543266540433089e-10,
	-.4646818653825730144081661058933e-11,
	+.7973350192007419656460767175359e-12,
	-.1368078209830916025799499172309e-12,
	+.2347319486563800657233471771688e-13,
	-.4027432614949066932766570534699e-14,
	+.6910051747372100912138336975257e-15,
	-.1185584500221992907052387126192e-15,
	+.2034148542496373955201026051932e-16,
	-.3490054341717405849274012949108e-17,
	+.5987993856485305567135051066026e-18,
	-.1027378057872228074490069778431e-18,
	+.1762702816060529824942759660748e-19,
	-.3024320653735306260958772112042e-20,
	+.5188914660218397839717833550506e-21,
	-.8902770842456576692449251601066e-22,
	+.1527474068493342602274596891306e-22,
	-.2620731256187362900257328332799e-23,
	+.4496464047830538670331046570666e-24,
	-.7714712731336877911703901525333e-25,
	+.1323635453126044036486572714666e-25,
	-.2270999412942928816702313813333e-26,
	+.3896418998003991449320816639999e-27,
	-.6685198115125953327792127999999e-28,
	+.1146998663140024384347613866666e-28,
	-.1967938586345134677295103999999e-29,
	+.3376448816585338090334890666666e-30,
	-.5793070335782135784625493333333e-31
    };

    std::vector<Float> gamcs_vec;
    for(int i = 0; i < 42; ++i)
      gamcs_vec.push_back(gamcs[i]);

#ifdef NOMORE_FOR_THREADS
    static int ngam = 0;
    static Float xmin = 0, xmax = 0., xsml = 0., dxrel = 0.;

    /* Initialize machine dependent constants, the first time gamma() is called.
	FIXME for threads ! */
    if (ngam == 0) {
	ngam = chebyshev_init(gamcs, 42, DBL_EPSILON/20);/*was .1*d1mach(3)*/
	gammalims(&xmin, &xmax);/*-> ./gammalims.c */
	xsml = exp(fmax2(log(DBL_MIN), -log(DBL_MAX)) + 0.01);
	/*   = exp(.01)*DBL_MIN = 2.247e-308 for IEEE */
	dxrel = sqrt(DBL_EPSILON);/*was sqrt(d1mach(4)) */
    }
#else
/* For IEEE double precision DBL_EPSILON = 2^-52 = 2.220446049250313e-16 :
 * (xmin, xmax) are non-trivial, see ./gammalims.c
 * xsml = exp(.01)*DBL_MIN
 * dxrel = sqrt(DBL_EPSILON) = 2 ^ -26
*/
# define ngam 22
# define xmin -170.5674972726612
#ifdef xmax
#undef xmax
#endif
# define xmax  171.61447887182298
# define xsml 2.2474362225598545e-308
# define dxrel 1.490116119384765696e-8
#endif

    if(ISNAN(x)) return x;

    /* If the argument is exactly zero or a negative integer
     * then return NaN. */
    if (x == 0 || (x < 0 && x == round(x))) {
	// ML_WARNING(ME_DOMAIN, "gammafn");
	return ML_NAN;
    }

    int i;
    Float y = fabs(x), value;

    if (y <= 10) {

	/* Compute gamma(x) for -10 <= x <= 10
	 * Reduce the interval and find gamma(1 + y) for 0 <= y < 1
	 * first of all. */

      int n = (int) trunc(x);
	if(x < 0) --n;
	y = x - n;/* n = floor(x)  ==>	y in [ 0, 1 ) */
	--n;
	value = chebyshev_eval((Float)(y * 2 - 1), gamcs_vec, ngam) + .9375;
	if (n == 0)
	    return value;/* x = 1.dddd = 1+y */

	if (n < 0) {
	    /* compute gamma(x) for -10 <= x < 1 */

	    /* exact 0 or "-n" checked already above */

	    /* The answer is less than half precision */
	    /* because x too near a negative integer. */
	    if (x < -0.5 && fabs(x - (int)trunc(x - 0.5) / x) < dxrel) {
		// ML_WARNING(ME_PRECISION, "gammafn");
	    }

	    /* The argument is so close to 0 that the result would overflow. */
	    if (y < xsml) {
		// ML_WARNING(ME_RANGE, "gammafn");
		if(x > 0) return ML_POSINF;
		else return ML_NEGINF;
	    }

	    n = -n;

	    for (i = 0; i < n; i++) {
		value /= (x + i);
	    }
	    return value;
	}
	else {
	    /* gamma(x) for 2 <= x <= 10 */

	    for (i = 1; i <= n; i++) {
		value *= (y + i);
	    }
	    return value;
	}
    }
    else {
	/* gamma(x) for	 y = |x| > 10. */

	if (x > xmax) {			/* Overflow */
	    // No warning: +Inf is the best answer
	    return ML_POSINF;
	}

	if (x < xmin) {			/* Underflow */
	    // No warning: 0 is the best answer
	    return 0.;
	}

	if(y <= 50 && y == (int)trunc(y)) { /* compute (n - 1)! */
	    value = 1.;
	    for (i = 2; i < y; i++) value *= i;
	}
	else { /* normal case */
	    value = exp((y - 0.5) * log(y) - y + M_LN_SQRT_2PI +
			((2*y == (int)trunc(2*y)) ? atomic::gamma_utils::stirlerr(y) : lgammacor(y)));
	}

	if (x > 0)
	    return value;
	// else:  x < 0, not an integer :

	if (fabs((x - (int)trunc(x - 0.5))/x) < dxrel) {
	    /* The answer is less than half precision because */
	    /* the argument is too near a negative integer. */

	    // ML_WARNING(ME_PRECISION, "gammafn");
	}

	Float sinpiy = sinpi(y);
	if (sinpiy == 0) {		/* Negative integer arg - overflow */
	    // ML_WARNING(ME_RANGE, "gammafn");
	    return ML_POSINF;
	}

	return -M_PI / (y * sinpiy * value);
    }
}
    


    /*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 2000-2020 The R Core Team
 *  Copyright (C) 1998 Ross Ihaka
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 *
 *  SYNOPSIS
 *
 *    #include <Rmath.h>
 *    double lgammafn_sign(double x, int *sgn);
 *    double lgammafn(double x);
 *
 *  DESCRIPTION
 *
 *    The function lgammafn computes log|gamma(x)|.  The function
 *    lgammafn_sign in addition assigns the sign of the gamma function
 *    to the address in the second argument if this is not NULL.
 *
 *  NOTES
 *
 *    This routine is a translation into C of a Fortran subroutine
 *    by W. Fullerton of Los Alamos Scientific Laboratory.
 *
 *    The accuracy of this routine compares (very) favourably
 *    with those of the Sun Microsystems portable mathematical
 *    library.
 *
 *  ./toms708.c  has  gamln()
 */

  /*
   * Modified 2024 Christoffer Moesgaard Albertsen, Technical University of Denmark
   */

    
 template<class Float>
    Float lgammafn_sign(Float x, int *sgn)
{
    Float ans, y, sinpiy;

#ifdef NOMORE_FOR_THREADS
    static Float xmax = 0.;
    static Float dxrel = 0.;

    if (xmax == 0) {/* initialize machine dependent constants _ONCE_ */
	xmax = d1mach(2)/log(d1mach(2));/* = 2.533 e305	 for IEEE double */
	dxrel = sqrt (d1mach(4));/* sqrt(Eps) ~ 1.49 e-8  for IEEE double */
    }
#else
/* For IEEE double precision DBL_EPSILON = 2^-52 = 2.220446049250313e-16 :
   xmax  = DBL_MAX / log(DBL_MAX) = 2^1024 / (1024 * log(2)) = 2^1014 / log(2)
   dxrel = sqrt(DBL_EPSILON) = 2^-26 = 5^26 * 1e-26 (is *exact* below !)
 */
    #ifdef xmax
    #undef xmax
    #endif
#define xmax  2.5327372760800758e+305
    #ifdef dxrel
    #undef dxrel
    #endif
#define dxrel 1.490116119384765625e-8
#endif

    if (sgn != NULL) *sgn = 1;

#ifdef IEEE_754
    if(ISNAN(x)) return x;
#endif

    if (sgn != NULL && x < 0 && myfmod((Float)floor(-x), (Float)2.) == 0)
	*sgn = -1;

    if (x <= 0 && x == trunc(x)) { /* Negative integer argument */
	// No warning: this is the best answer; was  ML_WARNING(ME_RANGE, "lgamma");
	return ML_POSINF;/* +Inf, since lgamma(x) = log|gamma(x)| */
    }

    y = fabs(x);

    if (y < 1e-306) return -log(y); // denormalized range, R change
    if (y <= 10) return log(fabs(gammafn(x)));
    /*
      ELSE  y = |x| > 10 ---------------------- */

    if (y > xmax) {
	// No warning: +Inf is the best answer
	return ML_POSINF;
    }

    if (x > 0) { /* i.e. y = x > 10 */
#ifdef IEEE_754
	if(x > 1e17)
	    return(x*(log(x) - 1.));
	else if(x > 4934720.)
	    return(M_LN_SQRT_2PI + (x - 0.5) * log(x) - x);
	else
#endif
	    return M_LN_SQRT_2PI + (x - 0.5) * log(x) - x + atomic::gamma_utils::lgammacor(x);
    }
    /* else: x < -10; y = -x */
    sinpiy = fabs(sinpi(y));

    if (sinpiy == 0) { /* Negative integer argument ===
			  Now UNNECESSARY: caught above */
	MATHLIB_WARNING(" ** should NEVER happen! *** [lgamma.c: Neg.int, y=%g]\n",y);
	ML_WARN_return_NAN;
    }

    ans = M_LN_SQRT_PId2 + (x - 0.5) * log(y) - x - log(sinpiy) - atomic::gamma_utils::lgammacor(y);

    if(fabs((x - trunc(x - 0.5)) * ans / x) < dxrel) {

	/* The answer is less than half precision because
	 * the argument is too near a negative integer; e.g. for  lgamma(1e-7 - 11) */

	ML_WARNING(ME_PRECISION, "lgamma");
    }

    return ans;
}

template<class Float>
 Float lgammafn(Float x)
{
    return lgammafn_sign(x, NULL);
}


/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 1998 Ross Ihaka
 *  Copyright (C) 2000-12 The R Core Team
 *  Copyright (C) 2003 The R Foundation
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 *
 *  SYNOPSIS
 *
 *    #include <Rmath.h>
 *    double lbeta(double a, double b);
 *
 *  DESCRIPTION
 *
 *    This function returns the value of the log beta function
 *
 *	log B(a,b) = log G(a) + log G(b) - log G(a+b)
 *
 *  NOTES
 *
 *    This routine is a translation into C of a Fortran subroutine
 *    by W. Fullerton of Los Alamos Scientific Laboratory.
 */

template<class Float>
Float lbeta(Float a, Float b)
{
    Float corr, p, q;

#ifdef IEEE_754
    if(ISNAN(a) || ISNAN(b))
	return a + b;
#endif
    p = q = a;
    if(b < p) p = b;/* := min(a,b) */
    if(b > q) q = b;/* := max(a,b) */

    /* both arguments must be >= 0 */
    if (p < 0)
      ML_WARN_return_NAN;
    else if (p == 0) {
	return ML_POSINF;
    }
    else if (!R_FINITE(q)) { /* q == +Inf */
	return ML_NEGINF;
    }

    if (p >= 10) {
	/* p and q are big. */
	corr = lgammacor(p) + lgammacor(q) - lgammacor(p + q);
	return log(q) * -0.5 + M_LN_SQRT_2PI + corr
		+ (p - 0.5) * log(p / (p + q)) + q * log1p(-p / (p + q));
    }
    else if (q >= 10) {
	/* p is small, but q is big. */
	corr = lgammacor(q) - lgammacor(p + q);
	return lgammafn(p) + corr + p - p * log(p + q)
		+ (q - 0.5) * log1p(-p / (p + q));
    }
    else {
	/* p and q are small: p <= q < 10. */
	/* R change for very small args */
	if (p < 1e-306) return lgammafn(p) + (lgammafn(q) - lgammafn(p+q));
	else return log(gammafn(p) * (gammafn(q) / gammafn(p + q)));
    }
}




/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
 *  Copyright (C) 2000-2007   The R Core Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 */

template<class Float>
Float pt(Float x, Float n, int lower_tail, int log_p)
{
/* return  P[ T <= x ]	where
 * T ~ t_{n}  (t distrib. with n degrees of freedom).

 *	--> ./pnt.c for NON-central
 */
    Float val, nx;
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(n))
	return x + n;
#endif
    if (n <= 0.0) ML_WARN_return_NAN;

    if(!R_FINITE(x))
	return (x < 0) ? R_DT_0 : R_DT_1;
    if(!R_FINITE(n))
      return pnorm5_raw((Float)x, (Float)0.0, (Float)1.0, (Float)lower_tail, (Float)log_p);

#ifdef R_version_le_260
    if (n > 4e5) { /*-- Fixme(?): test should depend on `n' AND `x' ! */
	/* Approx. from	 Abramowitz & Stegun 26.7.8 (p.949) */
	val = 1./(4.*n);
	return pnorm5_raw((Float)(x*(1. - val)/sqrt(1. + x*x*2.*val)), (Float)0.0, (Float)1.0,
			  (Float)lower_tail, (Float)log_p);
    }
#endif

    nx = 1 + (x/n)*x;
    /* FIXME: This test is probably losing rather than gaining precision,
     * now that pbeta(*, log_p = TRUE) is much better.
     * Note however that a version of this test *is* needed for x*x > D_MAX */
    if(nx > 1e100) { /* <==>  x*x > 1e100 * n  */
	/* Danger of underflow. So use Abramowitz & Stegun 26.5.4
	   pbeta(z, a, b) ~ z^a(1-z)^b / aB(a,b) ~ z^a / aB(a,b),
	   with z = 1/nx,  a = n/2,  b= 1/2 :
	*/
	Float lval;
	lval = -0.5*n*(2*log(fabs(x)) - log(n))
	  - lbeta((Float)(0.5*n), (Float)0.5) - log(0.5*n);
	val = log_p ? lval : exp(lval);
    } else {
	val = (n > x * x)
	  ? atomic::toms708::pbeta ((Float)(x * x / (n + x * x)), (Float)0.5, (Float)(n / 2.), /*lower_tail*/(int)0, (int)log_p)
	  : atomic::toms708::pbeta ((Float)(1. / nx), (Float)(n / 2.), (Float)0.5, /*lower_tail*/(int)1, (int)log_p);
    }

    /* Use "1 - v"  if	lower_tail  and	 x > 0 (but not both):*/
    if(x <= 0.)
	lower_tail = !lower_tail;

    if(log_p) {
	if(lower_tail) return log1p(-0.5*exp(val));
	else return val - M_LN2; /* = log(.5* pbeta(....)) */
    }
    else {
	val /= 2.;
	return R_D_Cval(val);
    }
}

template<class Float>
Float Rtanpi(Float x)
{
#ifdef IEEE_754
    if (ISNAN(x)) return x;
#endif
    if(!R_FINITE(x)) ML_WARN_return_NAN;

    x = myfmod((Float)x, (Float)1.); // tan(pi(x + k)) == tan(pi x)  for all integer k
    // map (-1,1] --> (-1/2, 1/2] :
    if(x <= -0.5){
      x = x + 1.0;
    }else if(x > 0.5){
      x = x - 1.0;
    }
    return (x == 0.) ? 0. :
	((x ==  0.5 ) ? ML_NAN :
	((x ==  0.25) ?  1. :
	((x == -0.25) ? -1. :
			tan(M_PI * x)
	    )));
}


template<class Float>
Float qt(Float p, Float ndf, int lower_tail, int log_p)
{
  const static Float eps = 1.e-12;

  Float P, q;

#ifdef IEEE_754
  if (ISNAN(p) || ISNAN(ndf))
    return p + ndf;
#endif

  R_Q_P01_boundaries(p, ML_NEGINF, ML_POSINF);

  if (ndf <= 0) ML_WARN_return_NAN;

  if (ndf < 1) { /* based on qnt */
    const static Float accu = 1e-13;
    const static Float Eps = 1e-11; /* must be > accu */

    Float ux, lx, nx, pp;

    int iter = 0;

    p = R_DT_qIv(p);

    /* Invert pt(.) :
     * 1. finding an upper and lower bound */
    if(p > 1 - DBL_EPSILON) return ML_POSINF;
    pp = myfmin2((Float)(1 - DBL_EPSILON), (Float)(p * (1 + Eps)));
    //pp = 0.5 * ((1 - DBL_EPSILON) + (p * (1 + Eps))) - 0.5 * fabs(((1 - DBL_EPSILON) - (p * (1 + Eps))));
    for(ux = 1.; ux < DBL_MAX && pt(ux, ndf, TRUE, FALSE) < pp; ux *= 2);
    pp = p * (1 - Eps);
    for(lx =-1.; lx > -DBL_MAX && pt(lx, ndf, TRUE, FALSE) > pp; lx *= 2);

    /* 2. interval (lx,ux)  halving
       regula falsi failed on qt(0.1, 0.1)
    */
    do {
      nx = 0.5 * (lx + ux);
      if (pt(nx, ndf, TRUE, FALSE) > p) ux = nx; else lx = nx;
    } while ((ux - lx) / fabs(nx) > accu && ++iter < 1000);

    if(iter >= 1000) ML_WARNING(ME_PRECISION, "qt");

    return 0.5 * (lx + ux);
  }

  /* Old comment:
   * FIXME: "This test should depend on  ndf  AND p  !!
   * -----  and in fact should be replaced by
   * something like Abramowitz & Stegun 26.7.5 (p.949)"
   *
   * That would say that if the qnorm value is x then
   * the result is about x + (x^3+x)/4df + (5x^5+16x^3+3x)/96df^2
   * The differences are tiny even if x ~ 1e5, and qnorm is not
   * that accurate in the extreme tails.
   */
  if (ndf > 1e20) return qnorm5(p, (Float)0., (Float)1., lower_tail, log_p);

  P = R_D_qIv(p); /* if exp(p) underflows, we fix below */

  bool neg = (!lower_tail || P < 0.5) && (lower_tail || P > 0.5);
  bool is_neg_lower = (lower_tail == neg); /* both TRUE or FALSE == !xor */
  if(neg)
    P = 2 * (log_p ? (lower_tail ? P : -expm1(p)) : R_D_Lval(p));
  else
    P = 2 * (log_p ? (lower_tail ? -expm1(p) : P) : R_D_Cval(p));
  /* 0 <= P <= 1 ; P = 2*min(P', 1 - P')  in all cases */

  if (fabs(ndf - 2) < eps) {	/* df ~= 2 */
    if(P > DBL_MIN) {
      if(3* P < DBL_EPSILON) /* P ~= 0 */
	q = 1 / sqrt(P);
      else if (P > 0.9)	   /* P ~= 1 */
	q = (1 - P) * sqrt(2 /(P * (2 - P)));
      else /* eps/3 <= P <= 0.9 */
	q = sqrt(2 / (P * (2 - P)) - 2);
    }
    else { /* P << 1, q = 1/sqrt(P) = ... */
      if(log_p)
	q = is_neg_lower ? exp(- p/2) / M_SQRT2 : 1/sqrt(-expm1(p));
      else
	q = ML_POSINF;
    }
  }
  else if (ndf < 1 + eps) { /* df ~= 1  (df < 1 excluded above): Cauchy */
    if(P == 1.) q = 0; // some versions of tanpi give Inf, some NaN
    else if(P > 0)
      q = ((Float)1)/Rtanpi(P/(Float)2.);/* == - tan((P+1) * M_PI_2) -- suffers for P ~= 0 */

    else { /* P = 0, but maybe = 2*exp(p) ! */
      if(log_p) /* 1/tan(e) ~ 1/e */
	q = is_neg_lower ? M_1_PI * exp(-p) : -1./(M_PI * expm1(p));
      else
	q = ML_POSINF;
    }
  }
  else {		/*-- usual case;  including, e.g.,  df = 1.1 */
    Float x = 0., y, log_P2 = 0./* -Wall */,
      a = 1 / (ndf - 0.5),
      b = 48 / (a * a),
      c = ((20700 * a / b - 98) * a - 16) * a + 96.36,
      d = ((94.5 / (b + c) - 3) / b + 1) * sqrt(a * M_PI_2) * ndf;
    bool
      P_ok1 = P > DBL_MIN || !log_p,
      P_ok  = P_ok1; // when true (after check below), use "normal scale": log_p=FALSE
    if(P_ok1) {
      y = pow(d * P, 2.0 / ndf);
      P_ok = (y >= DBL_EPSILON);
    }
    if(!P_ok) {// log.p && P very.small  ||  (d*P)^(2/df) =: y < eps_c
      log_P2 = is_neg_lower ? R_D_log(p) : R_D_LExp(p); /* == log(P / 2) */
      x = (log(d) + M_LN2 + log_P2) / ndf;
      y = exp(2 * x);
    }

    if ((ndf < 2.1 && P > 0.5) || y > 0.05 + a) { /* P > P0(df) */
      /* Asymptotic inverse expansion about normal */
      if(P_ok)
	x = qnorm5((Float)(0.5 * P), (Float)0., (Float)1., /*lower_tail*/TRUE,  /*log_p*/FALSE);
      else /* log_p && P underflowed */
	x = qnorm5((Float)log_P2,  (Float)0., (Float)1., lower_tail,	        /*log_p*/ TRUE);

      y = x * x;
      if (ndf < 5)
	c += 0.3 * (ndf - 4.5) * (x + 0.6);
      c = (((0.05 * d * x - 5) * x - 7) * x - 2) * x + b + c;
      y = (((((0.4 * y + 6.3) * y + 36) * y + 94.5) / c
	    - y - 3) / b + 1) * x;
      y = expm1(a * y * y);
      q = sqrt(ndf * y);
    } else if(!P_ok && x < - M_LN2 * DBL_MANT_DIG) {/* 0.5* log(DBL_EPSILON) */
      /* y above might have underflown */
      q = sqrt(ndf) * exp(-x);
    }
    else { /* re-use 'y' from above */
      y = ((1 / (((ndf + 6) / (ndf * y) - 0.089 * d - 0.822)
		 * (ndf + 2) * 3) + 0.5 / (ndf + 4))
	   * y - 1) * (ndf + 1) / (ndf + 2) + 1 / y;
      q = sqrt(ndf * y);
    }


    /* Now apply 2-term Taylor expansion improvement (1-term = Newton):
     * as by Hill (1981) [ref.above] */

    /* FIXME: This can be far from optimal when log_p = TRUE
     *      but is still needed, e.g. for qt(-2, df=1.01, log=TRUE).
     *	Probably also improvable when  lower_tail = FALSE */

    if(P_ok1) {
      Float M = fabs(sqrt(DBL_MAX/2.) - ndf);
      int it=0;
      while(it++ < 10 && (y = dt(q, ndf, FALSE)) > 0 &&
	    R_FINITE(x = (pt(q, ndf, FALSE, FALSE) - P/2) / y) &&
	    fabs(x) > 1e-14*fabs(q)) {
	/* Newton (=Taylor 1 term):
	 *  q += x;
	 * Taylor 2-term : */
	Float F = (fabs(q) < M) ?
	  q * (ndf + 1) / (2 * (q * q + ndf)) :
	  (ndf + 1) / (2 * (q     + ndf/q)),
	  del_q = x * (1. + x * F);
	if(R_FINITE(del_q) && R_FINITE(q + del_q))
	  q += del_q;
	else if(R_FINITE(x) && R_FINITE(q + x))
	  q += x;
	else // FIXME??  if  q+x = +/-Inf is *better* than q should use it
	  break; // cannot improve  q  with a Newton/Taylor step
      }
    }
  }
  return neg ? -q : q;
}


/*
 *  Algorithm AS 275 Appl.Statist. (1992), vol.41, no.2
 *  original  (C) 1992	     Royal Statistical Society
 *
 *  Computes the noncentral chi-squared distribution function with
 *  positive real degrees of freedom df and nonnegative noncentrality
 *  parameter ncp.  pnchisq_raw is based on
 *
 *    Ding, C. G. (1992)
 *    Algorithm AS275: Computing the non-central chi-squared
 *    distribution function. Appl.Statist., 41, 478-482.

 *  Other parts
 *  Copyright (C) 2000-2019  The R Core Team
 *  Copyright (C) 2003-2015  The R Foundation
 */

/*
 * Modified 2024-09-10 Christoffer Moesgaard Albertsen, Technical University of Denmark (cmoe@aqua.dtu.dk)
 */

static const double _dbl_min_exp = M_LN2 * DBL_MIN_EXP;

// template<class Float>
// Float pgamma(Float x, Float alph, Float scale, int lower_tail, int log_p);

    #include "pgamma.hpp"


/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 1998	Ross Ihaka
 *  Copyright (C) 2000	The R Core Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 *
 *  DESCRIPTION
 *
 *     The distribution function of the chi-squared distribution.
 */
template<class Float>
Float pchisq(Float x, Float df, int lower_tail, int log_p)
{
  Float v = pgamma((Float)x, (Float)(df/2.), (Float)2., (int)lower_tail, (int)log_p);
  return v;
}


template<class Float>
Float attribute_hidden
pnchisq_raw(Float x, Float f, Float theta /* = ncp */,
	    Float errmax, Float reltol, int itrmax,
	    int lower_tail, int log_p)
{
    Float lam, x2, f2, term, bound, f_x_2n, f_2n;
    Float l_lam = -1., l_x = -1.; /* initialized for -Wall */
    int n;
    bool lamSml, tSml, is_r, is_b;
    Float ans, u, v, t, lt, lu =-1;

    if (x <= 0.) {
	if(x == 0. && f == 0.) { // chi^2_0(.) has point mass at zero
#define _L  (-0.5 * theta) // = -lambda
	    return lower_tail ? R_D_exp(_L) : (log_p ? atomic::robust_utils::R_Log1_Exp(_L) : -expm1(_L));
	}
	/* x < 0  or {x==0, f > 0} */
	return R_DT_0;
    }
    if(!R_FINITE(x))	return R_DT_1;
    
    /* This is principally for use from qnchisq */
#ifndef MATHLIB_STANDALONE
    R_CheckUserInterrupt();
#endif

    if(theta < 80) { /* use 110 for Inf, as ppois(110, 80/2, lower.tail=FALSE) is 2e-20 */
	Float ans;
	int i;
	// Have  pgamma(x,s) < x^s / Gamma(s+1) (< and ~= for small x)
	// ==> pchisq(x, f) = pgamma(x, f/2, 2) = pgamma(x/2, f/2)
	//                  <  (x/2)^(f/2) / Gamma(f/2+1) < eps
	// <==>  f/2 * log(x/2) - log(Gamma(f/2+1)) < log(eps) ( ~= -708.3964 )
	// <==>        log(x/2) < 2/f*(log(Gamma(f/2+1)) + log(eps))
	// <==> log(x) < log(2) + 2/f*(log(Gamma(f/2+1)) + log(eps))
	if(lower_tail && f > 0. &&
	   log(x) < M_LN2 + 2/f*(lgammafn(f/2. + 1) + _dbl_min_exp)) {
	    // all  pchisq(x, f+2*i, lower_tail, FALSE), i=0,...,110 would underflow to 0.
	    // ==> work in log scale
	    Float lambda = 0.5 * theta; // < 40
	    Float sum, sum2, pr = -lambda, log_lam = log(lambda);
	    sum = sum2 = ML_NEGINF;
	    /* we need to renormalize here: the result could be very close to 1 */
	    for(i = 0; i < 110;  pr += log_lam - log(++i)) {
		sum2 = logspace_add(sum2, pr);
		sum  = logspace_add((Float)sum , (Float)(pr + pchisq(x, (Float)(f+2*i), (int)lower_tail, (int)TRUE)));
		if (sum2 >= -1e-15) /*<=> EXP(sum2) >= 1-1e-15 */ break;
	    }
	    ans = sum - sum2;
#ifdef DEBUG_pnch
	    REprintf("pnchisq(x=%g, f=%g, th.=%g); th. < 80, logspace: i=%d, ans=(sum=%g)-(sum2=%g)\n",
		     x,f,theta, i, (Float)sum, (Float)sum2);
#endif
	    return (Float) (log_p ? ans : exp(ans));
	}
	else {
	    Float lambda = 0.5 * theta; // < 40
	    Float sum = 0, sum2 = 0, pr = exp(-lambda); // does this need a feature test?
	    /* we need to renormalize here: the result could be very close to 1 */
	    for(i = 0; i < 110;  pr *= lambda/++i) {
		// pr == exp(-lambda) lambda^i / i!  ==  dpois(i, lambda)
		sum2 += pr;
		// pchisq(*, i, *) is  strictly decreasing to 0 for lower_tail=TRUE
		//                 and strictly increasing to 1 for lower_tail=FALSE
		sum += pr * pchisq(x, (Float)(f+2*i), (int)lower_tail, (int)FALSE);
		if (sum2 >= 1-1e-15) break;
	    }
	    ans = sum/sum2;
#ifdef DEBUG_pnch
	    REprintf("pnchisq(x=%g, f=%g, theta=%g); theta < 80: i=%d, sum=%g, sum2=%g\n",
		     x,f,theta, i, (Float)sum, (Float)sum2);
#endif
	    return (Float) (log_p ? log(ans) : ans);
	}
    } // if(theta < 80)

    // else: theta == ncp >= 80 --------------------------------------------
#ifdef DEBUG_pnch
    REprintf("pnchisq(x=%g, f=%g, theta=%g >= 80): ",x,f,theta);
#endif
    // Series expansion ------- FIXME: log_p=TRUE, lower_tail=FALSE only applied at end ==> underflow

    lam = .5 * theta; // = lambda = ncp/2
    lamSml = (-lam < _dbl_min_exp);
    if(lamSml) {
	// originally error: "non centrality parameter too large for current algorithm"
        u = 0;
        lu = -lam;/* == ln(u) */
        l_lam = log(lam);
    } else {
	u = exp(-lam);
    }

    /* evaluate the first term */
    v = u;
    x2 = .5 * x;
    f2 = .5 * f;
    f_x_2n = f - x;

#ifdef DEBUG_pnch
    REprintf("-- v=exp(-th/2)=%g, x/2= %g, f/2= %g\n",v,x2,f2);
#endif

    if(f2 * DBL_EPSILON > 0.125 && /* very large f and x ~= f: probably needs */
       fabs(t = x2 - f2) <         /* another algorithm anyway */
       sqrt(DBL_EPSILON) * f2) {
	/* evade cancellation error */
	/* t = exp((1 - t)*(2 - t/(f2 + 1))) / sqrt(2*M_PI*(f2 + 1));*/
        lt = (1 - t)*(2 - t/(f2 + 1)) - M_LN_SQRT_2PI - 0.5 * log(f2 + 1);
#ifdef DEBUG_pnch
	REprintf(" (case I) ==> ");
#endif
    }
    else {
	/* Usual case 2: careful not to overflow .. : */
	lt = f2*log(x2) -x2 - lgammafn(f2 + 1);
    }
#ifdef DEBUG_pnch
    REprintf(" lt= %g", lt);
#endif

    tSml = (lt < _dbl_min_exp);
    if(tSml) {
#ifdef DEBUG_pnch
	REprintf(" is very small\n");
#endif
	if (x > f + theta +  5* sqrt( 2*(f + 2*theta))) {
	    /* x > E[X] + 5* sigma(X) */
	    return R_DT_1; /* FIXME: could be more accurate than 0. */
	} /* else */
	l_x = log(x);
	ans = term = 0.; t = 0;
    }
    else {
	t = exp(lt);
#ifdef DEBUG_pnch
 	REprintf(", t=exp(lt)= %g\n", t);
#endif
	ans = term = (Float) (v * t);
    }

    for (n = 1, f_2n = f + 2., f_x_2n += 2.; n <= itrmax ; n++, f_2n += 2, f_x_2n += 2) {
#ifdef DEBUG_pnch_n
	if(n % 1000 == 0)
	    REprintf("\n _OL_: n=%d,  f_x_2n = %g", n);
	else REprintf(n % 100 == 0 ? ".\n" : ".");
#endif
#ifndef MATHLIB_STANDALONE
	if(n % 1000 == 0) R_CheckUserInterrupt();
#endif
	/* f_2n    === f + 2*n
	 * f_x_2n  === f - x + 2*n   > 0  <==> (f+2n)  >   x */
	if (f_x_2n > 0) {
	    /* find the error bound and check for convergence */

	    bound = (Float) (t * x / f_x_2n);
#ifdef DEBUG_pnch_n
	    if(n % 1000 == 0)
		REprintf("\n L10: n=%d; term, ans = %g, %g; bound= %g",
			 n, term, ans, bound);
#endif
	    is_r = FALSE;
	    /* convergence only if BOTH absolute and relative error < 'bnd' */
	    if (((is_b = (bound <= errmax)) &&
                 (is_r = (term <= reltol * ans))))
            {
#ifdef DEBUG_pnch
		REprintf("BREAK from for(n=1 ..): n=%d; bound= %g %s; term=%g, rel.err= %g %s\n",
			 n,
			 bound, (is_b ? "<= errmax" : ""), term,
			 term/ans, (is_r ? "<= reltol" : ""));
#endif
		break; /* out completely */
            }
	}

	/* evaluate the next term of the */
	/* expansion and then the partial sum */

        if(lamSml) {
            lu += l_lam - log(n); /* u = u* lam / n */
            if(lu >= _dbl_min_exp) {
		/* no underflow anymore ==> change regime */
#ifdef DEBUG_pnch_n
                REprintf(" n=%d; nomore underflow in u = exp(lu) ==> change\n",
			 n);
#endif
                v = u = exp(lu); /* the first non-0 'u' */
                lamSml = FALSE;
            }
        } else {
	    u *= lam / n;
	    v += u;
	}
	if(tSml) {
            lt += l_x - log(f_2n);/* t <- t * (x / f2n) */
            if(lt >= _dbl_min_exp) {
		/* no underflow anymore ==> change regime */
#ifdef DEBUG_pnch
                REprintf("  n=%d; nomore underflow in t = exp(lt) ==> change\n", n);
#endif
                t = exp(lt); /* the first non-0 't' */
                tSml = FALSE;
            }
        } else {
	    t *= x / f_2n;
	}
        if(!lamSml && !tSml) {
	    term = (Float) (v * t);
	    ans += term;
	}

    } /* for(n ...) */

    if (n > itrmax) {
	MATHLIB_WARNING4(_("pnchisq(x=%g, f=%g, theta=%g, ..): not converged in %d iter."),
			 x, f, theta, itrmax);
    }
#ifdef DEBUG_pnch
    REprintf("\n == L_End: n=%d; term= %g; bound=%g: ans=%Lg\n",
	     n, term, bound, ans);
#endif
    Float dans = (Float) ans;
    return R_DT_val(dans);
}



template<class Float>
Float pnchisq(Float x, Float df, Float ncp, int lower_tail, int log_p)
{
    Float ans;
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(df) || ISNAN(ncp))
	return x + df + ncp;
    if (!R_FINITE(df) || !R_FINITE(ncp))
	ML_WARN_return_NAN;
#endif

    if (df < 0. || ncp < 0.) ML_WARN_return_NAN;

    ans = pnchisq_raw((Float)x, (Float)df, (Float)ncp, (Float)1e-12, (Float)(8*DBL_EPSILON), 1000000, (int)lower_tail, (int)log_p);

    if (x <= 0. || x == ML_POSINF)
	return ans; // because it's perfect

    if(ncp >= 80) {
	if(lower_tail) {
	  ans = myfmin2(ans, (Float)R_D__1);  /* e.g., pchisq(555, 1.01, ncp = 80) */
	} else { /* !lower_tail */
	    /* since we computed the other tail cancellation is likely */
	    // FIXME: There are cases where  ans == 0. if(!log_p) is perfect
	    if(ans < (log_p ? (-10. * M_LN10) : 1e-10)) ML_WARNING(ME_PRECISION, "pnchisq");
	    if(!log_p && ans < 0.) ans = 0.;  /* Precaution PR#7099 */
	}
    }
    /* MM: the following "hack" from c51179 (<--> PR#14216, by Jerry Lewis)
     * -- is "kind of ok" ... but potentially suboptimal: we do  log1p(- p(*, <other tail>, log=FALSE)),
     *    but that  p(*, log=FALSE) may already be an exp(.) or even expm1(..)
     *   <---> "in principle"  this check should happen there, not here  */
    if (!log_p || ans < -1e-8)
	return ans;
    else { // log_p (==> ans <= 0) &&  -1e-8 <= ans <= 0
	// prob. = exp(ans) is near one: we can do better using the other tail
#ifdef DEBUG_pnch
	REprintf("   pnchisq_raw(*, log_p): ans=%g => 2nd call, other tail\n", ans);
#endif
	ans = pnchisq_raw((Float)x, (Float)df, (Float)ncp, (Float)1e-12, (Float)(8*DBL_EPSILON), 1000000, (int)!lower_tail, (int)FALSE);
	return log1p(-ans);
    }
}

  } // End fromR namespace
