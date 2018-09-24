#define MATHLIB_STANDALONE

#include <jni.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <Rmath.h>

JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_norm_1rand
  (JNIEnv *env, jclass obj)
{
	return norm_rand();
}

JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_unif_1rand
  (JNIEnv *env, jclass obj)
{
	return unif_rand();
}

JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_exp_1rand
  (JNIEnv *env, jclass obj)
{
	return exp_rand();
}

JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_dnorm
  (JNIEnv *env, jclass obj, jdouble x, jdouble mu, jdouble sigma, jint give_log)
{
	return dnorm (x, mu, sigma, give_log);
}

JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_pnorm
  (JNIEnv *env, jclass obj, jdouble x, jdouble mu, jdouble sigma, jint lower_tail, jint log_p)
{
	return pnorm (x, mu, sigma, lower_tail, log_p);
}

JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_qnorm
  (JNIEnv *env, jclass obj, jdouble p, jdouble mu, jdouble sigma, jint lower_tail, jint log_p)
{
	return qnorm (p, mu, sigma, lower_tail, log_p);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_rnorm
  (JNIEnv *env, jclass obj, jdouble mu, jdouble sigma)
{
	return rnorm (mu, sigma);
}


JNIEXPORT void JNICALL Java_java_1Rmath_jniRmath_pnorm_1both
  (JNIEnv *env, jclass obj, jdouble x, jdoubleArray jcum, jdoubleArray jccum, jint i_tail, jint log_p)
{
  	double *cum = (*env)-> GetDoubleArrayElements (env, jcum, NULL);
	double *ccum = (*env)-> GetDoubleArrayElements (env, jccum, NULL);

	pnorm_both (x, cum, ccum, i_tail, log_p);

  	(*env)-> ReleaseDoubleArrayElements (env, jcum, cum, 0);
  	(*env)-> ReleaseDoubleArrayElements (env, jccum, ccum, 0);
}

JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_dunif
  (JNIEnv *env, jclass obj, jdouble x, jdouble a, jdouble b, jint give_log)
{
	return dunif (x, a, b, give_log);
}
JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_punif
  (JNIEnv *env, jclass obj, jdouble x, jdouble a, jdouble b, jint lower_tail, jint log_p)
{
	return punif (x, a, b, lower_tail, log_p);
}

JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_qunif
  (JNIEnv *env, jclass obj, jdouble p, jdouble a, jdouble b, jint lower_tail, jint log_p)
{
	return qunif (p, a, b, lower_tail, log_p);
}
JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_runif
  (JNIEnv *env, jclass obj, jdouble a, jdouble b)
{
	return runif (a, b);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_dgamma
  (JNIEnv *env, jclass obj, jdouble x, jdouble shape, jdouble scale, jint give_log)
{
	return dgamma (x, shape, scale, give_log);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_pgamma
  (JNIEnv *env, jclass obj, jdouble x, jdouble alpha, jdouble scale, jint lower_tail, jint log_p)
{
	return pgamma (x, alpha, scale, lower_tail, log_p);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_qgamma
  (JNIEnv *env, jclass obj, jdouble p, jdouble alpha, jdouble scale, jint lower_tail, jint log_p)
{
	return qgamma (p, alpha, scale, lower_tail, log_p);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_rgamma
  (JNIEnv *env, jclass obj, jdouble a, jdouble scale)
{
	return rgamma (a, scale);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_dbeta
  (JNIEnv *env, jclass obj, jdouble x, jdouble a, jdouble b, jint give_log)
{
	return dbeta (x, a, b, give_log);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_pbeta
  (JNIEnv *env, jclass obj, jdouble x, jdouble a, jdouble b, jint lower_tail, jint log_p)
{
	return pbeta (x, a, b, lower_tail, log_p);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_qbeta
  (JNIEnv *env, jclass obj, jdouble alpha, jdouble p, jdouble q, jint lower_tail, jint log_p)
{
	return qbeta (alpha, p, q, lower_tail, log_p);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_rbeta
  (JNIEnv *env, jclass obj, jdouble aa, jdouble bb)
{
	return rbeta (aa, bb);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_dlnorm
  (JNIEnv *env, jclass obj, jdouble x, jdouble meanlog, jdouble sdlog, jint give_log)
{
	return dlnorm (x, meanlog, sdlog, give_log);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_plnorm
  (JNIEnv *env, jclass obj, jdouble x, jdouble meanlog, jdouble sdlog, jint lower_tail, jint log_p)
{
	return plnorm (x, meanlog, sdlog, lower_tail, log_p);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_qlnorm
  (JNIEnv *env, jclass obj, jdouble p, jdouble meanlog, jdouble sdlog, jint lower_tail, jint log_p)
{
	return qlnorm (p, meanlog, sdlog, lower_tail, log_p);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_rlnorm
  (JNIEnv *env, jclass obj, jdouble meanlog, jdouble sdlog)
{
	return rlnorm (meanlog, sdlog);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_dchisq
  (JNIEnv *env, jclass obj, jdouble x, jdouble df, jint give_log)
{
	return dchisq (x, df, give_log);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_pchisq
  (JNIEnv *env, jclass obj, jdouble x, jdouble df, jint lower_tail, jint log_p)
{
	return pchisq (x, df, lower_tail, log_p);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_qchisq
  (JNIEnv *env, jclass obj, jdouble p, jdouble df, jint lower_tail, jint log_p)
{
	return qchisq (p, df, lower_tail, log_p);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_rchisq
  (JNIEnv *env, jclass obj, jdouble df)
{
	return rchisq (df);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_dnchisq
  (JNIEnv *env, jclass obj, jdouble x, jdouble df, jdouble ncp, jint give_log)
{
	return dnchisq (x, df, ncp, give_log);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_pnchisq
  (JNIEnv *env, jclass obj, jdouble x, jdouble df, jdouble ncp, jint lower_tail, jint log_p)
{
	return pnchisq (x, df, ncp, lower_tail, log_p);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_qnchisq
  (JNIEnv *env, jclass obj, jdouble p, jdouble df, jdouble ncp, jint lower_tail, jint log_p)
{
	return qnchisq (p, df, ncp, lower_tail, log_p);
}	


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_rnchisq
  (JNIEnv *env, jclass obj, jdouble df, jdouble lambda)
{
	return rnchisq (df, lambda);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_df
  (JNIEnv *env, jclass obj, jdouble x, jdouble m, jdouble n, jint give_log)
{
	return df (x, m, n, give_log);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_pf
  (JNIEnv *env, jclass obj, jdouble d, jdouble df1, jdouble df2, jint lower_tail, jint log_p)
{
	return pf (d, df1, df2, lower_tail, log_p);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_qf
  (JNIEnv *env, jclass obj, jdouble p, jdouble df1, jdouble df2, jint lower_tail, jint log_p)
{
	return qf (p, df1, df2, lower_tail, log_p);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_rf
  (JNIEnv *env, jclass obj, jdouble n1, jdouble n2)
{
	return rf (n1, n2);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_dt
  (JNIEnv *env, jclass obj, jdouble x, jdouble n, jint give_log)
{
	return dt (x, n, give_log);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_pt
  (JNIEnv *env, jclass obj, jdouble x, jdouble n, jint lower_tail, jint log_p)
{
	return pt (x, n, lower_tail, log_p);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_qt
  (JNIEnv *env, jclass obj, jdouble p, jdouble ndf, jint lower_tail, jint log_p)
{
	return qt (p, ndf, lower_tail, log_p);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_rt
  (JNIEnv *env, jclass obj, jdouble df)
{
	return rt (df);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_dbinom
  (JNIEnv *env, jclass obj, jdouble x, jdouble n, jdouble p, jint give_log)
{
	return dbinom (x, n, p, give_log);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_pbinom
  (JNIEnv *env, jclass obj, jdouble x, jdouble n, jdouble p, jint lower_tail, jint log_p)
{
	return pbinom (x, n, p, lower_tail, log_p);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_qbinom
  (JNIEnv *env, jclass obj, jdouble p, jdouble n, jdouble pr, jint lower_tail, jint log_p)
{
	return qbinom (p, n, pr, lower_tail, log_p);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_rbinom
  (JNIEnv *env, jclass obj, jdouble nin, jdouble pp)
{	
	return rbinom (nin, pp);
}


JNIEXPORT void JNICALL Java_java_1Rmath_jniRmath_rmultinom
  (JNIEnv *env, jclass obj, jint n, jdoubleArray jprob, jint k, jintArray jrN)
{
  	double *prob = (*env)-> GetDoubleArrayElements (env, jprob, NULL);
	int *rN = (*env)-> GetIntArrayElements (env, jrN, NULL);

	rmultinom (n, prob, k, rN);

  	(*env)-> ReleaseDoubleArrayElements (env, jprob, prob, 0);
  	(*env)-> ReleaseIntArrayElements (env, jrN, rN, 0);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_dcauchy
  (JNIEnv *env, jclass obj, jdouble x, jdouble location, jdouble scale, jint give_log)
{
	return dcauchy (x, location, scale, give_log);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_pcauchy
  (JNIEnv *env, jclass obj, jdouble x, jdouble location, jdouble scale, jint lower_tail, jint log_p)
{
	return pcauchy (x, location, scale, lower_tail, log_p);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_qcauchy
  (JNIEnv *env, jclass obj, jdouble p, jdouble location, jdouble scale, jint lower_tail, jint log_p)
{
	return qcauchy (p, location, scale, lower_tail, log_p);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_rcauchy
  (JNIEnv *env, jclass obj, jdouble location, jdouble scale)
{
	return rcauchy (location, scale);
}	


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_dexp
  (JNIEnv *env, jclass obj, jdouble x, jdouble scale, jint give_log)
{
	return dexp (x, scale, give_log);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_pexp
  (JNIEnv *env, jclass obj, jdouble x, jdouble scale, jint lower_tail, jint log_p)
{
	return pexp (x, scale, lower_tail, log_p);
}	


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_qexp
  (JNIEnv *env, jclass obj, jdouble p, jdouble scale, jint lower_tail, jint log_p)
{
	return qexp (p, scale, lower_tail, log_p);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_rexp
  (JNIEnv *env, jclass obj, jdouble scale)
{
	return rexp (scale);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_dgeom
  (JNIEnv *env, jclass obj, jdouble x, jdouble p, jint give_log)
{
	return dgeom (x, p, give_log);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_pgeom
  (JNIEnv *env, jclass obj, jdouble x, jdouble p, jint lower_tail, jint log_p)
{
	return pgeom (x, p, lower_tail, log_p);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_qgeom
  (JNIEnv *env, jclass obj, jdouble p, jdouble prob, jint lower_tail, jint log_p)
{
	return qgeom (p, prob, lower_tail, log_p);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_rgeom
  (JNIEnv *env, jclass obj, jdouble p)
{
	return rgeom (p);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_dhyper
  (JNIEnv *env, jclass obj, jdouble x, jdouble r, jdouble b, jdouble n, jint give_log)
{
	return dhyper (x, r, b, n, give_log);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_phyper
  (JNIEnv *env, jclass obj, jdouble x, jdouble NR, jdouble NB, jdouble n, jint lower_tail, jint log_p)
{
	return phyper (x, NR, NB, n, lower_tail, log_p);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_qhyper
  (JNIEnv *env, jclass obj, jdouble p, jdouble NR, jdouble NB, jdouble n, jint lower_tail, jint log_p)
{
	return qhyper (p, NR, NB, n, lower_tail, log_p);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_rhyper
  (JNIEnv *env, jclass obj, jdouble nn1in, jdouble nn2in, jdouble kkin)
{
	return rhyper (nn1in, nn2in, kkin);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_dnbinom
  (JNIEnv *env, jclass obj, jdouble x, jdouble size, jdouble prob, jint give_log)
{
	return dnbinom (x, size, prob, give_log);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_pnbinom
  (JNIEnv *env, jclass obj, jdouble x, jdouble size, jdouble prob, jint lower_tail, jint log_p)
{
	return pnbinom (x, size, prob, lower_tail, log_p);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_qnbinom
  (JNIEnv *env, jclass obj, jdouble p, jdouble size, jdouble prob, jint lower_tail, jint log_p)
{
	return qnbinom (p, size, prob, lower_tail, log_p);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_rnbinom
  (JNIEnv *env, jclass obj, jdouble size, jdouble prob)
{
	return rnbinom (size, prob);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_dpois
  (JNIEnv *env, jclass obj, jdouble x, jdouble lambda, jint give_log)
{
	return dpois (x, lambda, give_log);
}	


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_ppois
  (JNIEnv *env, jclass obj, jdouble x, jdouble lambda, jint lower_tail, jint log_p)
{
	return ppois (x, lambda, lower_tail, log_p);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_qpois
  (JNIEnv *env, jclass obj, jdouble p, jdouble lambda, jint lower_tail, jint log_p)
{
	return qpois (p, lambda, lower_tail, log_p);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_rpois
  (JNIEnv *env, jclass obj, jdouble mu)
{	
	return rpois (mu);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_dweibull
  (JNIEnv *env, jclass obj, jdouble x, jdouble shape, jdouble scale, jint give_log)
{	
	return dweibull (x, shape, scale, give_log);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_pweibull
  (JNIEnv *env, jclass obj, jdouble x, jdouble shape, jdouble scale, jint lower_tail, jint log_p)
{	
	return pweibull (x, shape, scale, lower_tail, log_p);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_qweibull
  (JNIEnv *env, jclass obj, jdouble p, jdouble shape, jdouble scale, jint lower_tail, jint log_p)
{	
	return qweibull (p, shape, scale, lower_tail, log_p);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_rweibull
  (JNIEnv *env, jclass obj, jdouble shape, jdouble scale)
{
	return rweibull (shape, scale);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_dlogis
  (JNIEnv *env, jclass obj, jdouble x, jdouble location, jdouble scale, jint give_log)
{
	return dlogis (x, location, scale, give_log);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_plogis
  (JNIEnv *env, jclass obj, jdouble x, jdouble location, jdouble scale, jint lower_tail, jint log_p)
{
	return plogis (x, location, scale, lower_tail, log_p);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_qlogis
  (JNIEnv *env, jclass obj, jdouble p, jdouble location, jdouble scale, jint lower_tail, jint log_p)
{
	return qlogis (p, location, scale, lower_tail, log_p);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_rlogis
  (JNIEnv *env, jclass obj, jdouble location, jdouble scale)
{
	return rlogis (location, scale);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_dnbeta
  (JNIEnv *env, jclass obj, jdouble x, jdouble a, jdouble b, jdouble ncp, jint give_log)
{
	return dnbeta (x, a, b, ncp, give_log);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_pnbeta
  (JNIEnv *env, jclass obj, jdouble x, jdouble a, jdouble b, jdouble ncp, jint lower_tail, jint log_p)
{
	return pnbeta (x, a, b, ncp, lower_tail, log_p);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_qnbeta
  (JNIEnv *env, jclass obj, jdouble p, jdouble a, jdouble b, jdouble ncp, jint lower_tail, jint log_p)
{	
	return qnbeta (p, a, b, ncp, lower_tail, log_p);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_rnbeta
  (JNIEnv *env, jclass obj, jdouble a, jdouble b, jdouble ncp)
{
	return rnbeta (a, b, ncp);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_dnf
  (JNIEnv *env, jclass obj, jdouble x, jdouble df1, jdouble df2, jdouble ncp, jint give_log)
{	
	return dnf (x, df1, df2, ncp, give_log);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_pnf
  (JNIEnv *env, jclass obj, jdouble x, jdouble df1, jdouble df2, jdouble ncp, jint lower_tail, jint log_p)
{
	return pnf (x, df2, df2, ncp, lower_tail, log_p);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_qnf
  (JNIEnv *env, jclass obj, jdouble p, jdouble df1, jdouble df2, jdouble ncp, jint lower_tail, jint log_p)
{
	return qnf (p, df1, df2, ncp, lower_tail, log_p);
}	


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_dnt
  (JNIEnv *env, jclass obj, jdouble x, jdouble df, jdouble ncp, jint give_log)
{	
	return dnt (x, df, ncp, give_log);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_pnt
  (JNIEnv *env, jclass obj, jdouble t, jdouble df, jdouble ncp, jint lower_tail, jint log_p)
{	
	return pnt (t, df, ncp, lower_tail, log_p);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_qnt
  (JNIEnv *env, jclass obj, jdouble p, jdouble df, jdouble ncp, jint lower_tail, jint log_p)
{
	return qnt (p, df, ncp, lower_tail, log_p);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_ptukey
  (JNIEnv *env, jclass obj, jdouble q, jdouble rr, jdouble cc, jdouble df, jint lower_tail, jint log_p)
{
	return ptukey (q, rr, cc, df, lower_tail, log_p);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_qtukey
  (JNIEnv *env, jclass obj, jdouble p, jdouble rr, jdouble cc, jdouble df, jint lower_tail, jint log_p)
{
	return qtukey (p, rr, cc, df, lower_tail, log_p);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_dwilcox
  (JNIEnv *env, jclass obj, jdouble x, jdouble m, jdouble n, jint give_log)
{
	return dwilcox (x, m, n, give_log);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_pwilcox
  (JNIEnv *env, jclass obj, jdouble q, jdouble m, jdouble n, jint lower_tail, jint log_p)
{
	return pwilcox (q, m, n, lower_tail, log_p);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_qwilcox
  (JNIEnv *env, jclass obj, jdouble x, jdouble m, jdouble n, jint lower_tail, jint log_p)
{
	return qwilcox (x, m, n, lower_tail, log_p);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_rwilcox
  (JNIEnv *env, jclass obj, jdouble m, jdouble n)
{
	return rwilcox (m, n);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_dsignrank
  (JNIEnv *env, jclass obj, jdouble x, jdouble n, jint give_log)
{	
	return dsignrank (x, n, give_log);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_psignrank
  (JNIEnv *env, jclass obj, jdouble x, jdouble n, jint lower_tail, jint log_p)
{
	return psignrank (x, n, lower_tail, log_p);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_qsignrank
  (JNIEnv *env, jclass obj, jdouble x, jdouble n, jint lower_tail, jint log_p)
{
	return qsignrank (x, n, lower_tail, log_p);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_rsignrank
  (JNIEnv *env, jclass obj, jdouble n)
{
	return rsignrank (n);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_gammafn
  (JNIEnv *env, jclass obj, jdouble x)
{
	return gammafn (x);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_lgammafn
  (JNIEnv *env, jclass obj, jdouble x)
{
	return lgammafn (x);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_lgammafn_1sign
  (JNIEnv *env, jclass obj, jdouble x, jintArray jsgn)
{
	double result;
	int *sgn = (*env)-> GetIntArrayElements (env, jsgn, NULL);

	result = lgammafn_sign (x, sgn);
	
	(*env)-> ReleaseIntArrayElements (env, jsgn, sgn, 0);

	return result;	
}


JNIEXPORT void JNICALL Java_java_1Rmath_jniRmath_dpsifn
  (JNIEnv *env, jclass obj, jdouble x, jint n, jint kode, jint m, jdoubleArray jans, jintArray jnz, jintArray jierr)
{
	double *ans = (*env)-> GetDoubleArrayElements (env, jans, NULL);	
	int *nz = (*env)-> GetIntArrayElements (env, jnz, NULL);
	int *ierr = (*env)-> GetIntArrayElements (env, jierr, NULL);

	dpsifn (x, n, kode, m, ans, nz, ierr);

	(*env)-> ReleaseDoubleArrayElements (env, jans, ans, 0);
	(*env)-> ReleaseIntArrayElements (env, jnz, nz, 0);
	(*env)-> ReleaseIntArrayElements (env, jierr, ierr, 0);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_psigamma
  (JNIEnv *env, jclass obj, jdouble x, jdouble deriv)
{
	return psigamma (x, deriv);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_digamma
  (JNIEnv *env, jclass obj, jdouble x)
{
	return digamma (x);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_trigamma
  (JNIEnv *env, jclass obj, jdouble x)
{
	return trigamma (x);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_tetragamma
  (JNIEnv *env, jclass obj, jdouble x)
{
	return tetragamma (x);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_pentagamma
  (JNIEnv *env, jclass obj, jdouble x)
{
	return pentagamma (x);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_beta
  (JNIEnv *env, jclass obj, jdouble a, jdouble b)
{
	return beta (a, b);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_lbeta
  (JNIEnv *env, jclass obj, jdouble a, jdouble b)
{	
	return lbeta (a, b);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_choose
  (JNIEnv *env, jclass obj, jdouble n, jdouble k)
{
	return choose (n, k);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_lchoose
  (JNIEnv *env, jclass obj, jdouble n, jdouble k)
{
	return lchoose (n, k);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_bessel_1i
  (JNIEnv *env, jclass obj, jdouble x, jdouble alpha, jdouble expo)
{
	return bessel_i (x, alpha, expo);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_bessel_1j
  (JNIEnv *env, jclass obj, jdouble x, jdouble alpha)
{
	return bessel_j (x, alpha);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_bessel_1k
  (JNIEnv *env, jclass obj, jdouble x, jdouble alpha, jdouble expo)
{
	return bessel_k (x, alpha, expo);
}


JNIEXPORT jdouble JNICALL Java_java_1Rmath_jniRmath_bessel_1y
  (JNIEnv *env, jclass obj, jdouble x, jdouble alpha)
{
	return bessel_y (x, alpha);
}


