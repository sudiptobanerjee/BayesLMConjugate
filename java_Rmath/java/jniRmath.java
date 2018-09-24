package java_Rmath;

/**CBLAS.java*/

public class jniRmath {
 private jniRmath() {}
 static {
     
    /* load library (which will contain wrapper for Rmath function.)*/
    System.loadLibrary("jni_java_Rmath");
 
 }


	/* Random Number Generators */

public static native double     norm_rand();
public static native double	unif_rand();
public static native double	exp_rand();

	/* Normal Distribution */

public static native double	dnorm(double x, double mu, double sigma, int give_log);
public static native double	pnorm(double x, double mu, double sigma, int lower_tail, int log_p);
public static native double	qnorm(double p, double mu, double sigma, int lower_tail, int log_p);
public static native double	rnorm(double mu, double sigma);
public static native void	pnorm_both(double x, double[] cum, double[] ccum, int i_tail, int log_p);/* both tails */

	/* Uniform Distribution */

public static native double	dunif(double x, double a, double b, int give_log);
public static native double	punif(double x, double a, double b, int lower_tail, int log_p);
public static native double	qunif(double p, double a, double b, int lower_tail, int log_p);
public static native double	runif(double a, double b);

	/* Gamma Distribution */

public static native double	dgamma(double x, double shape, double scale, int give_log);
public static native double	pgamma(double x, double shape, double scale, int lower_tail, int log_p);
public static native double	qgamma(double p, double shape, double scale, int lower_tail, int log_p);
public static native double	rgamma(double shape, double scale);

	/* Beta Distribution */

public static native double	dbeta(double x, double a, double b, int give_log);
public static native double	pbeta(double x, double a, double b, int lower_tail, int log_p);
public static native double	qbeta(double p, double a, double b, int lower_tail, int log_p);
public static native double	rbeta(double a, double b);

	/* Lognormal Distribution */

public static native double	dlnorm(double x, double meanlog, double sdlog, int give_log);
public static native double	plnorm(double x, double meanlog, double sdlog, int lower_tail, int log_p);
public static native double	qlnorm(double p, double meanlog, double sdlog, int lower_tail, int log_p);
public static native double	rlnorm(double meanlog, double sdlog);

	/* Chi-squared Distribution */

public static native double	dchisq(double x, double df, int give_log);
public static native double	pchisq(double x, double df, int lower_tail, int log_p);
public static native double	qchisq(double p, double df, int lower_tail, int log_p);
public static native double	rchisq(double df);

	/* Non-central Chi-squared Distribution */

public static native double	dnchisq(double x, double df, double ncp, int give_log);
public static native double	pnchisq(double x, double df, double ncp, int lower_tail, int log_p);
public static native double	qnchisq(double p, double df, double ncp, int lower_tail, int log_p);
public static native double	rnchisq(double df, double ncp);

	/* F Distibution */

public static native double	df(double x, double df1, double df2, int give_log);
public static native double	pf(double x, double df1, double df2, int lower_tail, int log_p);
public static native double	qf(double p, double df1, double df2, int lower_tail, int log_p);
public static native double	rf(double df1, double df2);

	/* Student t Distibution */

public static native double	dt(double x, double df, int give_log);
public static native double	pt(double x, double df, int lower_tail, int log_p);
public static native double	qt(double p, double df, int lower_tail, int log_p);
public static native double	rt(double df);

	/* Binomial Distribution */

public static native double	dbinom(double x, double n, double p, int give_log);
public static native double	pbinom(double x, double n, double p, int lower_tail, int log_p);
public static native double	qbinom(double p, double n, double pr, int lower_tail, int log_p);
public static native double	rbinom(double n, double p);

	/* Multnomial Distribution */

public static native void	rmultinom(int n, double[] prob, int k, int[] rN);

	/* Cauchy Distribution */

public static native double	dcauchy(double x, double location, double scale, int give_log);
public static native double	pcauchy(double x, double location, double scale, int lower_tail, int log_p);
public static native double	qcauchy(double p, double location, double scale, int lower_tail, int log_p);
public static native double	rcauchy(double location, double scale);

	/* Exponential Distribution */

public static native double	dexp(double x, double scale, int give_log);
public static native double	pexp(double x, double scale, int lower_tail, int log_p);
public static native double	qexp(double p, double scale, int lower_tail, int log_p);
public static native double	rexp(double scale);

	/* Geometric Distribution */

public static native double	dgeom(double x, double p, int give_log);
public static native double	pgeom(double x, double p, int lower_tail, int log_p);
public static native double	qgeom(double p, double pro, int lower_tail, int log_p);
public static native double	rgeom(double p);

	/* Hypergeometric Distibution */

public static native double	dhyper(double x, double r, double b, double n, int give_log);
public static native double	phyper(double x, double r, double b, double n, int lower_tail, int log_p);
public static native double	qhyper(double p, double r, double b, double n, int lower_tail, int log_p);
public static native double	rhyper(double r, double b, double n);

	/* Negative Binomial Distribution */

public static native double	dnbinom(double x, double size, double pro, int give_log);
public static native double	pnbinom(double x, double size, double pro, int lower_tail, int log_p);
public static native double	qnbinom(double p, double size, double pro, int lower_tail, int log_p);
public static native double	rnbinom(double size, double pro);


	/* Poisson Distribution */

public static native double	dpois(double x, double lambda, int give_log);
public static native double	ppois(double x, double lambda, int lower_tail, int log_p);
public static native double	qpois(double p, double lambda, int lower_tail, int log_p);
public static native double	rpois(double mu);

	/* Weibull Distribution */

public static native double	dweibull(double x, double shape, double scale, int give_log);
public static native double	pweibull(double x, double shape, double scale, int lower_tail, int log_p);
public static native double	qweibull(double p, double shape, double scale, int lower_tail, int log_p);
public static native double	rweibull(double shape, double scale);

	/* Logistic Distribution */

public static native double	dlogis(double x, double location, double scale, int give_log);
public static native double	plogis(double x, double location, double scale, int lower_tail, int log_p);
public static native double	qlogis(double p, double location, double scale, int lower_tail, int log_p);
public static native double	rlogis(double location, double scale);

	/* Non-central Beta Distribution */

public static native double	dnbeta(double x, double a, double b, double ncp, int give_log);
public static native double	pnbeta(double x, double a, double b, double ncp, int lower_tail, int log_p);
public static native double	qnbeta(double p, double a, double b, double ncp, int lower_tail, int log_p);
public static native double	rnbeta(double a, double b, double ncp);

	/* Non-central F Distribution */

public static native double  dnf(double x, double df1, double df2, double ncp, int give_log);
public static native double	pnf(double x, double df1, double df2, double ncp, int lower_tail, int log_p);
public static native double	qnf(double p, double df1, double df2, double ncp, int lower_tail, int log_p);

	/* Non-central Student t Distribution */

public static native double	dnt(double x, double df, double ncp, int give_log);
public static native double	pnt(double t, double df, double ncp, int lower_tail, int log_p);
public static native double	qnt(double p, double df, double ncp, int lower_tail, int log_p);

	/* Studentized Range Distribution */

public static native double	ptukey(double x, double rr, double cc, double df, int lower_tail, int log_p);
public static native double	qtukey(double p, double rr, double cc, double df, int lower_tail, int log_p);

	/* Wilcoxon Rank Sum Distribution */

public static native double dwilcox(double x, double m, double n, int give_log);
public static native double pwilcox(double q, double m, double n, int lower_tail, int log_p);
public static native double qwilcox(double x, double m, double n, int lower_tail, int log_p);
public static native double rwilcox(double m, double n);

	/* Wilcoxon Signed Rank Distribution */

public static native double dsignrank(double x, double n, int give_log);
public static native double psignrank(double x, double n, int lower_tail, int log_p);
public static native double qsignrank(double x, double n, int lower_tail, int log_p);
public static native double rsignrank(double n);

	/* Gamma and Related Functions */
public static native double	gammafn(double x);
public static native double	lgammafn(double x);
public static native double	lgammafn_sign(double x, int[] sgn);
public static native void    dpsifn(double x, int n, int kode, int m, double[] ans, int[] nz, int[] ierr);
public static native double	psigamma(double x, double deriv);
public static native double	digamma(double x);
public static native double	trigamma(double x);
public static native double	tetragamma(double x);
public static native double	pentagamma(double x);

public static native double	beta(double a, double b);
public static native double	lbeta(double a, double b);

public static native double	choose(double n, double k);
public static native double	lchoose(double n, double k);

	/* Bessel Functions */

public static native double	bessel_i(double x, double alpha, double expo);
public static native double	bessel_j(double x, double alpha);
public static native double	bessel_k(double x, double alpha, double expo);
public static native double	bessel_y(double x, double alpha);

    
}







