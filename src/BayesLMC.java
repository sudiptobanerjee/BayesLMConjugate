package BayesLMConjugate;

import JAMAJniLite.*;
import java_Rmath.*;
import java.io.*;
import java.util.zip.GZIPInputStream;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Locale;
import java.lang.*;


public class BayesLMC implements java.io.Serializable{
    /* ------------------------
    * Class variables
    * ------------------------ */
    double a, b, as, bs;//scalar,a*, b*
    int n,p,nsam;//nsam means the number of random samples to generate
    double[] sigma, sigma2, sigma2recip;
    double[][] beta,temp;
    Matrix X, Y, L, M, m, mu, mub, Vbinv, Minv, Temp;
    CholeskyDecomposition Chol;
    final int matrix_layout = Matrix.LAYOUT.RowMajor;
    final int trans = Matrix.TRANSPOSE.Trans;
    final int notrans = Matrix.TRANSPOSE.NoTrans;
    /* ------------------------
    * Constructor
    * ------------------------ */
    public BayesLMC(int nsamv, Matrix Yv, Matrix Xv, Matrix mubv, Matrix Vbinvv, double av, double bv) {
    	/*Initiation*/
    	nsam = nsamv;
    	Y = Yv;
    	X = Xv;
    	mub = mubv;
    	Vbinv = Vbinvv;
    	a = av;
    	b = bv;
    	n = Y.getRowDimension();
    	p = X.getColumnDimension();
    	sigma = new double[nsam];
    	sigma2 = new double[nsam];
    	sigma2recip = new double[nsam];
    	beta = new double[nsam][p];

    	/*Calculation*/
    	m = Vbinv.times(mub);
        m = m.plus(X.times(Y,trans,notrans));//m = Vbinv.times(mub); + Xtrans.times(Vyinv.times(Y));
        Minv = Vbinv.plus(X.times(X,trans,notrans));// Minv = Vbinv + Xtrans.times(Vyinv.times(X));
        M = Minv.inverse();
        Chol = M.chol(); 
        L = Chol.getL();
        mu = M.times(m); 
        as = a + (double)1/2 * (double)n;
        Temp = mub.times(Vbinv.times(mub),trans,notrans);
        Temp = Temp.plus(Y.times(Y,trans,notrans));
        Temp = Temp.minus(m.times(M.times(m),trans,notrans));
        temp = Temp.getArray();
        bs = b + (double)1/2 * temp[0][0];//bs=b+1/2*(mubtrans.times(Vbinv.times(mub))+Ytrans.times(Vyinv.times(Y))-mtrans.times(M.times(m)))
      	}
    /* ------------------------
    * Public Methods
    * ------------------------ */
    public double[][] getResult(){
    	for(int i=0; i<nsam; i++){
    		sigma2recip[i] = jniRmath.rgamma(as, 1/bs);
    		sigma2[i] = 1/sigma2recip[i];
    		sigma[i] = Math.sqrt(sigma2[i]);               
    		beta[i] = rmvnorm(nsam, p, mu, L, sigma[i]);  
    	}
    	Matrix Sigma2 = new Matrix(sigma2,nsam);
    	Matrix Beta = new Matrix(beta);
    	Matrix Result = new Matrix(nsam,(p+1));
    	Result.setMatrix(0,(nsam-1),0,(p-1),Beta);
    	Result.setMatrix(0,(nsam-1),p,p,Sigma2);
    	return Result.getArray();
    }
    //
    /* Multivariate random number generator */
    //
    private static double[] rmvnorm(int nsam, int p, Matrix mu, Matrix L, double sigma) {
    	double[] z = new double[p];
    	for(int i=0; i<p; i++){
    		z[i] = jniRmath.rnorm(0,1);
    	}
    	Matrix Z = new Matrix(z,p);
    	mu = mu.plus(L.times(sigma).times(Z));
	    return mu.getRowPackedCopy();
   	}
    //
    //* Print the matrix X */
    //
   	private static void printMatrix(String prompt, int layout, double[][] X, int I, int J) {
    System.out.println(prompt);
    if (layout == Matrix.LAYOUT.ColMajor) {
     	for (int i=0; i<I; i++) {
      	for (int j=0; j<J; j++)
       	System.out.print("\t" + string(X[j][i]));
     	System.out.println();
   		}
	}
 	else if (layout == Matrix.LAYOUT.RowMajor){
   		for (int i=0; i<I; i++) {
    	for (int j=0; j<J; j++)
     	System.out.print("\t" + string(X[i][j]));
   		System.out.println();
 		}
	}
	else{System.out.println("** Illegal layout setting");}
	}
  	//
	/* Print the array X */
  	//
	private static void printIntArray(String prompt, int[] X, int L) {
  		System.out.println(prompt);
  		for (int i=0; i<L; i++) {
   		System.out.print("\t" + string(X[i]));
 	}
 	System.out.println();
	}
  	//
	/* Shorter string for real number */
  	//
	private static String string(double re) {
  		String s="";
  		if (re == (long)re)
   		s += (long)re;
 	else
   		s += re;
 		return s;
	}

	private static final long serialVersionUID = 1;

}
