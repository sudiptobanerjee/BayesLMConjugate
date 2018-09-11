import JAMAJniLite.*;
import BayesLMConjugate.*;
import java.io.*;
import java.util.zip.GZIPInputStream;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Locale;


public final class Test {
	private Test() {}
	public static void main(String[] args) {
        //
        //Load data
        //
        int matrix_layout = Matrix.LAYOUT.RowMajor;
        final int nrow = 1342;
        final int ncol = 7;
        Matrix Data = new Matrix(nrow, ncol);
        String filename = "../data/FORMGMT.dat.tab";
        File file = new File(filename); 
        fileinput fi = new fileinput();
        fi.readmatrix(file, Data);
        //
        // Prepare the matrices and other parameters
        //
        double a = 2, b = 2;
        final int nsam = 5000;
        Matrix Y = Data.getMatrix(0,(nrow-1),0,0);
        Matrix X = Data.getMatrix(0,(nrow-1),1,6);
        int p = X.getColumnDimension();
        Matrix mub = new Matrix(p,1);//μβ, p*1
        Matrix Vbinv = new Matrix(p,p);//p*p
        double[][] result;
        //
        //Doing Bayesian Conjugate Linear Regression
        //
        BayesLMC Bayes = new BayesLMC(nsam, Y, X, mub, Vbinv, a, b);
        result = Bayes.getResult();
        Matrix Result = new Matrix(result);
        printMatrix("Beta and sigma square is = ", matrix_layout, Result.getArray(), Result.getRowDimension(), Result.getColumnDimension());
    }

    //
    /* Print the matrix X */
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
}
