package edu.mit.csail.psrg.georg.StatUtil;

import java.util.ArrayList;
import cern.jet.random.Gamma;
import cern.jet.random.Beta;
import cern.jet.random.Uniform;
import cern.jet.random.HyperGeometric;
import java.util.Arrays;

public class StatUtil {
	static cern.jet.random.engine.RandomEngine engine = new cern.jet.random.engine.MersenneTwister();
	
	public static double uniform_rnd() {
		return Uniform.staticNextDouble();
	}
	
	public static double beta_rnd(double a, double b) {
		return Beta.staticNextDouble(a,b);
	}
	
	public static void dirichlet_rnd(double[] x, double[] alpha, int vSize) {
		double sum = 0.0;
		double xt = 0.0;
		int i = 0;
		
		for(i=0;i<vSize;i++) {
			if (alpha[i] <= 0.0) {
				xt = 0.0;
			} else {
				xt = Gamma.staticNextDouble(alpha[i],1.0);
			}
			if (Double.isNaN(xt)) {
				xt = 0.0;
			}
			x[i] = xt;
			sum = sum + x[i];
		}
		
		for(i=0;i<vSize;i++) {
			x[i] = x[i]/sum;
		}
	}
	
	public static double logDirichletPDF(double[] x, double[] alpha, int vSize) {
		double alpha_sum = 0.0;
		int i = 0;
		double lik = 0.0;
		
		for (i=0;i<vSize;i++) {
			if (alpha[i] > 0.0) {
				if (x[i] > 0.0) {
					lik = lik + (alpha[i] - 1.0)*Math.log(x[i]);
				}
				lik = lik - cern.jet.stat.Gamma.logGamma(alpha[i]);
				alpha_sum = alpha_sum + alpha[i];
			}
		}
		
		lik = lik+cern.jet.stat.Gamma.logGamma(alpha_sum);
		
		return lik;
	}
	
	public static int multinomial_rnd(double[] p, int vSize) {
		double sum = 0.0;
		double mass = 0.0;
		int cc = 0;
		int i = 0;
		
		for (i=0;i<vSize;i++) {
			sum += p[i];
		}
		
		mass = Uniform.staticNextDouble() * sum;
		
		i = 0;
		while (true) {
			mass -= p[i];
			if (mass <= 0.0) break;
			i++;
			cc++;
		}
		return cc;
	}
	
	public static int find(ArrayList<Integer> x,int f) {
		int found = -1;
		int i = 0;
		
		while((i < x.size()) & (found == -1)) {
			if (x.get(i) == f) {
				found = i;
			}
			i++;
		}
		
		return found;
	}
	
	public static double median(double[] x) {
		if (x.length == 1) {
			return x[0];
		}
		
		double[] y = new double[x.length];
		System.arraycopy(x,0,y,0,x.length);
		Arrays.sort(y);
		int n = x.length;
		double m = ((double) n)/2.0;
		if (Math.floor(m) == m) {
			n = n / 2;
			return (y[n] + y[n-1])/2.0;
		} else {
			n = (n-1)/2;
			return y[n];
		}
	}
	
	public static double logNormPDF(double[] mu,double s,double[] x) {
		int i = 0;
		double l = 0.0;
		double t = 0.0;
		for (i=0;i<x.length;i++) {
			t = x[i];
			if (!Double.isNaN(t)) {
				l = l + -0.5*Math.log(2.0*Math.PI)-Math.log(s)-Math.pow(t-mu[i],2.0)/(2.0*Math.pow(s,2.0));
			}
		}
		return l;
	}
	
	// x = # observed successes in the sample
	// N = population size
	// s = # of successes in population (e.g., size of positive set in the population)
	// n = sample size
	public static double hyperGeometricCDF(int x,int N,int s,int n) {
	//	HyperGeometric h = new HyperGeometric(N,s,n,engine);
		double v = 0.0;
		int k = 0;
		for (k=0;k<=x;k++) {
		//	v += h.pdf(k);
			v += hyperGeometricPDF(k,N,s,n);
		}
		
		if (v > 1.0)
			v = 1.0;
		return v;
	}
	
	// x = # observed successes in the sample
	// N = population size
	// s = # of successes in population (e.g., size of positive set in the population)
	// n = sample size
	// works with the gamma function and in log space in an attempt to avoid
	// numerical problems
	public static double hyperGeometricPDF(int x,int N,int s,int n) {
		double xd = (double) x;
		double Nd = (double) N;
		double sd = (double) s;
		double nd = (double) n;
		double kx = cern.jet.stat.Gamma.logGamma(sd+1)-cern.jet.stat.Gamma.logGamma(x+1)-cern.jet.stat.Gamma.logGamma(sd-x+1);
		double mn = cern.jet.stat.Gamma.logGamma(Nd+1)-cern.jet.stat.Gamma.logGamma(nd+1)-cern.jet.stat.Gamma.logGamma(Nd-nd+1);
		double mknx = cern.jet.stat.Gamma.logGamma(Nd-sd+1)-cern.jet.stat.Gamma.logGamma(nd-x+1)-cern.jet.stat.Gamma.logGamma(Nd-sd-(nd-x)+1);
		return Math.exp(kx + mknx - mn);
	}
	
	public static int[] randPerm(int n) {
		int[] order = new int[n];
		ArrayList<Integer> nm = new ArrayList<Integer>();
		int i = 0;
		for (i=0;i<n;i++) {
			nm.add(i);
		}
		
		int c = 0;
		for (i=n-1;i>=0;i--) {
			if (i > 0) {
				c = Uniform.staticNextIntFromTo(0,i);
			} else {
				c = 0;
			}
			order[i] = nm.get(c);
			nm.remove(c);
		}
		return order;
		//Uniform.staticNextIntFromTo(0,numClusters-1);
	}
}
