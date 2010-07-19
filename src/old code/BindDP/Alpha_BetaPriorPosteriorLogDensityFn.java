package edu.mit.csail.psrg.georg.BindDP;

import edu.mit.csail.psrg.georg.SampleUtil.LogDensityFn;
import cern.jet.stat.Gamma;

public class Alpha_BetaPriorPosteriorLogDensityFn implements LogDensityFn {
	double s1 = 0.0;
	double w = 0.0;
	double s2 = 0.0;
	double M = 0.0;
	
	public void setParams(double w_beta,double[] beta) {
		M = (double) beta.length;
		w = w_beta;
		
		int j = 0;
		s1 = 0.0;
		s2 = 0.0;
		for (j=0;j<M;j++) {
			s1 = s1 + Math.log(beta[j]);
			s2 = s2 + beta[j];
		}
		s2 = s2 * w;
	}
	
	public double eval(double[] x) {
		double alpha = x[0];
		double h = 0;
		h = -(3.0/2.0)*Math.log(alpha) - 1.0/(2.0*alpha) + (M*alpha/2.0)*Math.log(alpha*w/2.0) - M*Gamma.logGamma(alpha/2.0);
		h = h + (alpha/2.0)*(s1 - s2);
		return h;
	}
	
	

}
