package edu.mit.csail.psrg.georg.IGMM;

import edu.mit.csail.psrg.georg.SampleUtil.LogDensityFn;
import cern.jet.stat.Gamma;

// returns the log of the beta posterior density
// (used for slice sampling)
public class BetaPosteriorLogDensityFn implements LogDensityFn {
	double Q = 0.0;
	double hc1 = 0.0;
	double hc2 = 0.0;
	double w = 0.0;
	
	public void init(double[] s,double wi,int numClusters) {
		w = wi;
		hc1 = 0.0;
		hc2 = 0.0;
		int i = 0;
		double s2 = 0.0;
		for (i=1;i<numClusters;i++) {
			s2 = 1.0/Math.pow(s[i],2.0);
			hc1 = hc1 + Math.log(w*s2)/2.0;
			hc2 = hc2 -(w/2.0)*s2;
		}
		Q = (double) numClusters;
	}
	
	public double eval(double[] x) {
		double beta = x[0];
		
		if (beta <= 0.0) {
			return -1.0/0.0;
		}
		
		double h = -Q*Gamma.logGamma(beta/2.0) + ((Q*beta-3.0)/2.0)*Math.log(beta/2.0) - (1.0/(2.0*beta));
		h = h + beta*(hc1 + hc2);
		
		return h;
	}

}
