package edu.mit.csail.psrg.georg.IGMM;

import java.io.Serializable;

import edu.mit.csail.psrg.georg.SampleUtil.LogDensityFnARSCR;
import edu.mit.csail.psrg.georg.SpecialMathFn.Digamma;

public class BetaPosteriorSampleFn2 implements LogDensityFnARSCR, Serializable {
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
		for (i=0;i<numClusters;i++) {
			s2 = 1.0/Math.pow(s[i],2.0);
			hc1 = hc1 + Math.log(w*s2)/2.0;
			hc2 = hc2 -(w/2.0)*s2;
		}
		Q = (double) numClusters;
	}
	
	public double eval(double z, double[] dh) {
		double ht = 0.0;
		double h = 0.0;
		
		ht = z - Q*cern.jet.stat.Gamma.logGamma(Math.exp(z)/2.0);
		ht = ht + ((Q*Math.exp(z) - 3.0)/2.0)*(z - Math.log(2)) - (1.0/(2.0*Math.exp(z)));
		ht = ht + Math.exp(z)*(hc1+hc2);
		h = ht;
			
		ht = 1.0 - Q*(Math.exp(z)/2.0)*Digamma.eval(Math.exp(z)/2.0);
		ht = ht + Q*(Math.exp(z)/2.0)*(z-Math.log(2)) + (Q*Math.exp(z) - 3.0)/2.0 + 1.0/(2.0*Math.exp(z));
		ht = ht + Math.exp(z)*(hc1+hc2);
		dh[0] = ht;

		return h;
	}

}
