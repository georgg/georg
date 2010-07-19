package edu.mit.csail.psrg.georg.IGMM;

import java.util.Vector;

import edu.mit.csail.psrg.georg.SampleUtil.AdaptSampleFn;
import edu.mit.csail.psrg.georg.SpecialMathFn.Digamma;

// the given distribution is actually p(log(beta) | .), where z = log(beta)
// so, to transform to samples of beta, you'll need to transform via
// exp(sample)
public class BetaPosteriorSampleFn extends AdaptSampleFn {
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
	
	public void eval(Vector<Double> z, Vector<Double> h, Vector<Double> dh) {
		double ht = 0.0;
		int i = 0;
		
		for (i=0;i<z.size();i++) {
			ht = z.get(i) - Q*cern.jet.stat.Gamma.logGamma(Math.exp(z.get(i))/2.0);
			ht = ht + ((Q*Math.exp(z.get(i)) - 3.0)/2.0)*(z.get(i) - Math.log(2)) - (1.0/(2.0*Math.exp(z.get(i))));
			ht = ht + Math.exp(z.get(i))*(hc1+hc2);
			h.set(i,ht);
			
			ht = 1.0 - Q*(Math.exp(z.get(i))/2.0)*Digamma.eval(Math.exp(z.get(i))/2.0);
			ht = ht + Q*(Math.exp(z.get(i))/2.0)*(z.get(i)-Math.log(2)) + (Q*Math.exp(z.get(i)) - 3.0)/2.0 + 1.0/(2.0*Math.exp(z.get(i)));
			ht = ht + Math.exp(z.get(i))*(hc1+hc2);
			dh.set(i,ht);
		}
	}
}
