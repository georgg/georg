package edu.mit.csail.psrg.georg.SampleUtil;

import java.util.Vector;

public class GaussianAdaptSampleFn extends AdaptSampleFn {
	double mu = 0.0;
	double s = 0.0;
	
	GaussianAdaptSampleFn(double m,double ss) {
		mu = m;
		s = ss;
	}
	
	public void eval(Vector<Double> x, Vector<Double> h, Vector<Double> dh) {
		int i = 0;
		
		for (i=0;i<x.size();i++) {
			h.set(i,-0.5*Math.log(2*Math.PI)-Math.log(s)-Math.pow(x.get(i)-mu,2.0)/(2.0*Math.pow(s,2.0)));
			dh.set(i,-(x.get(i)-mu)/Math.pow(s,2.0));
		}
	}
}
