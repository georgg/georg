package edu.mit.csail.psrg.georg.SampleUtil;

public class GaussAdaptSampleFn2 implements LogDensityFnARSCR {
	double mu = 0.0;
	double s = 1.0;
	
	public GaussAdaptSampleFn2(double mm,double ss) {
		mu = mm;
		s = ss;		
	}
	
	public double eval(double x, double[] dh) {
		double v = -0.5*Math.log(2*Math.PI)-Math.log(s)-Math.pow(x-mu,2.0)/(2.0*Math.pow(s,2.0));
		double dv = -(x-mu)/Math.pow(s,2.0);
		dh[0] = dv;
		return v;
	}
}
