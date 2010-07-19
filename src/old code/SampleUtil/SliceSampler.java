package edu.mit.csail.psrg.georg.SampleUtil;
import cern.jet.random.Uniform;

//Based on a Matlab implementation by Iain Murray May 2004 from Pseudo-code in David MacKay's text
// book p375
// draw N samples from distribution
// after discarding burn-in period using slice sampling. Starts at X. Step
// size(s) w. Dist need not be normalised.  S will return the samples
public class SliceSampler {
	public static void sample(int N,int burn,LogDensityFn logDist,double[] x,double[] w,double[][] S) {
		int D = x.length;
		int d = 0;
		double r = 0.0;
		double lPx = 0.0;
		double luprime = 0.0;
		double[] x_l = new double[D];
		double[] x_r = new double[D];
		double[] xprime = new double[D];
		lPx = logDist.eval(x);
		
		int i = 0;
		int z = 0;
		for(i=0;i<N+burn;i++) {
			luprime = Math.log(Uniform.staticNextDouble()) + lPx;
						
			// Setup stuff for sampling along random axis
			// (could do random direction, but this is simpler for now)
			d = (int) Math.ceil(Uniform.staticNextDouble()*((double) D))-1;
			System.arraycopy(x,0,x_l,0,D);
			System.arraycopy(x,0,x_r,0,D);
			System.arraycopy(x,0,xprime,0,D);
			
			// Create a horizontal interval (x_l,x_r) enclosing x
			r = Uniform.staticNextDouble();
			x_l[d] = x[d]-r*w[d];
			x_r[d] = x[d]+(1-r)*w[d];
			
			while(logDist.eval(x_l) > luprime) {
				x_l[d] = x_l[d] - w[d];
			}
			while (logDist.eval(x_r) > luprime) {
				x_r[d] = x_r[d] + w[d];
			}
			
			// Inner loop
			z = 0;
			while(true) {
				z++;
				xprime[d] = Uniform.staticNextDouble()*(x_r[d]-x_l[d]) + x_l[d];
				lPx = logDist.eval(xprime);
				if (lPx > luprime) {
					break; // this is the only way to leave the loop
				} else {
					// Shrink in
					if (xprime[d] > x[d]) {
						x_r[d] = xprime[d];
					} else {
						x_l[d] = xprime[d];
					}
				}
			}
			System.arraycopy(xprime,0,x,0,D);
			
			// record samples
			if (i >= burn) {
				System.arraycopy(x,0,S[i-burn],0,D);
			}
		}
	}
}
