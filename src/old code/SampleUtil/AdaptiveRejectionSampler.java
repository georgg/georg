package edu.mit.csail.psrg.georg.SampleUtil;

import java.util.Vector;
import edu.mit.csail.psrg.georg.StatUtil.*;
import cern.jet.random.Uniform;
import java.io.FileWriter;
import java.io.IOException;

//Implements the adaptive rejection sampling method of Gilks and Wild (1992)
//g_fn is a function that returns ln g(x) and d(ln g(x))/dx for a
//vector x, where g(x) is a density function (or proportional to a density
//function - i.e., an unnormalized density can be used)
//g_p are the parameters of the distribution function
//x_init is a vector of values that spans the distribution (can just be two values)
//z_min and z_max are the bounds on the distribution (-Inf and Inf if unbounded)
public class AdaptiveRejectionSampler {
	AdaptSampleFn g_fn;
	public double[] samples;
	double z_max;
	double z_min;
	Vector<Double> x = new Vector<Double>();
	Vector<Double> h = new Vector<Double>();
	Vector<Double> dh = new Vector<Double>();
	Vector<Double> z = new Vector<Double>();
	Vector<Double> cd = new Vector<Double>();
	Vector<Double> c = new Vector<Double>();
	Vector<Double> cn = new Vector<Double>();	
	double offset;
	double norm;
	int n_samples = 0;
	double MIN_DROP = 0.1;  /* minimum drop in log probability for initial points */
	double TINY = 1e-9;         /* minimum tolerance for inserting new points */
	
	public AdaptiveRejectionSampler(AdaptSampleFn g,Vector<Double> x_init,double imin,double imax,int ns) {
		g_fn = g;
		z_max = imax;
		z_min = imin;
		n_samples = ns;
		init(x_init);
		sample();
	}
	
	public AdaptiveRejectionSampler(AdaptSampleFn g,double x_guess,double scale,double imin,double imax,int ns) {
		g_fn = g;
		z_max = imax;
		z_min = imin;
		n_samples = ns;
		Vector<Double> x_init = findInit(x_guess,scale);
		init(x_init);
		sample();
	}
	
	public Vector<Double> findInit(double x_guess,double scale) {
		Vector<Double> x_init = new Vector<Double>();
		
		x_init.add(x_guess);
		x.add(x_guess);
		h.add(0.0);
		dh.add(0.0);
		
		double xt = x_guess;
		double ft = 0;
		double dht = 0.0;
		g_fn.eval(x,h,dh);
		ft = h.get(0);
		dht = dh.get(0);
		double max = ft;
		
		/* find point to the left of guess */
		while(dht < TINY || ft > max-MIN_DROP) {
			xt -= scale+(x_guess-xt);
		//	xt -= scale;
			x.clear();
			x.add(xt);
			g_fn.eval(x,h,dh);
			ft = h.get(0);
			dht = dh.get(0);
			if (ft>max) {
				max = ft;
			}
		}
		x_init.insertElementAt(xt,0);
		
		xt = x_guess;
		x.clear();
		x.add(xt);
		g_fn.eval(x,h,dh);
		ft = h.get(0);
		dht = dh.get(0);
		
		while(dht>-TINY || ft > max-MIN_DROP) {
		//	xt += scale+(xt-x_guess);
			xt += scale;
			x.clear();
			x.add(xt);
			g_fn.eval(x,h,dh);
			ft = h.get(0);
			dht = dh.get(0);
			if (ft>max) {
				max = ft;
			}
		}
		x_init.add(xt);
		
		return x_init;
	}
	
	public void init(Vector<Double> x_init) {
		int n = x_init.size();
		int j = 0;
		samples = new double[n_samples];
		
		x.clear();
		x.addAll(x_init);
		VectorUtil.zeros(h,n);
		VectorUtil.zeros(dh,n);
		
		g_fn.eval(x,h,dh);
		
	//	offset = VectorUtil.mean(h);
		offset = VectorUtil.max(h);
		VectorUtil.add(h,-offset);
		
		VectorUtil.zeros(z,n-1);
		
		for (j=0;j<n-1;j++) {
			z.set(j,(h.get(j+1)-h.get(j)-x.get(j+1)*dh.get(j+1)+x.get(j)*dh.get(j))/(dh.get(j)-dh.get(j+1)));
		}
		
		z.add(z_max);
		z.insertElementAt(z_min,0);
		
		double cnt = 0.0;
		
		for (j=1;j<z.size();j++) {
			cnt = Math.exp(h.get(j-1)-x.get(j-1)*dh.get(j-1))/dh.get(j-1);
			cn.add(cnt);
			cnt = cnt*(Math.exp(z.get(j)*dh.get(j-1)) - Math.exp(z.get(j-1)*dh.get(j-1)));
			c.add(cnt);
		}
		
		cd.clear();
		cd.addAll(c);
		
		for (j=1;j<c.size();j++) {
			cd.set(j,cd.get(j-1)+cd.get(j));
		}
		
		norm = VectorUtil.sum(c);
		VectorUtil.mult(cd,1.0/norm);
	}
	
	public void recomputeAll() {
		int n = x.size();
		int j = 0;
		
		
		double offsetNew = VectorUtil.max(h);
		VectorUtil.add(h,offset);
		offset = offsetNew;
		VectorUtil.add(h,-offset);
		
		VectorUtil.zeros(z,n-1);
		for (j=0;j<n-1;j++) {
			z.set(j,(h.get(j+1)-h.get(j)-x.get(j+1)*dh.get(j+1)+x.get(j)*dh.get(j))/(dh.get(j)-dh.get(j+1)));
		}
		z.add(z_max);
		z.insertElementAt(z_min,0);
		
		double cnt = 0.0;
		
		for (j=1;j<z.size();j++) {
			cnt = Math.exp(h.get(j-1)-x.get(j-1)*dh.get(j-1))/dh.get(j-1);
			cn.set(j-1,cnt);
			cnt = cnt*(Math.exp(z.get(j)*dh.get(j-1)) - Math.exp(z.get(j-1)*dh.get(j-1)));
			c.set(j-1,cnt);
		}
		
		cd.clear();
		cd.addAll(c);
		
		for (j=1;j<c.size();j++) {
			cd.set(j,cd.get(j-1)+cd.get(j));
		}
		
		norm = VectorUtil.sum(c);
		VectorUtil.mult(cd,1.0/norm);
	}
	
	void recompute(int j_min,int j_max) {
		int j = 0;
		double cnt = 0.0;
		
		for (j=j_min;j<=j_max;j++) {
			if (j<=x.size() - 2)
				z.set(j+1,(h.get(j+1)-h.get(j)-x.get(j+1)*dh.get(j+1)+x.get(j)*dh.get(j))/(dh.get(j)-dh.get(j+1)));
		}
		
		for (j=j_min;j<=j_max;j++) {
			cnt = Math.exp(h.get(j)-x.get(j)*dh.get(j))/dh.get(j);
	    	cn.set(j,cnt);
	    	cnt = cnt*(Math.exp(z.get(j+1)*dh.get(j)) - Math.exp(z.get(j)*dh.get(j)));
	    	c.set(j,cnt);
		}
		
		cd.clear();
		cd.addAll(c);
		
		for (j=1;j<c.size();j++) {
			cd.set(j,cd.get(j-1) + cd.get(j));
		}
		
		norm = VectorUtil.sum(c);
		VectorUtil.mult(cd,1.0/norm);
	}
	
	void insertX(double x_new,double h_new,double dh_new) {
		int fi = x.size() - 1;
		int fp = 0;
		Boolean sat = false;
		
		while(fi >= 0 && !sat) {
			// don't insert it if it's already in there
			if (x.get(fi) == x_new)
				return;
			
			if (x_new > x.get(fi)) {
				sat = true;
			}
			else {
				fi--;
			}
		}
				
		if (!sat)
			fi = -1;
		
		// insert at the end	
		if (fi == (x.size() - 1)) {
			x.add(x_new);
			h.add(h_new);
			dh.add(dh_new);
			cd.add(0.0);
			cn.add(0.0);
			c.add(0.0);
			z.insertElementAt(0.0,fi);
			fp = fi;
		//	recompute(fp-1,fp+1);
			recomputeAll();
			return;
		}
		
		// insert at the beginning
		if (fi == -1) {
			x.insertElementAt(x_new,0);
			h.insertElementAt(h_new,0);
			dh.insertElementAt(dh_new,0);
			cd.insertElementAt(0.0,0);
			cn.insertElementAt(0.0,0);
			c.insertElementAt(0.0,0);
			z.insertElementAt(0.0,1);
			fp = 0;
		//	recompute(fp,fp+1);
			recomputeAll();
			return;
		}
		
		// insert in the middle
		if (fi >= 0 && fi < (x.size() - 1)) {
			fp = fi+1;
			fi++;
			x.insertElementAt(x_new,fi);
			h.insertElementAt(h_new,fi);
			dh.insertElementAt(dh_new,fi);
			cd.insertElementAt(0.0,fi);
			cn.insertElementAt(0.0,fi);
			c.insertElementAt(0.0,fi);
			z.insertElementAt(0.0,fi);
		//	recompute(fp-1,fp+1);
			recomputeAll();
			return;
		}
	}
	
	void sampleS(double[] sm,int[] fi) {
		double w = Uniform.staticNextDouble();
		
		fi[0] = 0;
		double t = 0.0;
		Boolean sat = false;
		
		while(fi[0] < cd.size() && !sat) {
			if (w <= cd.get(fi[0])) {
				sat = true;
			}
			else {
				fi[0]++;
			}
		}
				
		if (!sat)
			fi[0] = -1;
			
		if (fi[0] > 0) {
			t = w-cd.get(fi[0]-1);
		}
		else {
			t = w;
		}
		
		if (fi[0] == -1) {
			sm[0] = 0.0/0.0;
			return;
		}
		
		t = t*norm/cn.get(fi[0]);
		t = t + Math.exp(dh.get(fi[0])*z.get(fi[0]));
		sm[0] = Math.log(t)/dh.get(fi[0]);
	}
	
	Boolean sample() {
		int i = 0;
		int reject = 1;
		double sm = 0.0;
		double[] sm_a = new double[1];
		double w = 0.0;
		int fi = 0;
		int[] fi_a = new int[1];
		double uk = 0.0;
		double lk = 0.0;
		Boolean sat = false;
		Vector<Double> new_h = new Vector<Double>();
		Vector<Double> new_dh = new Vector<Double>();
		Vector<Double> vsm = new Vector<Double>();
		VectorUtil.zeros(new_h,1);
		VectorUtil.zeros(new_dh,1);
		VectorUtil.zeros(vsm,1);
		
		for (i = 0; i < n_samples; i++) {
			reject = 1;
			while(reject == 1) {
				sampleS(sm_a,fi_a);
				sm = sm_a[0];
				fi = fi_a[0];
				if (Double.isNaN(sm)) {
					samples[i] = sm;
					return false;
				}
				w = Uniform.staticNextDouble();
				uk = h.get(fi) + (sm-x.get(fi))*dh.get(fi);
				
				sat = false;
				fi = x.size() - 1;
				while(fi >= 0 && !sat) {
					if (sm >= x.get(fi)) {
						sat = true;
					}
					else {
						fi--;
					}
				}
			
				if ((!sat) || fi >= (x.size() - 1)) {
					lk = -1.0/0.0;
				}
				else {
					lk = ((x.get(fi+1)-sm)*h.get(fi) + (sm-x.get(fi))*h.get(fi+1))/(x.get(fi+1)-x.get(fi));
				}
				
				if (w <= Math.exp(lk-uk)) {
					reject = 0;
					samples[i] = sm;
				}
				else {
					vsm.set(0,sm);
					g_fn.eval(vsm,new_h,new_dh);
					new_h.set(0,new_h.get(0) - offset);
					if (w <= Math.exp(new_h.get(0)-uk)) {
						reject = 0;
						samples[i] = sm;
					}
					else {
						insertX(sm,new_h.get(0),new_dh.get(0));
					}
				}
			}
		}
		
		return true;
	}
	
	public void writeSamplesToFile(String fname) throws IOException {
		int i = 0;
		String s = "";
		FileWriter outFile = new FileWriter(fname);
		for(i=0;i<samples.length;i++) {
			s = (new Double(samples[i])).toString();
			outFile.write(s);
			outFile.write("\n");
		}
		outFile.close();
	}
}
