package edu.mit.csail.psrg.georg.SampleUtil;

import java.io.Serializable;

import cern.jet.random.Uniform;

/* ARS.C - Procedure for performing Adaptive Rejection Sampling. */

/* Copyright (c) 1996 by Carl Edward Rasmussen and Radford M. Neal
 *
 * Permission is granted for anyone to copy, use, or modify this program 
 * for purposes of research or education, provided this copyright notice 
 * is retained, and note is made of any changes that have been made. 
 *
 * This program is distributed without any warranty, express or implied.
 * As this program was written for research purposes only, it has not been
 * tested to the degree that would be advisable in any important application.
 * All use of this program is entirely at the user's own risk.
 */


/* This module implements the "Adaptive Rejection Sampling" scheme due to
 * Gilks and Wiles.  See "Adaptive Rejection Sampling for Gibbs Sampling",
 * Applied Statistics, vol. 41, no. 2, pp. 337-348 (1992).  This is not
 * the most sophisticated possible implementation of the method, however. */

// Java adaptation by Georg Gerber (georg@mit.edu) 2006

public class AdaptiveRejectionSampleCR implements Serializable {
	static double MIN_DROP = 0.1;  /* minimum drop in log probability for initial points */
	static double TINY = 1e-9;         /* minimum tolerance for inserting new points */
	static int MAX_LIST = 100;     /* max number of segments in piecewise exponential */
	
	Segtype root = null;
	
	/* ADAPTIVE REJECTION SAMPLING PROCEDURE.
	 *
	 * The ars function returns a single point sampled at random from the univariate
	 * distribution specified by the function logp, provided by the user.  The
	 * logp function takes a point as its first argument, and returns the log 
	 * of the probability density at that point, plus any arbitrary constant (not
	 * depending on the point).  The logp function takes two additional arguments,
	 * the first a pointer to a place where it must store the derivative of the
	 * log probability density, the second a pointer to additional information
	 * describing the distribution, which is passed on unchanged from the last 
	 * argument of ars.
	 *
	 * The logp function passed MUST be log-concave. It is assumed that any real 
	 * number is legal as input to logp. 
	 *
	 * The user must also supply an initial guess, "init", and a typical "scale" 
	 * of variation.  It is not essential that these values be very accurate, but 
	 * performance will generally depend on their accuracy.
	 *
	 * The ars function first tries to locate points on either side of the mode;
	 * the derivative must have the right sign and be non-negligible, and the drop 
	 * from the max. probability seen must be of at least moderate size to qualify.
	 * Then a piece-wise exponential distribution is iteratively improved using 
	 * knowledge from rejected points, until a sample is accepted. At most MAX_LIST
	 * pieces are introduced in the approximation. Before pieces are inserted 
	 * various checks are made in order to prevent numerical problems. If new points
	 * don't qualify, the peice-wise exponential is simply not updated. A warning 
	 * will be issued (one time only) when 10000 rejections are exceeded */

	public double ars(
			  double init,                                              /* initial guess */
			  double scale,                         /* guess for scale of variation in x */
			  LogDensityFnARSCR logp       /* function to sample from */
			)
	{
		Segtype[] seg = new Segtype[MAX_LIST];
		int i = 0;
		for (i=0;i<MAX_LIST;i++) {
			seg[i] = new Segtype();
		}
		
		Segtype prv, cur, nxt;
		Segtype[] q = new Segtype[1];
		int no_seg = 2;
		boolean warning = false;
		double x, max, f, df;
		double[] dfe = new double[1];

		root = seg[0]; nxt = seg[1];
		root.prv = null; root.nxt = nxt; nxt.prv = root; nxt.nxt = null;
		
		x = init; f = logp.eval(x, dfe);       /* find point to the left of mode */
		df = dfe[0];
		max = f;
		while (df<TINY || f>max-MIN_DROP) {
			x -= scale+(init-x); 
			f = logp.eval(x, dfe);
			df = dfe[0];
		    if (f>max) max = f;
		  }
		root.x = x; root.a = df; root.b = f-x*df;
		
		x = init; f = logp.eval(x, dfe);      /* find point to the right of mode */
		df = dfe[0];
		while (df>-TINY || f>max-MIN_DROP) {
			x += scale+(x-init); 
		    f = logp.eval(x, dfe);
		    df = dfe[0];
		    if (f>max) max = f;
		}
		nxt.x = x; nxt.a = df; nxt.b = f-x*df;
		root.xmax = (nxt.b-root.b)/(root.a-nxt.a);
		
		for (i=0; ; i++) {                     /* repeat until a point is accepted */

		    if (i==10000 && !(warning)) {
		      System.out.println("WARNING: More than 10000 rejections in ars\n");
		      warning = true;
		    } 

		    cur = root;      /* find max y-value; needed to avoid numerical problems */
		    max = cur.a*cur.xmax+cur.b;
		    
		    cur = cur.nxt;
		    while(cur.nxt != null) {
		    	x = cur.a*cur.xmax+cur.b;
		    	if (x > max) max = x;
		    	cur = cur.nxt;
		    }

		    cur = root;                                            /* compute masses */
		    cur.mass = Math.exp(cur.a*cur.xmax+cur.b-max)/cur.a;
		    prv = cur;
		    cur = cur.nxt;
		    while (cur.nxt != null) {
		      cur.mass = (Math.exp(cur.a*cur.xmax+cur.b-max)-
		                   Math.exp(cur.a*prv.xmax+cur.b-max))/cur.a;
		      prv = cur;
		      cur = cur.nxt;
		    }
		    cur.mass = -Math.exp(cur.a*prv.xmax+cur.b-max)/cur.a;

		    x = rpwed(q);                               /* this is the new sample */
		    cur = q[0];
		    f = logp.eval(x, dfe);
		    df = dfe[0];
		    if (Uniform.staticNextDouble() <= Math.exp(f-cur.a*x-cur.b)) return x;      /* success! */
		    
		    /* Now, insert a new piece in the piece-wise approximation if the situation is
		     * appropriate. Eg, if we have enough memory, if the slope at the new x isn't
		     * too small (the exponential distribution will degenerate), and if the slope
		     * isn't too close to the slope of currrent, previous and next segment (or
		     * which ever of these may exist), since this may cause numerical problems. */

		    if (no_seg < MAX_LIST && Math.abs(df) > TINY && Math.abs(df-cur.a) > TINY &&
		          ((cur.prv == null) || Math.abs(df-cur.prv.a) > TINY) &&
		          ((cur.nxt == null) || Math.abs(df-cur.nxt.a) > TINY)) {      

		          if (x<cur.x) cur = cur.prv;               /* now, insert *after* cur */
		          prv = cur; cur = seg[no_seg++]; cur.prv = prv;
		          if (prv != null)
		            { cur.nxt = prv.nxt; prv.nxt = cur; }
		          else
		            { cur.nxt = root; root = cur; }
		          nxt = cur.nxt; 
		          if (nxt != null) nxt.prv = cur;
		          cur.x = x; cur.a = df; cur.b = f-x*df;

		          if (prv != null) prv.xmax = (cur.b-prv.b)/(prv.a-cur.a);
		          if (nxt != null) cur.xmax = (nxt.b-cur.b)/(cur.a-nxt.a);
		   }

		}

	}
	
	/* Private function to sample from piece-wise exponential distribution. First a
	 * piece under the piece-wise distribution is sampled at random. Then a random
	 * sample is drawn from this piece.  A pointer to the segment, or piece which
	 * was used is returned in q, and the function returns the random sample. Care
	 * is taken to avoid numerical over and underflow.  */
	double rpwed(Segtype[] q)
	{
	  double mass = 0.0, t, u;
	  
	  q[0] = root; while (q[0] != null) { mass += q[0].mass;  q[0] = q[0].nxt; }
	  t = mass*Uniform.staticNextDouble();
	  q[0] = root; while ((q[0].nxt != null) && ((t -= q[0].mass) >= 0.0)) { q[0] = q[0].nxt; } 

	  u = Uniform.staticNextDouble();
	  if (q[0].prv == null)
	    return q[0].xmax+Math.log(u)/q[0].a;
	  if (q[0].nxt == null) 
	    return q[0].prv.xmax+Math.log(u)/q[0].a;
	  t = Math.log(u+(1.0-u)*Math.exp(Math.abs(q[0].a)*(q[0].prv.xmax-q[0].xmax)))/q[0].a;
	  return (q[0].a > 0) ? q[0].xmax+t : q[0].prv.xmax+t;
	}
}
