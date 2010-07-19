package edu.mit.csail.psrg.georg.SampleUtil;

import java.util.Vector;

// base class for function to be evaluated by
// the adaptive rejection sampler
public abstract class AdaptSampleFn {
	public abstract void eval(Vector<Double> x, Vector<Double> h, Vector<Double> dh);
}
