package edu.mit.csail.psrg.georg.SampleUtil;

public interface LogDensityFnARSCR {
	// returns the log of a density function (for adapative rejection sampling)
	// x is the value to evaluate the function at, and dh is for return of the
	// derivative
	double eval(double x, double[] dh);
}
