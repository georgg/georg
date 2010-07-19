package edu.mit.csail.psrg.georg.SpecialMathFn;

//adapted from Tony Minka's Lightspeed 2.0 Matlab library
/*
 * TRIGAMMA   Trigamma function.
 * TRIGAMMA(X) returns trigamma(x) = d**2 log(gamma(x)) / dx**2
 * If X is a matrix, returns the trigamma function evaluated at each element.
 *
 * Reference:
 *
 *    B Schneider,
 *    Trigamma Function,
 *    Algorithm AS 121,
 *    Applied Statistics, 
 *    Volume 27, Number 1, page 97-99, 1978.
 * From http://www.psc.edu/~burkardt/src/dirichlet/dirichlet.f
 */
public class Trigamma {
	static double small = 1e-4;
	static double large = 8;
	static double c = Math.pow(Math.PI,2.0)/6.0;
	static double c1 = -2.404113806319188570799476;
	static double b2 =  1.0/6.0;
	static double b4 = -1.0/30.0;
	static double b6 =  1.0/42.0;
	static double b8 = -1.0/30.0;
	static double b10 = 5.0/66.0;
	
	public static double eval(double x) {
		double y = 0.0;
		double z = 0.0;
		
		if (Double.isInfinite(x) || Double.isNaN(x)) {
			return 0.0/0.0;
		}
		
		// zero or negative integer
		if (x <= 0.0 && Math.floor(x) == x) {
			return 1.0/0.0;
		}
		
		// Negative non-integer
		if (x < 0 && Math.floor(x) != x) {
			// Use the derivative of the digamma reflection formula:
			// -trigamma(-x) = trigamma(x+1) - (pi*csc(pi*x))^2
			y = -Trigamma.eval(-x+1.0) + Math.pow(Math.PI*(1.0/Math.sin(-Math.PI*x)),2.0);
			return y;
		}
		
		// Small value approximation
		if (x<= small) {
			y = 1.0/(x*x) + c + c1*x;
			return y;
		}

		// Reduce to trigamma(x+n) where ( X + N ) >= large.
		while(true) {
			if (x>small && x<large) {
				y = y + 1.0/(x*x);
				x = x + 1.0;
			} else {
				break;
			}
		}
		
		if (x>=large) {
			z = 1.0/(x*x);
			y = y + 0.5*z + (1.0 + z*(b2 + z*(b4 + z*(b6 + z*(b8 + z*b10)))))/x;
		}
		
		return y;
	}
}
