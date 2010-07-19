package edu.mit.csail.psrg.georg.SampleUtil;

import java.io.Serializable;

// utility class for AdaptiveRejectionSamplerCS
public class Segtype implements Serializable {
	double a;           /* in the log domain, the piece wise... */
    double b;     /* exponential is a piece wise linear: y=ax+b */
    double x;                        /* interior point in piece */
    double xmax;                   /* upper limit of this piece */
    double mass;          /* the probability mass of this piece */
    Segtype prv = null;               /* ptr to previous piece */
    Segtype nxt = null; 
}
