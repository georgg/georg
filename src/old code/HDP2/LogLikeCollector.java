package edu.mit.csail.psrg.georg.HDP2;

import java.util.ArrayList;

public class LogLikeCollector extends HDPStatCollector {
	/**
	 * 
	 */
	private static final long serialVersionUID = -2833356977616700177L;
	//	 sampled log likelihood
	public ArrayList<Double> logLikeSamples = new ArrayList<Double>();
	
	
	public LogLikeCollector(HDP h) {
		super(h);
	}

	public void updateStats() {
		logLikeSamples.add(myHDP.calcLogLike());
	}

}
