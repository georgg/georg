package edu.mit.csail.psrg.georg.HDPTopicModel;

import java.util.ArrayList;

public class LogLikeCollector extends HDPStatCollector {
	//	 sampled log likelihood
	ArrayList<Double> logLikeSamples = new ArrayList<Double>();
	
	
	public LogLikeCollector(HierDirichletProcess h) {
		super(h);
	}

	public void updateStats() {
		logLikeSamples.add(myHDP.calcLogLike());
	}

}
