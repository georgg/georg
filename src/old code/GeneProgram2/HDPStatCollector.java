package edu.mit.csail.psrg.georg.GeneProgram2;

import java.io.Serializable;

public abstract class HDPStatCollector implements Serializable {
	HDP myHDP = null;
	public int numSamples = 0;
	
	public HDPStatCollector(HDP h) {
		myHDP = h;
	}
	
	public abstract void updateStats();
}
