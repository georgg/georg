package edu.mit.csail.psrg.georg.HDPTopicModel;

import java.io.Serializable;

public abstract class HDPStatCollector implements Serializable {
	HierDirichletProcess myHDP = null;
	
	public HDPStatCollector(HierDirichletProcess h) {
		myHDP = h;
	}
	
	public abstract void updateStats();
}
