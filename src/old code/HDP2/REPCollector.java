package edu.mit.csail.psrg.georg.HDP2;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;

public class REPCollector extends HDPStatCollector {
	/**
	 * 
	 */
	private static final long serialVersionUID = 6277386663346045605L;
	// thresh-hold for p-value to determine that a document controls a module
	double tissueControlPvalThresh = 0.05;
	// minimum number of genes in intersection
	int minGenesInREP = 5;
	
	// maps sets of tissues to recurrent expression programs
	LinkedHashMap<HashMap<ControllingTissue,ControllingTissue>,REP> REPMap = new LinkedHashMap<HashMap<ControllingTissue,ControllingTissue>,REP>();
	// maximum number of strong modules that will be maintained
	int maxNumREPs = 2500;
	// interval to "clean" the number of strong modules
	int cleanREPInterval = 1000;
	int cleanREPCounter = 0;
	
	public REPCollector(HDP h) {
		super(h);
	}

	public void updateStats() {
		int i = 0;
  		int j = 0;
  		double p = 0;
  		DirichletProcess dp = null;
  		REP[] REPs = new REP[myHDP.numTopics];
  		numSamples++;
  		
  		for (i=0;i<myHDP.numTopics;i++) {
  			REPs[i] = new GenericREP(myHDP.totalGenes,myHDP.topics.get(i).posCounts);
  		}
  		
  		for (j=0;j<myHDP.DP.length;j++) {
  			dp = myHDP.DP[j];
  			if (dp.state != DirichletProcess.HELDOUT & dp.getClass() == TissueDP.class) {
				((TissueDP) dp).allocateGenesToTopicMap();
				for (i=0;i<myHDP.numTopics;i++) {
					p = REPs[i].controlPVal(myHDP.topics.get(i).topicMap,(TissueDP) dp,myHDP.totalGenes,minGenesInREP);
	  				if (p <= tissueControlPvalThresh) {
				//	if (p <= 0.005) {
	  					REPs[i].addControllingTissue((TissueDP) dp,myHDP.topics.get(i).topicMap,myHDP.topics.get(i).topicArray,((TissueDP) dp).upDownTopicChoice[i]);
	  				}
				}
			}
		}
  		
  		for (i=0;i<myHDP.numTopics;i++) {
  			if (REPs[i].controllingTissues.size() >= 1 & REPs[i].totalGenesInREP >= minGenesInREP) {
  				REPs[i].addToREPMap(REPMap,myHDP.totalGenes);
  			}
  		}
  		
  		cleanREPCounter++;
  		if (cleanREPCounter >= cleanREPInterval) {
  			cleanREPCounter = 0;
  			if (REPMap.size() > maxNumREPs) {
  				REP module = null;
  				REP[] marray = new REP[REPMap.size()];
				marray = REPMap.values().toArray(marray);
				List<REP> mlist = Arrays.asList(marray);
				Collections.sort(mlist);
  				int ss = REPMap.size();
  				Iterator<REP> iter = mlist.iterator();
  				while(ss > maxNumREPs & iter.hasNext()) {
  					module = iter.next();
  					module = REPMap.remove(module.controllingTissues);
  					ss--;
  				}
  			}
  		}
	}
}
