package edu.mit.csail.psrg.georg.HDPTopicModel;

import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;

public class StrongModuleCollector extends HDPStatCollector {
	// thresh-hold for p-value to determine that a document controls a module
	double documentControlPvalThresh = 0.05;
	// minimum number of genes in intersection
	int minGenesInModule = 5;
	
	// determines whether to collect module statistics
	boolean collectStrongModules = true;
	// maps sets of documents to modules
	HashMap<HashSet<Document>,StrongModule> strongModuleMap = new HashMap<HashSet<Document>,StrongModule>();
	// maximum number of strong modules that will be maintained
	int maxNumStrongModules = 2500;
	// interval to "clean" the number of strong modules
	int cleanStrongModulesInterval = 1000;
	int cleanStrongModulesCounter = 0;
	int numSamples = 0;
	
	public StrongModuleCollector(HierDirichletProcess h) {
		super(h);
	}

	public void updateStats() {
		int i = 0;
  		int j = 0;
  		int k = 0;
  		double p = 0;
  		DirichletProcess dp = null;
  		Document[] docs = null;
  		Document doc = null;
  		StrongModule[] modules = new StrongModule[myHDP.numClusters];
  		numSamples++;
  		
  		for (i=0;i<myHDP.numClusters;i++) {
  			modules[i] = new CrossSpeciesStrongModule(myHDP.totalGenes,myHDP.posCounts[i]);
  		//	modules[i] = new SingleStrongModule(totalGenes,posCounts[i]);
  		}
  		
  		for (j=0;j<myHDP.DP.length;j++) {
  			dp = myHDP.DP[j];
  			if (dp.state != DirichletProcess.HELDOUT) {
				docs = dp.documents;
				if (docs != null) {
					for (k=0;k<docs.length;k++) {
						doc = docs[k];
						doc.allocateGenesToClusterMap(myHDP.tempClusterMap,myHDP.tempClusterArray,myHDP.numClusters);
						for (i=0;i<myHDP.numClusters;i++) {
							p = modules[i].controlPVal(myHDP.tempClusterMap[i],doc,myHDP.totalGenes,minGenesInModule);
	  						if (p <= documentControlPvalThresh) {
	  							modules[i].addControllingDocument(doc,myHDP.tempClusterMap[i],myHDP.tempClusterArray[i]);
	  						}
						}
					}
				}
  			}
		}
  		
  		for (i=0;i<myHDP.numClusters;i++) {
  			if (modules[i].controllingDocuments.size() >= 1 & modules[i].totalGenesInModule >= minGenesInModule) {
  				modules[i].addToModuleMap(strongModuleMap,myHDP.totalGenes);
  			}
  		}
  		
  		cleanStrongModulesCounter++;
  		if (cleanStrongModulesCounter >= cleanStrongModulesInterval) {
  			cleanStrongModulesCounter = 0;
  			if (strongModuleMap.size() > maxNumStrongModules) {
  				StrongModule module = null;
  				StrongModule[] marray = new StrongModule[strongModuleMap.size()];
				marray = strongModuleMap.values().toArray(marray);
				List<StrongModule> mlist = Arrays.asList(marray);
				Collections.sort(mlist);
  				int ss = strongModuleMap.size();
  				Iterator<StrongModule> iter = mlist.iterator();
  				while(ss > maxNumStrongModules & iter.hasNext()) {
  					module = iter.next();
  					module = strongModuleMap.remove(module.controllingDocuments);
  					ss--;
  				}
  			}
  		}
	}
}
