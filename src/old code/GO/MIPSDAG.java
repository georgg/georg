package edu.mit.csail.psrg.georg.GO;

import java.util.Iterator;

public class MIPSDAG extends OntDAG {
	void propagateGeneCounts() {
		OntTerm[] list = new OntTerm[DAG.size()];
		DAG.values().toArray(list);
		int i = 0;
		int j = 0;
		Iterator<OntGene> iter = null;
		
		for (i=list.length-1;i>=0;i--) {
			if (list[i].parents != null && list[i].genes != null) {
				for (j=0;j<list[i].parents.length;j++) {
					list[i].parents[j].addGenes(list[i].genes);
				}
			}
		}
		
		for (i=list.length-1;i>=0;i--) {
			if (list[i].genes != null) {
				list[i].numGenes = list[i].genes.size();
			} else {
				list[i].numGenes = 0;
			}
		}
 	}
}
