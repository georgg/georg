package edu.mit.csail.psrg.georg.GO;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;

public class OntDAG {
	LinkedHashMap<String,OntTerm> DAG = new LinkedHashMap<String,OntTerm>();
	
	public void addTerm(OntTerm term) {
		DAG.put(term.ID,term);
	}
	
	public OntTerm getTerm(String ID) {
		return DAG.get(ID);
	}
	
	public void associateParents() {
		OntTerm[] mv = new OntTerm[DAG.size()];
		DAG.values().toArray(mv);
		associateParentsUtil(mv);
	}
	
	void associateParentsUtil(OntTerm[] mv) {
		int i = 0;
		for(i=0;i<mv.length;i++) {
			mv[i].associateParents(this);
		}
	}
	
	void topologicalSort() throws DAGException {
		ArrayList<OntTerm> list = new ArrayList<OntTerm>();
		Iterator<OntTerm> iter = DAG.values().iterator();
		OntTerm term = null;
		while(iter.hasNext()) {
			term = iter.next();
			if (term.color == OntTerm.White) {
				DFSVisit(term,list);
			}
		}
		DAG.clear();
		iter = list.iterator();
		while(iter.hasNext()) {
			term = iter.next();
			addTerm(term);
		//	System.out.println(term.ID);
		}
	}
	
	void DFSVisit(OntTerm term,ArrayList<OntTerm> list) throws DAGException {
		if (term.color != OntTerm.White) {
			System.out.println("DAG error detected in " + term.ID);
			throw new DAGException();
		}
		OntTerm child = null;
		term.color = OntTerm.Gray;
		if (term.children != null) {
			Iterator<OntTerm> iter = term.children.iterator();
			while(iter.hasNext()) {
				child = iter.next();
				if (child.color == OntTerm.White) {
					DFSVisit(child,list);
				}
			}
		}
		term.color = OntTerm.Black;
		list.add(0,term);
	}
	
	void propagateGeneCounts() {
		OntTerm[] list = new OntTerm[DAG.size()];
		DAG.values().toArray(list);
		int i = 0;
		int j = 0;
		
		for (i=list.length-1;i>=0;i--) {
			if (list[i].genes != null) {
				list[i].numGenes = list[i].numGenes + list[i].genes.size();
			}
			if (list[i].parents != null) {
				for (j=0;j<list[i].parents.length;j++) {
					list[i].parents[j].numGenes = list[i].parents[j].numGenes + list[i].numGenes;
				}
			}
		}
 	}
}
