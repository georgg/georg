package com.newmacondo.georg.Dynamo;

import java.util.ArrayList;
import java.util.Iterator;

public class DBN {
	ArrayList<League> nodes = new ArrayList<League>();
	ArrayList<OTU> OTUs = new ArrayList<OTU>();
	
	void addNode(League l) {
		nodes.add(l);
	}
	
	void removeNode(League l) {
		l.removeAllLinks();
		l.myDBN = null;
		nodes.remove(l);
	}
	
	public boolean isAcyclic() {
		boolean testAcyclic = true;
		Iterator<League> nodeIter = nodes.iterator();
		League l = null;
		
		while(nodeIter.hasNext()) {
			l = nodeIter.next();
			l.nodeColor = 0;
		}
		
		nodeIter = nodes.iterator();
		
		while(nodeIter.hasNext() && testAcyclic) {
			l = nodeIter.next();
			if (l.nodeColor == 0) {
				testAcyclic = DFS(l);
			}
		}
		
		return testAcyclic;
	}
	
	public boolean DFS(League l) {
		l.nodeColor = 1;
		Iterator<League> nodeIter = l.intraChildren.iterator();
		
		boolean testAcyclic = true;
		League l2 = null;
		
		while(nodeIter.hasNext() && testAcyclic) {
			l2 = nodeIter.next();
			if (l2.nodeColor == 0) {
				testAcyclic = DFS(l2);
			} else {
				if (l2.nodeColor == 1) {
					testAcyclic = false;
				}
			}
		}
	l.nodeColor = 2;
	return testAcyclic;
	}
}
