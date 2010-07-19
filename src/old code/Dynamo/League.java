package com.newmacondo.georg.Dynamo;

import java.util.HashSet;
import java.util.Iterator;

public class League {
// A League is a collection of OTUs with very similar dynamics/interactions
// interslice dependencies are assuming to be at lag 1
	HashSet<League> intraParents = new HashSet<League>();
	HashSet<League> intraChildren = new HashSet<League>();
	HashSet<League> interParents = new HashSet<League>();
	HashSet<League> interChildren = new HashSet<League>();
	DBN myDBN = null;
	// used for determining if graph is cyclic
	int nodeColor = 0;
	
	public League(DBN d) {
		myDBN = d;
	}
	
	public void removeAllLinks() {
		Iterator<League> iter = intraParents.iterator();
		League l = null;
		while (iter.hasNext()) {
			l = iter.next();
			l.removeIntraChild(this);
		}
		iter = interParents.iterator();
		while (iter.hasNext()) {
			l = iter.next();
			l.removeInterChild(this);
		}
		iter = intraChildren.iterator();
		while (iter.hasNext()) {
			l = iter.next();
			removeIntraChild(l);
		}
		iter = interChildren.iterator();
		while (iter.hasNext()) {
			l = iter.next();
			removeInterChild(l);
		}
	}
	
	public void removeIntraChild(League l) {
		intraChildren.remove(l);
		l.intraParents.remove(this);
	}
	
	public boolean addIntraChild(League l) {
		if (l == this) {
			return false;
		}
		if (intraChildren.contains(l) || intraParents.contains(l)) {
			return false;
		}
		intraChildren.add(l);
		l.intraParents.add(this);
		return true;
	}
	
	public void removeInterChild(League l) {
		interChildren.remove(l);
		l.interParents.remove(l);
	}
	
	public boolean addInterChild(League l) {
		if (intraChildren.contains(l)) {
			return false;
		}
		interChildren.add(l);
		l.interParents.add(this);
		return true;
	}
}
