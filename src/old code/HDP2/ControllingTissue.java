package edu.mit.csail.psrg.georg.HDP2;

import java.io.Serializable;

public class ControllingTissue implements Serializable, Comparable {
	/**
	 * 
	 */
	private static final long serialVersionUID = 2800414813191274668L;
	public TissueDP dp = null;
	// the percent load of the controlling tissue
	public double load = 0.0;
	// whether the tissue uses the topic for up or down expression
	public boolean UDuse = true;
	String sID;
	int code = 0;
	
	public ControllingTissue(TissueDP d,double l,boolean u) {
		dp = d;
		load = l;
		UDuse = u;
		
		sID = d.label;
		if (UDuse) {
			sID = sID + " UP";
		} else {
			sID = sID + " DN";
		}
		code = sID.hashCode();
	}
	
	public int hashCode() {
		return code;
	}
	
	public void merge(ControllingTissue t2) {
		load += t2.load;
	}
	
	public int compareTo(Object o2) {
		ControllingTissue ct2 = (ControllingTissue) o2;
		
		if (load > ct2.load)
			return -1;
		
		if (load == ct2.load)
			return 0;
		
		return 1;
	}
	
	public boolean equals(Object o) {
		if (((ControllingTissue) o).code == code) {
			return true;
		}
		
		return false;
	}
}
