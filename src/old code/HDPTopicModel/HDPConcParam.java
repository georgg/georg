package edu.mit.csail.psrg.georg.HDPTopicModel;

import java.io.Serializable;
import java.util.ArrayList;
import cern.jet.random.Beta;
import cern.jet.random.Uniform;
import cern.jet.random.Gamma;

public class HDPConcParam implements Serializable {
	double alpha = 1.0;
	double alphaa = 0.5;
	double alphab = 0.5;
	
	DirichletProcess[] DPs;
	int[] totalNumTables;
	int[] totalNumData;
	
	HDPConcParam(double a,double b) {
		alphaa = a;
		alphab = b;
		alpha = alphaa/alphab;
	}
	
	void addDP(DirichletProcess dd) {
		ArrayList<DirichletProcess> dd2 = new ArrayList<DirichletProcess>();
		dd2.add(dd);
		addDPs(dd2);
	}
	
	void addDPs(ArrayList<DirichletProcess> dd) {
		int i = 0;
		DPs = new DirichletProcess[dd.size()];
		totalNumTables = new int[dd.size()];
		totalNumData = new int[dd.size()];
		
		for (i=0;i<DPs.length;i++) {
			DPs[i] = dd.get(i);
			totalNumTables[i] = 0;
			totalNumData[i] = 0;
			DPs[i].concParam = this;
		}
	}
	
	int findDPNum(DirichletProcess dp) {
		int i = 0;
		while(i<DPs.length) {
			if (DPs[i] == dp) {
				return i;
			}
			i++;
		}
		return -1;
	}
	
	void removeDP(DirichletProcess dp) {
		int dpNum = findDPNum(dp);
		if (dpNum == -1) {
			return;
		}
		DirichletProcess[] DPs2 = new DirichletProcess[DPs.length-1];
		int[] totalNumTables2 = new int[DPs.length-1];
		int[] totalNumData2 = new int[DPs.length-1];
		if (dpNum > 0) {
			System.arraycopy(DPs,0,DPs2,0,dpNum);
			System.arraycopy(totalNumTables,0,totalNumTables2,0,dpNum);
			System.arraycopy(totalNumData,0,totalNumData2,0,dpNum);
		}
		System.arraycopy(DPs,dpNum+1,DPs2,dpNum,DPs.length-dpNum-1);
		System.arraycopy(totalNumTables,dpNum+1,totalNumTables2,dpNum,DPs.length-dpNum-1);
		System.arraycopy(totalNumData,dpNum+1,totalNumData2,dpNum,DPs.length-dpNum-1);
		
		DPs = DPs2;
		totalNumTables = totalNumTables2;
		totalNumData = totalNumData2;
	}
	
	void appendDP(DirichletProcess dp) {
		DirichletProcess[] DPs2 = new DirichletProcess[DPs.length+1];
		int[] totalNumTables2 = new int[DPs.length+1];
		int[] totalNumData2 = new int[DPs.length+1];
		System.arraycopy(DPs,0,DPs2,0,DPs.length);
		System.arraycopy(totalNumTables,0,totalNumTables2,0,DPs.length);
		System.arraycopy(totalNumData,0,totalNumData2,0,DPs.length);
		DPs = DPs2;
		totalNumTables = totalNumTables2;
		totalNumData = totalNumData2;
		
		DPs[DPs.length-1] = dp;
	}
	
	void calcTotals() {
		int jj = 0;
		int cc = 0;
		
		for (jj=0;jj<DPs.length;jj++) {
			totalNumData[jj] = 0;
			totalNumTables[jj] = 0;
			for (cc=0;cc<DPs[jj].clusterNumData.length;cc++) {
				if (DPs[jj].state != DirichletProcess.HELDOUT) {
					totalNumData[jj] = totalNumData[jj] + DPs[jj].clusterNumData[cc];
					totalNumTables[jj] = totalNumTables[jj] + DPs[jj].clusterNumTables[cc];
				}
			}
		}
	}
	
	void updateDPs() {
		int jj = 0;
			
		for (jj=0;jj<DPs.length;jj++) {
			DPs[jj].alpha = alpha;
		}
	}
	
	void resample(int numiter) {
		int iter, jj, zz;
	  	double aa, bb, xx, nd;
	  	for (iter = 0 ;iter < numiter;iter++ ) {
	    	aa = alphaa;
	    	bb = alphab;
	    	for ( jj = 0 ; jj < totalNumData.length ; jj++ ) {
	    		if (DPs[jj].state != DirichletProcess.HELDOUT) {
		      		nd = (double) totalNumData[jj];
		      		if (nd <= 0) {
		      			xx = 1.0;
		      		} else {
		      			xx = Beta.staticNextDouble(alpha+1.0, nd);
		      		}
		      		if ( Uniform.staticNextDouble() * (alpha + nd) < nd ) {
		      			zz = 1;
		      		}
		      		else {
		      			zz = 0;
		      		}
		      
		      		aa += totalNumTables[jj] - zz;
		      		bb -= Math.log(xx);
	    		}
	    	}
	    	if (aa >= 1.0) {
	    		alpha = Gamma.staticNextDouble(aa,1.0) / bb;
	    	}
	  	}
	  	updateDPs();
	}
}
