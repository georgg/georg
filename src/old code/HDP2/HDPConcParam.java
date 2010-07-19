package edu.mit.csail.psrg.georg.HDP2;

import java.io.Serializable;
import java.util.ArrayList;
import cern.jet.random.Beta;
import cern.jet.random.Uniform;
import cern.jet.random.Gamma;
import edu.mit.csail.psrg.georg.StatUtil.StatUtil;

public class HDPConcParam implements Serializable {
	/**
	 * 
	 */
	private static final long serialVersionUID = 4143272203282313196L;
	public double alpha = 1.0;
	public double alphaa = 0.5;
	public double alphab = 0.5;
	
	public DirichletProcess[] DPs;
	public int[] totalNumTables;
	public int[] totalNumData;
	
	public HDPConcParam(double a,double b) {
		alphaa = a;
		alphab = b;
		alpha = alphaa/alphab;
	}
	
	public void addDP(DirichletProcess dd) {
		ArrayList<DirichletProcess> dd2 = new ArrayList<DirichletProcess>();
		dd2.add(dd);
		addDPs(dd2);
	}
	
	public void addDPs(ArrayList<DirichletProcess> dd) {
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
	
	public int findDPNum(DirichletProcess dp) {
		int i = 0;
		while(i<DPs.length) {
			if (DPs[i] == dp) {
				return i;
			}
			i++;
		}
		return -1;
	}
	
	public void removeDP(DirichletProcess dp) {
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
	
	public void appendDP(DirichletProcess dp) {
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
	
	public void calcTotals(int numTopics) {
		int jj = 0;
		int cc = 0;
		
		for (jj=0;jj<DPs.length;jj++) {
			totalNumData[jj] = 0;
			totalNumTables[jj] = 0;
			for (cc=0;cc<numTopics;cc++) {
				if (DPs[jj].state != DirichletProcess.HELDOUT) {
					totalNumData[jj] = totalNumData[jj] + DPs[jj].clusterNumData[cc];
					totalNumTables[jj] = totalNumTables[jj] + DPs[jj].clusterNumTables[cc];
				}
			}
		}
	}
	
	public void updateDPs() {
		int jj = 0;
			
		for (jj=0;jj<DPs.length;jj++) {
			DPs[jj].alpha = alpha;
		}
	}
	
	// deal with resampling if there is only one DP
	// using the concentration parameter
	public void resampleSingle(int numIters) {
			int i = 0;
			double eta = 0.0;
			int zz = 0;
			double[] w = new double[2];
			int numClusters = DPs[0].myModel.numTopics;
			if (numClusters <= 1) {
				return;
			}
			
			int numPresent = totalNumData[0];
			for (i=0;i<numIters;i++) {
				eta = Beta.staticNextDouble(alpha+1.0,(double) numPresent); 
				if (eta > 0) {
					w[0] = (alphaa + (double) (numClusters - 1))/(alphab - Math.log(eta));
					w[1] = (double) numPresent;
					zz = StatUtil.multinomial_rnd(w,2);
					w[0] = alphab - Math.log(eta);
					if (w[0] >= 0.0) {
						if (zz == 0) {
							alpha = Gamma.staticNextDouble(alphaa+((double) numClusters),w[0]);
						} else {
							alpha = Gamma.staticNextDouble(alphaa+((double) numClusters-1),w[0]);
						}
					}
				}
			}
	}
	
	public void resample(int numiter) {
		int iter, jj, zz;
	  	double aa, bb, xx, nd;
	  	
	  	if (DPs.length > 1) {
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
	  	} else {
	  		resampleSingle(numiter);
	  	}
	  	updateDPs();
	}
}
