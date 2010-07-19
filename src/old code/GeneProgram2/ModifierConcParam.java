package edu.mit.csail.psrg.georg.GeneProgram2;

import java.io.Serializable;
import java.util.ArrayList;

import cern.jet.random.Beta;
import cern.jet.random.Gamma;
import cern.jet.random.Uniform;

public class ModifierConcParam implements Serializable {
	private static final long serialVersionUID = -4401568941593761364L;
	public double[] alpha_theta = null;
	public double alpha_theta_a = 0.5;
	public double alpha_theta_b = 0.5;
	
	public ParentDP[] DPs;
	
	public ModifierConcParam(double a,double b,int numModifiers) {
		alpha_theta_a = a;
		alpha_theta_b = b;
		alpha_theta = new double[numModifiers];
		int i = 0;
		for (i=0;i<numModifiers;i++) {
			alpha_theta[i] = alpha_theta_a/alpha_theta_b;
		}
	}
	
	public void addDPs(ArrayList<ParentDP> dd) {
		int i = 0;
		DPs = new ParentDP[dd.size()];
		
		for (i=0;i<DPs.length;i++) {
			DPs[i] = dd.get(i);
			DPs[i].myModifierConcParam = this;
		}
	}
	
	public void appendDP(ParentDP dp) {
		ParentDP[] DPs2 = new ParentDP[DPs.length+1];
		System.arraycopy(DPs,0,DPs2,0,DPs.length);
		DPs = DPs2;
		DPs[DPs.length-1] = dp;
	}
	
	public int findDPNum(ParentDP dp) {
		int i = 0;
		while(i<DPs.length) {
			if (DPs[i] == dp) {
				return i;
			}
			i++;
		}
		return -1;
	}
	
	public void removeDP(ParentDP dp) {
		int dpNum = findDPNum(dp);
		if (dpNum == -1) {
			return;
		}
		ParentDP[] DPs2 = new ParentDP[DPs.length-1];
		if (dpNum > 0) {
			System.arraycopy(DPs,0,DPs2,0,dpNum);
		}
		System.arraycopy(DPs,dpNum+1,DPs2,dpNum,DPs.length-dpNum-1);
		
		DPs = DPs2;
	}
	
	public void updateDPs() {
		int jj = 0;
		int i = 0;
		for (jj=0;jj<DPs.length;jj++) {
			for (i=0;i<alpha_theta.length;i++) {
				DPs[jj].alpha_thetaMod[i] = alpha_theta[i];
			}
		}
	}
	
	public void resample(int numiter, int numTopics) {
		int iter, jj, zz, cc;
	  	double aa, bb, xx, nd;
	  	int totalNumData = 0;
	  	int totalNumTables = 0;
	  	int i = 0;
	  	int j = 0;
	  	for (i=0;i<alpha_theta.length;i++) { 
		  	for (iter = 0 ;iter < numiter;iter++ ) {
		    	aa = alpha_theta_a;
		    	bb = alpha_theta_b;
		    	for ( jj = 0 ; jj < DPs.length ; jj++ ) {
		    		if (DPs[jj].state != DirichletProcess.HELDOUT) {
		    			for (cc=0;cc<numTopics;cc++) {
		    				totalNumData = 0;
		    			  	totalNumTables = 0;
		    				if (DPs[jj].clusterNumData[cc] > 0) {
		    					for (j=0;j<DPs[jj].numModifierData[i].length;j++) {
		    						totalNumData += DPs[jj].numModifierData[i][j][cc];
		    						totalNumTables += DPs[jj].numModifierTables[i][j][cc];
		    					}
		    				}
		    				if (totalNumData > 0) {
			    				nd = (double) totalNumData;
			    				if (nd <= 0) {
			    					xx = 1.0;
			    				} else {
			    					xx = Beta.staticNextDouble(alpha_theta[i]+1.0, nd);
			    				}
			    				if ( Uniform.staticNextDouble() * (alpha_theta[i] + nd) < nd ) {
			    					zz = 1;
			    				}
			    				else {
			    					zz = 0;
			    				}
				      
			    				aa += totalNumTables - zz;
			    				bb -= Math.log(xx);
		    				}
		    			}
		    		}
		    	}
		    	if (aa >= 1.0) {
		    		alpha_theta[i] = Gamma.staticNextDouble(aa,1.0) / bb;
		    	}
		  	}
	  	}
	  	updateDPs();
	}
}
