package edu.mit.csail.psrg.georg.GeneProgram2;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedHashSet;

import cern.jet.random.Gamma;
import cern.jet.random.Uniform;
import edu.mit.csail.psrg.georg.StatUtil.StatUtil;
import java.io.Serializable;

public abstract class DirichletProcess implements Serializable {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1165762274729634726L;
	public double alpha;
	public double[] beta;
	public int[] clusterNumData;
	public int[] clusterNumTables;
	public ParentDP parent = null;
	public HDPConcParam concParam = null;
	public String label;
	
	public int dpID = 0;
	public static int ACTIVE = 0;
	public static int HELDOUT = 1;
	public static int FROZEN = 2;
	public int state = 0;
	public boolean generateNewClusters = true;
	public int code = 0;
	public HDP myModel = null;
  	
	public DirichletProcess(HDP mm, ParentDP p, String l) {
  		myModel = mm;
  		parent = p;
  		label = l;
  		alpha = 1.0;
  		if (p != null) {
  			p.addChild(this);
  		}
  	}
  	
  	public int getCode() { return code; }
  	public void setCode(int c) { code = c; }
  	
  	public void addExpressionProgram(int numExpressionPrograms) {
  		if (state != DirichletProcess.HELDOUT) {
			beta[numExpressionPrograms+1] = 0;
			clusterNumData[numExpressionPrograms+1] = 0;
			clusterNumTables[numExpressionPrograms+1] = 0;
  		}
  	}
  	
  	public void deleteExpressionProgram(int c, int numClusters) {
  		System.arraycopy(beta,c+1,beta,c,numClusters-c);
  		System.arraycopy(clusterNumData,c+1,clusterNumData,c,numClusters-c);
  		System.arraycopy(clusterNumTables,c+1,clusterNumTables,c,numClusters-c);
  	}
  	
  	public void growExpressionPrograms() {
  		beta = myModel.growDoubleArray(beta);
  		clusterNumData = myModel.growIntArray(clusterNumData);
  		clusterNumTables = myModel.growIntArray(clusterNumTables);
  	}
  	
  	public void sampleNewParameters(int numClusters) {
  		sampleNewBeta(numClusters);
  	}
  	
  	public void sampleNewBeta(int numClusters) {
  		double b1 = 0;
  		double b2 = 0;
  		double bb = 0;
  		
		if (parent == null) {
			b1 = Gamma.staticNextDouble(1.0,1.0);
			if (alpha <=0.0) {
				b2 = 0.0;
			} else {
				b2 = Gamma.staticNextDouble(alpha,1.0);
			}
		}
		else {
			b1 = alpha*parent.beta[numClusters];
			if (b1 <= 0.0) {
				b1 = 0.0;
			} else {
				b1 = Gamma.staticNextDouble(b1,1.0);
			}
			b2 = alpha*parent.beta[numClusters+1];
			if (b2 <= 0.0) {
				b2 = 0.0;
			} else {
				b2 = Gamma.staticNextDouble(b2,1.0);
			}
		}
		
		if (Double.isNaN(b1)) {
			b1 = 0.0;
		}
		
		if (Double.isNaN(b2)) {
			b2 = 0.0;
		}
		
		bb = beta[numClusters] / (b1 + b2);
		
		if (Double.isNaN(bb)) {
			bb = 0.0;
		}
		
		beta[numClusters] = bb * b1;
      	beta[numClusters+1] = bb * b2;
      	
      	if (Double.isNaN(beta[numClusters]) || Double.isInfinite(beta[numClusters])) {
      		beta[numClusters] = 0.0;
      	}
      	
      	if (Double.isNaN(beta[numClusters+1]) || Double.isInfinite(beta[numClusters+1])) {
      		beta[numClusters+1] = 0.0;
      	}
  	}
  	
  	public void resampleBeta(int numClusters,double[] condLike,double[] beta2) {
  		int cc = 0;
  		
  		if (parent == null) {
  			for(cc=0;cc<numClusters;cc++) {
  				condLike[cc] = clusterNumData[cc];
  			}
  			condLike[numClusters] = clusterNumData[numClusters] + alpha;
  		}
  		else {
  			for(cc=0;cc<numClusters+1;cc++) {
  				condLike[cc] = clusterNumData[cc] + alpha*parent.beta[cc];
  			}
  		}
  		
  		StatUtil.dirichlet_rnd(beta2,condLike,numClusters+1);
  	}
  	
  	public void resampleParameters(int numClusters,double[] condLike) {
  		resampleBeta(numClusters,condLike);
  	}
  	
  	public void resampleBeta(int numClusters,double[] condLike) {
  		resampleBeta(numClusters,condLike,beta);
  	}
  	
  	public int numTables_rnd(double alpha_s, int numData) {
  		int ii, numtable;

  	  	if ( numData == 0 ) {
  	    	return 0;
  	  	} else {
  	    	numtable = 1;
  	    	for ( ii = 1 ; ii < numData ; ii++ ) {
  	      		if ( Uniform.staticNextDouble() < alpha_s / (((double) ii)+alpha_s) ) numtable++;
  	    	}
  	    	return numtable;
  	  	}
  	}
  	
  	public void sampleClusterNumTables(int numClusters) {
  		int cc = 0;
  		
  		if (parent == null) {
  			for (cc=0;cc<numClusters;cc++) {
  				if (clusterNumData[cc] > 0) {
  					clusterNumTables[cc] = 1;
  				}
  			}
  		}
  		else {
  			for (cc=0;cc<numClusters;cc++) {
  				parent.clusterNumData[cc] = parent.clusterNumData[cc] - clusterNumTables[cc];
  				clusterNumTables[cc] = numTables_rnd(alpha*parent.beta[cc],clusterNumData[cc]);
  				parent.clusterNumData[cc] = parent.clusterNumData[cc] + clusterNumTables[cc];
  			}
  		}
  	}
  	
  	public void activate() {
  		int cc = 0;
  		beta = new double[myModel.allocExpressionPrograms];
  		clusterNumData = new int[myModel.allocExpressionPrograms];
  		clusterNumTables = new int[myModel.allocExpressionPrograms];
  		
  		for (cc=0;cc<clusterNumData.length;cc++) {
  			clusterNumTables[cc] = clusterNumData[cc];
  		}
  		
  		if (parent != null) {
  			for (cc=0;cc<clusterNumData.length;cc++) {
  				parent.clusterNumData[cc] = parent.clusterNumData[cc] + clusterNumTables[cc];
  			}
  		}
  		
  		state = DirichletProcess.ACTIVE;
  		alpha = concParam.alpha;
  	}
  	
  	public void holdout() {
  		int cc = 0;
  		DirichletProcess ancestor = null;
  		if (state != DirichletProcess.HELDOUT) {
  			if (parent != null) {
  				for (cc=0;cc<clusterNumData.length;cc++) {
  					parent.clusterNumData[cc] = parent.clusterNumData[cc] - clusterNumTables[cc];
  				}
  				ancestor = parent;
  				while(ancestor != null) {
  					ancestor.sampleClusterNumTables(myModel.numExpressionPrograms);
  					ancestor = ancestor.parent;
  				}
  			}
  		}
  	}
  	
  	public void freeze() {
  		state = DirichletProcess.FROZEN;
  	}
  	
  	public String toString() {
  		String s = "DP " + label + ": [ ";
  		int cc = 0;
  		for (cc=0;cc<beta.length;cc++) {
  			s = s + beta[cc] + " ";
  		}
  		s = s + "]";
  		return s;
  	}
  	
  	public String getLabel() {
  		if (label.length() > 0) {
  			return label;
  		}
  		
  		String label2 = "none";
  		return label2;
  	}
}
