package edu.mit.csail.psrg.georg.HDP2;

import java.util.ArrayList;

import edu.mit.csail.psrg.georg.StatUtil.StatUtil;

public class ParentDP extends DirichletProcess {
	/**
	 * 
	 */
	private static final long serialVersionUID = 7752976779536105909L;

	public ArrayList<DirichletProcess> children = null;
	
	// used for sampling theta
	public int[] numUpData = null;
	public int[] numUpTables = null;
	public int[] numDownData = null;
	public int[] numDownTables = null;
	
	// the probabilities for choosing whether
	// topics represent up or down regulated genes
	public double[] theta = null;
	// parameter for pushing down influence of parents
	// on theta distribution
	public double alpha_theta = 1.0;
	public UpDownConcParam myUpDownConcParam = null;
	// temporary variables used for sampling theta
	public double[] ud = new double[2];
	public double[] udt = new double[2];
	
	public ParentDP(HDP mm, ParentDP p, String l) {
		super(mm, p, l);
	}

	public void addTopic(int numTopics) {
		super.addTopic(numTopics);
  		if (state != DirichletProcess.HELDOUT) {
  			numUpData[numTopics+1] = 0;
  			numUpTables[numTopics+1] = 0;
  			numDownData[numTopics+1] = 0;
  			numDownTables[numTopics+1] = 0;
  			theta[numTopics+1] = 0;
  		}
	}
	
	public void deleteTopic(int c, int numClusters) {
		super.deleteTopic(c,numClusters);
  		System.arraycopy(theta,c+1,theta,c,numClusters-c);
  		System.arraycopy(numUpData,c+1,numUpData,c,numClusters-c);
  		System.arraycopy(numUpTables,c+1,numUpTables,c,numClusters-c);
  		System.arraycopy(numDownData,c+1,numDownData,c,numClusters-c);
  		System.arraycopy(numDownTables,c+1,numDownTables,c,numClusters-c);
	}
	
	public void growTopics() {
		super.growTopics();
  		theta = myModel.growDoubleArray(theta);
  		numUpData = myModel.growIntArray(numUpData);
  		numUpTables = myModel.growIntArray(numUpTables);
  		numDownData = myModel.growIntArray(numDownData);
  		numDownTables = myModel.growIntArray(numDownTables);
  	}
	
	public void addChild(DirichletProcess c) {
  		if (children == null) {
  			children = new ArrayList<DirichletProcess>();
  		}
  		children.add(c);
  	}
  	
  	public int numChildren() {
  		if (children == null) {
  			return 0;
  		} else {
  			return children.size();
  		}
  	}
  	
  	public void removeChild(DirichletProcess c) {
  		children.remove(c);
  	}
  	
  	public String getLabel() {
  		if (label.length() > 0) {
  			return label;
  		}
  		
  		String label2 = "GRP";
  		int i = 0;
  		if (numChildren() > 0) {
  			for (i=0;i<numChildren();i++) {
  				label2 = label2 + "_" + children.get(i).getLabel();
  			}
  		}
  		return label2;
  	}
  	
  	public void activate() {
  		super.activate();
  		theta = new double[myModel.allocTopics];
  		numUpData = new int[myModel.allocTopics];
  		numUpTables = new int[myModel.allocTopics];
  		numDownData = new int[myModel.allocTopics];
  		numDownTables = new int[myModel.allocTopics];
  	}
  	
  	public void resampleTheta(int numClusters,double[] theta2) {
  		int cc = 0;
  		for (cc=0;cc<numClusters+1;cc++) {
  			theta2[cc] = sampleSingleTheta(cc);
  		}
  	}
  	
  	public void resampleTheta(int numClusters) {
  		resampleTheta(numClusters,theta);
  		
  	/*	int tt = 0;
  		for (tt = 0;tt<numClusters+1;tt++) {
  			System.out.print(theta[tt]);
  			System.out.print(" ");
  		}
  		System.out.println(""); */
  	}
  	
  	public double sampleSingleTheta(int cc) {
  		ud[0] = (double) numUpData[cc];
		ud[1] = (double) numDownData[cc];
		if (parent == null) {
			ud[0] += myModel.UD_a;
			ud[1] += myModel.UD_b;
		} else {
			ud[0] += alpha_theta * parent.theta[cc];
			ud[1] += alpha_theta * (1-parent.theta[cc]);
		}
		StatUtil.dirichlet_rnd(udt,ud,2);
		return udt[0];
  	}
  	
  	public void resampleParameters(int numClusters,double[] condLike) {
  		super.resampleParameters(numClusters,condLike);
  		resampleTheta(numClusters);
  	}
  	
  	public void sampleNewTheta(int numClusters) {
  		theta[numClusters] = sampleSingleTheta(numClusters);
  	}
  	
  	public void sampleNewParameters(int numClusters) {
  		super.sampleNewParameters(numClusters);
  		sampleNewTheta(numClusters);
  	}
  	
  	public void sampleClusterNumTables(int numClusters) {
  		super.sampleClusterNumTables(numClusters);
  		sampleUDNumTables(numClusters);
  	}
  	
  	public void sampleUDNumTables(int numClusters) {
  		int cc = 0;
  		
  		if (parent == null) {
  			for (cc=0;cc<numClusters;cc++) {
  				if (numUpData[cc] > 0) {
  					numUpTables[cc] = 1;
  				}
  				if (numDownData[cc] > 0) {
  					numDownTables[cc] = 1;
  				}
  			}
  		}
  		else {
  			for (cc=0;cc<numClusters;cc++) {
  				parent.numUpData[cc] = parent.numUpData[cc] - numUpTables[cc];
  				parent.numDownData[cc] = parent.numDownData[cc] - numDownTables[cc];
  				
  				numUpTables[cc] = numTables_rnd(alpha_theta*parent.theta[cc],numUpData[cc]);
  				parent.numUpData[cc] += numUpTables[cc];  
  				
  				numDownTables[cc] = numTables_rnd(alpha_theta*(1.0-parent.theta[cc]),numDownData[cc]);
  				parent.numDownData[cc] += numDownTables[cc];  
  			}
  		}
  	}
}
