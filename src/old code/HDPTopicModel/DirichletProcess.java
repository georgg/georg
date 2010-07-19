package edu.mit.csail.psrg.georg.HDPTopicModel;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedHashSet;

import cern.jet.random.Gamma;
import cern.jet.random.Uniform;
import edu.mit.csail.psrg.georg.StatUtil.StatUtil;
import java.io.Serializable;

public class DirichletProcess implements Serializable {
	double alpha;
	double[] beta;
  	int[] clusterNumData;
  	int[] clusterNumTables;
  	Document[] documents = null;
  	DirichletProcess parent = null;
  	HDPConcParam concParam = null;
  	String label;
  	ArrayList<DirichletProcess> children = null;
  	int dpID = 0;
  	static int ACTIVE = 0;
  	static int HELDOUT = 1;
  	static int FROZEN = 2;
  	int state = 0;
  	boolean generateNewClusters = true;
  	int code = 0;
  	
  	DirichletProcess(DirichletProcess p, String l) {
  		parent = p;
  		label = l;
  		alpha = 1.0;
  		if (p != null) {
  			p.addChild(this);
  		}
  	}
  	
  	public int getCode() { return code; }
  	public void setCode(int c) { code = c; }
  	
  	void addChild(DirichletProcess c) {
  		if (children == null) {
  			children = new ArrayList<DirichletProcess>();
  		}
  		children.add(c);
  	}
  	
  	int numChildren() {
  		if (children == null) {
  			return 0;
  		} else {
  			return children.size();
  		}
  	}
  	
  	void removeChild(DirichletProcess c) {
  		children.remove(c);
  	}
  	
  	void addDocument(Document d) {
  		ArrayList<Document> dd = new ArrayList<Document>();
  		dd.add(d);
  		addDocuments(dd);
  		d.parent = this;
  	}
  	
  	int numDocuments() {
  		if (documents == null) {
  			return 0;
  		} else {
  			return documents.length;
  		}
  	}
  	
  	void removeDocument(Document doc) {
  		int docNum = findDocNum(doc);
		if (docNum == -1) {
			return;
		}
		
		if (documents.length - 1 == 0) {
			documents = null;
			return;
		}
		
		Document[] documents2 = new Document[documents.length-1];
		int i = 0;
		for(i=0;i<docNum;i++) {
			documents2[i] = documents[i];
		}
		for(i=docNum+1;i<documents.length;i++) {
			documents2[i-1] = documents[i];
		}
		
	/*	if (docNum > 0) {
			System.arraycopy(documents,0,documents2,0,docNum);
		}
		System.arraycopy(documents,docNum+1,documents2,docNum,documents.length-docNum-1); */
		documents = documents2;
  	}
  	
  	void appendDocument(Document doc) {
  		if (documents == null) {
  			documents = new Document[1];
  			documents[0] = doc;
  			return;
  		}
  		
  		Document[] documents2 = new Document[documents.length+1];
  		System.arraycopy(documents,0,documents2,0,documents.length);
  		documents = documents2;
  		documents[documents.length-1] = doc;
  	}
  	
  	int findDocNum(Document doc) {
  		if (documents == null) {
  			return -1;
  		}
  		
		int i = 0;
		while(i<documents.length) {
			if (documents[i] == doc) {
				return i;
			}
			i++;
		}
		return -1;
	}
  	
  	void addDocuments(ArrayList<Document> dd) {
  		int i = 0;
  		documents = new Document[dd.size()];
  		
  		for (i=0;i<documents.length;i++) {
  			documents[i] = dd.get(i);
  			documents[i].parent = this;
  		}
  	}
  	
  	void deleteCluster(int c, int numClusters) {
  		int i = 0;
  		
  		if (documents != null) {
  			for (i=0;i<documents.length;i++) {
  				documents[i].deleteCluster(c,numClusters);
  			}
  		}
  		
  		System.arraycopy(beta,c+1,beta,c,numClusters-c);
  		System.arraycopy(clusterNumData,c+1,clusterNumData,c,numClusters-c);
  		System.arraycopy(clusterNumTables,c+1,clusterNumTables,c,numClusters-c);
  	}
  	
  	void growClusters(int growClustersAlloc) {
  		int currSize = beta.length;
  		
  		double[] beta2 = new double[currSize + growClustersAlloc];
  		System.arraycopy(beta,0,beta2,0,currSize);
  		beta = beta2;
  		
  		int[] clusterNumData2 = new int[currSize + growClustersAlloc];
  		System.arraycopy(clusterNumData,0,clusterNumData2,0,currSize);
  		clusterNumData = clusterNumData2;
  		
  		int[] clusterNumTables2 = new int[currSize + growClustersAlloc];
  		System.arraycopy(clusterNumTables,0,clusterNumTables2,0,currSize);
  		clusterNumTables = clusterNumTables2;
  		
  		if (documents != null) {
  			int dd = 0;
  			for (dd=0;dd<documents.length;dd++) {
  				documents[dd].growClusters(growClustersAlloc);
  			}
  		}
  	}
  	
  	void sampleNewBeta(int numClusters) {
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
  	
  	void resampleBeta(int numClusters,double[] condLike,double[] beta2) {
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
  	
  	void resampleBeta(int numClusters,double[] condLike) {
  		resampleBeta(numClusters,condLike,beta);
  	}
  	
  	int numTables_rnd(double alpha, int numData) {
  		int ii, numtable;

  	  	if ( numData == 0 ) {
  	    	return 0;
  	  	} else {
  	    	numtable = 1;
  	    	for ( ii = 1 ; ii < numData ; ii++ ) {
  	      		if ( Uniform.staticNextDouble() < alpha / (((double) ii)+alpha) ) numtable++;
  	    	}
  	    	return numtable;
  	  	}
  	}
  	
  	void sampleClusterNumTables(int numClusters) {
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
  	
  	void activate(HierDirichletProcess myModel) {
  		int numClusters = myModel.numClusters;
  		int i = 0;
  		int gg = 0;
  		int gene = 0;
  		int cc = 0;
  		int geneID = 0;
  		Document doc = null;
  		beta = new double[myModel.allocClusters];
  		clusterNumData = new int[myModel.allocClusters];
  		clusterNumTables = new int[myModel.allocClusters];
  		
  		if (documents != null) {
	  		for (i=0;i<documents.length;i++) {
	  			doc = documents[i];
	  			doc.clusterNumData = new int[myModel.allocClusters];
	  			for (gg=0;gg<(doc.genes.length);gg++) {
	  				gene = doc.genes[gg];
	  				cc = Uniform.staticNextIntFromTo(0,numClusters-1);
	  				doc.geneClusterAssigns[gg] = cc;
	  				myModel.posCounts[cc][gene]++;
	  				myModel.clusterCounts[cc]++;
	  			//	myModel.clusterMap[cc].add(new Integer(geneID));
	  				clusterNumData[cc]++;
	  				doc.clusterNumData[cc]++;
	  			}
	  		}
  		}
  		
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
  	
  	void holdout(HierDirichletProcess myModel) {
  		int cc = 0;
  		DirichletProcess ancestor = null;
  		if (state != DirichletProcess.HELDOUT) {
  			if (parent != null) {
  				for (cc=0;cc<clusterNumData.length;cc++) {
  					parent.clusterNumData[cc] = parent.clusterNumData[cc] - clusterNumTables[cc];
  				}
  				ancestor = parent;
  				while(ancestor != null) {
  					ancestor.sampleClusterNumTables(myModel.numClusters);
  					ancestor = ancestor.parent;
  				}
  			}
  		}
  		
  		int i = 0;
  		int gg = 0;
  		int gene = 0;
  		Document doc = null;
  		
  		if (documents != null) {
	  		for (i=0;i<documents.length;i++) {
	  			doc = documents[i];
	  			for (gg=0;gg<(doc.genes.length);gg++) {
	  				gene = doc.genes[gg];
	  				cc = doc.geneClusterAssigns[gg];
	  				myModel.posCounts[cc][gene]--;
	  				myModel.clusterCounts[cc]--;
	  				clusterNumData[cc]--;
	  				doc.clusterNumData[cc]--;
	  			}
	  		}
  		}
  		
  		state = DirichletProcess.HELDOUT;
  	}
  	
  	void freeze() {
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
  		
  		String label2 = "GRP";
  		int i = 0;
  		if (numChildren() > 0) {
  			for (i=0;i<numChildren();i++) {
  				label2 = label2 + "_" + children.get(i).getLabel();
  			}
  		} else {
  			for (i=0;i<numDocuments();i++) {
  				label2 = label2 + "_" + documents[i].annotation;
  			}
  		}
  		return label2;
  	}
}
