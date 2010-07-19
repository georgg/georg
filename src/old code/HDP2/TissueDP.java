package edu.mit.csail.psrg.georg.HDP2;

import java.util.ArrayList;
import java.util.HashSet;

import cern.jet.random.Uniform;
import edu.mit.csail.psrg.georg.StatUtil.StatUtil;

public class TissueDP extends DirichletProcess {
	/**
	 * 
	 */
	private static final long serialVersionUID = 3062155651551882080L;
	// map of expression units to genes
	public int[] genes;
	// whether each unit of expression is up or down for the topic
	public boolean[] expressionUp = null;
	// assigns each expression unit to a topic
	public int[] geneTopicAssigns;
	// number of unique genes with some expression in the topic
	public int numUniqueGenes = 0;
	// decision as to whether each topic contains up or down regulated expression
	public boolean[] upDownTopicChoice = null;
	// a count of how many genes are up-regulated in the given topic
	public int[] upTopicCount = null;
	
	// temporary variables used for sampling
	double pU = 0.0;
	double pD = 0.0;
	double minLik = 0.0;
	double u = 0.0;
	
	public void deleteTopic(int c, int numClusters) {
		super.deleteTopic(c,numClusters);
		int i = 0;
		
		for (i=0;i<geneTopicAssigns.length;i++) {
			if (geneTopicAssigns[i] > c) {
				geneTopicAssigns[i] = geneTopicAssigns[i] - 1;
			}
		}
		System.arraycopy(upDownTopicChoice,c+1,upDownTopicChoice,c,numClusters-c);
		System.arraycopy(upTopicCount,c+1,upTopicCount,c,numClusters-c);
	}
	
	public void growTopics() {
		super.growTopics();
  		upDownTopicChoice = myModel.growBooleanArray(upDownTopicChoice);
  		upTopicCount = myModel.growIntArray(upTopicCount);
  	}
	
	public void addTopic(int numTopics) {
		super.addTopic(numTopics);
  		if (state != DirichletProcess.HELDOUT) {
  			upDownTopicChoice[numTopics+1] = true;
  			upTopicCount[numTopics+1] = 0;
  		}
	}
	
	public TissueDP(HDP mm, ParentDP p, String l) {
		super(mm, p, l);
	}

	public TissueDP(HDP mm, ParentDP p, String l,ArrayList<Integer> glist, ArrayList<Boolean> udlist) {
		super(mm, p,l);
		
		int i = 0;
		
		genes = new int[glist.size()];
		expressionUp = new boolean[udlist.size()];
		geneTopicAssigns = new int[glist.size()];
		
		for (i=0;i<glist.size();i++) {
			genes[i] = glist.get(i);
			expressionUp[i] = udlist.get(i);
			geneTopicAssigns[i] = 0;
		}
		
		HashSet<Integer> uniqueGenes = new HashSet<Integer>();
		for (i=0;i<genes.length;i++) {
			uniqueGenes.add(genes[i]);
		}
		numUniqueGenes = uniqueGenes.size();
	}
	
	public void allocateGenesToTopicMap() {
		int i = 0;
		int gene = 0;
		
		ArrayList<Topic> topics = myModel.topics;
		
		for (i=0;i<topics.size();i++) {
			topics.get(i).topicMap.clear();
			topics.get(i).topicArray.clear();
		}
		
		for (i=0;i<genes.length;i++) {
			gene = genes[i];
			topics.get(geneTopicAssigns[i]).topicMap.add(new Integer(gene));
			topics.get(geneTopicAssigns[i]).topicArray.add(new Integer(gene));
		}
	}
	
	public void activate() {
  		int numTopics = myModel.numTopics;
  		int i = 0;
  		int gg = 0;
  		int gene = 0;
  		int cc = 0;
  		int geneID = 0;
  
  		beta = new double[myModel.allocTopics];
  		clusterNumData = new int[myModel.allocTopics];
  		clusterNumTables = new int[myModel.allocTopics];
  		upDownTopicChoice = new boolean[myModel.allocTopics];
  		upTopicCount = new int[myModel.allocTopics];
  		Topic myTopic = null;
  		
	  	for (gg=0;gg<genes.length;gg++) {
	  		gene = genes[gg];
	  		cc = Uniform.staticNextIntFromTo(0,numTopics-1);
	  		geneTopicAssigns[gg] = cc;
	  		myTopic = myModel.topics.get(cc);
	  		myTopic.addGeneCount(gene);
	  		clusterNumData[cc]++;
	  		if (expressionUp[gg]) {
	  			upTopicCount[cc]++;
	  		}
	  	}
  		
  		for (cc=0;cc<numTopics;cc++) {
  			clusterNumTables[cc] = clusterNumData[cc];
  			u = Uniform.staticNextDouble();
  			if (u <= (myModel.UD_a/(myModel.UD_a + myModel.UD_b))) {
  				upDownTopicChoice[cc] = true;
  				parent.numUpData[cc]++;
  			} else {
  				upDownTopicChoice[cc] = false;
  				parent.numDownData[cc]++;
  			}
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
		super.holdout();
		
  		int cc = 0;
  		int gg = 0;
  		int gene = 0;
  		Topic myTopic = null;
  
	  	for (gg=0;gg<genes.length;gg++) {
	  		gene = genes[gg];
	  		cc = geneTopicAssigns[gg];
	  		myTopic = myModel.topics.get(cc);
	  		myTopic.removeGeneCount(gene);
	  		clusterNumData[cc]--;
	  	}
  	}
	
	public void sampleDataTopicAssigns() {
  		int gg = 0;
  		int cc = 0;
  		int oldCC;
  		int gene = 0;
  		double mlt = 0.0;
  		double[] condLike = myModel.condLike;
  		
  		Topic oldTopic = null;
  		Topic newTopic = null;
  		
  		double matchProb = 0.0;
  		double integrateUDProb = 0.0;
  		double upChoicePart = 0.0;
  		double downChoicePart = 0.0;
  		double useTheta = 0.0;
  //	double upr = myModel.UD_a/(myModel.UD_a + myModel.UD_b);
  		
		for (gg=0;gg<genes.length;gg++) {
			oldCC = geneTopicAssigns[gg];
			gene = genes[gg];
			
			oldTopic = myModel.topics.get(oldCC);
			
			// remove the gene from consideration
			oldTopic.removeGeneCount(gene);
			clusterNumData[oldCC]--;
			if (expressionUp[gg]) {
				upTopicCount[oldCC]--;
			}
			
			for (cc=0;cc<=myModel.numTopics;cc++) {
				// this is the marginal conditional likelihood that the gene
				// comes from the cluster (we marginalize out the parameters under
				// a beta prior)
				newTopic = myModel.topics.get(cc);
				matchProb = myModel.upProb;
				
				if (cc < myModel.numTopics) {
					useTheta = parent.theta[cc];
				} else {
					useTheta = StatUtil.beta_rnd(myModel.UD_a,myModel.UD_b);
				}	
				
			/*	upChoicePart = ((double) upTopicCount[cc])*Math.log(matchProb) + ((double) clusterNumData[cc] - upTopicCount[cc])*Math.log(1.0-matchProb);
				downChoicePart = ((double) upTopicCount[cc])*Math.log(1.0-matchProb) + ((double) clusterNumData[cc] - upTopicCount[cc])*Math.log(matchProb);
				if (expressionUp[gg]) {
					upChoicePart += Math.log(matchProb);
					downChoicePart += Math.log(1.0-matchProb);
				} else {
					upChoicePart += Math.log(1.0-matchProb);
					downChoicePart += Math.log(matchProb);
				}
				integrateUDProb = useTheta*Math.exp(upChoicePart) + (1.0-useTheta)*Math.exp(downChoicePart); */
				
				condLike[cc] = (((double) newTopic.posCounts[gene]) + myModel.priorPos)/(((double) newTopic.totalCount) + myModel.priorPos*myModel.totalGenes);
				
				if (cc < myModel.numTopics) {
					if ((expressionUp[gg] & !upDownTopicChoice[cc]) | (!expressionUp[gg] & upDownTopicChoice[cc])) {
						matchProb = 1.0 - matchProb;
					}
				} else {
					if (expressionUp[gg]) {
						matchProb = useTheta*matchProb + (1.0-useTheta)*(1.0-matchProb);
					} else {
						matchProb = useTheta*(1.0-matchProb) + (1.0-useTheta)*(matchProb);
					}
				}
				
				condLike[cc] = condLike[cc] * matchProb;
				
			//	condLike[cc] = condLike[cc] * integrateUDProb;
				
				mlt = clusterNumData[cc];
				if (parent == null) {
					if (cc == myModel.numTopics) {
						mlt = mlt + alpha;
					}
				}
				else {
					mlt = mlt + alpha*parent.beta[cc];
				}
				condLike[cc] = condLike[cc]*mlt;
			}
			
			cc = StatUtil.multinomial_rnd(condLike,myModel.numTopics+1);
			
			// add the gene to the cluster it just joined
			newTopic = myModel.topics.get(cc);
			newTopic.addGeneCount(gene);
			clusterNumData[cc]++;
			geneTopicAssigns[gg] = cc;
			if (expressionUp[gg]) {
				upTopicCount[cc]++;
			}
			
			// if it's a new cluster, we'll have to expand the # of clusters
			if (cc==myModel.numTopics) {
				myModel.addTopic();
				condLike = myModel.condLike;
			}
		}
  	}
	
/*	public void sampleUpDownTopicPrior(int numTopics) {
		double u = 0.0;
		u = Uniform.staticNextDouble();
		if (u <= parent.theta[numTopics]) {
			upDownTopicChoice[numTopics] = true;
			parent.numUpData[numTopics]++;
		} else {
			upDownTopicChoice[numTopics] = false;
			parent.numDownData[numTopics]++;
		}
	} */
	
	public void resampleParameters(int numClusters,double[] condLike) {
  		super.resampleParameters(numClusters,condLike);
  		sampleUpDownTopicPosterior();
	}
	
	public void sampleNewParameters(int numClusters) {
  		super.sampleNewParameters(numClusters);
  		sampleSingleUpDownTopicPosterior(numClusters);
  		if (upDownTopicChoice[numClusters]) {
			parent.numUpData[numClusters]++;
		} else {
			parent.numDownData[numClusters]++;
		}
	}
	
	public void sampleUpDownTopicPosterior() {
		int cc = 0;
		for (cc=0;cc<myModel.numTopics+1;cc++) {
			sampleSingleUpDownTopicPosterior(cc);
		}
	}
	
	public void sampleSingleUpDownTopicPosterior(int cc) {
		if (cc < myModel.numTopics) {
			if (upDownTopicChoice[cc]) {
				parent.numUpData[cc]--;
			} else {
				parent.numDownData[cc]--;
			}
		}
		
		if (cc < myModel.numTopics) {
			pU = ((double) upTopicCount[cc])*Math.log(myModel.upProb) + ((double) (clusterNumData[cc]-upTopicCount[cc]))*Math.log(1.0-myModel.upProb);
			pD = ((double) upTopicCount[cc])*Math.log(1.0-myModel.upProb) + ((double) (clusterNumData[cc]-upTopicCount[cc]))*Math.log(myModel.upProb);
			if (pU > pD) {
				minLik = pU;
			} else { minLik = pD; }
			pU = pU - minLik;
			pD = pD - minLik;
		} else {
			pU = 0;
			pD = 0;
		}
		
		pU = parent.theta[cc]*Math.exp(pU);
		pD = (1.0-parent.theta[cc])*Math.exp(pD);
		
		if (Double.isInfinite(-pD) | Double.isNaN(pD)) {
			pD = 0;
		} else {
			if (Double.isInfinite(pD)) {
				pD = 1.0;
			}
		}
		
		if (Double.isInfinite(-pU) | Double.isNaN(pU)) {
			pU = 0;
		} else {
			if (Double.isInfinite(pU)) {
				pU = 1.0;
			}
		}
		
	/*	System.out.println("pU = " + pU);
		System.out.println("pD = " + pD); */
		
		pU = pU/(pU + pD);
		
		if (Double.isNaN(pU)) {
			pU = 0.5;
		}
		
		u = Uniform.staticNextDouble();
		if (u <= pU) {
			upDownTopicChoice[cc] = true;
		} else {
			upDownTopicChoice[cc] = false;
		}
		
		if (cc < myModel.numTopics) {
			if (upDownTopicChoice[cc]) {
				parent.numUpData[cc]++;
			} else {
				parent.numDownData[cc]++;
			}
		}
	}
}