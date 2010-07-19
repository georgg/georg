package edu.mit.csail.psrg.georg.GeneProgram2;

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
	//	whether each unit of expression is up or down for the program
	public boolean[] expressionUp = null;
	// the modifier(s) for each unit of expression
	public int[][] expressionModify = null;
	// assigns each expression unit to an expression program
	public int[] geneExpressionProgramAssigns;
	// number of unique genes with some expression in the program
	public int numUniqueGenes = 0;
	// the modifier chosen for the expression program (dim. 1 = modifiers; dim. 2 = programs)
	public int[][] modifierExpressionProgramChoice = null;
	// a count of how many genes have a particular modifier in the given program
	// (dim. 1 = modifiers; dim. 2 = modifier values; dim. 3 = programs) 
	public int[][][] modifierExpressionProgramCount = null;
	public boolean[] upDownExpressionProgramChoice = null;
	// a count of how many genes are up-regulated in the given topic
	public int[] upExpressionProgramCount = null;
	
	// temporary variables used for sampling
	double pU = 0.0;
	double pD = 0.0;
	double minLik = 0.0;
	double u = 0.0;
	
	public void deleteExpressionProgram(int c, int numClusters) {
		super.deleteExpressionProgram(c,numClusters);
		int i = 0;
		int j = 0;
		
		for (i=0;i<geneExpressionProgramAssigns.length;i++) {
			if (geneExpressionProgramAssigns[i] > c) {
				geneExpressionProgramAssigns[i] = geneExpressionProgramAssigns[i] - 1;
			}
		}
		
		if (myModel.useUpDown) {
			System.arraycopy(upDownExpressionProgramChoice,c+1,upDownExpressionProgramChoice,c,numClusters-c);
			System.arraycopy(upExpressionProgramCount,c+1,upExpressionProgramCount,c,numClusters-c);
		}
		
		if (modifierExpressionProgramChoice != null) {
			for (i=0;i<modifierExpressionProgramChoice.length;i++) {
				System.arraycopy(modifierExpressionProgramChoice[i],c+1,modifierExpressionProgramChoice[i],c,numClusters-c);
				for (j=0;j<modifierExpressionProgramCount[i].length;j++) {
					System.arraycopy(modifierExpressionProgramCount[i][j],c+1,modifierExpressionProgramCount[i][j],c,numClusters-c);
				}
			}
		}
	}
	
	public void growExpressionPrograms() {
		super.growExpressionPrograms();
		
		if (myModel.useUpDown) {
			upDownExpressionProgramChoice = myModel.growBooleanArray(upDownExpressionProgramChoice);
			upExpressionProgramCount = myModel.growIntArray(upExpressionProgramCount);
		}
  		
  		if (modifierExpressionProgramChoice != null) {
  			int i = 0;
  	  		int j = 0;
  			for (i=0;i<modifierExpressionProgramChoice.length;i++) {
  				modifierExpressionProgramChoice[i] = myModel.growIntArray(modifierExpressionProgramChoice[i]);
  				for (j=0;j<modifierExpressionProgramCount[i].length;j++) {
  					modifierExpressionProgramCount[i][j] = myModel.growIntArray(modifierExpressionProgramCount[i][j]);
  				}
  			}
  		}
  	}
	
	public void addExpressionProgram(int numTopics) {
		super.addExpressionProgram(numTopics);
  		if (state != DirichletProcess.HELDOUT) {
  			
  			if (myModel.useUpDown) {
  				upDownExpressionProgramChoice[numTopics+1] = true;
  				upExpressionProgramCount[numTopics+1] = 0;
  			}
  			
  			if (modifierExpressionProgramChoice != null) {
  				int i = 0;
  	  	  		int j = 0;
  	  			for (i=0;i<modifierExpressionProgramChoice.length;i++) {
  	  				modifierExpressionProgramChoice[i][numTopics+1] = 0;
  	  				for (j=0;j<modifierExpressionProgramCount[i].length;j++) {
  	  					modifierExpressionProgramCount[i][j][numTopics+1] = 0;
  	  				}
  	  			}
  			}
  		}
	}
	
	public TissueDP(HDP mm, ParentDP p, String l) {
		super(mm, p, l);
	}

	public TissueDP(HDP mm, ParentDP p, String l,ArrayList<Integer> glist, ArrayList<Boolean> udlist, ArrayList<ArrayList<Integer> > modifierList) {
		super(mm, p,l);
		
		int i = 0;
		int j = 0;
		
		genes = new int[glist.size()];
		
		if (udlist != null) {
			expressionUp = new boolean[udlist.size()];
		}
		
		if (modifierList != null) {
			expressionModify = new int[modifierList.size()][glist.size()];
		}
		geneExpressionProgramAssigns = new int[glist.size()];
		
		for (i=0;i<glist.size();i++) {
			genes[i] = glist.get(i);
			geneExpressionProgramAssigns[i] = 0;
			
			if (udlist != null) {
				expressionUp[i] = udlist.get(i);
			}	
		}
		
		if (modifierList != null) {
			ArrayList<Integer> tlist = null;
			for (j=0;j<modifierList.size();j++) {
				tlist = modifierList.get(j);
				for (i=0;i<tlist.size();i++) {
					expressionModify[j][i] = tlist.get(i);
				}
			}
		}
		
		HashSet<Integer> uniqueGenes = new HashSet<Integer>();
		for (i=0;i<genes.length;i++) {
			uniqueGenes.add(genes[i]);
		}
		numUniqueGenes = uniqueGenes.size();
	}
	
	public void allocateGenesToExpressionProgramMap() {
		int i = 0;
		int gene = 0;
		
		ArrayList<ExpressionProgram> programs = myModel.expressionPrograms;
		
		for (i=0;i<programs.size();i++) {
			programs.get(i).expressionProgramMap.clear();
			programs.get(i).expressionProgramArray.clear();
		}
		
		for (i=0;i<genes.length;i++) {
			gene = genes[i];
			programs.get(geneExpressionProgramAssigns[i]).expressionProgramMap.add(new Integer(gene));
			programs.get(geneExpressionProgramAssigns[i]).expressionProgramArray.add(new Integer(gene));
		}
	}
	
	public void activate() {
  		int numTopics = myModel.numExpressionPrograms;
  		int i = 0;
  		int gg = 0;
  		int gene = 0;
  		int cc = 0;
  		int geneID = 0;
  
  		beta = new double[myModel.allocExpressionPrograms];
  		clusterNumData = new int[myModel.allocExpressionPrograms];
  		clusterNumTables = new int[myModel.allocExpressionPrograms];
  		
  		if (myModel.useUpDown) {
  			upDownExpressionProgramChoice = new boolean[myModel.allocExpressionPrograms];
  			upExpressionProgramCount = new int[myModel.allocExpressionPrograms];
  		}
  		
  		ExpressionProgram myTopic = null;
  		
  		if (expressionModify != null) {
  			modifierExpressionProgramChoice = new int[expressionModify.length][myModel.allocExpressionPrograms];
  			modifierExpressionProgramCount = new int[expressionModify.length][][];
  			for (i=0;i<expressionModify.length;i++) {
  				modifierExpressionProgramCount[i] = new int[myModel.modifierLevels[i]][myModel.allocExpressionPrograms];
  			}
  		}
  		
	  	for (gg=0;gg<genes.length;gg++) {
	  		gene = genes[gg];
	  		cc = Uniform.staticNextIntFromTo(0,numTopics-1);
	  		geneExpressionProgramAssigns[gg] = cc;
	  		myTopic = myModel.expressionPrograms.get(cc);
	  		myTopic.addGeneCount(gene);
	  		clusterNumData[cc]++;
	  		
	  		if (myModel.useUpDown) {
	  			if (expressionUp[gg]) {
	  				upExpressionProgramCount[cc]++;
	  			}
	  		}
	  		
	  		if (expressionModify != null) {
	  			for (i=0;i<expressionModify.length;i++) {
	  				modifierExpressionProgramCount[i][expressionModify[i][gg]][cc]++;
	  			}
	  		}
	  	}
  		
  		for (cc=0;cc<numTopics;cc++) {
  			clusterNumTables[cc] = clusterNumData[cc];
  			u = Uniform.staticNextDouble();
  			
  			if (myModel.useUpDown) {
  				if (u <= (myModel.UD_a/(myModel.UD_a + myModel.UD_b))) {
  					upDownExpressionProgramChoice[cc] = true;
  					parent.numUpData[cc]++;
  				} else {
  					upDownExpressionProgramChoice[cc] = false;
  					parent.numDownData[cc]++;
  				}
  			}
  			
  			if (expressionModify != null) {
  				for (i=0;i<expressionModify.length;i++) {
  					modifierExpressionProgramChoice[i][cc] = StatUtil.multinomial_rnd(myModel.modifierPrior[i],myModel.modifierPrior[i].length);
  					parent.numModifierData[i][modifierExpressionProgramChoice[i][cc]][cc]++;
  				}
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
  		ExpressionProgram myTopic = null;
  
	  	for (gg=0;gg<genes.length;gg++) {
	  		gene = genes[gg];
	  		cc = geneExpressionProgramAssigns[gg];
	  		myTopic = myModel.expressionPrograms.get(cc);
	  		myTopic.removeGeneCount(gene);
	  		clusterNumData[cc]--;
	  	}
  	}
	
	public void sampleDataExpressionProgramAssigns() {
  		int gg = 0;
  		int cc = 0;
  		int oldCC;
  		int gene = 0;
  		double mlt = 0.0;
  		double[] condLike = myModel.condLike;
  		
  		ExpressionProgram oldTopic = null;
  		ExpressionProgram newTopic = null;
  		
  		double matchProb = 0.0;
  		double useTheta = 0.0;
  
  		double[][] useThetaMod = null;
  		double[] matchProbMod = null;
  		int i = 0;
  		int j = 0;
  		if (expressionModify != null) {
  			useThetaMod = new double[expressionModify.length][];
  			matchProbMod = new double[expressionModify.length];
  			for (i=0;i<expressionModify.length;i++) {
  				useThetaMod[i] = new double[myModel.modifierLevels[i]];
  			}
  		}
  		
		for (gg=0;gg<genes.length;gg++) {
			oldCC = geneExpressionProgramAssigns[gg];
			gene = genes[gg];
			
			oldTopic = myModel.expressionPrograms.get(oldCC);
			
			// remove the gene from consideration
			oldTopic.removeGeneCount(gene);
			clusterNumData[oldCC]--;
			
			if (myModel.useUpDown) {
				if (expressionUp[gg]) {
					upExpressionProgramCount[oldCC]--;
				}
			}
			
			if (expressionModify != null) {
				for (i=0;i<expressionModify.length;i++) {
					modifierExpressionProgramCount[i][expressionModify[i][gg]][oldCC]--;
				}
			}
			
			for (cc=0;cc<=myModel.numExpressionPrograms;cc++) {
				// this is the marginal conditional likelihood that the gene
				// comes from the cluster (we marginalize out the parameters under
				// a beta prior)
				newTopic = myModel.expressionPrograms.get(cc);
				matchProb = myModel.upProb;
				
				if (cc < myModel.numExpressionPrograms) {
					if (myModel.useUpDown) {
						useTheta = parent.theta[cc];
					}
					
					if (expressionModify != null) {
						for (i=0;i<expressionModify.length;i++) {
							for (j=0;j<myModel.modifierLevels[i];j++) {
								useThetaMod[i][j] = parent.thetaMod[i][j][cc];
							}
						}
					}
				} else {
					if (myModel.useUpDown) {
						useTheta = StatUtil.beta_rnd(myModel.UD_a,myModel.UD_b);
					}
					
					if (expressionModify != null) {
						for (i=0;i<expressionModify.length;i++) {
							StatUtil.dirichlet_rnd(useThetaMod[i],myModel.modifierPrior[i],myModel.modifierLevels[i]);
						}
					}
				}	
				
				condLike[cc] = (((double) newTopic.posCounts[gene]) + myModel.priorPos)/(((double) newTopic.totalCount) + myModel.priorPos*myModel.totalGenes);
				
				if (cc < myModel.numExpressionPrograms) {
					
					if (myModel.useUpDown) {
						if ((expressionUp[gg] & !upDownExpressionProgramChoice[cc]) | (!expressionUp[gg] & upDownExpressionProgramChoice[cc])) {
							matchProb = 1.0 - matchProb;
						}
					}
					
					if (expressionModify != null) {
						for (i=0;i<expressionModify.length;i++) {
							if (modifierExpressionProgramChoice[i][cc] == expressionModify[i][gg]) {
								matchProbMod[i] = myModel.matchProbMod[i];
							} else {
								matchProbMod[i] = myModel.unmatchProbMod[i];
							}
						}
					}
					
				} else {
					if (myModel.useUpDown) {
						if (expressionUp[gg]) {
							matchProb = useTheta*matchProb + (1.0-useTheta)*(1.0-matchProb);
						} else {
							matchProb = useTheta*(1.0-matchProb) + (1.0-useTheta)*(matchProb);
						}
					}
					
					if (expressionModify != null) {
						for (i=0;i<expressionModify.length;i++) {
							matchProbMod[i] = 0.0;
							for (j=0;j<myModel.modifierLevels[i];j++) {
								if (expressionModify[i][gg] == j) {
									matchProbMod[i] += useThetaMod[i][j]*myModel.matchProbMod[i];
								} else {
									matchProbMod[i] += useThetaMod[i][j]*myModel.unmatchProbMod[i];
								}
							}
						}
					}
				}
				
				if (myModel.useUpDown) {
					condLike[cc] = condLike[cc] * matchProb;
				}
				
				if (expressionModify != null) {
					for (i=0;i<expressionModify.length;i++) {
						condLike[cc] = condLike[cc]*matchProbMod[i];
					}
				}
				
				mlt = clusterNumData[cc];
				if (parent == null) {
					if (cc == myModel.numExpressionPrograms) {
						mlt = mlt + alpha;
					}
				}
				else {
					mlt = mlt + alpha*parent.beta[cc];
				}
				condLike[cc] = condLike[cc]*mlt;
			}
			
			cc = StatUtil.multinomial_rnd(condLike,myModel.numExpressionPrograms+1);
			
			// add the gene to the cluster it just joined
			newTopic = myModel.expressionPrograms.get(cc);
			newTopic.addGeneCount(gene);
			clusterNumData[cc]++;
			geneExpressionProgramAssigns[gg] = cc;
			
			if (myModel.useUpDown) {
				if (expressionUp[gg]) {
					upExpressionProgramCount[cc]++;
			
				}
			}
			
			if (expressionModify != null) {
				for (i=0;i<expressionModify.length;i++) {
					modifierExpressionProgramCount[i][expressionModify[i][gg]][cc]++;
				}
			}
			
			// if it's a new cluster, we'll have to expand the # of clusters
			if (cc==myModel.numExpressionPrograms) {
				myModel.addExpressionProgram();
				condLike = myModel.condLike;
			}
		}
  	}
	
	public void resampleParameters(int numClusters,double[] condLike) {
  		super.resampleParameters(numClusters,condLike);
  		sampleModifierExpressionProgramPosterior();
	}
	
	public void sampleNewParameters(int numClusters) {
  		super.sampleNewParameters(numClusters);
  		sampleSingleModifiersExpressionProgramPosterior(numClusters);
  		
  		if (myModel.useUpDown) {
  			if (upDownExpressionProgramChoice[numClusters]) {
  				parent.numUpData[numClusters]++;
  			} else {
  				parent.numDownData[numClusters]++;
  			}
  		}
  		
  		if (expressionModify != null) {
  			int i = 0;
  			for (i=0;i<expressionModify.length;i++) {
  				parent.numModifierData[i][modifierExpressionProgramChoice[i][numClusters]][numClusters]++;
  			}
  		}
	}
	
	public void sampleModifierExpressionProgramPosterior() {
		int cc = 0;
		for (cc=0;cc<myModel.numExpressionPrograms+1;cc++) {
			sampleSingleModifiersExpressionProgramPosterior(cc);
		}
	}
	
	public void sampleSingleModifiersExpressionProgramPosterior(int cc) {
		int i = 0;
		int j = 0;
		double[] pC = null;
		
		if (cc < myModel.numExpressionPrograms) {
			
			if (myModel.useUpDown) {
				if (upDownExpressionProgramChoice[cc]) {
					parent.numUpData[cc]--;
				} else {
					parent.numDownData[cc]--;
				}
			}
			
			if (expressionModify != null) {
	  			for (i=0;i<expressionModify.length;i++) {
	  				parent.numModifierData[i][modifierExpressionProgramChoice[i][cc]][cc]--;
	  			}
	  		}
		}
		
		if (myModel.useUpDown) {
			if (cc < myModel.numExpressionPrograms) {
				pU = ((double) upExpressionProgramCount[cc])*Math.log(myModel.upProb) + ((double) (clusterNumData[cc]-upExpressionProgramCount[cc]))*Math.log(1.0-myModel.upProb);
				pD = ((double) upExpressionProgramCount[cc])*Math.log(1.0-myModel.upProb) + ((double) (clusterNumData[cc]-upExpressionProgramCount[cc]))*Math.log(myModel.upProb);
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
				upDownExpressionProgramChoice[cc] = true;
			} else {
				upDownExpressionProgramChoice[cc] = false;
			}
		}
		
		if (expressionModify != null) {
			for (i=0;i<expressionModify.length;i++) {
				pC = new double[myModel.modifierLevels[i]];
				if (cc < myModel.numExpressionPrograms) {
					for (j=0;j<myModel.modifierLevels[i];j++) {
						pC[j] = ((double) modifierExpressionProgramCount[i][j][cc])*Math.log(myModel.matchProbMod[i]) + ((double) (clusterNumData[cc] - modifierExpressionProgramCount[i][j][cc]))*Math.log(myModel.unmatchProbMod[i]);
					}
					minLik = pC[0];
					for (j=0;j<pC.length;j++) {
						if (minLik < pC[j]) {
							minLik = pC[j];
						}
					}
					for (j=0;j<pC.length;j++) {
						pC[j] = pC[j] - minLik;
					}
				}
				
				for (j=0;j<pC.length;j++) {
					pC[j] = parent.thetaMod[i][j][cc]*Math.exp(pC[j]);
				}
				
				modifierExpressionProgramChoice[i][cc] = StatUtil.multinomial_rnd(pC,myModel.modifierLevels[i]);
			}
		}
		
		if (cc < myModel.numExpressionPrograms) {
			if (myModel.useUpDown) {
				if (upDownExpressionProgramChoice[cc]) {
					parent.numUpData[cc]++;
				} else {
					parent.numDownData[cc]++;
				}
			}
			
			if (expressionModify != null) {
	  			for (i=0;i<expressionModify.length;i++) {
	  				parent.numModifierData[i][modifierExpressionProgramChoice[i][cc]][cc]++;
	  			}
	  		}
		}
	}
}