package edu.mit.csail.psrg.georg.GeneProgram2;

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
	
	// used for sampling thetaMod
	// dim. 1 = modifiers, dim. 2 = modifier level, dim. 3 = expression programs
	public int[][][] numModifierData = null;
	public int[][][] numModifierTables = null;
	
	// the probabilities for choosing whether
	// expression programs represent up or down regulated genes
	public double[] theta = null;
	// corresponding probabilities for modifiers
	// (dim. 1 = modifiers, dim. 2 = modifier levels, dim. 3 = programs)
	public double[][][] thetaMod = null;
	// parameter for pushing down influence of parents
	// on theta distribution
	public double alpha_theta = 1.0;
	public UpDownConcParam myUpDownConcParam = null;
	public double[] alpha_thetaMod = null;
	public ModifierConcParam myModifierConcParam = null;
	// temporary variables used for sampling theta
	public double[] ud = new double[2];
	public double[] udt = new double[2];
	public double[][] md = null;
	public double[][] mdt = null;
	
	public ParentDP(HDP mm, ParentDP p, String l) {
		super(mm, p, l);
	}

	public void addExpressionProgram(int numTopics) {
		super.addExpressionProgram(numTopics);
  		if (state != DirichletProcess.HELDOUT) {
  			
  			if (myModel.useUpDown) {
	  			numUpData[numTopics+1] = 0;
	  			numUpTables[numTopics+1] = 0;
	  			numDownData[numTopics+1] = 0;
	  			numDownTables[numTopics+1] = 0;
	  			theta[numTopics+1] = 0.0;
  			}
  			
  			if (myModel.modifierLevels != null) {
  				int i = 0;
  				int j = 0;
  				for (i=0;i<myModel.modifierLevels.length;i++) {
  					for (j=0;j<myModel.modifierLevels[i];j++) {
  						numModifierData[i][j][numTopics+1] = 0;
  						numModifierTables[i][j][numTopics+1] = 0;
  						thetaMod[i][j][numTopics+1] = 0.0;
  					}
  				}
  			}
  		}
	}
	
	public void deleteExpressionProgram(int c, int numClusters) {
		super.deleteExpressionProgram(c,numClusters);
		
		if (myModel.useUpDown) {
	  		System.arraycopy(theta,c+1,theta,c,numClusters-c);
	  		System.arraycopy(numUpData,c+1,numUpData,c,numClusters-c);
	  		System.arraycopy(numUpTables,c+1,numUpTables,c,numClusters-c);
	  		System.arraycopy(numDownData,c+1,numDownData,c,numClusters-c);
	  		System.arraycopy(numDownTables,c+1,numDownTables,c,numClusters-c);
		}
  		
  		if (myModel.modifierLevels != null) {
			int i = 0;
			int j = 0;
			for (i=0;i<myModel.modifierLevels.length;i++) {
				for (j=0;j<myModel.modifierLevels[i];j++) {
					System.arraycopy(numModifierData[i][j],c+1,numModifierData[i][j],c,numClusters-c);
					System.arraycopy(numModifierTables[i][j],c+1,numModifierTables[i][j],c,numClusters-c);
					System.arraycopy(thetaMod[i][j],c+1,thetaMod[i][j],c,numClusters-c);
				}
			}
		}
	}
	
	public void growExpressionPrograms() {
		super.growExpressionPrograms();
		
		if (myModel.useUpDown) {
	  		theta = myModel.growDoubleArray(theta);
	  		numUpData = myModel.growIntArray(numUpData);
	  		numUpTables = myModel.growIntArray(numUpTables);
	  		numDownData = myModel.growIntArray(numDownData);
	  		numDownTables = myModel.growIntArray(numDownTables);
		}
  		
  		if (myModel.modifierLevels != null) {
			int i = 0;
			int j = 0;
			for (i=0;i<myModel.modifierLevels.length;i++) {
				for (j=0;j<myModel.modifierLevels[i];j++) {
					numModifierData[i][j] = myModel.growIntArray(numModifierData[i][j]);
					numModifierTables[i][j] = myModel.growIntArray(numModifierTables[i][j]);
					thetaMod[i][j] = myModel.growDoubleArray(thetaMod[i][j]);
				}
			}
		}
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
  		
  		if (myModel.useUpDown) {
	  		theta = new double[myModel.allocExpressionPrograms];
	  		numUpData = new int[myModel.allocExpressionPrograms];
	  		numUpTables = new int[myModel.allocExpressionPrograms];
	  		numDownData = new int[myModel.allocExpressionPrograms];
	  		numDownTables = new int[myModel.allocExpressionPrograms];
  		}
  		
  		if (myModel.modifierLevels != null) {
  			alpha_thetaMod = new double[myModel.modifierLevels.length];
  			numModifierData = new int[myModel.modifierLevels.length][][];
  			numModifierTables = new int[myModel.modifierLevels.length][][];
  			thetaMod = new double[myModel.modifierLevels.length][][];
  			int i = 0;
  			md = new double[myModel.modifierLevels.length][];
  			mdt = new double[myModel.modifierLevels.length][];
  			for (i=0;i<numModifierData.length;i++) {
  				md[i] = new double[myModel.modifierLevels[i]];
  				mdt[i] = new double[myModel.modifierLevels[i]];
  				numModifierData[i] = new int[myModel.modifierLevels[i]][myModel.allocExpressionPrograms];
  				numModifierTables[i] = new int[myModel.modifierLevels[i]][myModel.allocExpressionPrograms];
  				thetaMod[i] = new double[myModel.modifierLevels[i]][myModel.allocExpressionPrograms];
  				alpha_thetaMod[i] = 1.0;
  			}
  		}
  	}
  	
  	public void resampleTheta(int numClusters,double[] theta2,double[][][] thetaMod2) {
  		int cc = 0;
  		int i = 0;
  		int j = 0;
  		for (cc=0;cc<numClusters+1;cc++) {
  			
  			if (myModel.useUpDown) {
  				theta2[cc] = sampleSingleTheta(cc);
  			}
  		
  			if (myModel.modifierLevels != null) {
  				sampleSingleThetaMod(cc);
  				for (i=0;i<mdt.length;i++) {
  					for (j=0;j<mdt[i].length;j++) {
  						thetaMod2[i][j][cc] = mdt[i][j];
  					}
  				}
  			}
  		}
  	}
  	
  	public void resampleTheta(int numClusters) {
  		resampleTheta(numClusters,theta,thetaMod);
  	}
  	
  	public double sampleSingleTheta(int cc) {
  		
  		if (myModel.useUpDown) {
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
  		
  		return -1.0;
  	}
  	
  	public void sampleSingleThetaMod(int cc) {
  		if (myModel.modifierLevels != null) {
  			int i = 0;
  			int j = 0;
  			for (i=0;i<md.length;i++) {
  				for (j=0;j<myModel.modifierLevels[i];j++) {
  					md[i][j] = (double) numModifierData[i][j][cc];
  					if (parent == null) {
  						md[i][j] += myModel.modifierPrior[i][j];
  					} else {
  						md[i][j] += alpha_thetaMod[i] * parent.thetaMod[i][j][cc];
  					}
  				}
  				StatUtil.dirichlet_rnd(mdt[i],md[i],myModel.modifierLevels[i]);
  			}
  		}
  	}
  	
  	public void resampleParameters(int numClusters,double[] condLike) {
  		super.resampleParameters(numClusters,condLike);
  		resampleTheta(numClusters);
  	}
  	
  	public void sampleNewTheta(int numClusters) {
  		
  		if (myModel.useUpDown) {
  			theta[numClusters] = sampleSingleTheta(numClusters);
  		}
  		
  		if (myModel.modifierLevels != null) {
  			int i = 0;
			int j = 0;
			sampleSingleThetaMod(numClusters);
			for (i=0;i<mdt.length;i++) {
				for (j=0;j<mdt[i].length;j++) {
					thetaMod[i][j][numClusters] = mdt[i][j];
				}
			}
  		}
  	}
  	
  	public void sampleNewParameters(int numClusters) {
  		super.sampleNewParameters(numClusters);
  		sampleNewTheta(numClusters);
  	}
  	
  	public void sampleClusterNumTables(int numClusters) {
  		super.sampleClusterNumTables(numClusters);
  		sampleModifierNumTables(numClusters);
  	}
  	
  	public void sampleModifierNumTables(int numClusters) {
  		int cc = 0;
  		int i = 0;
  		int j = 0;
  		
  		if (parent == null) {
  			for (cc=0;cc<numClusters;cc++) {
  				
  				if (myModel.useUpDown) {
	  				if (numUpData[cc] > 0) {
	  					numUpTables[cc] = 1;
	  				}
	  				if (numDownData[cc] > 0) {
	  					numDownTables[cc] = 1;
	  				}
  				}
  			
	  			if (myModel.modifierLevels != null) {
	  				for (i=0;i<myModel.modifierLevels.length;i++) {
	  					for (j=0;j<myModel.modifierLevels[i];j++) {
	  						for (cc=0;cc<numClusters;cc++) {
	  							if (numModifierData[i][j][cc] > 0) {
	  								numModifierTables[i][j][cc] = 1;
	  							}
	  						}
	  					}
	  				}
	  			}
  			}
  		}
  		else {
  			for (cc=0;cc<numClusters;cc++) {
  				
  				if (myModel.useUpDown) {
	  				parent.numUpData[cc] = parent.numUpData[cc] - numUpTables[cc];
	  				parent.numDownData[cc] = parent.numDownData[cc] - numDownTables[cc];
	  				
	  				numUpTables[cc] = numTables_rnd(alpha_theta*parent.theta[cc],numUpData[cc]);
	  				parent.numUpData[cc] += numUpTables[cc];  
	  				
	  				numDownTables[cc] = numTables_rnd(alpha_theta*(1.0-parent.theta[cc]),numDownData[cc]);
	  				parent.numDownData[cc] += numDownTables[cc];  
  				}
  			
	  			if (myModel.modifierLevels != null) {
	  				for (i=0;i<myModel.modifierLevels.length;i++) {
	  					for (j=0;j<myModel.modifierLevels[i];j++) {
	  						parent.numModifierData[i][j][cc] = parent.numModifierData[i][j][cc] - numModifierTables[i][j][cc];
	  						
	  						numModifierTables[i][j][cc] = numTables_rnd(alpha_thetaMod[i]*parent.thetaMod[i][j][cc],numModifierData[i][j][cc]);
	  						parent.numModifierData[i][j][cc] += numModifierTables[i][j][cc];
	  					}
	  				}
	  			}
  			}
  		}
  	}
}
