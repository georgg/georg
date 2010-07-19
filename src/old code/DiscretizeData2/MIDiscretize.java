package edu.mit.csail.psrg.georg.DiscretizeData2;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;

import edu.mit.csail.psrg.georg.StatUtil.InfoTheoryUtil;

public class MIDiscretize {
	public double[][] expression = null;
	public int[][] dExpression = null;
	public int[][] pdExpression = null;
	public double[][] levelsUP = null;
	public double[][] levelsDN = null;
	double[] bestMI = null;
	int[][] bestMerge = null;  
	int init_numLevels = 10;
	int numLevels = 0;
	int numGenes = 0;
	int numExperiments = 0;
	
	public MIDiscretize(double[][] E,int nl) {
		expression = E;
		init_numLevels = nl;
		numGenes = expression.length;
		numExperiments = expression[0].length;
	}
	
	// performs an agglomerative discretization, computing
	// maximal mutual information between experiments
	// fName is a file for output of the mutual information and
	// discretization merges
	public void discretizeDownTree(String fName) throws IOException {
		numLevels = init_numLevels;
		
		bestMI = new double[numLevels];
		bestMerge = new int[numLevels][2];
		int[] pX = new int[2*numLevels+1];
		int[] pY = new int[2*numLevels+1];
		int[][] pXY = new int[2*numLevels+1][2*numLevels+1];
		int i = 0;
		int j = 0;
		double MI = 0.0;
		
		initialDiscretizationRank();
		bestMerge[numLevels-1][0] = 0;
		bestMerge[numLevels-1][1] = 0;
		updatepdExpression();
		bestMI[numLevels-1] = computeMI(pX,pY,pXY);
		int[][] D2 = new int[numExperiments][numGenes];
		int[][] bestD = new int[numExperiments][numGenes];
		copyArray(dExpression,bestD);
		
		while(numLevels > 2) {
			numLevels = numLevels-1;
			bestMI[numLevels-1] = 0.0;
			bestMerge[numLevels-1][0] = 0;
			bestMerge[numLevels-1][1] = 0;
			pX = new int[2*numLevels+1];
			pY = new int[2*numLevels+1];
			pXY = new int[2*numLevels+1][2*numLevels+1];
			for (i=0;i<numLevels;i++) {
				for (j=0;j<numLevels;j++) {
					copyArray(bestD,dExpression);
					mergeDown(i+1,j+1);
					updatepdExpression();
					MI = computeMI(pX,pY,pXY);
					System.out.println((new Integer(numLevels)).toString() + " " + (new Integer(i)).toString() + " " + (new Integer(j)).toString() + " " + (new Double(MI)).toString());
					if (MI >= bestMI[numLevels-1]) {
						bestMI[numLevels-1] = MI;
						bestMerge[numLevels-1][0] = i+1;
						bestMerge[numLevels-1][1] = j+1;
						copyArray(dExpression,D2);
					}
				}
			}
			copyArray(D2,bestD);
		}
		
		FileWriter outFile = new FileWriter(fName);
		
		for (i=0;i<init_numLevels;i++) {
			outFile.write((new Double(bestMI[i])).toString());
			outFile.write("\t");
			outFile.write((new Integer(bestMerge[i][0])).toString());
			outFile.write("\t");
			outFile.write((new Integer(bestMerge[i][1])).toString());
			outFile.write("\n");
		}
		
		outFile.close();
	}
	
	void mergeDownLevels(int stopMerge,String fName) throws IOException {
		int i = 0;
		int j = 0;
		int kUP = 0;
		int kDN = 0;
		numLevels = init_numLevels;

		initialDiscretizationRank();
		for (j=0;j<numExperiments;j++) {
			for(i=init_numLevels-2;i>=stopMerge-1;i--) {
				kUP = bestMerge[i][0]-1;
				kDN = bestMerge[i][0]-1;
				System.arraycopy(levelsUP[j],kUP+1,levelsUP[j],kUP,levelsUP[j].length-kUP-1);
				System.arraycopy(levelsDN[j],kDN+1,levelsDN[j],kDN,levelsDN[j].length-kDN-1);
			}
		}
		FileWriter outFileUP = new FileWriter(fName+"_UP.txt");
		FileWriter outFileDN = new FileWriter(fName+"_DN.txt");
		for(j=0;j<numExperiments;j++) {
			for (i=0;i<stopMerge;i++) {
				outFileUP.write((new Double(levelsUP[j][i])).toString());
				outFileUP.write("\t");
				outFileDN.write((new Double(-levelsDN[j][i])).toString());
				outFileDN.write("\t");
			}
			outFileUP.write("\n");
			outFileDN.write("\n");
		}
		outFileUP.close();
		outFileDN.close();
		
		for(i=init_numLevels-2;i>=stopMerge-1;i--) {
			kUP = bestMerge[i][0];
			kDN = bestMerge[i][1];
			mergeDown(kUP,kDN);
		}
	}
	
	// shifts the discretized expression matrix so values are >= 0
	// for efficient MI calculation
	void updatepdExpression() {
		int i = 0;
		int j = 0;
		for (i=0;i<dExpression.length;i++) {
			for (j=0;j<dExpression[0].length;j++) {
				pdExpression[i][j] = dExpression[i][j]+numLevels;
			}
		}
	}
	
	// compute pairwise mutual information between all pairs of experiments
	double computeMI(int[] pX,int[] pY,int[][] pXY) {
		double m = 0.0;
		int i = 0;
		int j = 0;
		
		for (i=0;i<numExperiments-1;i++) {
			for (j=i+1;j<numExperiments;j++) {
				m = m + InfoTheoryUtil.discreteMutualInformation(pdExpression[i],pdExpression[j],pX,pY,pXY);
			}
		}
		
		return m;
	}
	
	void initialDiscretizationRank() {
		dExpression = new int[numExperiments][numGenes];
		pdExpression = new int[numExperiments][numGenes];
		levelsUP = new double[numExperiments][numLevels];
		levelsDN = new double[numExperiments][numLevels];
		double e = 0.0;
		int i = 0;
		int j = 0;
		int k = 0;
		
		for (j=0;j<numExperiments;j++) {
			setRankColumnLevels(j,numLevels);
		}
		
		for (j=0;j<numExperiments;j++) {
			for (k=0;k<numGenes;k++) {
				e = expression[k][j];
				if (!Double.isNaN(e)) {
					i = 0;
					if (e >= 0.0) {
						while(e > levelsUP[j][i]) {
							i++;
						}
						dExpression[j][k] = i+1;
					} else {
						e = -e;
						while(e > levelsDN[j][i]) {
							i++;
						}
						dExpression[j][k] = -(i+1);
					}
				}
			}
		}
	}
	
	void setRankColumnLevels(int c,int nl) {
		int i = 0;
		int numUP = 0;
		int numDN = 0;
		
		for (i=0;i<numGenes;i++) {
			if (expression[i][c] < 0 & !Double.isNaN(expression[i][c])) {
				numDN++;
			}
			if (expression[i][c] >= 0 & !Double.isNaN(expression[i][c])) {
				numUP++;
			}
		}
		
		Double[] sUP = new Double[numUP];
		Double[] sDN = new Double[numDN];
		int jUP = 0;
		int jDN = 0;
		for (i=0;i<numGenes;i++) {
			if (expression[i][c] < 0 & !Double.isNaN(expression[i][c])) {
				sDN[jDN] = -expression[i][c];
				jDN++;
			}
			if (expression[i][c] >= 0 & !Double.isNaN(expression[i][c])) {
				sUP[jUP] = expression[i][c];
				jUP++;
			}
		}
		Arrays.sort(sUP);
		Arrays.sort(sDN);
		
		double incUP = numUP/nl;
		double incDN = numDN/nl;
		
		double jd = 0.0;
		int j = 0;
		jd = 0.0;
		
		for (i=0;i<nl;i++) {
			jd += incUP;
			j = (int) Math.round(jd);
			if (j > numUP-1) {
				j = numUP-1;
			} else {
				levelsUP[c][i] = sUP[j];
			}
		}
		if (numUP >= 1) {
			levelsUP[c][nl-1] = sUP[numUP-1]+10e-8;
		}
		
		j = 0;
		jd = 0.0;
		for (i=0;i<nl;i++) {
			jd += incDN;
			j = (int) Math.round(jd);
			if (j > numDN-1) {
				j = numDN-1;
			} else {
				levelsDN[c][i] = sDN[j];
			}
		}
		if (numDN >= 1) {
			levelsDN[c][nl-1] = sDN[numDN-1]+10e-8;
		}
	}
	
	void mergeDown(int lUP,int lDN) {
		int i = 0;
		int j = 0;
		
		for (i=0;i<numExperiments;i++) {
			for (j=0;j<numGenes;j++) {
				if (dExpression[i][j] > 0 & dExpression[i][j] > lUP) {
					dExpression[i][j]--;
				}
				if (dExpression[i][j] < 0 & -dExpression[i][j] > lDN) {
					dExpression[i][j]++;
				}
			}
		}
	}
	
	void copyArray(int[][] src,int[][] dest) {
		int i = 0;
		for (i=0;i<src.length;i++) {
			System.arraycopy(src[i],0,dest[i],0,src[i].length);
		}
	}
	
	void transposeDExpression() {
		int[][] D2 = new int[dExpression[0].length][dExpression.length];
		int i = 0;
		int j = 0;
		for (i=0;i<dExpression.length;i++) {
			for (j=0;j<dExpression[0].length;j++) {
				D2[j][i] = dExpression[i][j];
			}
		}
		dExpression = D2;
		D2 = null;
	}
}
