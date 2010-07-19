package edu.mit.csail.psrg.georg.DiscretizeData;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;

import edu.mit.csail.psrg.georg.StatUtil.InfoTheoryUtil;

public class DiscretizeExpression {
	public double[][] expression = null;
	public int[][] dExpression = null;
	public double[][] levels = null;
	double[] bestMI = null;
	int[] bestMerge = null; 
	int init_numLevels = 10;
	int numLevels = 0;
	int numGenes = 0;
	int numExperiments = 0;
	
	public DiscretizeExpression(double[][] E,int nl) {
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
		bestMerge = new int[numLevels];
		int[] pX = new int[numLevels+1];
		int[] pY = new int[numLevels+1];
		int[][] pXY = new int[numLevels+1][numLevels+1];
		int i = 0;
		double MI = 0.0;
		
	//	initialDiscretization();
		initialDiscretizationRank();
		bestMerge[numLevels-1] = 0;
		bestMI[numLevels-1] = computeMI(pX,pY,pXY);
		int[][] D2 = new int[numExperiments][numGenes];
		int[][] bestD = new int[numExperiments][numGenes];
		copyArray(dExpression,bestD);
		
		while(numLevels > 2) {
			numLevels = numLevels-1;
			bestMI[numLevels-1] = 0.0;
			bestMerge[numLevels-1] = 0;
			pX = new int[numLevels+1];
			pY = new int[numLevels+1];
			pXY = new int[numLevels+1][numLevels+1];
			for (i=0;i<numLevels;i++) {
				copyArray(bestD,dExpression);
				mergeDown(i+1);
				MI = computeMI(pX,pY,pXY);
				System.out.println((new Integer(numLevels)).toString() + " " + (new Integer(i)).toString() + " " + (new Double(MI)).toString());
				if (MI >= bestMI[numLevels-1]) {
					bestMI[numLevels-1] = MI;
					bestMerge[numLevels-1] = i+1;
					copyArray(dExpression,D2);
				}
			}
			copyArray(D2,bestD);
		}
		
		FileWriter outFile = new FileWriter(fName);
		
		for (i=0;i<init_numLevels;i++) {
			outFile.write((new Double(bestMI[i])).toString());
			outFile.write("\t");
			outFile.write((new Double(bestMerge[i])).toString());
			outFile.write("\n");
		}
		
		outFile.close();
	}
	
	void mergeDownLevels(int stopMerge,String fName) throws IOException {
		int i = 0;
		int j = 0;
		int k = 0;
		numLevels = init_numLevels;
	//	initialDiscretization();
		initialDiscretizationRank();
		for (j=0;j<numExperiments;j++) {
			for(i=init_numLevels-2;i>=stopMerge-1;i--) {
				k = bestMerge[i]-1;
				System.arraycopy(levels[j],k+1,levels[j],k,levels[j].length-k-1);
			}
		}
		FileWriter outFile = new FileWriter(fName);
		for(j=0;j<numExperiments;j++) {
			for (i=0;i<stopMerge;i++) {
				outFile.write((new Double(levels[j][i])).toString());
				outFile.write("\t");
			}
			outFile.write("\n");
		}
		outFile.close();
		
		for(i=init_numLevels-2;i>=stopMerge-1;i--) {
			k = bestMerge[i];
			mergeDown(k);
		}
	}
	
	// compute pairwise mutual information between all pairs of experiments
	double computeMI(int[] pX,int[] pY,int[][] pXY) {
		double m = 0.0;
		int i = 0;
		int j = 0;
		
		for (i=0;i<numExperiments-1;i++) {
			for (j=i+1;j<numExperiments;j++) {
				m = m + InfoTheoryUtil.discreteMutualInformation(dExpression[i],dExpression[j],pX,pY,pXY);
			}
		}
		
		return m;
	}
	
	void initialDiscretization() {
		// dExpression is transposed to make MI computations easier
		dExpression = new int[numExperiments][numGenes];
		levels = new double[numExperiments][numLevels];
		double maxE = 0.0;
		double minE = 0.0;
		double e = 0.0;
		double inc = 0.0;
		int i = 0;
		int j = 0;
		int k = 0;
		
	//	minE = 0.0+10e-8;
	//	maxE = 1.0+10e-8;
		
		for (j=0;j<numExperiments;j++) {
			maxE = maxColumn(j) + 10e-8;
			minE = minColumn(j) - 10e-8;
			inc = (maxE-minE)/((double) numLevels);
			for (i=0;i<numLevels;i++) {
				levels[j][i] = ((double) i+1)*inc + minE;
			}
		}
		
		for (j=0;j<numExperiments;j++) {
			for (k=0;k<numGenes;k++) {
				e = expression[k][j];
				if (!Double.isNaN(e)) {
					i = 0;
					while(e > levels[j][i]) {
						i++;
					}
					dExpression[j][k] = i+1;
				}
			}
		}
	}
	
	void initialDiscretizationRank() {
		dExpression = new int[numExperiments][numGenes];
		levels = new double[numExperiments][numLevels];
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
					while(e > levels[j][i]) {
						i++;
					}
					dExpression[j][k] = i+1;
				}
			}
		}
	}
	
	void setRankColumnLevels(int c,int nl) {
		int i = 0;
		int numNotNaN = 0;
		for (i=0;i<numGenes;i++) {
			if (!Double.isNaN(expression[i][c])) {
				numNotNaN++;
			}
		}
		
		Double[] s = new Double[numNotNaN];
		int j = 0;
		for (i=0;i<numGenes;i++) {
			if (!Double.isNaN(expression[i][c])) {
				s[j] = expression[i][c];
				j++;
			}
		}
		Arrays.sort(s);
		
		double inc = numNotNaN/nl;
	/*	if (inc < 1) {
			inc = 1;
		} */
		
		double jd = 0.0;
		j = 0;
		jd = 0.0;
		for (i=0;i<nl;i++) {
			jd += inc;
			j = (int) Math.round(jd);
			if (j > numNotNaN-1) {
				j = numNotNaN-1;
			} else {
				levels[c][i] = s[j];
			}
		}
		if (numNotNaN >= 1) {
			levels[c][nl-1] = s[numNotNaN-1]+10e-8;
		}
	}
	
	double maxColumn(int c) {
		double m = -1.0/0.0;
		int i = 0;
		for (i=0;i<numGenes;i++) {
			if (expression[i][c] > m) {
				m = expression[i][c];
			}
		}
		return m;
	}
	
	double minColumn(int c) {
		double m = 1.0/0.0;
		int i = 0;
		for (i=0;i<numGenes;i++) {
			if (expression[i][c] < m) {
				m = expression[i][c];
			}
		}
		return m;
	}
	
	void mergeDown(int l) {
		int i = 0;
		int j = 0;
		
		for (i=0;i<numExperiments;i++) {
			for (j=0;j<numGenes;j++) {
				if (dExpression[i][j] > l) {
					dExpression[i][j]--;
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
