package edu.mit.csail.psrg.georg.GeneProgram2;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashSet;

public class ExpressionProgram implements Serializable {
	/**
	 * 
	 */
	private static final long serialVersionUID = -2156340130982474428L;
	// how many positive counts (for each gene) are associated with each expression program
	public int[] posCounts = null;
	public int totalCount = 0;
	// temporary map used for evaluating tissue loading on expression programs
	public HashSet<Integer> expressionProgramMap = new HashSet<Integer>();
	// temporary array for evaluating tissue loading on expression programs
	public ArrayList<Integer> expressionProgramArray = new ArrayList<Integer>();
	
	public ExpressionProgram(int numGenes) {
		posCounts = new int[numGenes];
	}
	
	public void zeroCounts() {
		int i = 0;
		for (i=0;i<posCounts.length;i++) {
			posCounts[i] = 0;
		}
		totalCount = 0;
	}
	
	public void removeGeneCount(int gene) {
		posCounts[gene]--;
		totalCount--;
	}
	
	public void addGeneCount(int gene) {
		posCounts[gene]++;
		totalCount++;
	}
}
