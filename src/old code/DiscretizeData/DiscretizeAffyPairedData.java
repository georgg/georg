package edu.mit.csail.psrg.georg.DiscretizeData;

import edu.mit.csail.psrg.georg.StatUtil.InfoTheoryUtil;

public class DiscretizeAffyPairedData extends DiscretizeAffyData {
	int[][] matchPairs = null;
	
	public DiscretizeAffyPairedData(double[][] E, int[][] P, int[][] mp,int nl) {
		super(E, P, nl);
		matchPairs = mp;
	}
	
/*	double computeMI(int[] pX,int[] pY,int[][] pXY) {
		double m = 0.0;
		int k = 0;
		int i = 0;
		int j = 0;
		
		for (k=0;k<matchPairs.length;k++) {
			i = matchPairs[k][0];
			j = matchPairs[k][1];
			m = m + InfoTheoryUtil.discreteMutualInformation(dExpression[i],dExpression[j],pX,pY,pXY);
		}
		
		return m;
	} */
}
