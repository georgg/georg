package edu.mit.csail.psrg.georg.DiscretizeData;

import java.util.ArrayList;

import edu.mit.csail.psrg.georg.StatUtil.StatUtil;

public class DiscretizeAffyData extends DiscretizeExpression {
	public int[][] PMA = null;
	
	public DiscretizeAffyData(double[][] E, int[][] P,int nl) {
		super(E, nl);
		PMA = P;
	}

	// this will take the log of the intensities
	// and subtract the median log of absent/marginal values
	// for each gene
	public void adjustIntensities() {
		double baseline = 0.0;
		int i = 0;
		int j = 0;
		ArrayList<Double> absentValues = new ArrayList<Double>();
		double[] v = null;
		
		for (i=0;i<numGenes;i++) {
			absentValues.clear();
			for (j=0;j<numExperiments;j++) {
				if (expression[i][j] <= 0.0) {
					expression[i][j] = 0.0;
				} else {
					expression[i][j] = Math.log(expression[i][j]);
				}
				if (PMA[i][j] == 0) {
					absentValues.add(expression[i][j]);
				}
			}
			baseline = 0.0;
			if (!absentValues.isEmpty()) {
				v = new double[absentValues.size()];
				for (j=0;j<v.length;j++) {
					v[j] = absentValues.get(j);
				}
				baseline = StatUtil.median(v);
			}
			for (j=0;j<numExperiments;j++) {
				if (PMA[i][j] > 0) {
					expression[i][j] = expression[i][j]-baseline;
					// remove expression that is below the baseline
					if (expression[i][j] <= 0.0) {
						expression[i][j] = 0.0/0.0;
					}
				} else {
					expression[i][j] = 0.0/0.0;
				}
			}
			
		}
	}
}
