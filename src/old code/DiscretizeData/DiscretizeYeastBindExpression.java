package edu.mit.csail.psrg.georg.DiscretizeData;

import java.io.IOException;

import edu.mit.csail.psrg.georg.DataAccess.MicroArrayData;

public class DiscretizeYeastBindExpression {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		processData();
	}

	public static void processData() {
		String dirPath = "c:\\research_data\\combine_yeast_bind_expression\\";
		String yeastDataFile = "yeast_combined_data.txt";
		String discretizedMI_fName = "yeast_combined_discretized_MI.txt";
		String discretizedLevels_fName = "yeast_combined_discretized_levels.txt";
		String discretizedData_fName = "yeast_combined_discretized_data.txt";
		
		MicroArrayData yeastData = new MicroArrayData();
		
		try {
			yeastData.readFile(dirPath + yeastDataFile);
		} catch(IOException e) {
			System.out.println(e);
		}
		
		DiscretizeExpression discretizer = new DiscretizeExpression(yeastData.values,10);
		
		try {
			discretizer.discretizeDownTree(dirPath+discretizedMI_fName);
			discretizer.mergeDownLevels(3,dirPath+discretizedLevels_fName);
			discretizer.transposeDExpression();
			yeastData.setDiscrete();
			yeastData.values = null;
			yeastData.dvalues = discretizer.dExpression;
			yeastData.writeFile(dirPath+discretizedData_fName);
		} catch(IOException e) {
			System.out.println(e);
		}
	}
}
