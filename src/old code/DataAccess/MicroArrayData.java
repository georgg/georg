package edu.mit.csail.psrg.georg.DataAccess;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.StringTokenizer;

public class MicroArrayData {
	public double[][] values = null;
	public int[][] dvalues = null;
	public boolean discreteData = false;
	public int numRows;
	public int numCols;
	public int numNaN;
	public String[] experimentNames = null;
	public String[] geneNames = null;
	public int[] experimentCode = null;
	public int[] experimentCode2 = null;
	
	public void setDiscrete() {
		discreteData = true;
	}
	
	public void setContinuous() {
		discreteData = false;
	}
	
	public void readFileBlind(String fileName) throws IOException {
		if (discreteData) {
			readFileBlindDiscrete(fileName);
		} else {
			readFileBlindContinuous(fileName);
		}
	}
	
	public void readFileBlindContinuous(String fileName) throws IOException {
		BufferedReader is = new BufferedReader(new FileReader(fileName));
		
		String line = is.readLine();
		String[] splitLine = null;
		int i = 0;
		int j = 0;
		
		splitLine = line.split("\t");
		numCols = splitLine.length-1;
		experimentNames = new String[numCols];
		for (i=1;i<splitLine.length;i++) {
			experimentNames[i-1] = splitLine[i];
		}
		
		ArrayList<Double[]> valueArray = new ArrayList<Double[]>();
		ArrayList<String> names = new ArrayList<String>();
		line = is.readLine();
		Double[] valueLine = null;
		
		while(line != null) {
			splitLine = line.split("\t");
			names.add(splitLine[0]);
			valueLine = new Double[splitLine.length-1];
			for (i=1;i<splitLine.length;i++) {
				if (splitLine[i].equals("null")) {
					valueLine[i-1] = 0.0/0.0;
				} else {
					valueLine[i-1] = new Double(splitLine[i]);
				}
			}
			valueArray.add(valueLine);
			line = is.readLine();
		}
		numRows = names.size();
		geneNames = new String[numRows];
		for (i=0;i<names.size();i++) {
			geneNames[i] = names.get(i);
		}
		
		values = new double[numRows][numCols];
		
		for (i=0;i<valueArray.size();i++) {
			valueLine = valueArray.get(i);
			for (j=0;j<valueLine.length;j++) {
				values[i][j] = valueLine[j];
			}
		}
	}
	
	public void readFileBlindDiscrete(String fileName) throws IOException {
		BufferedReader is = new BufferedReader(new FileReader(fileName));
		
		String line = is.readLine();
		String[] splitLine = null;
		int i = 0;
		int j = 0;
		
		splitLine = line.split("\t");
		numCols = splitLine.length-1;
		experimentNames = new String[numCols];
		for (i=1;i<splitLine.length;i++) {
			experimentNames[i-1] = splitLine[i];
		}
		
		ArrayList<Integer[]> valueArray = new ArrayList<Integer[]>();
		ArrayList<String> names = new ArrayList<String>();
		line = is.readLine();
		Integer[] valueLine = null;
		
		while(line != null) {
			splitLine = line.split("\t");
			names.add(splitLine[0]);
			valueLine = new Integer[splitLine.length-1];
			for (i=1;i<splitLine.length;i++) {
				if (splitLine[i].equals("null")) {
					valueLine[i-1] = 0;
				} else {
					valueLine[i-1] = new Integer(splitLine[i]);
				}
			}
			valueArray.add(valueLine);
			line = is.readLine();
		}
		numRows = names.size();
		geneNames = new String[numRows];
		for (i=0;i<names.size();i++) {
			geneNames[i] = names.get(i);
		}
		
		dvalues = new int[numRows][numCols];
		
		for (i=0;i<valueArray.size();i++) {
			valueLine = valueArray.get(i);
			for (j=0;j<valueLine.length;j++) {
				dvalues[i][j] = valueLine[j];
			}
		}
		
		setDiscrete();
	}
	
	public void readFile(String fileName) throws IOException {
		BufferedReader is = new BufferedReader(new FileReader(fileName));
		String line = "";
		
		line = is.readLine();
		StringTokenizer st = new StringTokenizer(line, "\t");
		
		int i = 0;
		int j = 0;
		int hasGeneNames = 0;
		int hasExperimentNames = 0;
		int hasExperimentCode = 0;
		int hasExperimentCode2 = 0;
		
		numRows = new Integer(st.nextToken());
		numCols = new Integer(st.nextToken());
		hasGeneNames = new Integer(st.nextToken());
		hasExperimentNames = new Integer(st.nextToken());
		hasExperimentCode = new Integer(st.nextToken());
		hasExperimentCode2 = new Integer(st.nextToken());
		
		if (!discreteData) {
			values = new double[numRows][numCols];
		} else {
			dvalues = new int[numRows][numCols];
		}
		
		if (hasExperimentNames == 1) {
			experimentNames = new String[numCols];
			line = is.readLine();
			st = new StringTokenizer(line, "\t");
			// skip first column
			if (hasGeneNames == 1)
				st.nextToken();
			
			for (i=0;i<numCols;i++) {
				experimentNames[i] = st.nextToken();
			}
		}
		
		if (hasExperimentCode == 1) {
			experimentCode = new int[numCols];
			line = is.readLine();
			st = new StringTokenizer(line, "\t");
			// skip first column
			if (hasGeneNames == 1)
				st.nextToken();
			
			for (i=0;i<numCols;i++) {
				experimentCode[i] = new Integer(st.nextToken());
			}
		}
		
		if (hasExperimentCode2 == 1) {
			experimentCode2 = new int[numCols];
			line = is.readLine();
			st = new StringTokenizer(line, "\t");
			// skip first column
			if (hasGeneNames == 1)
				st.nextToken();
			
			for (i=0;i<numCols;i++) {
				experimentCode2[i] = new Integer(st.nextToken());
			}
		}
		
		if (hasGeneNames == 1) {
			geneNames = new String[numRows];
		}
		
		String tt = null;
		
		for (i=0;i<numRows;i++) {
			line = is.readLine();
			st = new StringTokenizer(line, "\t");
			if (hasGeneNames == 1) {
				geneNames[i] = st.nextToken();
			}
			for (j=0;j<numCols;j++) {
				tt = st.nextToken();
				if (!discreteData) {
					if (tt.equals("null")) {
						values[i][j] = 0.0/0.0;
					} else {
						values[i][j] = new Double(tt);
					}
				} else {
					dvalues[i][j] = new Integer(tt);
				}
			}
		}
		is.close();
	}
	
	public void writeFile(String fname) throws IOException {
		int i = 0;
		int j = 0;
		FileWriter outFile = new FileWriter(fname);
		
		String t = "";
		t = (new Integer(numRows)).toString();
		t = t + "\t";
		t = t + (new Integer(numCols)).toString() + "\t";
		
		if (geneNames != null) {
			t = t + "1";
		} else {
			t = t + "0";
		}
		t = t + "\t";
		
		if (experimentNames != null) {
			t = t + "1";
		} else {
			t = t + "0";
		}
		t = t + "\t";
		
		if (experimentCode != null) {
			t = t + "1";
		} else {
			t = t + "0";
		}
		t = t + "\t";
		
		if (experimentCode2 != null) {
			t = t + "1";
		} else {
			t = t + "0";
		}
		
		t = t + "\n";
		
		outFile.write(t);
		
		if (experimentNames != null) {
			outFile.write("genes");
			for(i=0;i<numCols;i++) {
				outFile.write("\t");
				outFile.write(experimentNames[i]);
			}
			outFile.write("\n");
		}
		
		if (experimentCode != null) {
			outFile.write("genes");
			for(i=0;i<numCols;i++) {
				outFile.write("\t");
				t = (new Integer(experimentCode[i])).toString();
				outFile.write(t);
			}
			outFile.write("\n");
		}
		
		if (experimentCode2 != null) {
			outFile.write("genes");
			for(i=0;i<numCols;i++) {
				outFile.write("\t");
				t = (new Integer(experimentCode2[i])).toString();
				outFile.write(t);
			}
			outFile.write("\n");
		}
		
		for(i=0;i<numRows;i++) {
			if (geneNames != null) {
				outFile.write(geneNames[i]);
				outFile.write("\t");
			}
			for(j=0;j<numCols;j++) {
				if (!discreteData) {
					t = (new Double(values[i][j])).toString();
				//	t = String.format("%.4f",values[i][j]);
				} else {
					t = (new Integer(dvalues[i][j])).toString();
				}
				outFile.write(t);
				if (j<numCols-1) {
					outFile.write("\t");
				} else {
					outFile.write("\n");
				}
			}
		}
		
		outFile.close();
	}
	
	public void normCols() {
		double[] mu = new double[numCols];
		double[] s = new double[numCols];
		int[] n = new int[numCols];
		double t = 0.0;
		int i = 0;
		int j = 0;
		
		if (numRows == 0 && numCols == 0)
			return;
			
		for (j=0; j < numCols; j++) {
			mu[j] = 0;
			n[j] = 0;
			for (i = 0; i < numRows; i++) {
				if (!Double.isNaN(values[i][j])) {
					mu[j] = mu[j] + values[i][j];
					n[j]++;
				}
			}
			mu[j] = mu[j]/((double) n[j]);
		}
		
		for (j=0; j < numCols; j++) {
			s[j] = 0;
			for (i = 0; i < numRows; i++) {
				if (!Double.isNaN(values[i][j])) {
					t = (values[i][j] - mu[j]);
					s[j] = s[j] + t*t;
				}
			}
			s[j] = Math.sqrt(s[j]/((double) n[j]));
		}
		
		for (i=0;i<numRows;i++) {
			for (j=0;j<numCols;j++) {
				values[i][j] = (values[i][j]-mu[j])/s[j];
			}
		}
	}
	
	public void normRows() {
		double[] mu = new double[numRows];
		double[] s = new double[numRows];
		int[] n = new int[numRows];
		double t = 0.0;
		int i = 0;
		int j = 0;
		
		if (numRows == 0 && numCols == 0)
			return;
			
		for (i=0; i < numRows; i++) {
			mu[i] = 0;
			n[i] = 0;
			for (j = 0; j < numCols; j++) {
				if (!Double.isNaN(values[i][j])) {
					mu[i] = mu[i] + values[i][j];
					n[i]++;
				}
			}
			mu[i] = mu[i]/((double) n[i]);
		}
		
		for (i=0; i < numRows; i++) {
			s[i] = 0;
			for (j = 0; j < numCols; j++) {
				if (!Double.isNaN(values[i][j])) {
					t = (values[i][j] - mu[i]);
					s[i] = s[i] + t*t;
				}
			}
			s[i] = Math.sqrt(s[i]/((double) n[i]));
		}
		
		for (i=0;i<numRows;i++) {
			for (j=0;j<numCols;j++) {
				values[i][j] = (values[i][j]-mu[i])/s[i];
			}
		}
	}
	
	public double[] meanCols() {
		double[] mu = new double[numCols];
		int i = 0;
		int j = 0;
		double muj = 0.0;
		double nj = 0.0;
		
		for (j=0;j<numCols;j++) {
			muj = 0.0;
			nj = 0.0;
			for (i=0;i<numRows;i++) {
				if (!Double.isNaN(values[i][j])) {
					muj = muj + values[i][j];
					nj = nj + 1.0;
				}
			}
			mu[j] = muj/nj;
		}
		return mu;
	}
	
	public double totalSTD(double[] mu) {
		if (mu.length != numCols) {
			return -1.0;
		}
		
		double tv = 0.0;
		double s = 0.0;
		int i = 0;
		int j = 0;
		
		for (i=0;i<numRows;i++) {
			for (j=0;j<numCols;j++) {
				if (!Double.isNaN(values[i][j])) {
					s = s + Math.pow(values[i][j]-mu[j],2.0);
					tv = tv + 1.0;
				}
			}
		}
		
		s = Math.sqrt(s/(tv-1.0));
		return s;
	}
	
	public static MicroArrayData mergeCols(MicroArrayData d1,MicroArrayData d2) {
		MicroArrayData d3 = new MicroArrayData();
		int i = 0;
		
		d3.numRows = d1.numRows;
		d3.numCols = d1.numCols + d2.numCols;
		
		if (d1.discreteData) {
			d3.dvalues = new int[d1.numRows][d1.numCols+d2.numCols];
			d3.discreteData = true;
			for (i=0;i<d1.numRows;i++) {
				System.arraycopy(d1.dvalues[i],0,d3.dvalues[i],0,d1.numCols);
				System.arraycopy(d2.dvalues[i],0,d3.dvalues[i],d1.numCols,d2.numCols);
			}
		} else {
			d3.values = new double[d1.numRows][d1.numCols+d2.numCols];
			for (i=0;i<d1.numRows;i++) {
				System.arraycopy(d1.values[i],0,d3.values[i],0,d1.numCols);
				System.arraycopy(d2.values[i],0,d3.values[i],d1.numCols,d2.numCols);
			}
		}
		
		if (d1.experimentNames != null) {
			d3.experimentNames = new String[d1.numCols+d2.numCols];
			System.arraycopy(d1.experimentNames,0,d3.experimentNames,0,d1.numCols);
			System.arraycopy(d2.experimentNames,0,d3.experimentNames,d1.numCols,d2.numCols);
		}
		
		if (d1.experimentCode != null) {
			d3.experimentCode = new int[d1.numCols+d2.numCols];
			System.arraycopy(d1.experimentCode,0,d3.experimentCode,0,d1.numCols);
			System.arraycopy(d2.experimentCode,0,d3.experimentCode,d1.numCols,d2.numCols);
		}
		
		if (d1.experimentCode2 != null) {
			d3.experimentCode2 = new int[d1.numCols+d2.numCols];
			System.arraycopy(d1.experimentCode2,0,d3.experimentCode2,0,d1.numCols);
			System.arraycopy(d2.experimentCode2,0,d3.experimentCode2,d1.numCols,d2.numCols);
		}
		
		if (d1.geneNames != null) {
			d3.geneNames = new String[d1.geneNames.length];
			System.arraycopy(d1.geneNames,0,d3.geneNames,0,d3.numRows);
		}
		
		return d3;
	}
	
	public void deleteCols(HashSet<Integer> dcols) {
		
		if (dcols.isEmpty()) {
			return;
		}
		
		int newNumCols = numCols - dcols.size();
		double[][] newValues = null;
		int[][] newValuesD = null;
		String[] newExperimentNames = null;
		int[] newExperimentCode = null;
		int[] newExperimentCode2 = null;
		
		if (discreteData) {
			newValuesD = new int[numRows][newNumCols];
		} else {
			newValues = new double[numRows][newNumCols];
		}
		if (experimentNames != null) {
			newExperimentNames = new String[newNumCols];
		}
		if (experimentCode != null) {
			newExperimentCode = new int[newNumCols];
		}
		if (experimentCode2 != null) {
			newExperimentCode2 = new int[newNumCols];
		}
		
		int i = 0;
		int j = 0;
		int k = 0;
		
		for (j=0;j<numCols;j++) {
			if (!dcols.contains(j)) {
				if (experimentNames != null) {
					newExperimentNames[k] = experimentNames[j];
				}
				if (experimentCode != null) {
					newExperimentCode[k] = experimentCode[j];
				}
				if (experimentCode2 != null) {
					newExperimentCode2[k] = experimentCode2[j];
				}
				for (i=0;i<numRows;i++) {
					if (discreteData) {
						newValuesD[i][k] = dvalues[i][j];
					} else {
						newValues[i][k] = values[i][j];
					}
				}
				k++;
			}
		}
		
		if (discreteData) {
			dvalues = newValuesD;
		} else {
			values = newValues;
		}
		experimentNames = newExperimentNames;
		experimentCode = newExperimentCode;
		experimentCode2 = newExperimentCode2;
		numCols = newNumCols;
	}
	
	public void deleteRows(HashSet<Integer> drows) {
		
		if (drows.isEmpty()) {
			return;
		}
		
		int newNumRows = numRows - drows.size();
		double[][] newValues = null;
		int[][] newValuesD = null;
		String[] newGeneNames = null;
		
		if (discreteData) {
			newValuesD = new int[newNumRows][numCols];
		} else {
			newValues = new double[newNumRows][numCols];
		}
		if (geneNames != null) {
			newGeneNames = new String[newNumRows];
		}
		
		int i = 0;
		int j = 0;
		int k = 0;
		
		for (i=0;i<numRows;i++) {
			if (!drows.contains(i)) {
				if (geneNames != null) {
					newGeneNames[k] = geneNames[i];
				}
				for (j=0;j<numCols;j++) {
					if (discreteData) {
						newValuesD[k][j] = dvalues[i][j];
					} else {
						newValues[k][j] = values[i][j];
					}
				}
				k++;
			}
		}
		
		if (discreteData) {
			dvalues = newValuesD;
		} else {
			values = newValues;
		}
		geneNames = newGeneNames;
		numRows = newNumRows;
	}
	
	public void log2Data() {
		int i = 0;
		int j = 0;
		
		for (i=0;i<numRows;i++) {
			for (j=0;j<numCols;j++) {
				values[i][j] = Math.log(values[i][j])/Math.log(2);
			}
		}
	}
	
	public double[][] copyContinousData() {
		double[][] values2 = new double[numRows][numCols];
		int i = 0;
		
		for (i=0;i<numRows;i++) {
			System.arraycopy(values[i],0,values2[i],0,numCols);
		}
		
		return values2;
	}
}

