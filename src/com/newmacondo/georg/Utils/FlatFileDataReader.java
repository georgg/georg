package com.newmacondo.georg.Utils;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;

public class FlatFileDataReader {
	public double[][] values = null;
	public int[][] dvalues = null;
	public boolean discreteData = false;
	public int numRows;
	public int numCols;
	public int numNaN;
	public String[] colNames = null;
	public String[] rowNames = null;
	
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
			readFileContinuous(fileName);
		}
	}
	
	public void readFileContinuous(String fileName) throws IOException {
		BufferedReader is = new BufferedReader(new FileReader(fileName));
		
		String line = is.readLine();
		String[] splitLine = null;
		int i = 0;
		int j = 0;
		
		splitLine = line.split("\t");
		numCols = splitLine.length-1;
		colNames = new String[numCols];
		for (i=1;i<splitLine.length;i++) {
			colNames[i-1] = splitLine[i];
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
		rowNames = new String[numRows];
		for (i=0;i<names.size();i++) {
			rowNames[i] = names.get(i);
		}
		
		values = new double[numRows][numCols];
		
		for (i=0;i<valueArray.size();i++) {
			valueLine = valueArray.get(i);
			for (j=0;j<valueLine.length;j++) {
				values[i][j] = valueLine[j];
			}
		}
	}
	
	public void readFileDiscrete(String fileName) throws IOException {
		BufferedReader is = new BufferedReader(new FileReader(fileName));
		
		String line = is.readLine();
		String[] splitLine = null;
		int i = 0;
		int j = 0;
		
		splitLine = line.split("\t");
		numCols = splitLine.length-1;
		colNames = new String[numCols];
		for (i=1;i<splitLine.length;i++) {
			colNames[i-1] = splitLine[i];
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
		rowNames = new String[numRows];
		for (i=0;i<names.size();i++) {
			rowNames[i] = names.get(i);
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
	
	public void deleteRows(HashSet<Integer> drows) {
		
		if (drows.isEmpty()) {
			return;
		}
		
		int newNumRows = numRows - drows.size();
		double[][] newValues = null;
		int[][] newValuesD = null;
		String[] newRowNames = null;
		
		if (discreteData) {
			newValuesD = new int[newNumRows][numCols];
		} else {
			newValues = new double[newNumRows][numCols];
		}
		if (rowNames != null) {
			newRowNames = new String[newNumRows];
		}
		
		int i = 0;
		int j = 0;
		int k = 0;
		
		for (i=0;i<numRows;i++) {
			if (!drows.contains(i)) {
				if (rowNames != null) {
					newRowNames[k] = rowNames[i];
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
		rowNames = newRowNames;
		numRows = newNumRows;
	}
}

