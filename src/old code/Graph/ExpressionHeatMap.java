package edu.mit.csail.psrg.georg.Graph;

import java.io.IOException;

public class ExpressionHeatMap {
	double[][] expression = null;
	String[] geneNames = null;
	String[] experimentNames = null;
	
	public ExpressionHeatMap(double[][] e,String[] gn,String[] en) {
		expression = e;
		geneNames = gn;
		experimentNames = en;
	}
	
	public void outputHeatMap(String fName) {
		int i = 0;
		int j = 0;
		TabularComponent sheet = new TabularComponent(expression.length+1,expression[0].length+1,30);
		int[] colWidths = new int[expression[0].length+2];
		for (i=0;i<colWidths.length;i++) {
			colWidths[i] = 60;
		}
		colWidths[0] = 20;
		sheet.setColWidths(colWidths);
		sheet.setGrid(true);
		
		double maxExpression = Math.abs(expression[0][0]);
		double minExpression = Math.abs(expression[0][0]);
		for (i=0;i<expression.length;i++) {
			for (j=0;j<expression[0].length;j++) {
				if (maxExpression < Math.abs(expression[i][j])) {
					maxExpression = Math.abs(expression[i][j]);
				}
				if (minExpression > Math.abs(expression[i][j])) {
					minExpression = Math.abs(expression[i][j]);
				}
			}
		}
		
		TabText myTab = null;
		float color = 0.0f;
		double minIntensity = 0.0;
		double maxIntensity = 1.0;
		double v = 0.0;
		double a = (maxIntensity-minIntensity)/(maxExpression-minExpression);
		double b = minIntensity - minExpression*a;
		
		for (j=0;j<experimentNames.length;j++) {
			myTab = new TabText(experimentNames[j]);
			sheet.setCell(0,j+1,myTab);
		}
		
		for (i=0;i<expression.length;i++) {
			myTab = new TabText(geneNames[i]);
			sheet.setCell(i+1,0,myTab);
			for (j=0;j<expression[0].length;j++) {
				myTab = new TabText("");
				v = Math.abs(expression[i][j])*a + b;
				color = (float) v;
				if (color > 1.0f) {
					color = 1.0f;
				}
				if (color < 0.0f) {
					color = 0.0f;
				}
				myTab = new TabText("");
				if (expression[i][j] <= 0) {
					myTab.setFillColor(0.0f,color,0.0f);
				} else {
					myTab.setFillColor(color,0.0f,0.0f);
				}
				myTab.setFill(true);
				sheet.setCell(i+1,j+1,myTab);
			}
		}
		
		try {
			sheet.outputSVG(fName);
		} catch(IOException e) {
			System.out.println(e);
		}
	}
}
