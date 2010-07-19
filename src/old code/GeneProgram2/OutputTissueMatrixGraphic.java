package edu.mit.csail.psrg.georg.GeneProgram2;

import java.awt.Color;
import java.awt.Container;
import java.awt.EventQueue;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.Vector;

import javax.swing.JFrame;

import edu.mit.csail.psrg.georg.DataAccess.ClusterReader;
import edu.mit.csail.psrg.georg.DataAccess.MicroArrayData;
import edu.mit.csail.psrg.georg.Graph3.Table;
import edu.mit.csail.psrg.georg.Graph3.TableCell;
import edu.mit.csail.psrg.georg.StatUtil.VectorUtil;

public class OutputTissueMatrixGraphic extends JFrame {
	public OutputTissueMatrixGraphic() {
		String baseFileName = "c:\\research_data\\Jenner\\GeneProgram\\output\\infection_first_2l_group\\infection";
		String modNamesFileName = "c:\\research_data\\Jenner\\GeneProgram\\first_2l\\infection_first_2l_mod_names.txt";
	//	String baseFileName = "c:\\research_data\\infection\\infection_first_2l_group\\infection_first_2l_group";
		
//	String baseFileName = "c:\\research_data\\zon\\2_16_07\\lenient\\zon_lenient";
//	String modNamesFileName = "c:\\research_data\\zon\\2_16_07\\expression_direction_mod_names.txt";
		
		String outputFileName = baseFileName + "_REP_tissueMatrix_graphic.svg";
		
		Container cp = getContentPane();
		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		int startEP = 76;
		int endEP = 104;
		
		Table table = createTable(baseFileName,modNamesFileName,startEP,endEP);
		cp.add(table);
		
		pack();
		
		try {
			table.outputSVG(outputFileName);
		} catch(IOException e) {
			System.out.println(e);
		}
	}
	
	public static void main(String[] args) {
		final OutputTissueMatrixGraphic frame = new OutputTissueMatrixGraphic();
		
		EventQueue.invokeLater(new Runnable() {
			public void run() { frame.setVisible(true); }
		});

	}
	
	public Table createTable(String baseFileName,String modNamesFileName,int startEP, int endEP) {
		 MicroArrayData data = new MicroArrayData();
		 Vector<Double> generality = null;
		 int i = 0;
		 int j = 0;
		 int nameCellWidth = 180;
		 int REPCellWidth = 19;
		 int generalityTickFrequency = 8;
		 int numEPs = endEP - startEP + 1;
		 
		 String tissueLoadFileName = null;
		 String groupsFileName = null;
		 String modFileName = null;
		 String generalityFileName = null;
		 
		 tissueLoadFileName = baseFileName + "REPs_tissueMatrix_revise.txt";
		 groupsFileName = baseFileName + "_groups.txt";
		 generalityFileName = baseFileName + "REPs_generality.txt";
		 modFileName = baseFileName + "REPs_modMatrix_revise.txt";
		 
		 MicroArrayData modData = new MicroArrayData();
		 ClusterReader groupsData = new ClusterReader();
		 ClusterReader generalityData = new ClusterReader();
		 ClusterReader modNames = new ClusterReader();
		 
		 try {
			 data.readFileBlindContinuous(tissueLoadFileName);
			 modData.readFileBlindDiscrete(modFileName);
			 groupsData.readFile(groupsFileName);
			 generalityData.readFile(generalityFileName);
			 modNames.readFile(modNamesFileName);
		 } catch(IOException e) {
			 System.out.println(e);
		 }
		 
		 generality = new Vector<Double>();
		 for (i=0;i<generalityData.clusters.get(0).size();i++) {
			 generality.add(new Double(generalityData.clusters.get(0).get(i)));
		 }
		 
			
		int headerRows = 1;
		Table table = new Table(data.numRows+headerRows+1,numEPs+2);
		TableCell cell = null;
		
		for (i=0;i<table.cells.length;i++) {
			for (j=0;j<table.cells[i].length;j++) {
				cell = new TableCell();
				table.setCell(i,j,cell);
			}
		}
		
		int numGroups = groupsData.clusters.size();
		
		// set up modifier colors
		int numModLevels = -2;
		for (i=0;i<modData.dvalues.length;i++) {
			for (j=0;j<modData.dvalues[i].length;j++) {
				if (modData.dvalues[i][j] > numModLevels) {
					numModLevels = modData.dvalues[i][j];
				}
			}
		}
		numModLevels++;
			
		Color[] modColors = createGroupColors(numModLevels);
		LinkedHashMap<String,Integer> tissueMap = new LinkedHashMap<String,Integer>();
		ArrayList<String> orderedTissues = new ArrayList<String>();
		int k = 0;
		for (i=0;i<groupsData.clusters.size();i++) {
			for (j=0;j<groupsData.clusters.get(i).size();j++) {
				orderedTissues.add(groupsData.clusters.get(i).get(j));
				cell = table.cells[k+1][0];
				cell.width = REPCellWidth;
				cell.text = " ";
				cell.fill = false;
				cell.bottomBorder = false;
				if (j > 0) {
					cell.topBorder = false;
					if (j == groupsData.clusters.get(i).size() - 1) {
						cell.bottomBorder = true;
					}
				} else {
				//	cell.text = (new Integer(i+1)).toString();
					cell.text = String.valueOf((char) (i+65));
				}
				k++;
			}
		}
		
		for (j=0;j<data.numRows;j++) {
			tissueMap.put(data.geneNames[j],j);
		}
			
		cell = table.cells[0][0];
		cell.width = REPCellWidth;
		
		Iterator<String> sIter = orderedTissues.iterator();
			
		// set up tissue cells
		i = 0;
		while(sIter.hasNext()) {
			cell = table.cells[i+1][1];
			cell.width = nameCellWidth;
			cell.text = sIter.next();
			i++;
		}
		cell = table.cells[0][1];
		cell.width = nameCellWidth;
			
	//	int[] sortOrder = VectorUtil.sortOrder(generality);
			
		for (j=startEP-1;j<endEP;j++) {
			cell = table.cells[0][j-(startEP-1)+2];
			cell.width = REPCellWidth;
		//	cell.text = ((new Integer(sortOrder[j] + 1))).toString();
			cell.text = ((new Integer(j + 1))).toString();
			cell.textAlign = TableCell.CENTER;
		}
			
		double maxLoad = 0.0;
		double minLoad = 1000.0;
		for (i=0;i<data.numRows;i++) {
			for (j=0;j<data.numCols;j++) {
				if (data.values[i][j] > maxLoad) {
					maxLoad = data.values[i][j];
				}
				if (data.values[i][j] < minLoad & data.values[i][j] > 0.0) {
					minLoad = data.values[i][j];
				}
			}
		}
			
		float color = 0.0f;
		double minIntensity = 0.05;
		double maxIntensity = 1.0;
		double v = 0.0;
		double a = (maxIntensity-minIntensity)/(maxLoad-minLoad);
		double b = minIntensity - minLoad*a;
		int modV = 0;
		float[] RGBc = new float[3]; 
		
		sIter = orderedTissues.iterator();
		int ri = 0;
		while(sIter.hasNext()) {
			i = tissueMap.get(sIter.next());
			for (j=startEP-1;j<endEP;j++) {
			//	k = sortOrder[j];
				k = j;
				
				modV = modData.dvalues[i][k];
				
				v = data.values[i][k]*a + b;
				color = (float) v;
				if (color > 1.0f) {
					color = 1.0f;
				}
				if (color < 0.0f) {
					color = 0.0f;
				}
				if (Math.abs(data.values[i][k]) > 0.0) {
					color = (float) Math.pow(color,1.0/2.0);
				}
				
				cell = table.cells[ri+1][j-(startEP-1)+2];
				
				if (modV >= 0) {
					modColors[modV].getRGBColorComponents(RGBc);
					cell.fillColor = new Color(RGBc[0],RGBc[1],RGBc[2],color);
					cell.textAlign = cell.CENTER;
					cell.text = modNames.clusters.get(0).get(modV);
				//	cell.textColor = new Color(0.0f,0.0f,0.0f,color);
				//	cell.text = (new Integer(modV+1)).toString();
				} else {
					cell.fillColor = Color.WHITE;
				}
				cell.fill = true;
				cell.width = REPCellWidth;
				cell.borderColor = Color.BLACK;
			}
			ri++;
		}
		
		cell = table.cells[table.cells.length -1][0];
		cell.topBorder = false;
		cell.bottomBorder = false;
		cell.rightBorder = false;
		cell.leftBorder = false;
		cell.width = REPCellWidth;
		
		cell = table.cells[table.cells.length -1][1];
		cell.topBorder = false;
		cell.bottomBorder = false;
		cell.rightBorder = false;
		cell.leftBorder = false;
		cell.width = nameCellWidth;
		
		double maxGenerality = -1.0;
		for (j=0;j<generality.size();j++) {
			if (generality.get(j) > maxGenerality) {
				maxGenerality = generality.get(j);
			}
		}
		
		i = generalityTickFrequency;
		
		for (j=startEP-1;j<endEP;j++) {
		//	k = sortOrder[j];
			k = j;
			
			cell = table.cells[table.cells.length -1][j-(startEP-1)+2];
			cell.topBorder = false;
			cell.bottomBorder = false;
			cell.rightBorder = false;
			cell.leftBorder = false;
			
			cell.width = REPCellWidth;
			
			if (i == generalityTickFrequency) {
				cell.leftBorder = true;
				cell.borderColor = Color.BLACK;
			//	cell.text = String.format("%.3f",generality.get(k)/maxGenerality);
				cell.text = String.format("%.3f",generality.get(k));
				cell.textAlign = TableCell.LEFT;
				i = 1;
			}
			
			i++;
		}
				
		return table;
	}
	
//	 disperse colors using HSL model
	 public static Color[] createGroupColors(int numGroups) {
		 double hInc = 360.0/((double) numGroups);
		 double sInc = 0.05/((double) numGroups);
		 double lInc = 0.05/((double) numGroups);
		 
		 int i = 0;
		 
		 double H = 0.0;
		 double S = 0.95;
		 double L = 0.60;
		 float[] hls = null;
		 float[] rgb = null;
		 
		 Color[] groupColors = new Color[numGroups];
		 
		 for (i=0;i<numGroups;i++) {
			 hls = new float[3];
			 hls[0] = (float) H;
			 hls[2] = (float) S;
			 hls[1] = (float) L;
			 rgb = HLStoRGB(hls);
			 groupColors[i] = new Color(rgb[0],rgb[1],rgb[2]);
			 H += hInc;
			 S += sInc;
			 L -= lInc;
		 }
		 
		 return groupColors;
	 }
	
	public static float RGBHelper(float q1,float q2,float hue2) {
		 float hue = hue2;
		 if (hue > 360) hue = hue - 360;
		 if (hue <0) hue = hue+360;
		 if (hue<60) {
			 return (q1+(q2-q1)*hue/60);
		 } else {
			 if (hue<180) {
				 return q2;
			 }
			 if (hue<240) {
				 return(q1+(q2-q1)*(240-hue)/60);
			 }
		 }
		 return q1;
	 }
	 
	 public static float[] HLStoRGB(float[] hls) {
		 float[] rgb = new float[3];
		 float p1,p2;
		 float H = hls[0];
		 float L = hls[1];
		 float S = hls[2];
		 
		 if (L <= 0.5f) {
			 p2 = L*(1.0f+S);
		 } else {
			 p2 = L + S-(L*S);
		 }
		 p1 = 2.0f*L-p2;
		 if (S==0) {
			 rgb[0] = L;	rgb[1] = L;	rgb[2] = L;
		 } else {
			 rgb[0] = RGBHelper(p1,p2,H+120);
			 rgb[1] = RGBHelper(p1,p2,H);
			 rgb[2] = RGBHelper(p1,p2,H-120);
		 }
		 return rgb;
	 }
}
