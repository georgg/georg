package edu.mit.csail.psrg.georg.Graph2;

import java.awt.Color;
import java.awt.Container;
import java.awt.EventQueue;
import java.io.IOException;
import java.util.Vector;

import javax.swing.JFrame;

import edu.mit.csail.psrg.georg.DataAccess.ClusterReader;
import edu.mit.csail.psrg.georg.DataAccess.MicroArrayData;
import edu.mit.csail.psrg.georg.StatUtil.VectorUtil;

public class OutputTopicMatrix extends JFrame {
	public OutputTopicMatrix() {
		String dirPath = "c:\\research_data\\mouse_human\\b47mm7\\geo\\EPGraph\\";
		String tissueLoadFileName = dirPath + "tissue_load.txt";
		String topicEntropyFileName = dirPath + "topic_entropy.txt";
		String outputFileName = dirPath + "topicMatrix.svg";
		
		Container cp = getContentPane();
		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		
		Table table = createTable(tissueLoadFileName,topicEntropyFileName);
		cp.add(table);
		
		pack();
	}
	
	public static void main(String[] args) {
		final OutputTopicMatrix frame = new OutputTopicMatrix();
		EventQueue.invokeLater(new Runnable() {
			public void run() { frame.setVisible(true); }
		});

	}
	
	public Table createTable(String tissueLoadFileName,String topicEntropyFileName) {
		 MicroArrayData data = new MicroArrayData();
		 Vector<Double> entropy = null;
		 int i = 0;
		 int j = 0;
		 int nameCellWidth = 75;
		 int topicCellWidth = 12;
		 
		 try {
			 data.readFile(tissueLoadFileName);
		 } catch(IOException e) {
			 System.out.println(e);
		 }
		 
		 if (topicEntropyFileName != null) {
			 ClusterReader entropyData = new ClusterReader();
			 try {
				 entropyData.readFile(topicEntropyFileName);
			 } catch(IOException e) {
				 System.out.println(e);
			 }
		 	
			 entropy = new Vector<Double>();
		
			 for (i=0;i<entropyData.clusters.size();i++) {
				 entropy.add(new Double(entropyData.clusters.get(i).get(0)));
			 }
		 }
		 
//		 normalize rows of the data matrix (loading on each tissue)
			double norm = 0.0;
			for (i=0;i<data.numRows;i++) {
				norm = 0.0;
				for (j=1;j<data.numCols;j++) {
					norm += data.values[i][j];
				}
				if (norm > 0.0) {
					for (j=1;j<data.numCols;j++) {
						data.values[i][j] = data.values[i][j]/norm;
					}
				}
			}
			
			int headerRows = 1;
			Table table = new Table(data.numRows+headerRows,data.numCols-1+2);
			TableCell cell = null;
			
			// generate group color bars
			int numGroups = 0;
			j = -2;
			for (i=0;i<data.numRows;i++) {
				if (((int) data.values[i][0]) != j) {
					j = ((int) data.values[i][0]);
					numGroups++;
				}
			}
			
			Color[] groupColors = createGroupColors(numGroups);
			
			j = -2;
			int grpIdx = -1;
			for (i=0;i<data.numRows;i++) {
				if (((int) data.values[i][0]) != j) {
					j = ((int) data.values[i][0]);
					grpIdx++;
				}
				cell = new TableCell();
				cell.width = topicCellWidth;
				cell.text = (new Integer(j)).toString();
				cell.fill = true;
				cell.fillColor = groupColors[grpIdx];
				table.setCell(i+1,0,cell);
			}
			cell = new TableCell();
			cell.width = topicCellWidth;
			table.setCell(0,0,cell);
			
			// set names
			for (i=0;i<data.numRows;i++) {
				cell = new TableCell();
				cell.width = nameCellWidth;
				cell.text = data.geneNames[i];
				table.setCell(i+1,1,cell);
			}
			cell = new TableCell();
			cell.width = nameCellWidth;
			table.setCell(0,1,cell);
			
			int[] sortOrder = VectorUtil.sortOrder(entropy);
			
			for (j=1;j<data.numCols;j++) {
				cell = new TableCell();
				cell.width = topicCellWidth;
				cell.text = ((new Integer(sortOrder[j-1] + 1))).toString();
				cell.textAlign = TableCell.CENTER;
				table.setCell(0,j+1,cell);
			}
			
			double maxLoad = 0.0;
			double minLoad = 1000.0;
			for (i=0;i<data.numRows;i++) {
				for (j=1;j<data.numCols;j++) {
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
			int k = 0;
			
			for (i=0;i<data.numRows;i++) {
				for (j=1;j<data.numCols;j++) {
					k = sortOrder[j-1];
					
					v = data.values[i][k+1]*a + b;
					color = (float) v;
					if (color > 1.0f) {
						color = 1.0f;
					}
					if (color < 0.0f) {
						color = 0.0f;
					}
					if (Math.abs(data.values[i][k+1]) > 0.0) {
						color = (float) Math.pow(color,1.0/3.0);
					}
					cell = new TableCell();
					cell.fillColor = new Color(color,0.0f,0.0f);
					cell.fill = true;
					cell.width = topicCellWidth;
					table.setCell(i+1,j+1,cell);	
				}
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
