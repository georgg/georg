package edu.mit.csail.psrg.georg.Graph;

import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.color.ColorSpace;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Rectangle2D;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Vector;

import javax.swing.JPanel;

import org.apache.batik.dom.GenericDOMImplementation;
import org.apache.batik.svggen.SVGGraphics2D;
import org.w3c.dom.DOMImplementation;
import org.w3c.dom.Element;

import edu.mit.csail.psrg.georg.DataAccess.ClusterReader;
import edu.mit.csail.psrg.georg.DataAccess.MicroArrayData;

public class ExpressionProgramGraph extends JPanel {
	public ArrayList<ExpressionProgram> programs = new ArrayList<ExpressionProgram>();
	
	// for each tissue, designate the group it belongs to
	public int[] tissueToGroups = null;
	
	public String[] groupNames = null;
	public Color[] groupColors = null;
	
	int numEntropyBins = 15;
	double[] entropyIntervals = null;
	ArrayList<Vector<ExpressionProgram>> entropyBins = new ArrayList<Vector<ExpressionProgram>>();
	
	public double pageHeight = 950.0;
	public double pageWidth = 1350.0;
	
	public double legendBoxHeight = 20.0;
	public double legendBoxWidth = 30.0;
	
	public boolean outputLegend = false;
	
	public double aInterp = 0.0;
	public double bInterp = 0.0;
	
	public void addExpressionProgram(ExpressionProgram p) {
		programs.add(p);
	}
	
	 public void paint(Graphics g) {
	        Graphics2D g2 = (Graphics2D) g;	        
	        g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
	        Font myFont = new Font("Arial",Font.BOLD,10);
	        g2.setFont(myFont);

	        int i = 0;
	        
	        if (!outputLegend) {
	        	for (i=0;i<programs.size();i++) {
	        		programs.get(i).draw(g2);
	        	}
	        } else {
	        	drawLegend(g2);
	        }
	 }
	 
	 // disperse colors using HSL model
	 void createGroupColors(int numGroups) {
		 double hInc = 360.0/((double) numGroups);
		 double sInc = 0.05/((double) numGroups);
		 double lInc = 0.10/((double) numGroups); 
		 
	/*	 double hInc = 360.0/((double) numGroups);
		 double sInc = 0.50/((double) numGroups);
		 double lInc = 0.10/((double) numGroups); */
		 
		 int i = 0;
		 
		 double H = 0.0;
		 double S = 0.95;
		 double L = 0.45; 
	/*	 double H = 0.0;
		 double S = 0.50;
		 double L = 0.50; */
		 float[] hls = null;
		 float[] rgb = null;
		 
		 groupColors = new Color[numGroups];
		 Color color = null;
		 
		 for (i=0;i<numGroups;i++) {
			 hls = new float[3];
			 hls[0] = (float) H;
			 hls[2] = (float) S;
			 hls[1] = (float) L;
			 rgb = HLStoRGB(hls);
			 color = new Color(rgb[0],rgb[1],rgb[2]);
			 if (((double) i)/2.0 - ((double) (i/2)) == 0.0) {
				 groupColors[i] = color;
			 } else {
				 groupColors[numGroups-i-1] = color;
			 }
			 H += hInc;
			 S += sInc;
			 L -= lInc;
		 }
	 }
	 
	 float RGBHelper(float q1,float q2,float hue2) {
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
	 
	 public float[] HLStoRGB(float[] hls) {
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
	 
	 public void loadData(String tissueLoadFileName,String topicEntropyFileName) {
		 MicroArrayData data = new MicroArrayData();
		 
		 try {
			 data.readFile(tissueLoadFileName);
		 } catch(IOException e) {
			 System.out.println(e);
		 }
		 
		 int i = 0;
		 int j = 0;
		 
		 int numGroups = 0;
		 tissueToGroups = new int[data.numRows];
		 for (i=0;i<data.numRows;i++) {
			 tissueToGroups[i] = (int) data.values[i][0];
			 if (tissueToGroups[i] > numGroups) {
				 numGroups = tissueToGroups[i];
			 }
		 }
		 numGroups++;
		 
		 createGroupColors(numGroups);
		 
		 groupNames = new String[numGroups];
		 for (i=0;i<numGroups;i++) {
			 groupNames[i] = String.valueOf((char) (i+65));
		 }
		 
		 double maxLoad = -50000;
		 double minLoad = 50000;
		 
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
		 
		 double minIntensity = 0.05;
		 double maxIntensity = 1.0;
		 
		 ArrayList<Double> tissueIntensity = null;
		 ArrayList<Integer> tissues = null;
		 ExpressionProgram program = null;
		 for (j=1;j<data.numCols;j++) {
	//	 	j = 24;
		 	tissueIntensity = new ArrayList<Double>();
			tissues = new ArrayList<Integer>();
			for (i=0;i<data.numRows;i++) {
				if (data.values[i][j] < minLoad) {
					minLoad = data.values[i][j];
				}
				if (data.values[i][j] > maxLoad) {
					maxLoad = data.values[i][j];
				}
				if (data.values[i][j] > 0.0) {
					tissues.add(i);
					tissueIntensity.add(data.values[i][j]);
				}
			}
			program = new ExpressionProgram(this,tissues,tissueIntensity);
			program.programLabel = (new Integer(j)).toString();
			addExpressionProgram(program);
		 }
		 	aInterp = (maxIntensity-minIntensity)/(maxLoad-minLoad);
		 	bInterp = minIntensity - minLoad*aInterp;
			loadTopicEntropy(topicEntropyFileName);
	 }
	 
	 public void loadTopicEntropy(String fname) {
		 ClusterReader data = new ClusterReader();
		 try {
			 data.readFile(fname);
		 } catch(IOException e) {
			 System.out.println(e);
		 }
		 
		 entropyIntervals = new double[numEntropyBins];
		 double minEntropy = 10e8;
		 double maxEntropy = -10e8;
		 double entropy = 0.0;
		 int i = 0;
		 int j = 0;
		 for (i=0;i<data.clusters.size();i++) {
			 entropy = new Double(data.clusters.get(i).get(0));
			 programs.get(i).entropy = entropy;
			 if (entropy > maxEntropy) {
				 maxEntropy = entropy;
			 }
			 if (entropy < minEntropy) {
				 minEntropy = entropy;
			 }
		 }
		 
		 double entropyInc = (maxEntropy-minEntropy + 10e-6)/((double) numEntropyBins);
		 entropy = minEntropy + entropyInc;
		 for (i=0;i<numEntropyBins;i++) {
			 entropyIntervals[i] = entropy;
			 entropy += entropyInc;
		 }
		 
		 for (i=0;i<numEntropyBins;i++) {
			 entropyBins.add(new Vector<ExpressionProgram>());
		 }
		 
		 Vector<ExpressionProgram> bin = null;
		 for (j=0;j<programs.size();j++) {
			 entropy = programs.get(j).entropy;
			 i = 0;
			 while (entropy > entropyIntervals[i]) {
				 i++;
			 }
			 bin = entropyBins.get(i);
			 bin.add(programs.get(j));
		 }
		 
		 System.out.println("Entropy bins");
		 arrangePrograms();
	 }
	 
	 public void arrangePrograms() {
		 int numNonZeroBins = 0;
		 int i = 0;
		 int j = 0;
		 
		 // for now, get rid of first bin
		 entropyBins.get(0).clear();
		 
		 for (i=0;i<entropyBins.size();i++) {
			 if (entropyBins.get(i).size() > 0) {
				 numNonZeroBins++;
			 }
		 }
		 
	//	 double xInc = pageWidth/((double) numNonZeroBins);
		 double xInc = 90.0;
		 double yInc = 90.0;
		 
		 double xPos = xInc/2.0;
		 double yPos = pageHeight;
		
		 Vector<ExpressionProgram> bin = null;
		 ExpressionProgram program = null;
		 for (i=0;i<entropyBins.size();i++) {
			 yPos = pageHeight - yInc/2.0;
			 bin = entropyBins.get(i);
			 Collections.sort(bin);
			 for (j=0;j<bin.size();j++) {
				 program = bin.get(j);
				 program.xCenter = xPos;
				 program.yCenter = yPos;
				 yPos -= yInc;
			 }
			 if (bin.size() > 0) {
				 xPos += xInc;
				 System.out.println((entropyIntervals[i-1]+entropyIntervals[i])/2.0);
			 }
		 }
	 }
	 
	 public void drawLegend(Graphics2D g) {
		 int i = 0;
		 double boxStart = 0.0;
		 double textX = 0.0;
		 double textY = 0.0;
		 
		 Rectangle2D rect = null;
		 FontMetrics fm = g.getFontMetrics(g.getFont());
		 textY = ((double)fm.getAscent())/2.0;
		 String label;
		 for (i=0;i<groupNames.length;i++) {
			 g.setColor(groupColors[i]);
			 rect = new Rectangle2D.Double(0.0,boxStart,legendBoxWidth,legendBoxHeight);
			 label = groupNames[i];
			 textX = ((double)fm.stringWidth(label))/2.0;
			 g.fill(rect);
			 g.setColor(Color.BLACK);
			 g.draw(rect);
		//	 g.drawString(label,(float) (11.0),(float) (boxStart + legendBoxHeight/2.0 + textY));
			 boxStart += legendBoxHeight;
		 }
	 }
	 
	 public void outputSVG(String fName) throws IOException {
//		 Get a DOMImplementation
        DOMImplementation domImpl =
            GenericDOMImplementation.getDOMImplementation();

        // Create an instance of org.w3c.dom.Document
        org.w3c.dom.Document document = domImpl.createDocument(null, "svg", null);

        // Create an instance of the SVG Generator
        SVGGraphics2D svgGenerator = new SVGGraphics2D(document);

        // Ask the test to render into the SVG Graphics2D implementation
        paint(svgGenerator);
        
//      get the root element (the svg element)
        Element svgRoot = document.getDocumentElement();

        // Finally, stream out SVG to the standard output using UTF-8
        // character to byte encoding
        boolean useCSS = true; // we want to use CSS style attribute
        FileOutputStream fos = new FileOutputStream(fName);
        OutputStreamWriter out = new OutputStreamWriter(fos, "UTF-8");
        svgGenerator.stream(out, useCSS);
	}

}
