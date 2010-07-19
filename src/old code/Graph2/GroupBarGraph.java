package edu.mit.csail.psrg.georg.Graph2;

import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;

import javax.swing.JPanel;

import org.apache.batik.dom.GenericDOMImplementation;
import org.apache.batik.svggen.SVGGraphics2D;
import org.w3c.dom.DOMImplementation;
import org.w3c.dom.Element;

public class GroupBarGraph extends JPanel {
	
	public int maxBarWidth = 60;
	public double barHeightPercent = 0.67;
	public int horizontalBarSpace = 3;
	public int groupVerticalSpace = 5;
	public int groupHorizontalSpace = 5;
	public int groupColorWidth = 10;
	public int fontSize = 12;
	public String fontFamily = "Arial";
	public Color barColor = Color.BLUE;
	
	public int[] groupAssignment = null;
	public String[] itemName = null;
	public double[] itemValue = null;
	public int[] itemType = null;
	
	public GroupBarGraph(int numItems) {
		groupAssignment = new int[numItems];
		itemName = new String[numItems];
		itemValue = new double[numItems];
		itemType = new int[numItems];
	}
	
	public void paint(Graphics g) {
	        Graphics2D g2 = (Graphics2D) g;	        
	        g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
	        int i = 0;
	        
	        Font boldFont = new Font(fontFamily,Font.BOLD,fontSize);
	        Font italicFont = new Font(fontFamily,Font.ITALIC,fontSize);
	        FontMetrics fm = g2.getFontMetrics(boldFont);
			int height = fm.getLeading() + fm.getMaxAscent() + fm.getMaxDescent();
	        int maxNameWidth = 0;
	        int nameWidth = 0;
	        double v = 0.0;
	        int barWidth = 0;
	        int barHeight = (int) (((double) height)*barHeightPercent);
	        double minV = Double.POSITIVE_INFINITY;
	        double maxV = Double.NEGATIVE_INFINITY;
	        
	        for (i=0;i<itemName.length;i++) {
	        	nameWidth = fm.stringWidth(itemName[i]);
	        	if (nameWidth > maxNameWidth) {
	        		maxNameWidth = nameWidth;
	        	}
	        	if (itemValue[i] < minV) {
	        		minV = itemValue[i];
	        	}
	        	if (itemValue[i] > maxV) {
	        		maxV = itemValue[i];
	        	}
	        }
	        
	        int lastGroup = groupAssignment[0];
	        int yPos = height*2;
	        int xPos = 0;
	        int textY = 0;
	        int numGroups = -1;
	        
	        for (i=0;i<groupAssignment.length;i++) {
	        	if (numGroups < groupAssignment[i]) {
	        		numGroups = groupAssignment[i];
	        	}
	        }
	        numGroups++;
	        
	        Color[] groupColors = ColorGenerator.createGroupColors(numGroups);
	        
	        for (i=0;i<groupAssignment.length;i++) {
	        	xPos = 0;
	        	if (lastGroup != groupAssignment[i]) {
	        		yPos += groupVerticalSpace;
	        		lastGroup = groupAssignment[i];
	        	}
	        	g2.setColor(groupColors[groupAssignment[i]]);
				g2.fillRect(xPos,yPos,groupColorWidth,height);
				
				xPos += (groupColorWidth+groupHorizontalSpace);
				
				if (itemType[i] == 0) {
					g2.setFont(boldFont);
				} else {
					g2.setFont(italicFont);
				}
				textY = (height + fm.getAscent())/2 - fm.getDescent();
				g2.setColor(Color.BLACK);
				g2.drawString(itemName[i],xPos,yPos+textY);
				
				xPos += (maxNameWidth + horizontalBarSpace);
				
				v = (itemValue[i] - minV)/(maxV-minV);
				barWidth = (int) (((double) maxBarWidth)*v);
				g2.setColor(barColor);
				g2.fillRect(xPos,yPos+barHeight/2,barWidth,barHeight);
				
				yPos += height;
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
      
//    get the root element (the svg element)
      Element svgRoot = document.getDocumentElement();

      // Finally, stream out SVG to the standard output using UTF-8
      // character to byte encoding
      boolean useCSS = true; // we want to use CSS style attribute
      FileOutputStream fos = new FileOutputStream(fName);
      OutputStreamWriter out = new OutputStreamWriter(fos, "UTF-8");
      svgGenerator.stream(out, useCSS);
	}
}
