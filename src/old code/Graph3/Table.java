package edu.mit.csail.psrg.georg.Graph3;

import java.awt.Color;
import java.awt.Dimension;
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

public class Table extends JPanel {
	int originX = 5;
	int originY = 5;
	int numRows = 0;
	int numCols = 0;
	int maxX = 100;
	int maxY = 100;
	public TableCell[][] cells = null;
	
	public Table(int r,int c) {
		numRows = r;
		numCols = c;
		int i = 0;
		
		cells = new TableCell[numRows][];
		for (i=0;i<numRows;i++) {
			cells[i] = new TableCell[numCols];
		}
	}
	
	public int getNumRows() {
		return numRows;
	}
	
	public int getNumCols() {
		return numCols;
	}
	
	public void setCell(int r,int c,TableCell cell) {
		cells[r][c] = cell;
		cell.tableCol = c;
		cell.tableRow = r;
		cell.myTable = this;
	}
	
	public void paint(Graphics g) {
		Graphics2D g2 = (Graphics2D) g;	 
		g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
		
		int i = 0;
		int j = 0;
		int x = originX;
		int y = originY;
		for (i=0;i<numRows;i++) {
			x = originX;
			for (j=0;j<numCols;j++) {
				if (cells[i][j] != null) {
					cells[i][j].xPos = x;
					cells[i][j].yPos = y;
					cells[i][j].height = cells[i][j].calculateCellHeight(g2);
					x += cells[i][j].width;
				}
			}
			if (cells[i][0] != null) {
				y += cells[i][0].height;
			}
		}
		
		for (i=0;i<numRows;i++) {
			for (j=0;j<numCols;j++) {
				cells[i][j].render(g2);
			}
		}
	}
	
	public Dimension getPreferredSize() {
		return new Dimension(maxX,maxY);
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
       
//     get the root element (the svg element)
       Element svgRoot = document.getDocumentElement();

       // Finally, stream out SVG to the standard output using UTF-8
       // character to byte encoding
       boolean useCSS = true; // we want to use CSS style attribute
       FileOutputStream fos = new FileOutputStream(fName);
       OutputStreamWriter out = new OutputStreamWriter(fos, "UTF-8");
       svgGenerator.stream(out, useCSS);
	}
}
