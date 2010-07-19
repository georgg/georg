package edu.mit.csail.psrg.georg.Graph;

import java.awt.*;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.io.OutputStreamWriter;
import org.apache.batik.svggen.SVGGraphics2D;
import org.apache.batik.dom.GenericDOMImplementation;
import org.w3c.dom.DOMImplementation;
import org.w3c.dom.Element;


public class TabularComponent {
	int originX = 5;
	int originY = 5;
	int numRows = 0;
	int numCols = 0;
	int[] colX = null;
	int cellHeight = 0;
	int cellWidth = 0;
	int maxX = 100;
	int maxY = 100;
	TabCell[][] cells = null;
	Font myFont = null;
	boolean grid = false;
	int gridWidth = 1;
	int alternateInterval = 5;
	Color gridColor = Color.BLACK;
	Color altGridColor = Color.RED;
	
	public TabularComponent(int r,int c,int cw) {
		numRows = r;
		numCols = c;
		int i = 0;
		colX = new int[numCols+1];
		cellWidth = cw;
		
		cells = new TabCell[numRows][];
		for (i=0;i<numRows;i++) {
			cells[i] = new TabCell[numCols];
		}
		
		maxX = originX;
		for (i=0;i<numCols+1;i++) {
			colX[i] = maxX;
			maxX += cellWidth;
		}
		myFont = new Font("Arial",Font.PLAIN,12);
	}
	
	public int getNumRows() {
		return numRows;
	}
	
	public int getNumCols() {
		return numCols;
	}
	
	public void setFont(Font f) {
		myFont = f;
	}
	
	public void setGrid(boolean g) {
		grid = g;
	}
	
	public void initSize(SVGGraphics2D g) {
		cellHeight = calculateCellHeight(g);
		int i = 0;
		maxY = originY;
		for (i=0;i<numRows;i++) {
			maxY += cellHeight;
		}
		g.setSVGCanvasSize(new Dimension(maxX,maxY));
	}
	
	public void setColWidths(int[] colWidths) {
		int i = 0;
		maxX = originX;
		colX[0] = originX;
		for (i=1;i<numCols+1;i++) {
			maxX += colWidths[i-1];
			colX[i] = maxX;
		}
	}
	
	public void setCell(int r,int c,TabCell cell) {
		cells[r][c] = cell;
	}
	
	public void paint(SVGGraphics2D g) {
		g.setFont(myFont);
		initSize(g);
		g.setColor(Color.WHITE);
		int i = 0;
		int j = 0;
		int y = originY;
		int x = 0;
		for (i=0;i<numRows;i++) {
			for (j=0;j<numCols;j++) {
				x = colX[j];
				if (cells[i][j] != null) {
					cells[i][j].render(g,this,x,y,colX[j+1]-colX[j],cellHeight);
				}
			}
			y += cellHeight;
		}
		if (grid) {
			drawGrid(g);
		}
	}
	
	public void drawGrid(SVGGraphics2D g) {
		int i = 0;
		int y = originY;
		int altCounter = 0;
		
		Color tempColor = g.getColor();
		g.setColor(gridColor);
		
		for (i=0;i<numRows;i++) {
			altCounter++;
			if (altCounter == alternateInterval) {
				g.setColor(altGridColor);
			}
			g.drawLine(originX,y,colX[colX.length-1],y);
			if (altCounter == alternateInterval) {
				g.setColor(gridColor);
				altCounter = 0;
			}
			y += cellHeight;
		}
		
		int j = 0;
		for (j=0;j<colX.length;j++) {
			altCounter++;
			if (altCounter == alternateInterval) {
				g.setColor(altGridColor);
			}
			g.drawLine(colX[j],originY,colX[j],y);
			if (altCounter == alternateInterval) {
				g.setColor(gridColor);
				altCounter = 0;
			}
		}
		
		g.setColor(tempColor);
	}
	
	public Dimension getPreferredSize() {
		return new Dimension(maxX,maxY);
	}
	
	int calculateCellHeight(Graphics2D g) {
		FontMetrics fm = g.getFontMetrics(myFont);
		return fm.getLeading() + fm.getMaxAscent() + fm.getMaxDescent();
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
