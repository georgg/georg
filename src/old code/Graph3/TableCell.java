package edu.mit.csail.psrg.georg.Graph3;

import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics2D;

public class TableCell {
	public int tableRow = 0;
	public int tableCol = 0;
	public Table myTable = null;
	public int xPos = 0;
	public int yPos = 0;
	public int width = 0;
	public int height = 0;
	public String text = "";
	public static int RIGHT = 1;
	public static int LEFT = 2;
	public static int CENTER = 3;
	public int textAlign = LEFT;
	public boolean fill = false;
	public Color fillColor = null;
	public Color textColor = Color.BLACK;
	public boolean topBorder = true;
	public boolean leftBorder = true;
	public boolean rightBorder = true;
	public boolean bottomBorder = true;
	public int borderWidth = 1;
	public Color borderColor = Color.BLACK;
	public Font myFont = new Font("Arial",Font.PLAIN,10);
	public int leftSpace = 2;
	public boolean doDraw = true;
	
	public void render(Graphics2D g) {
		if (!doDraw)
			return;
		
		int textX = 0;
		int textY = 0;
		
		Color tempColor = g.getColor();
		
		if (fill) {
			g.setColor(fillColor);
			g.fillRect(xPos,yPos,width,height);
			g.setColor(tempColor);
		}
		
		if (text == null) {
			return;
		}
		
		Font tempFont = g.getFont();
		g.setFont(myFont);
		
		g.setColor(textColor);
		FontMetrics fm = g.getFontMetrics(g.getFont());
		textY = (height + fm.getAscent())/2 - fm.getDescent();
		
		if (textAlign == TableCell.CENTER) {
			textX = (width-fm.stringWidth(text))/2;
		}
		
		if (textAlign == TableCell.RIGHT) {
			textX = (width-fm.stringWidth(text));
		}
		
		if (textAlign == TableCell.LEFT) {
			textX = leftSpace;
		}
		
		g.drawString(text,xPos+textX,yPos+textY);
		g.setColor(tempColor);
		g.setFont(myFont);
		
		if (checkTopOk()) {
			g.setColor(borderColor);
			g.drawLine(xPos,yPos,xPos+width,yPos);
			g.setColor(tempColor);
		}
		
		if (checkBottomOk()) {
			g.setColor(borderColor);
			g.drawLine(xPos,yPos+height,xPos+width,yPos+height);
			g.setColor(tempColor);
		}
		
		if (checkLeftOk()) {
			g.setColor(borderColor);
			g.drawLine(xPos,yPos,xPos,yPos+height);
			g.setColor(tempColor);
		}
		
		if (checkRightOk()) {
			g.setColor(borderColor);
			g.drawLine(xPos+width,yPos,xPos+width,yPos+height);
			g.setColor(tempColor);
		}
	}
	
	boolean checkTopOk() {
		return topBorder;
	}
	
	boolean checkBottomOk() {
		if (!bottomBorder)
			return false;
		
		if (tableRow == myTable.cells.length-1)
			return true;
		
		if (!myTable.cells[tableRow+1][tableCol].topBorder)
			return true;
		
		return false;
	}
	
	boolean checkLeftOk() {
		return leftBorder;
	}
	
	boolean checkRightOk() {
		if (!rightBorder)
			return false;
		
		if (tableCol == myTable.cells[tableRow].length -1)
			return true;
		
		if (!myTable.cells[tableRow][tableCol+1].leftBorder)
			return true;
		
		return false;
	}
	
	int calculateCellHeight(Graphics2D g) {
		FontMetrics fm = g.getFontMetrics(myFont);
		return (fm.getLeading() + fm.getMaxAscent() + fm.getMaxDescent() + 2*borderWidth);
	}
}
