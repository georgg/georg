package edu.mit.csail.psrg.georg.Graph2;

import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics2D;

public class TableCell {
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
	public Color borderColor = Color.GRAY;
	public Font myFont = new Font("Arial",Font.PLAIN,4);
	
	public void render(Graphics2D g) {
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
		
		g.drawString(text,xPos+textX,yPos+textY);
		g.setColor(tempColor);
		g.setFont(myFont);
		
		if (topBorder) {
			g.setColor(borderColor);
			g.drawLine(xPos,yPos,xPos+width,yPos);
			g.setColor(tempColor);
		}
		
		if (bottomBorder) {
			g.setColor(borderColor);
			g.drawLine(xPos,yPos+height,xPos+width,yPos+height);
			g.setColor(tempColor);
		}
		
		if (leftBorder) {
			g.setColor(borderColor);
			g.drawLine(xPos,yPos,xPos,yPos+height);
			g.setColor(tempColor);
		}
		
		if (rightBorder) {
			g.setColor(borderColor);
			g.drawLine(xPos+width,yPos,xPos+width,yPos+height);
			g.setColor(tempColor);
		}
	}
	
	int calculateCellHeight(Graphics2D g) {
		FontMetrics fm = g.getFontMetrics(myFont);
		return (fm.getLeading() + fm.getMaxAscent() + fm.getMaxDescent() + 2*borderWidth);
	}
}
