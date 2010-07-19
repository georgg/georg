package edu.mit.csail.psrg.georg.Graph;

import java.awt.Color;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Graphics2D;

public class TabText extends TabCell {
	String myText = "";
	public static int RIGHT = 1;
	public static int LEFT = 2;
	public static int CENTER = 3;
	int textAlign = LEFT;
	boolean fill = false;
	Color fillColor = null;
	Color textColor = Color.BLACK;
	
	public void setTextAlign(int s) {
		textAlign = s;
	}
	
	public void setFill(boolean f) {
		fill = f;
	}
	
	public void setFillColor(float r,float g,float b) {
		fillColor = new Color(r,g,b);
	}
	
	public void setFillColor(Color c) {
		fillColor = c;
	}
	
	public TabText(String t) {
		myText = t;
	}
	
	public void render(Graphics2D g,TabularComponent parent, int x, int y,int w,int h) {
		int textX = 0;
		int textY = 0;
		
		Color tempColor = null;
		
		if (fill) {
			tempColor = g.getColor();
			g.setColor(fillColor);
			g.fillRect(x,y,w,h);
			g.setColor(tempColor);
		}
		
		if (myText == null) {
			return;
		}
		
		tempColor = g.getColor();
		g.setColor(textColor);
		
		FontMetrics fm = g.getFontMetrics(g.getFont());
		textY = (h + fm.getAscent())/2 - fm.getDescent();
		
		if (textAlign == TabText.CENTER) {
			textX = (w-fm.stringWidth(myText))/2;
		}
		
		if (textAlign == TabText.RIGHT) {
			textX = (w-fm.stringWidth(myText));
		}
		
		g.drawString(myText,x+textX,y+textY);
		g.setColor(tempColor);
	}

}
