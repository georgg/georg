package edu.mit.csail.psrg.georg.Graph;

import java.awt.Graphics;
import java.awt.Graphics2D;

public abstract class TabCell {
	public abstract void render(Graphics2D g,TabularComponent parent, int x, int y,int w,int h);
}
