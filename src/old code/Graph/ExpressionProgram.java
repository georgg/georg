package edu.mit.csail.psrg.georg.Graph;

import java.awt.Color;
import java.awt.FontMetrics;
import java.awt.Graphics2D;
import java.awt.geom.AffineTransform;
import java.awt.geom.Arc2D;
import java.awt.geom.Line2D;
import java.awt.geom.Ellipse2D;
import java.util.ArrayList;

public class ExpressionProgram implements Comparable {
	public double xCenter = 100.0;
	public double yCenter = 100.0;
	public double labelRadius = 15.0;
	public double innerRadius = 25.0;
	public double groupRadius = innerRadius+7.0;
	public double outerRadius = groupRadius+13.0;
	
	public double textRadius = outerRadius+10.0;
	public double textOffset = 10.0;
	
	public double entropy = 0.0;
	
	ExpressionProgramGraph myEPG = null;
	ArrayList<Integer> tissues = null;
	ArrayList<Integer> groups = new ArrayList<Integer>();
	ArrayList<Integer> groupCounts = new ArrayList<Integer>();
	
	// specifies the intensities for tissues
	public ArrayList<Double> tissueIntensity = null;
	
	public String programLabel = "";
	
	public ExpressionProgram(ExpressionProgramGraph epg,ArrayList<Integer> t,ArrayList<Double> ti) {
		myEPG = epg;
		tissues = t;
		tissueIntensity = ti;
		
		int i = 0;
		int[] groupTally = new int[myEPG.groupNames.length];
		for (i=0;i<tissues.size();i++) {
			groupTally[myEPG.tissueToGroups[tissues.get(i)]]++;
		}
		for (i=0;i<groupTally.length;i++) {
			if (groupTally[i] > 0) {
				groups.add(i);
				groupCounts.add(groupTally[i]);
			}
		}
	}
	
	public void draw(Graphics2D g) {
		int i = 0;
		
		AffineTransform origTransform = g.getTransform();
		
		double angleInc = 360.0/((double) tissues.size());
		Arc2D arc = null;
		double angleStart = 0.0;
		double angle = 0.0;
		
		// draw labels and spikes
		angleStart = 0.0;
		double xEnd = 0.0;
		double yEnd = 0.0;
		FontMetrics fm = g.getFontMetrics(g.getFont());
	//	double textY = ((double)fm.getAscent())/2.0 - ((double)fm.getDescent());
		double textY = ((double)fm.getAscent())/2.0;
		double textX = 0.0;
		String label = "";
		for (i=0;i<tissueIntensity.size();i++) {
			angle = angleStart;
			if (tissueIntensity.size() > 1)
				angle += angleInc/2.0;
			angle = angle*Math.PI/(180.0);
			xEnd = textRadius*Math.cos(angle) + xCenter;
			yEnd = textRadius*Math.sin(angle) + yCenter;
			g.setColor(Color.BLACK);
		//	g.draw(new Line2D.Double(xCenter,yCenter,xEnd,yEnd));
			
			
			label = (new Integer(tissues.get(i)+1)).toString();
			textX = ((double)fm.stringWidth(label))/2.0;
			xEnd = (textRadius+textOffset)*Math.cos(angle);
			yEnd = (textRadius+textOffset)*Math.sin(angle);
			
			g.translate(xCenter,yCenter);
			g.rotate(angle);
			g.translate(textRadius+textOffset,textY);
			
			if (yEnd < 0.0) {
			//	yEnd -= textY;
			} else {
				yEnd += textY;
			}
			
			if (xEnd < 0.0) {
				xEnd -= textX;
			} else {
			//	xEnd += textX;
			}
			
		//	g.drawString(label,(float) (xEnd+xCenter),(float) (yEnd+yCenter));
			
		//	g.drawString(label,0,0);
			g.setTransform(origTransform);
			angleStart += angleInc;
		}
		
		// draw the outer circle
		Ellipse2D outerCircle = new Ellipse2D.Double(xCenter-outerRadius,yCenter-outerRadius,2*outerRadius,2*outerRadius);
		
		g.setColor(Color.WHITE);
		g.fill(outerCircle);
		g.setColor(Color.BLACK);
		g.draw(outerCircle);
		
		Ellipse2D groupCircle = new Ellipse2D.Double(xCenter-groupRadius,yCenter-groupRadius,2*groupRadius,2*groupRadius);
		g.setColor(Color.WHITE);
		g.fill(groupCircle);
		g.setColor(Color.BLACK);
		g.draw(groupCircle);
		
		angleStart = 360.0;
		for (i=0;i<groups.size();i++) {
			angle = angleInc * ((double) groupCounts.get(i));
			arc = new Arc2D.Double();
		//	arc.setArcByCenter(xCenter,yCenter,outerRadius,angleStart-angle,angle,Arc2D.PIE);
			arc.setArcByCenter(xCenter,yCenter,groupRadius,angleStart-angle,angle,Arc2D.PIE);
			angleStart -= angle;
			g.setColor(myEPG.groupColors[groups.get(i)]);
		//	g.setColor(Color.WHITE);
			if (groups.size() > 1) { 
				g.fill(arc);
			} else {
			//	g.fill(new Ellipse2D.Double(xCenter-outerRadius,yCenter-outerRadius,2*outerRadius,2*outerRadius));
				g.fill(new Ellipse2D.Double(xCenter-groupRadius,yCenter-groupRadius,2*groupRadius,2*groupRadius));
			}
			g.setColor(Color.BLACK);
			if (groups.size() > 1) {
				g.draw(arc);
			}
		}
		
		// draw group labels
		angleStart = 0.0;
		fm = g.getFontMetrics(g.getFont());
	//	textY = ((double)fm.getAscent())/2.0 - ((double)fm.getDescent());
		textY = ((double)fm.getAscent())/2.0;
		textX = 0.0;
		for (i=0;i<groups.size();i++) {
			angle = angleStart;
			if (groups.size() > 1)
				angle += angleInc * ((double) groupCounts.get(i))/2.0;
			angle = angle*Math.PI/(180.0);
			xEnd = (innerRadius+textOffset)*Math.cos(angle);
			yEnd = (innerRadius+textOffset)*Math.sin(angle);
			
			label = myEPG.groupNames[groups.get(i)];
			textX = ((double)fm.stringWidth(label))/2.0;
			
			if (yEnd < 0.0) {
			//	yEnd -= textY;
			} else {
				yEnd += textY;
			}
			
			if (xEnd < 0.0) {
				xEnd -= textX;
			} else {
			//	xEnd += textX;
			}
			
			g.translate(xCenter,yCenter);
			g.rotate(angle);
		//	g.translate(innerRadius+textX,textY);
			g.translate(innerRadius+textOffset,textY);
			
		//	g.drawString(label,(float) (xEnd+xCenter),(float) (yEnd+yCenter));
			g.setColor(Color.BLACK);
			g.drawString(label,0,0);
			g.setTransform(origTransform);
			angleStart += angleInc*((double) groupCounts.get(i));
			g.setColor(Color.BLACK);
		}
		
		Ellipse2D innerCircle = new Ellipse2D.Double(xCenter-innerRadius,yCenter-innerRadius,2*innerRadius,2*innerRadius);
		g.setColor(Color.WHITE);
		g.fill(innerCircle);
		g.setColor(Color.BLACK);
		g.draw(innerCircle);
		
		float color = 0.0f;
		double v = 0.0;
		
		// draw tissue intensities
		angleStart = 360.0;
		double intensity = 0.0;
		for (i=0;i<tissueIntensity.size();i++) {
			arc = new Arc2D.Double();
			arc.setArcByCenter(xCenter,yCenter,innerRadius,angleStart-angleInc,angleInc,Arc2D.PIE);
			angleStart -= angleInc;
			intensity = tissueIntensity.get(i);
			v = intensity*myEPG.aInterp + myEPG.bInterp;
			color = (float) v;
			if (color > 1.0f) {
				color = 1.0f;
			}
			if (color < 0.0f) {
				color = 0.0f;
			}
			if (Math.abs(intensity) > 0.0) {
				color = (float) Math.pow(color,1.0/3.0);
			}
		//	g.setColor(new Color(color,0.0f,0.0f));
			g.setColor(new Color(1.0f-0.75f*color,1.0f-0.75f*color,1.0f-0.75f*color));
			if (tissueIntensity.size() > 1) {
				g.fill(arc);
			} else {
				g.fill(new Ellipse2D.Double(xCenter-innerRadius,yCenter-innerRadius,2*innerRadius,2*innerRadius));
			}
		//	g.setColor(Color.GRAY);
			g.setColor(Color.BLACK);
			if (tissueIntensity.size() > 1) 
				g.draw(arc);
		}
		
		textX = ((double)fm.stringWidth(programLabel))/2.0;
		textY = ((double)fm.getAscent())/2.0;
		Ellipse2D labelCircle = new Ellipse2D.Double(xCenter-labelRadius,yCenter-labelRadius,2*labelRadius,2*labelRadius);
		g.setColor(Color.WHITE);
		g.fill(labelCircle);
		g.setColor(Color.BLACK);
		g.draw(labelCircle);
		
		g.setColor(Color.BLACK);
	//	g.drawString(programLabel,(float) (xCenter-textRadius-10.0-textX),(float) (yCenter));
		g.drawString(programLabel,(float) (xCenter-textX),(float) (yCenter+textY));
	}

	public int compareTo(Object o) {
		ExpressionProgram p2 = (ExpressionProgram) o;

		if (p2.entropy < entropy)
			return 1;
		
		if (p2.entropy > entropy)
			return -1;
		
		return 0;
	}
}
