package edu.mit.csail.psrg.georg.Graph2;

import java.awt.Color;

public class ColorGenerator {
//	 disperse colors using HSL model
	 public static Color[] createGroupColors(int numGroups) {
		 double hInc = 360.0/((double) numGroups);
		 double sInc = 0.05/((double) numGroups);
		 double lInc = 0.10/((double) numGroups); 
		 
		 int i = 0;
		 
		 double H = 0.0;
		 double S = 0.95;
		 double L = 0.45; 
		 float[] hls = null;
		 float[] rgb = null;
		 
		 Color[] groupColors = new Color[numGroups];
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
