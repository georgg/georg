package edu.mit.csail.psrg.georg.Graph;

import javax.swing.*;
import java.awt.*;
import java.io.IOException;

public class OutputYeastTFTabular extends JFrame {

	/**
	 * @param args
	 */
	public OutputYeastTFTabular() {
	//	Container cp = getContentPane();
	//	setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		
		TabularComponent sheet = new TabularComponent(2,2,100);
		sheet.setCell(0,0,new TabText("Hello"));
		sheet.setCell(0,1,new TabText("Buffo cowboy"));
		sheet.setCell(1,1,new TabText("Amy's daily"));
		TabText t11 = new TabText("");
		t11.setFill(true);
		t11.setFillColor(0.5f,0.0f,0.0f);
		sheet.setCell(1,0,t11);
		int[] colWidths = {50,100};
		sheet.setColWidths(colWidths);
		
	//	cp.add(sheet);
		try {
			sheet.outputSVG("c:\\research_data\\test.svg");
		} catch (IOException e) {
			System.out.println(e);
		}
		
		pack();
	}
	
	public static void main(String[] args) {
	/*	final OutputYeastTFTabular frame = new OutputYeastTFTabular();
		EventQueue.invokeLater(new Runnable() {
			public void run() { frame.setVisible(true); }
		}); */
		OutputYeastTFTabular frame = new OutputYeastTFTabular();
		frame.setVisible(true);
	}

}
