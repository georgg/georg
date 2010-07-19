package edu.mit.csail.psrg.georg.Graph;

import java.awt.Container;
import java.awt.EventQueue;
import java.io.IOException;

import javax.swing.JFrame;

public class RunProgramGraph extends JFrame {

	/**
	 * @param args
	 */
	public RunProgramGraph() {
		String dirPath = "c:\\research_data\\mouse_human\\b47mm7\\geo\\EPGraph\\";
		String tissueLoadFile = dirPath + "tissue_load.txt";
		String topicEntropyFile = dirPath + "topic_entropy.txt";
		String SVGoutFile = dirPath + "EPGraph_high_entropy3.svg";
		String LegendOutFile = dirPath + "Legend.svg";
		
		Container cp = getContentPane();
		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		
		ExpressionProgramGraph graph = new ExpressionProgramGraph();
		graph.loadData(tissueLoadFile,topicEntropyFile);
	//	graph.outputLegend = true;
		
		cp.add(graph);
		
		pack();
		try {
			graph.outputSVG(SVGoutFile);
		//	graph.outputSVG(LegendOutFile);
		} catch(IOException e) {
			System.out.println(e);
		}
	}
	
	public static void main(String[] args) {
		final RunProgramGraph frame = new RunProgramGraph();
		EventQueue.invokeLater(new Runnable() {
			public void run() { frame.setVisible(true); }
		});

	}

}
