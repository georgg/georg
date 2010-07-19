package edu.mit.csail.psrg.georg.Graph2;

import java.awt.Container;
import java.awt.EventQueue;
import java.io.IOException;
import javax.swing.JFrame;
import edu.mit.csail.psrg.georg.DataAccess.ClusterReader;

public class OutputEntropyBarGraph extends JFrame {
	
	public OutputEntropyBarGraph() {
		String dirPath = "c:\\research_data\\mouse_human\\b47mm7\\geo\\EPGraph\\";
		String avgEntropyFileName = dirPath + "avg_tissue_entropy.txt";
		String outputFileName = dirPath + "entropy_bar_graph.svg";
		
		Container cp = getContentPane();
		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		
		GroupBarGraph graph = loadData(avgEntropyFileName);
		cp.add(graph);
		pack();
		try {
			graph.outputSVG(outputFileName);
		} catch(IOException e) {
			System.out.println(e);
		}
	}
	
	public static void main(String[] args) {
		final OutputEntropyBarGraph frame = new OutputEntropyBarGraph();
		EventQueue.invokeLater(new Runnable() {
			public void run() { frame.setVisible(true); }
		});
	}
	
	public GroupBarGraph loadData(String avgEntropyFileName) {
		 
		ClusterReader entropyData = new ClusterReader();
		try {
			 entropyData.readFile(avgEntropyFileName);
		 } catch(IOException e) {
			 System.out.println(e);
		 }
		 
		 int numItems = entropyData.clusters.size();
		 GroupBarGraph graph = new GroupBarGraph(numItems);
		 
		 int i = 0;
		 String type = null;
		 
		 for (i=0;i<entropyData.clusters.size();i++) {
			 graph.groupAssignment[i] = new Integer(entropyData.clusters.get(i).get(0));
			 graph.itemName[i] = entropyData.clusters.get(i).get(1);
			 type = entropyData.clusters.get(i).get(2);
			 if (type.equals("H")) {
				 graph.itemType[i] = 0;
			 } else {
				 graph.itemType[i] = 1;
			 }
			 graph.itemValue[i] = new Double(entropyData.clusters.get(i).get(3));
		 }
		 
		 return graph;
	}
}
