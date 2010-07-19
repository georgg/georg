package edu.mit.csail.psrg.georg.MouseHuman2;

import java.util.HashMap;
import java.util.LinkedHashMap;

import edu.mit.csail.psrg.georg.HDP2.ControllingTissue;
import edu.mit.csail.psrg.georg.HDP2.DirichletProcess;
import edu.mit.csail.psrg.georg.HDP2.GenericREP;
import edu.mit.csail.psrg.georg.HDP2.HDP;
import edu.mit.csail.psrg.georg.HDP2.REP;
import edu.mit.csail.psrg.georg.HDP2.REPOverlapBuilder;
import edu.mit.csail.psrg.georg.HDP2.SummarizeInferredModel;
import edu.mit.csail.psrg.georg.HDP2.TissueDP;

public class SummarizeModel {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
	//	String dirPath = "C:\\research_data\\mouse_human\\b47mm7\\geo\\UD\\";
	//	String dirPath = "C:\\research_data\\ramaswamy\\GeneProgram\\results\\rama_PM_log2\\";
	//	String dirPath = "C:\\research_data\\Bild\\Results\\bild_ovarian_log2_hier\\";
	//	String dirPath = "C:\\research_data\\mouse_human\\simulated\\";
		String dirPath = "C:\\research_data\\zon\\zfish_grps_gcrma\\";
		String snapFileName = "iter_snap_99999.persist";
	//	String refDirPath = "C:\\research_data\\mouse_human\\b47mm7\\";
	//	String humanAssocFile = refDirPath + "Hs.genelist.go.plus";
	//	String refDirPath = "C:\\research_data\\ramaswamy\\GeneProgram\\";
	//	String humanAssocFile = refDirPath + "Rama.hg18.genelist.go.plus";
	//	String refDirPath = "C:\\research_data\\Bild\\";
		String refDirPath = "C:\\research_data\\zon\\";
		String humanAssocFile = refDirPath + "zfin.go.plus";
	//	String classifyFileName = refDirPath + "Rama_PM_hierarchy.txt";
		String classifyFileName = null;
		
		SummarizeInferredModel sum = new SummarizeInferredModel();
		
		sum.refDirPath = refDirPath;
		sum.humanAssocFile = humanAssocFile;
		
		sum.outputFiles(dirPath,snapFileName,classifyFileName);
		
	/*	REPOverlapBuilder builder = new REPOverlapBuilder();
		builder.refDirPath = refDirPath;
		builder.humanAssocFile = humanAssocFile;
	//	builder.processSamples("C:\\research_data\\mouse_human\\simulated\\nohier_overlap\\iter_burned_snap_","C:\\research_data\\mouse_human\\simulated\\nohier_overlap\\REPs2",99,100,100);  	
		builder.processSamples("C:\\research_data\\zon\\zfish_gcrma\\iter_burned_snap_","C:\\research_data\\zon\\zfish_gcrma\\REPs2",99,100,100);  	
	*/
		
	//	testREP();
	}
	
	public static void testREP() {
		String dirPath = "C:\\research_data\\Bild\\Results\\test\\";
		String snapFileName = "iter_snap_9999.persist";
		HDP myHDP = HDP.restoreFromFile(dirPath + snapFileName);
		int j = 0;
		int i = 0;
		DirichletProcess dp = null;
		REP[] REPs = new REP[myHDP.numTopics];
		REP[] REPs2 = new REP[myHDP.numTopics];
		double p = 0.0;
		
		LinkedHashMap<HashMap<ControllingTissue,ControllingTissue>,REP> REPMap = new LinkedHashMap<HashMap<ControllingTissue,ControllingTissue>,REP>();
		
		for (i=0;i<myHDP.numTopics;i++) {
  			REPs[i] = new GenericREP(myHDP.totalGenes,myHDP.topics.get(i).posCounts);
  			REPs2[i] = new GenericREP(myHDP.totalGenes,myHDP.topics.get(i).posCounts);
  		}
		
		for (j=0;j<myHDP.DP.length;j++) {
  			dp = myHDP.DP[j];
  			if (dp.state != DirichletProcess.HELDOUT & dp.getClass() == TissueDP.class) {
				((TissueDP) dp).allocateGenesToTopicMap();
				for (i=0;i<myHDP.numTopics;i++) {
					p = REPs[i].controlPVal(myHDP.topics.get(i).topicMap,(TissueDP) dp,myHDP.totalGenes,5);
	  				if (p <= 0.05) {
	  					REPs[i].addControllingTissue((TissueDP) dp,myHDP.topics.get(i).topicMap,myHDP.topics.get(i).topicArray,((TissueDP) dp).upDownTopicChoice[i]);
	  					REPs2[i].addControllingTissue((TissueDP) dp,myHDP.topics.get(i).topicMap,myHDP.topics.get(i).topicArray,((TissueDP) dp).upDownTopicChoice[i]);
	  				}
				}
			}
		}
		
		for (i=0;i<REPs.length;i++) {
			REPs[i].addToREPMap(REPMap,myHDP.totalGenes);
		//	REPs2[i].addToREPMap(REPMap,myHDP.totalGenes);
		}
		
		System.out.println("Check point");
	}


}
