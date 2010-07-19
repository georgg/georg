package edu.mit.csail.psrg.georg.DataAccess;

import java.io.IOException;

public class DataAccessTest {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		MicroArrayData myData = new MicroArrayData();
		try {
			myData.readFile("C:\\CPPDataFiles\\test_docs.txt");
			System.out.println("Done.");
		} catch(IOException e) {
			System.out.println(e);
		}
	}

}
