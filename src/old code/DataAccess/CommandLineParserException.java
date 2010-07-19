package edu.mit.csail.psrg.georg.DataAccess;

public class CommandLineParserException extends Exception {
	String description = "";
	public CommandLineParserException(String d) {
		description = d;
	}
	
	public String toString() {
		return "parameter type not recognized: " + description;
	}
}
