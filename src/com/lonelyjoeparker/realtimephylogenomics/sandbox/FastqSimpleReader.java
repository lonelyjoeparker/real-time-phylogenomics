/**
 * 
 */
package com.lonelyjoeparker.realtimephylogenomics.sandbox;

import java.io.File;
import java.util.ArrayList;

import uk.ac.qmul.sbcs.evolution.convergence.util.CapitalisedFileReader;

/**
 * Simple parser for Fastq files. Not warranted or tested in any way whatsoever, at all!!
 * @author <a href="http://github.com/lonelyjoeparker">@lonelyjoeparker</a>
 * @since 30 Oct 2015
 * @version 0.1
 */
public class FastqSimpleReader {
	ArrayList<String> inputData;
	ArrayList<String> parsedData;
	File inputFile;
	
	/**
	 * Constructor. Takes input as java.io.File. Assumes input file
	 * contains FASTQ file formatted as sequential lines of:
	 * <ul>
	 * 	<li>Sequence name / identifier, '@' at first char. will be parsed</li>
	 * 	<li>Sequence of nucleosides, sequential (one line only). will be parsed</li>
	 * 	<li>optional ID repeater line, '+' at first line. will not be parsed</li>
	 * 	<li>QC scores for each position, sequential (one line only). will not be parsed</li>
	 * </ul>
	 * @param inputFile
	 */
	public FastqSimpleReader(File inputFile){
		inputData = new CapitalisedFileReader().loadSequences(inputFile);
		// parse input
		ArrayList<String> parsedData = new ArrayList<String>();
		boolean inSeqLine = false;
		for(String line:inputData){
			if(inSeqLine){
				parsedData.add(line);
				inSeqLine = false;
			}else{
				char firstChar = line.charAt(0);
				switch(firstChar){
					case '@':{parsedData.add(line);inSeqLine = true;}
					case '+':{inSeqLine = false;}
					default:{
						//do nothing
						}
					}
				}	
			}
		}

	public ArrayList<String> getParsedData() {
		return parsedData;
	}
}
