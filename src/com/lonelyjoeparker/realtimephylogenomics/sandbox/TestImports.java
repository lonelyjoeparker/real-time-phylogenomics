package com.lonelyjoeparker.realtimephylogenomics.sandbox;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;

import uk.ac.qmul.sbcs.evolution.convergence.AlignedSequenceRepresentation;
import uk.ac.qmul.sbcs.evolution.convergence.util.BasicFileWriter;
import uk.ac.qmul.sbcs.evolution.convergence.util.VerboseSystemCommand;

public class TestImports {

	String [] args;
	WatchDir watcher;
	AlignedSequenceRepresentation asr;
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		new TestImports(args).go();
	}
	
	public TestImports(String[] arguments){
		this.args = arguments;
	}
	
	private void go(){
	       // parse arguments
        if (args.length == 0 || args.length > 2)
            WatchDir.usage();
        boolean recursive = false;
        int dirArg = 0;
        if (args[0].equals("-r")) {
            if (args.length < 2)
                WatchDir.usage();
            recursive = true;
            dirArg++;
        }

        // register directory and process its events
        Path dir = Paths.get(args[dirArg]);
        try {
			watcher = new WatchDir(dir, recursive);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
        /* 
         * process dir events 
         * run in separate thread
         * 
        watcher.processEvents();
         * */
        asr = new AlignedSequenceRepresentation();
        asr.getNumberOfSites();
        VerboseSystemCommand vsc = new VerboseSystemCommand("echo foo");
        System.out.println(vsc.output.toString());
        FastqSimpleReader fqr = new FastqSimpleReader(new File("./examples/example_fastq.fastq"));
        System.out.println(fqr.getParsedData());
        File outputFile = new File("./examples/example.output.converted.fasta");
        new BasicFileWriter(outputFile,fqr.getParsedData());
	}

}
