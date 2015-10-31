/**
 * 
 */
package com.lonelyjoeparker.realtimephylogenomics;

import static java.nio.file.LinkOption.NOFOLLOW_LINKS;
import static java.nio.file.StandardWatchEventKinds.ENTRY_CREATE;
import static java.nio.file.StandardWatchEventKinds.ENTRY_DELETE;
import static java.nio.file.StandardWatchEventKinds.ENTRY_MODIFY;
import static java.nio.file.StandardWatchEventKinds.OVERFLOW;

import java.io.File;
import java.io.IOException;
import java.nio.file.FileSystems;
import java.nio.file.FileVisitResult;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.SimpleFileVisitor;
import java.nio.file.WatchEvent;
import java.nio.file.WatchKey;
import java.nio.file.WatchService;
import java.nio.file.attribute.BasicFileAttributes;
import java.util.HashMap;
import java.util.Map;

import com.lonelyjoeparker.realtimephylogenomics.sandbox.FastqSimpleReader;

import uk.ac.qmul.sbcs.evolution.convergence.AlignedSequenceRepresentation;
import uk.ac.qmul.sbcs.evolution.convergence.handlers.RAxMLAnalysisSGE;
import uk.ac.qmul.sbcs.evolution.convergence.util.BasicFileWriter;
import uk.ac.qmul.sbcs.evolution.convergence.util.TaxaLimitException;
import uk.ac.qmul.sbcs.evolution.convergence.util.VerboseSystemCommand;

/**
 * A simple class to monitor a directory, and (when new .fastq files appear)
 * align, them and build a phylogeny
 * @author <a href="http://github.com/lonelyjoeparker">@lonelyjoeparker</a>
 * @since 30 Oct 2015
 * @version 0.1
 */
public class RealtimeReadAnalyzer {

	private final Map<WatchKey,Path> keys;
	private boolean trace = false;
	private final WatchService watcher;

    @SuppressWarnings("unchecked")
    static <T> WatchEvent<T> cast(WatchEvent<?> event) {
        return (WatchEvent<T>)event;
    }

	/**
     * Creates a WatchService and registers the given directory
     */
    public RealtimeReadAnalyzer(Path dir) throws IOException {
        this.watcher = FileSystems.getDefault().newWatchService();
        this.keys = new HashMap<WatchKey,Path>();
        register(dir);

        // enable trace after initial registration
        this.trace = true;
    }

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		if(args.length == 1){
			try {
				RealtimeReadAnalyzer analyser = new RealtimeReadAnalyzer(Paths.get(args[0]));
				//analyser.watchForFilesAndAnalyse();
				analyser.concatenateNewReadToAlignment();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}
	
	/**
	 * Process all events for keys queued to the watcher
	 */
	public void watchForFilesAndAnalyse() {
	    for (;;) {
	
	        // wait for key to be signalled
	        WatchKey key;
	        try {
	            key = watcher.take();
	        } catch (InterruptedException x) {
	            return;
	        }
	
	        Path dir = keys.get(key);
	        if (dir == null) {
	            System.err.println("WatchKey not recognized!!");
	            continue;
	        }
	
	        for (WatchEvent<?> event: key.pollEvents()) {
	            WatchEvent.Kind kind = event.kind();
	
	            // Context for directory entry event is the file name of entry
	            WatchEvent<Path> ev = cast(event);
	            Path name = ev.context();
	            Path child = dir.resolve(name);
	
	            if(event.kind().name().equals("ENTRY_CREATE")){
		            // print out event
		            System.out.format("%s: %s\n", event.kind().name(), child);
		            convertFast5();
		            concatenateNewReadToAlignment();
		            realign();
		            buildPhylogeny();
	            }
	
	        }
	
	        // reset key and remove from set if directory no longer accessible
	        boolean valid = key.reset();
	        if (!valid) {
	            keys.remove(key);
	
	            // all directories are inaccessible
	            if (keys.isEmpty()) {
	                break;
	            }
	        }
	    }
	}

	/**
	 * build a phylogeny from output.fa.stops.removed.phy with raxml
	 */
	private void buildPhylogeny() {
		// TODO Auto-generated method stub
		long seed = System.currentTimeMillis();
		File alignment = new File("/Users/joeparker/Documents/all_work/programming/metrichor/analyses/output.fa.stops.removed.phy");
		File workDir = new File("~/Documents/all_work/programming/metrichor/analyses/");
		String runID = "realtme"+seed;
		RAxMLAnalysisSGE ra = new RAxMLAnalysisSGE(alignment, workDir, null, runID, RAxMLAnalysisSGE.NTmodelOptions.GTRCAT, RAxMLAnalysisSGE.algorithmOptions.e);
		ra.setTreeConstraint(false);
		ra.setMultifuricatingConstraint(false);
		ra.setNoStartTree(true);
		ra.setBinaryDir(new File("/Applications/Phylogenetics/RAxML/RAxML-7.2.8-ALPHA/raxmlHPC"));
	//	ra.setWorkingDir(this.workDir);
		//ra.RunAnalysis();
		new VerboseSystemCommand("/Applications/Phylogenetics/RAxML/RAxML-7.2.8-ALPHA/raxmlHPC -m GTRCAT -n "+runID+" -s /Users/joeparker/Documents/all_work/programming/metrichor/analyses/output.fa.stops.removed.phy");
	}

	/**
	 * align input.fa with muscle  (to output.fa) and convert to .phy with PrepareFilesForPaml
	 */
	private void realign() {
		// TODO Auto-generated method stub
		// muscle -in analyses/input.fa -out analyses/output.fa
		new VerboseSystemCommand("muscle -in /Users/joeparker/Documents/all_work/programming/metrichor/analyses/input.fa -out /Users/joeparker/Documents/all_work/programming/metrichor/analyses/output.fa");
		// java -jar /Users/joeparker/Documents/all_work/convergence_pipeline_alphas/PrepareFilesForPaml.jar analyses/output.fa
		new VerboseSystemCommand("java -jar /Users/joeparker/Documents/all_work/convergence_pipeline_alphas/PrepareFilesForPaml.jar /Users/joeparker/Documents/all_work/programming/metrichor/analyses/output.fa");
		//AlignedSequenceRepresentation asr = new AlignedSequenceRepresentation();
	}

	/**
	 * read in output.fa and also the new read from FQR, concatenate them, and write out to input.fa
	 */
	private void concatenateNewReadToAlignment() {
		// TODO Auto-generated method stub
		// FQR
		// CapitalizedFileReader
		// BasicFileWriter
		//
        FastqSimpleReader fqr = new FastqSimpleReader(new File("/Users/joeparker/Documents/all_work/programming/metrichor/converted_reads/out.fasta"));
        System.out.println(fqr.getParsedDataAsString());
        File outputFile = new File("/Users/joeparker/Documents/all_work/programming/metrichor/analyses/input.fa");
        new BasicFileWriter(outputFile,fqr.getParsedDataAsString());
		AlignedSequenceRepresentation asr = new AlignedSequenceRepresentation();
		try {
			asr.loadSequences(fqr.getParsedData(), false);
			asr.writePhylipFile(outputFile, false);
		} catch (TaxaLimitException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}

	private void convertFast5(){
		new VerboseSystemCommand("printenv");		
		new VerboseSystemCommand("/Users/joeparker/Documents/all_work/programming/metrichor/runme.sh /Users/joeparker/Documents/all_work/programming/metrichor/converted_reads/out.fasta");		
	}
	
	/**
	 * Register the given directory with the WatchService
	 */
	private void register(Path dir) throws IOException {
	    WatchKey key = dir.register(watcher, ENTRY_CREATE, ENTRY_DELETE, ENTRY_MODIFY);
	    if (trace) {
	        Path prev = keys.get(key);
	        if (prev == null) {
	            System.out.format("register: %s\n", dir);
	        } else {
	            if (!dir.equals(prev)) {
	                System.out.format("update: %s -> %s\n", prev, dir);
	            }
	        }
	    }
	    keys.put(key, dir);
	}

	/**
	 * Register the given directory, and all its sub-directories, with the
	 * WatchService.
	 */
	private void registerAll(final Path start) throws IOException {
	    // register directory and sub-directories
	    Files.walkFileTree(start, new SimpleFileVisitor<Path>() {
	        @Override
	        public FileVisitResult preVisitDirectory(Path dir, BasicFileAttributes attrs)
	            throws IOException
	        {
	            register(dir);
	            return FileVisitResult.CONTINUE;
	        }
	    });
	}

}
