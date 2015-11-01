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
	private RealtimeReadAnalyzerMVC GUI;

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
        GUI = new RealtimeReadAnalyzerMVC();
    }

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		if(args.length == 1){
			try {
				RealtimeReadAnalyzer analyser = new RealtimeReadAnalyzer(Paths.get(args[0]));
				analyser.watchForFilesAndAnalyse();
				//analyser.analyse();
			} catch (IOException e) {
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
		            analyse();
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
	 * convert fast5 to fastq then fasta
	 * concatenate alignment
	 * re-align
	 * infer phylogeny
	 */
	private void analyse(){
		// extract fasta from fast5
        convertFast5();
        // may be redundant
        convertFastQ();
        // muscle
        realign();
        // RAxML
        File finishedPhylogeny = buildPhylogeny();
        // update GUI with newest phylogeny
        GUI.updatePhylogeny(finishedPhylogeny);
	}

	/**
	 * build a phylogeny from output.fa.stops.removed.phy with raxml
	 */
	private File buildPhylogeny() {
		long seed = System.currentTimeMillis();
		String runID = "realtme"+seed;
		new VerboseSystemCommand("/Applications/Phylogenetics/RAxML/RAxML-7.2.8-ALPHA/raxmlHPC -m GTRCAT -n "+runID+" -s /Users/joeparker/Documents/all_work/programming/metrichor/analyses/aligned.phy");
		return new File("/Users/joeparker/Documents/all_work/programming/metrichor/RAxML_bestTree."+runID);
	}

	/**
	 * align input.fa with muscle  (to output.fa) and convert to .phy with PrepareFilesForPaml
	 */
	private void realign() {
		// muscle -in analyses/input.fa -out analyses/output.fa
		new VerboseSystemCommand("muscle -in /Users/joeparker/Documents/all_work/programming/metrichor/analyses/input.fa -out /Users/joeparker/Documents/all_work/programming/metrichor/analyses/aligned.fa -maxiters 1 -diags");
		// convert output to phylip
		AlignedSequenceRepresentation asr = new AlignedSequenceRepresentation();
		try {
			asr.loadSequences(new File("/Users/joeparker/Documents/all_work/programming/metrichor/analyses/aligned.fa"), false);
			asr.writePhylipFile(new File("/Users/joeparker/Documents/all_work/programming/metrichor/analyses/aligned.phy"), false);
		} catch (TaxaLimitException e) {
			e.printStackTrace();
		}
	}

	/**
	 * read in new read from FQR, concatenate them, and write out to input.fa
	 */
	private void convertFastQ() {
       FastqSimpleReader fqr = new FastqSimpleReader(new File("/Users/joeparker/Documents/all_work/programming/metrichor/converted_reads/out.fastq"));
        System.out.println(fqr.getParsedDataAsString());
        //
        File outputFile = new File("/Users/joeparker/Documents/all_work/programming/metrichor/analyses/input.fa");
        new BasicFileWriter(outputFile,fqr.getParsedDataAsString());
	}

	/* runs the mahesh script which takes all fast5 in a dir and concatenates them to a single fasta file */
	private void convertFast5(){
		new VerboseSystemCommand("printenv");		
		new VerboseSystemCommand("/Users/joeparker/Documents/all_work/programming/metrichor/runme.sh /Users/joeparker/Documents/all_work/programming/metrichor/converted_reads/out.fastq");		
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
