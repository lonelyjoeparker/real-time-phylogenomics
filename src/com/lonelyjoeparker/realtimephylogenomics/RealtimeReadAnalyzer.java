/**
 * 
 */
package com.lonelyjoeparker.realtimephylogenomics;

import static java.nio.file.LinkOption.NOFOLLOW_LINKS;
import static java.nio.file.StandardWatchEventKinds.ENTRY_CREATE;
import static java.nio.file.StandardWatchEventKinds.ENTRY_DELETE;
import static java.nio.file.StandardWatchEventKinds.ENTRY_MODIFY;
import static java.nio.file.StandardWatchEventKinds.OVERFLOW;

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
				analyser.watchForFilesAndAnalyse();
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
	
	            // print out event
	            System.out.format("%s: %s\n", event.kind().name(), child);
	            concatenateNewReadToAlignment();
	            realign();
	            buildPhylogeny();
	
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
		new VerboseSystemCommand("echo 'raxml -m GTRCAT etc'");
	}

	/**
	 * align input.fa with muscle  (to output.fa) and convert to .phy with PrepareFilesForPaml
	 */
	private void realign() {
		// TODO Auto-generated method stub
		new VerboseSystemCommand("echo 'muscle -in -out etc'");		
	}

	/**
	 * read in output.fa and also the new read from FQR, concatenate them, and write out to input.fa
	 */
	private void concatenateNewReadToAlignment() {
		// TODO Auto-generated method stub
		// FQR
		// CapitalizedFileReader
		// BasicFileWriter
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
