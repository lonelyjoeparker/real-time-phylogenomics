package com.lonelyjoeparker.realtimephylogenomics;

import java.io.File;

import javax.swing.JComponent;
import javax.swing.JFrame;
import javax.swing.JTabbedPane;

import uk.ac.qmul.sbcs.evolution.convergence.gui.controllers.PhylogeniesController;
import uk.ac.qmul.sbcs.evolution.convergence.gui.models.PhylogeniesModel;
import uk.ac.qmul.sbcs.evolution.convergence.gui.views.PhylogeniesView;

/**
 * Small class containing model/view/controller for a wireframe real-time-read analyser GUI
 * @author <a href="http://github.com/lonelyjoeparker">@lonelyjoeparker</a>
 * @since 1 Nov 2015
 * @version 0.1
 */
public class RealtimeReadAnalyzerMVC {

	PhylogeniesModel model;
	PhylogeniesView view;
	PhylogeniesController controller;
	RealtimeReadAnalyserFrame frame;
	/**
	 * Default no-arg constructor
	 */
	public RealtimeReadAnalyzerMVC(){
		model = new PhylogeniesModel();
		view = new PhylogeniesView();
		controller = new PhylogeniesController(model,view);
		frame = new RealtimeReadAnalyserFrame(controller);
	}
			
	private class RealtimeReadAnalyserFrame extends JFrame{
		/**
		 * 
		 */
		private static final long serialVersionUID = 8534580046578443687L;
		JTabbedPane mainTabPane;
		/**
		 * Default no-arg, deprecated
		 */
		@Deprecated
		RealtimeReadAnalyserFrame(){}

		public RealtimeReadAnalyserFrame(PhylogeniesController controller) {
			// TODO Auto-generated constructor stub
			super("Real-time reads analysed into phylogenies");
			mainTabPane = new JTabbedPane();
			mainTabPane.setTabLayoutPolicy(JTabbedPane.SCROLL_TAB_LAYOUT);
			mainTabPane.setOpaque(true);
			add(mainTabPane);
			addTab(controller.getView(),"Phylogenies");
			setDefaultCloseOperation(EXIT_ON_CLOSE);
			setSize(1080,960);
			setVisible(true);
		}

		/**
		 * Method to add tabs to main tab pane
		 * @param comp
		 * @param tabName
		 */
		public void addTab(JComponent comp, String tabName){
			mainTabPane.add(comp, tabName);
		}
	}

	public void updatePhylogeny(File finishedPhylogeny) {
		// TODO Auto-generated method stub
		controller.getModel().addPhylogenyRow(finishedPhylogeny);
		controller.forceUpdateViewWithLastModelRow();

	}
}
