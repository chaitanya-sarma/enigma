/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.icrisat.genomicSelection.genoTypeDB;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import javax.swing.JDialog;

import org.icrisat.genomicSelection.helper.components.urlPanel.*;

import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;

/**
 *
 * @author Chaitanya Creating GUI dialog box
 */
public class ConnectToGenotypeDB extends JDialog implements ActionListener {

	private URLPanel gobiiURLPanel, germinateURLPanel, g4rURLPanel;

	public ConnectToGenotypeDB(java.awt.Frame parent, boolean modal) {
		// super(parent, "Connect to Phenotype Databases", modal);
		super(parent, modal);
		setSize(new Dimension(500, 250));
		setLocationRelativeTo(parent);
		gobiiURLPanel = new GobiiURLPanel(parent, "GOBII");
		germinateURLPanel = new GerminateURLPanel(parent, "Germinate");
		g4rURLPanel = new G4rURLPanel(parent, "G4R");
		setLayout(new GridBagLayout());
		GridBagConstraints gc = new GridBagConstraints();

		gc.weightx = 1;
		gc.weighty = 1;

		gc.gridx = 0;
		gc.gridy = 0;
		add(gobiiURLPanel, gc);

		gc.gridx = 0;
		gc.gridy = 1;
		add(germinateURLPanel, gc);

		gc.gridx = 0;
		gc.gridy = 2;
		add(g4rURLPanel, gc);

	}

	@Override
	public void actionPerformed(ActionEvent evt) {

	}
}
