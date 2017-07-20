/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.icrisat.genomicSelection.PhenoTypeDB;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import javax.swing.JDialog;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.Frame;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;

/**
 *
 * @author Chaitanya Creating GUI dialog box
 */
public class ConnectToPhenotypeDB extends JDialog implements ActionListener {

	private URLPanel bmsURLPanel, b4rURLPanel;

	public ConnectToPhenotypeDB(java.awt.Frame parent, boolean modal) {
		//super(parent, "Connect to Phenotype Databases", modal);
		super(parent, modal);
		setSize(new Dimension(500, 200));
		setLocationRelativeTo(parent);
		bmsURLPanel = new URLPanel("BMS URL");
		b4rURLPanel = new URLPanel("B4R URL");
		setLayout(new GridBagLayout());
		GridBagConstraints gc = new GridBagConstraints();

		gc.weightx = 1;
		gc.weighty = 1;

		gc.gridx = 0;
		gc.gridy = 0;
		add(bmsURLPanel, gc);

		gc.weightx = 1;
		gc.weighty = 1.5;

		gc.gridx = 0;
		gc.gridy = 1;
		add(b4rURLPanel, gc);
	}

	@Override
	public void actionPerformed(ActionEvent evt) {
		Object source = evt.getSource(); // gettig the source (name of the
											// button)
	}
}
