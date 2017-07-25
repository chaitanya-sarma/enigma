package org.icrisat.genomicSelection.phenoTypeDB;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import javax.swing.JDialog;

import org.icrisat.genomicSelection.helper.components.urlPanel.B4rURLPanel;
import org.icrisat.genomicSelection.helper.components.urlPanel.BmsURLPanel;
import org.icrisat.genomicSelection.helper.components.urlPanel.URLPanel;

import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;

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
		bmsURLPanel = new BmsURLPanel(parent,"BMS");
		b4rURLPanel = new B4rURLPanel(parent,"B4R");
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
	public void actionPerformed(ActionEvent arg0) {
		// TODO Auto-generated method stub
		
	}
}
