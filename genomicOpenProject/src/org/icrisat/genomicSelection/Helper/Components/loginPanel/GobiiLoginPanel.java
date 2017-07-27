package org.icrisat.genomicSelection.helper.components.loginPanel;

import java.awt.Frame;
import java.awt.event.ActionEvent;

import javax.swing.JOptionPane;

import org.icrisat.genomicSelection.helper.components.gobiiPanels.GobiiMainPanel;

public class GobiiLoginPanel extends LoginPanel {

	private Frame parent;
	private GobiiMainPanel mainPanel;
	public GobiiLoginPanel(Frame parent) {
		super(parent, true, "Gobii Login");
		this.parent = parent;
		mainPanel = new GobiiMainPanel(parent);
	}

	@Override
	public void actionPerformed(ActionEvent e) {
		validateUserFields();
		JOptionPane.showMessageDialog(parent, "\tSuccessful login.\n \tWelcome to GOBII.");
		this.setVisible(false);
	}
}
