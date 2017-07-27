package org.icrisat.genomicSelection.components.gobiiPanels;

import java.awt.Frame;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JDialog;
import javax.swing.JOptionPane;

public class GobiiMainPanel extends JDialog implements ActionListener{

	private Frame parent;
	
	public GobiiMainPanel(Frame parent) {
		this.parent = parent;
		System.out.println("Gobii Main Panel");
	}

	@Override
	public void actionPerformed(ActionEvent e) {
	}
}
