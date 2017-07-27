package org.icrisat.genomicSelection.components;

import java.awt.Frame;
import java.awt.event.ActionEvent;

import javax.swing.JOptionPane;

import org.icrisat.genomicSelection.components.loginPanel.LoginPanel;

public class GobiiURLPanel extends URLPanel {
	private Frame parent;

	public GobiiURLPanel(Frame parent, String title, LoginPanel loginPanel) {
		super(parent, title, loginPanel);
		System.out.println("bladddd");
		this.parent = parent;
	}

	@Override
	public void actionPerformed(ActionEvent e) {
		System.out.println(urlField.getText());
		if (urlField.getText().equals("URL")) {
			JOptionPane.showMessageDialog(parent, "Atleast enter something. Please dont leave me blank.");
		} else {
			loginPanel.setVisible(true);
		}

	}
}
