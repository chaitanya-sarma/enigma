package org.icrisat.genomicSelection.helper.components.urlPanel;

import java.awt.Frame;
import java.awt.event.ActionEvent;

import org.icrisat.genomicSelection.helper.components.loginPanel.G4rLoginPanel;

public class G4rURLPanel extends URLPanel {

	private Frame parent;
	private G4rLoginPanel g4rLoginPanel;
	public G4rURLPanel(Frame parent, String title) {
		super(title);
		this.parent = parent;
		g4rLoginPanel = new G4rLoginPanel(parent);
	}
	@Override
	public void actionPerformed(ActionEvent e) {
		System.out.println("Hello G4r");
		g4rLoginPanel.clearPasswordField();
		g4rLoginPanel.setVisible(true);
	}

}
