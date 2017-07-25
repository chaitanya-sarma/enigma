package org.icrisat.genomicSelection.helper.components.urlPanel;

import java.awt.Frame;
import java.awt.event.ActionEvent;

import org.icrisat.genomicSelection.helper.components.loginPanel.GobiiLoginPanel;

public class GobiiURLPanel extends URLPanel {

	private Frame parent;
	private GobiiLoginPanel gobiiLoginPanel;
	public GobiiURLPanel(Frame parent, String title) {
		super(title);
		this.parent = parent;
		gobiiLoginPanel = new GobiiLoginPanel(parent);
	}
	@Override
	public void actionPerformed(ActionEvent e) {
		System.out.println("Hello Gobii");
		gobiiLoginPanel.clearPasswordField();
		gobiiLoginPanel.setVisible(true);
	}

}
