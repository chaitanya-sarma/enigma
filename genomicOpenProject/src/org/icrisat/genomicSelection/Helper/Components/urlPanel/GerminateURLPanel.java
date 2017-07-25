package org.icrisat.genomicSelection.helper.components.urlPanel;

import java.awt.Frame;
import java.awt.event.ActionEvent;

import org.icrisat.genomicSelection.helper.components.loginPanel.GerminateLoginPanel;

public class GerminateURLPanel extends URLPanel {

	private Frame parent;
	private GerminateLoginPanel germinateLoginPanel;
	public GerminateURLPanel(Frame parent, String title) {
		super(title);
		this.parent = parent;
		germinateLoginPanel = new GerminateLoginPanel(parent);
	}

	@Override
	public void actionPerformed(ActionEvent e) {
		System.out.println("Hello Germinate");
		germinateLoginPanel.clearPasswordField();
		germinateLoginPanel.setVisible(true);
	}

}
