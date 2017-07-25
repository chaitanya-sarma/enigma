package org.icrisat.genomicSelection.helper.components.urlPanel;

import java.awt.Frame;
import java.awt.event.ActionEvent;

import org.icrisat.genomicSelection.helper.components.loginPanel.BMSLoginPanel;

public class BmsURLPanel extends URLPanel {

	private Frame parent;
	private BMSLoginPanel bmsLoginPanel;
	public BmsURLPanel(Frame parent, String title) {
		super(title);
		this.parent = parent;
		bmsLoginPanel = new BMSLoginPanel(parent);
	}

	@Override
	public void actionPerformed(ActionEvent e) {
		System.out.println("Hello BMS");
		bmsLoginPanel.clearPasswordField();
		bmsLoginPanel.setVisible(true);
	}

}
