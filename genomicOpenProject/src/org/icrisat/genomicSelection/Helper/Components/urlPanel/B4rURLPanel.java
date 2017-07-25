package org.icrisat.genomicSelection.helper.components.urlPanel;

import java.awt.Frame;
import java.awt.event.ActionEvent;

import org.icrisat.genomicSelection.helper.components.loginPanel.B4rLoginPanel;

public class B4rURLPanel extends URLPanel {

	private Frame parent;
	private B4rLoginPanel b4rLoginPanel;
	public B4rURLPanel(Frame parent, String title) {
		super(title);
		this.parent = parent;
		b4rLoginPanel = new B4rLoginPanel(parent);
	}

	@Override
	public void actionPerformed(ActionEvent e) {
		System.out.println("HELLO B4R");
		b4rLoginPanel.clearPasswordField();
		b4rLoginPanel.setVisible(true);
	}

}
