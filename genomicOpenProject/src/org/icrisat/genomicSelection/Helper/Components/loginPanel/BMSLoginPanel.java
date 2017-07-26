package org.icrisat.genomicSelection.helper.components.loginPanel;

import java.awt.Frame;
import java.awt.event.ActionEvent;
import java.awt.event.ComponentEvent;
import java.awt.event.HierarchyEvent;

public class BMSLoginPanel extends LoginPanel {

	private Frame parent;

	public BMSLoginPanel(Frame parent) {
		super(parent, true, "BMS Login");
		this.parent = parent;
	}

	@Override
	public void actionPerformed(ActionEvent e) {
		System.out.println("Welcome to BMS Login");
		String password = String.valueOf(passwordField.getPassword());
		System.out.println("UserName :" + usernameField.getText() + "\t password  :" + password);
	}
}
