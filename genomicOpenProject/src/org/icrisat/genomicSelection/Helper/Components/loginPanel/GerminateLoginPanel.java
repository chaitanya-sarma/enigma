package org.icrisat.genomicSelection.helper.components.loginPanel;

import java.awt.Frame;
import java.awt.event.ActionEvent;
import java.awt.event.ComponentEvent;

public class GerminateLoginPanel extends LoginPanel {

	private Frame parent;

	public GerminateLoginPanel(Frame parent) {
		super(parent, true, "Germinate Login");
		this.parent = parent;
	}

	@Override
	public void actionPerformed(ActionEvent e) {
		System.out.println("Welcome to Germinate Login");
		String password = String.valueOf(passwordField.getPassword());
		System.out.println("UserName :" + usernameField.getText() + "\t password  :" + password);

	}
}
