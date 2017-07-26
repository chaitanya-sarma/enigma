package org.icrisat.genomicSelection.helper.components.loginPanel;

import java.awt.Frame;
import java.awt.event.ActionEvent;
import java.awt.event.ComponentEvent;

public class G4rLoginPanel extends LoginPanel {

	private Frame parent;

	public G4rLoginPanel(Frame parent) {
		super(parent, true, "G4R Login");
		this.parent = parent;
	}

	@Override
	public void actionPerformed(ActionEvent e) {
		System.out.println("Welcome to G4r Login");
		String password = String.valueOf(passwordField.getPassword());
		System.out.println("UserName :" + usernameField.getText() + "\t password  :" + password);

	}

}
