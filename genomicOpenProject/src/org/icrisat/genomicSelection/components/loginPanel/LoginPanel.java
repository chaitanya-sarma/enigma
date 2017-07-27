package org.icrisat.genomicSelection.components.loginPanel;

import java.awt.Dimension;
import java.awt.Font;
import java.awt.Frame;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionListener;

import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPasswordField;
import javax.swing.JTextField;

import org.icrisat.genomicSelection.util.*;

public abstract class LoginPanel extends JDialog implements ActionListener{

	private JLabel usernameLbl, passwordLbl;
	protected JPasswordField passwordField;
	protected JTextField usernameField;
	private JButton connect;
	protected String url, userName, password;
	private Frame parent ;
	public LoginPanel(Frame parent, boolean modal, String title) {
		super(parent, modal);
		this.parent = parent;
		setSize(new Dimension(400, 200));
		setLocationRelativeTo(parent);
		setTitle(title);
		setLayout(new GridBagLayout());
		initialize();
	}

	public void setURL(String url){
		this.url = url;
	}
	
	public void clearURL(){
		this.url = "";
	}
	
	public void clearPasswordField(){
		passwordField.setText("");
	}
	
	private void initialize() {
		usernameLbl = new JLabel("USERNAME");
		usernameLbl.setIcon(Helper.createIcon("users-icon.png"));
		usernameLbl.setFont(new Font("DejaVu Sans", Font.BOLD, 15));
		usernameLbl.setHorizontalTextPosition(JLabel.LEFT);
		usernameLbl.setIconTextGap(120);

		usernameField = new JTextField(20);

		passwordLbl = new JLabel("PASSWORD");
		passwordLbl.setIcon(Helper.createIcon("password.png"));
		passwordLbl.setFont(new Font("DejaVu Sans", Font.BOLD, 15));
		passwordLbl.setHorizontalTextPosition(JLabel.LEFT);
		passwordLbl.setIconTextGap(110);

		passwordField = new JPasswordField(20);
		connect = new JButton("Login Now");
		connect.setFont(new Font("DejaVu Sans", Font.BOLD, 14));

		GridBagConstraints gc = new GridBagConstraints();

		gc.gridx = 0;
		gc.gridy = 0;
		gc.anchor = GridBagConstraints.FIRST_LINE_START;
		gc.insets = new Insets(0, 0, 10, 0);
		add(usernameLbl, gc);

		gc.gridx = 0;
		gc.gridy = 1;
		gc.insets = new Insets(0, 0, 10, 0);
		add(usernameField, gc);

		gc.gridx = 0;
		gc.gridy = 2;
		gc.insets = new Insets(0, 0, 10, 0);
		gc.anchor = GridBagConstraints.FIRST_LINE_START;
		add(passwordLbl, gc);

		gc.gridx = 0;
		gc.gridy = 3;
		gc.insets = new Insets(0, 0, 10, 0);
		add(passwordField, gc);

		gc.gridx = 0;
		gc.gridy = 5;
		gc.insets = new Insets(0, 0, 10, 0);
		gc.anchor = GridBagConstraints.LAST_LINE_END;
		add(connect, gc);

		connect.addActionListener(this);

	}

	// Validates username and password.
	protected void validateUserFields(){
		userName = usernameField.getText().trim();
		password = String.valueOf(passwordField.getPassword());
		password = password.trim();
		if (userName == null || password == null || userName.isEmpty() || password.isEmpty()) {
			JOptionPane.showMessageDialog(parent, "Please check the entered Username and Password.");
		}
	}
}
