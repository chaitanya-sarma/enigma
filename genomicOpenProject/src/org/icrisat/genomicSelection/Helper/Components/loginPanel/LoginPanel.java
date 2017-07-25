package org.icrisat.genomicSelection.helper.components.loginPanel;

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
import javax.swing.JPasswordField;
import javax.swing.JTextField;
import javax.swing.SwingConstants;

import org.icrisat.genomicSelection.helper.Helper;

import com.sun.glass.events.WindowEvent;
import com.sun.javafx.stage.WindowHelper.WindowAccessor;

public abstract class LoginPanel extends JDialog implements ActionListener{

	private JLabel username, password;
	protected JPasswordField passwordField;
	protected JTextField usernameField;
	private JButton connect;

	public LoginPanel(Frame parent, boolean modal) {
		super(parent, modal);
		setSize(new Dimension(400, 250));
		setLocationRelativeTo(parent);
		setTitle("Login");
		setLayout(new GridBagLayout());
		initialize();
	}

	public void clearPasswordField(){
		passwordField.setText("");
	}
	
	private void initialize() {
		username = new JLabel("USERNAME");
		username.setIcon(Helper.createIcon("users-icon.png"));
		username.setFont(new Font("DejaVu Sans", Font.BOLD, 15));
		username.setHorizontalTextPosition(JLabel.LEFT);
		username.setIconTextGap(120);

		usernameField = new JTextField(20);

		password = new JLabel("PASSWORD");
		password.setIcon(Helper.createIcon("password.png"));
		password.setFont(new Font("DejaVu Sans", Font.BOLD, 15));
		password.setHorizontalTextPosition(JLabel.LEFT);
		password.setIconTextGap(110);

		passwordField = new JPasswordField(20);
		connect = new JButton("Login Now");
		connect.setFont(new Font("DejaVu Sans", Font.BOLD, 14));

		GridBagConstraints gc = new GridBagConstraints();

		gc.gridx = 0;
		gc.gridy = 0;
		gc.anchor = GridBagConstraints.FIRST_LINE_START;
		gc.insets = new Insets(0, 0, 10, 0);
		add(username, gc);

		gc.gridx = 0;
		gc.gridy = 1;
		gc.insets = new Insets(0, 0, 10, 0);
		add(usernameField, gc);

		gc.gridx = 0;
		gc.gridy = 2;
		gc.insets = new Insets(0, 0, 10, 0);
		gc.anchor = GridBagConstraints.FIRST_LINE_START;
		add(password, gc);

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
}
