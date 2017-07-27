package org.icrisat.genomicSelection.helper.components;

import java.awt.Dimension;
import java.awt.Font;
import java.awt.Frame;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.border.Border;
import javax.swing.border.TitledBorder;
import javax.ws.rs.core.Response;

import org.icrisat.genomicSelection.helper.UtilWebService;
import org.icrisat.genomicSelection.helper.components.loginPanel.LoginPanel;

import com.sun.xml.internal.bind.v2.runtime.reflect.Lister.Pack;

public class URLPanel extends JPanel implements ActionListener{
	protected CustomTextField urlField;
	private JButton connect;
	protected LoginPanel loginPanel;
	private Frame parent;
	public URLPanel(Frame parent, String title, LoginPanel loginPanel) {

		this.parent = parent;
		this.loginPanel = loginPanel;
		urlField = new CustomTextField(30);
		urlField.setPlaceholder("URL");
		
		connect = new JButton("Connect");
		connect.setPreferredSize(new Dimension(100, 23));
		connect.setFont(new Font("DejaVu Sans", Font.BOLD, 15));
		TitledBorder inner = new TitledBorder(title);
		inner.setTitleFont(new Font("DejaVu Sans", Font.BOLD, 20));
		Border outer = BorderFactory.createEmptyBorder(5, 5, 5, 5);
		setBorder(BorderFactory.createCompoundBorder(outer, inner));
		setLayout(new GridBagLayout());
		GridBagConstraints gc = new GridBagConstraints();
		
		gc.weightx = 1;
		gc.weighty = 1;

		gc.gridx = 0;
		gc.gridy = 0;
		gc.insets = new Insets(0, 0, 0, 5);
		gc.fill = GridBagConstraints.HORIZONTAL;
		gc.anchor = GridBagConstraints.LINE_END;
		add(urlField, gc);

		gc.weightx = 1;
		gc.weighty = 1;

		gc.gridx = 1;
		gc.gridy = 0;
		gc.insets = new Insets(0, 0, 0, 0);
		gc.fill = GridBagConstraints.NONE;
		gc.anchor = GridBagConstraints.LINE_START;
		add(connect, gc);
		
		connect.addActionListener(this);
	}
	
	@Override
	public void actionPerformed(ActionEvent e) {
		Response response = null;
		loginPanel.clearURL();
		try{
			response = UtilWebService.authenticate(urlField.getText(), "userName", "password");
			if(response == null)
				JOptionPane.showMessageDialog(parent, "Please check the entered URL.");
			else{
				// URL is good. As username and password or wrong, we got this error.
				if(response.getStatus() == 405){
					loginPanel.setURL(urlField.getText());
					loginPanel.clearPasswordField();
					loginPanel.setVisible(true);
				}else{
					JOptionPane.showMessageDialog(parent, "Please check the entered URL.");
				}
			}
			
		}catch(Exception e1){
			JOptionPane.showMessageDialog(parent, "Please check the entered URL.");
		}
		
	}
}
