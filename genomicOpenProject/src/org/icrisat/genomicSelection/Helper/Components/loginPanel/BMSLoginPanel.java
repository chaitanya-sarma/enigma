package org.icrisat.genomicSelection.helper.components.loginPanel;

import java.awt.Frame;
import java.awt.event.ActionEvent;

import javax.swing.JOptionPane;
import javax.ws.rs.core.Response;

import org.apache.derby.tools.sysinfo;
import org.icrisat.genomicSelection.helper.AuthObject;
import org.icrisat.genomicSelection.helper.UtilWebService;

import com.google.gson.JsonObject;
import com.google.gson.JsonParser;

public class BMSLoginPanel extends LoginPanel {

	private Frame parent;
	private AuthObject auth ;
	public BMSLoginPanel(Frame parent) {
		super(parent, true, "BMS Login");
		this.parent = parent;
	}

	@Override
	public void actionPerformed(ActionEvent e) {
		validateUserFields();
		Response response = null;
		JsonParser parser = new JsonParser();
		try {
			url = url + "/authenticate";
			response = UtilWebService.authenticate(url, userName, password);
			if (response == null)
				JOptionPane.showMessageDialog(parent, "Please check the entered URL.");
			else{
				if(response.getStatus() != 200){
					JOptionPane.showMessageDialog(parent, "Please check the entered login credentials.");
				}else{
					JOptionPane.showMessageDialog(parent, "\tSuccessful login.\n \tWelcome to BMS.");
					auth = new AuthObject();
					String res = response.readEntity(String.class);
					JsonObject json = (JsonObject) parser.parse(res);
					auth.setToken(json.get("token").getAsString());
					auth.setExpires(json.get("expires").getAsString());
					this.setVisible(false);
				}
			}
		} catch (Exception e1) {
			JOptionPane.showMessageDialog(parent, "Please check the entered URL and try again.");
		}
		//System.out.println(auth.getToken());
	}
}
