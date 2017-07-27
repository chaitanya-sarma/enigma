package org.icrisat.genomicSelection.util;

import javax.ws.rs.client.Client;
import javax.ws.rs.client.ClientBuilder;
import javax.ws.rs.client.Entity;
import javax.ws.rs.client.WebTarget;
import javax.ws.rs.core.Form;
import javax.ws.rs.core.Response;

public class UtilWebService {

	private static Client client = ClientBuilder.newClient();
	private static WebTarget webTarget;
	private static Response response;

	public static Response authenticate(String baseURI, String userName, String password) {
		webTarget = client.target(baseURI);
		Form formData = new Form();
		formData.param("username", userName);
		formData.param("password", password);
		response = webTarget.request().post(Entity.form(formData));
		return response;

	}
}
