package org.icrisat.genomicSelection.Helper;

import java.util.ArrayList;

import javax.ws.rs.client.Client;
import javax.ws.rs.client.ClientBuilder;
import javax.ws.rs.client.Entity;
import javax.ws.rs.client.WebTarget;
import javax.ws.rs.core.Form;
import javax.ws.rs.core.Response;

import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;

public class Util {

	private static Client client = ClientBuilder.newClient();
	private static AuthObject auth = new AuthObject();
	private static WebTarget webTarget;
	private static Response response;
	private static JsonParser parser = new JsonParser();
	private static final String baseURI = "http://bms.icrisat.ac.in:48080/bmsapi/";

	public static int authenticate(String userName, String password) {
		int status = 0;
		webTarget = client.target(baseURI + "authenticate");
		Form formData = new Form();
		formData.param("username", "arathore");
		formData.param("password", "abcd@1234");
		response = webTarget.request().post(Entity.form(formData));
		status = response.getStatus();
		if (status == 200) {
			String res = response.readEntity(String.class);
			JsonObject json = (JsonObject) parser.parse(res);
			auth.setToken(json.get("token").getAsString());
			auth.setExpires(json.get("expires").getAsString());
		}
		return status;

	}

	public static int getPrograms(ArrayList<ProgramList> programList) {
		int status = 0;
		webTarget = client.target(baseURI + "program/list");
		response = webTarget.request().header("X-Auth-Token", auth.getToken()).get();
		status = response.getStatus();

		if (status == 200) {
			String res = response.readEntity(String.class);
			JsonElement jsonTree = parser.parse(res);
			JsonArray jsonArr = null;
			if (jsonTree.isJsonArray()) {
				jsonArr = jsonTree.getAsJsonArray();
			}
			for (int i = 0; i < jsonArr.size(); i++) {
				JsonObject ele = jsonArr.get(i).getAsJsonObject();
				ProgramList pl = new ProgramList();
				pl.setId(ele.get("id").getAsString());
				pl.setUniqueId(ele.get("uniqueID").getAsString());
				pl.setName(ele.get("name").getAsString());
				pl.setCreatedBy(ele.get("createdBy").getAsString());
				pl.setMembers(ele.get("members").toString());
				pl.setCrop(ele.get("crop").toString());
				pl.setStartDate(ele.get("startDate").toString());
				programList.add(pl);
			}
		}
		return status;
	}

	public static int getCrops(ArrayList<String> cropList) {
		int status = 0;
		webTarget = client.target(baseURI + "crop/list");
		response = webTarget.request().header("X-Auth-Token", auth.getToken()).get();
		status = response.getStatus();
		if (status == 200) {
			String res = response.readEntity(String.class);
			JsonElement jsonTree = parser.parse(res);
			JsonArray jsonArr = null;
			if (jsonTree.isJsonArray()) {
				jsonArr = jsonTree.getAsJsonArray();
			}
			for (int i = 0; i < jsonArr.size(); i++) {
				cropList.add(jsonArr.get(i).getAsString());
			}
		}

		return status;
	}

	public static int getUsers(ArrayList<UserList> userList) {
		int status = 0;
		webTarget = client.target(baseURI + "brapi/v1/users");
		response = webTarget.request().header("X-Auth-Token", auth.getToken()).get();
		status = response.getStatus();

		if (status == 200) {
			String res = response.readEntity(String.class);
			JsonElement jsonTree = parser.parse(res);
			JsonArray jsonArr = null;
			if (jsonTree.isJsonArray()) {
				jsonArr = jsonTree.getAsJsonArray();
			}
			for (int i = 0; i < jsonArr.size(); i++) {
				JsonObject ele = jsonArr.get(i).getAsJsonObject();
				UserList pl = new UserList();
				pl.setId(ele.get("id").getAsString());
				pl.setUsername(ele.get("username").getAsString());
				pl.setRole(ele.get("role").getAsString());
				pl.setEmail(ele.get("email").getAsString());
				userList.add(pl);
			}
		}
		return status;
	}

	public static int getCropLocation(String crop, ArrayList<Location> locList) {
		int status = 0;
		webTarget = client.target(baseURI + "location/" + crop + "/types");
		response = webTarget.request().header("X-Auth-Token", auth.getToken()).get();
		status = response.getStatus();
		
		if (status == 200) {
			String res = response.readEntity(String.class);
			System.out.println(res);
			JsonElement jsonTree = parser.parse(res);
			JsonArray jsonArr = null;
			if (jsonTree.isJsonArray()) {
				jsonArr = jsonTree.getAsJsonArray();
			}
			for (int i = 0; i < jsonArr.size(); i++) {
				JsonObject ele = jsonArr.get(i).getAsJsonObject();
				Location loc = new Location();
				loc.setName(ele.get("name").toString());
				loc.setId(ele.get("id").toString());
				loc.setDescription(ele.get("description").toString());
				locList.add(loc);
			}
		}
		return status;
	}

	public static int getTrialSummary(String crop, ArrayList<TrialSummary> tsList) {
		int status = 0;
		webTarget = client.target(baseURI + crop + "/brapi/v1/trials");
		response = webTarget.request().header("X-Auth-Token", auth.getToken()).get();
		status = response.getStatus();
		
		if (status == 200) {
			String res = response.readEntity(String.class);
			JsonObject json = (JsonObject) parser.parse(res);
			JsonArray result = json.get("result").getAsJsonObject().get("data").getAsJsonArray();
			for (int i = 0; i < result.size(); i++) {
				TrialSummary ts = new TrialSummary();
				JsonObject ele = result.get(i).getAsJsonObject();
				ts.setActive(ele.get("active").toString());
				ts.setEndDate(ele.get("endDate").toString());
				ts.setStartDate(ele.get("startDate").toString());
				ts.setTrialDbId(ele.get("trialDbId").toString());
				ts.setTrialName(ele.get("trialName").toString());
				tsList.add(ts);
			}
		}
		return status;
	}
/**	
	public static void main(String args[]) {
		authenticate("", "");
		getCropLocation("bean", new ArrayList<Location>());
	}
*/}
