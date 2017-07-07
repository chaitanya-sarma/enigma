package org.icrisat.genomicSelection.Helper;

/**
 * Model Schema
 *  {
    "id": "",
    "uniqueID": "",
    "name": "",
    "createdBy": "",
    "members": [
      ""
    ],
    "crop": "",
    "startDate": ""
  }
 */



import java.util.ArrayList;

public class ProgramList {

	private String id;
	private String uniqueId;
	private String name;
	private String createdBy;
	private ArrayList<String> members;
	private String crop;
	private String startDate;
	public ProgramList() {
		super();
		this.id = "";
		this.uniqueId = "";
		this.name = "";
		this.createdBy = "";
		this.members = new ArrayList<String>();
		this.crop = "";
		this.startDate = "";
	}
	public String getId() {
		return id;
	}
	public void setId(String id) {
		this.id = id;
	}
	public String getUniqueId() {
		return uniqueId;
	}
	public void setUniqueId(String uniqueId) {
		this.uniqueId = uniqueId;
	}
	public String getName() {
		return name;
	}
	public void setName(String name) {
		this.name = name;
	}
	public String getCreatedBy() {
		return createdBy;
	}
	public void setCreatedBy(String createdBy) {
		this.createdBy = createdBy;
	}
	public ArrayList<String> getMembers() {
		return members;
	}
	public void setMembers(String members) {
		//this.members = members;
	}
	public String getCrop() {
		return crop;
	}
	public void setCrop(String crop) {
		this.crop = crop;
	}
	public String getStartDate() {
		return startDate;
	}
	public void setStartDate(String startDate) {
		this.startDate = startDate;
	}
	
}
