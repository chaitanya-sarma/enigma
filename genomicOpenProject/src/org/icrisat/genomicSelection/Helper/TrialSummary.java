package org.icrisat.genomicSelection.helper;
/**
 *  "data": [
      {
        "trialDbId": 3505,
        "trialName": "Bean Nursery template",
        "startDate": "20150507",
        "endDate": "20151231",
        "active": true,
        "studies": [
          
        ],
        "additionalInfo": {
          "STUDY_OBJECTIVE": "Nursery template Bean"
        }
      }
 */
import java.util.List;

public class TrialSummary {
	private String trialDbId;
	private String trialName;
	private String startDate;
	private String endDate;
	private String active;
	// TODO: Below two fields are modified as required. Declared just for fun.
	List<String> studies;
	List<String> additionalInfo;
	public String getTrialDbId() {
		return trialDbId;
	}
	public void setTrialDbId(String trialDbId) {
		this.trialDbId = trialDbId;
	}
	public String getTrialName() {
		return trialName;
	}
	public void setTrialName(String trialName) {
		this.trialName = trialName;
	}
	public String getStartDate() {
		return startDate;
	}
	public void setStartDate(String startDate) {
		this.startDate = startDate;
	}
	public String getEndDate() {
		return endDate;
	}
	public void setEndDate(String endDate) {
		this.endDate = endDate;
	}
	public String getActive() {
		return active;
	}
	public void setActive(String active) {
		this.active = active;
	}

}
