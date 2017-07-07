package org.icrisat.genomicSelection.Helper;


public class AuthObject {
	String token;
	String expires;
	
	public AuthObject() {
		this.setToken("");
		this.setExpires("");
	}
	
	public String getToken() {
		return token;
	}
	public void setToken(String token) {
		this.token = token;
	}
	public String getExpires() {
		return expires;
	}
	public void setExpires(String expires) {
		this.expires = expires;
	}
	
}
