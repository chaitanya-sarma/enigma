package org.icrisat.genomicSelection.Helper;

/**
 * Schema
 * [
  {
    "id": 0,
    "username": "",
    "firstName": "",
    "lastName": "",
    "role": "",
    "status": "",
    "email": ""
  }
]
 * @author CSarma
 *
 */

public class UserList {
	private String id, username, role, email;

	public String getId() {
		return id;
	}

	public void setId(String id) {
		this.id = id;
	}

	public String getUsername() {
		return username;
	}

	public void setUsername(String username) {
		this.username = username;
	}

	public String getRole() {
		return role;
	}

	public void setRole(String role) {
		this.role = role;
	}

	public String getEmail() {
		return email;
	}

	public void setEmail(String email) {
		this.email = email;
	}
	
}
