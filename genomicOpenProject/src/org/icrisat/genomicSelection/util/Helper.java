package org.icrisat.genomicSelection.util;

import java.awt.Image;
import java.net.URL;

import javax.swing.ImageIcon;

import org.icrisat.genomicSelection.OpenDial;

public class Helper {
	public static ImageIcon createIcon(String name) {
		URL url = OpenDial.class.getResource("./image/"+name);
		if (url == null) {
			url = OpenDial.class.getResource("./image/"+"default.png");
		}
		ImageIcon icon = new ImageIcon(new ImageIcon(url).getImage().getScaledInstance(20, 20, Image.SCALE_DEFAULT));
		return icon;
	}
	
	
	
}
