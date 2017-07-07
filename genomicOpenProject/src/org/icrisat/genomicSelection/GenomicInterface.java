/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.icrisat.genomicSelection;

import chrriis.dj.nativeswing.swtimpl.components.JWebBrowser;
import javax.swing.JPanel;
import javax.swing.JScrollPane;

/**
 *
 * @author Mohan
 * all the method needed in the gsMethod class 
 */
public interface GenomicInterface {
    //display the html file
    public JWebBrowser browser(String path); 
    
    //display pdf file 
    public JPanel pdfviewer(String path); 
    
    //display csv file
    public JPanel csvreader(String path);

    //display text file
    public JScrollPane textFileReader(String path);

    //write the content to text file
    public void writeToFile(java.util.List<String> porSavelist, String path);
}
