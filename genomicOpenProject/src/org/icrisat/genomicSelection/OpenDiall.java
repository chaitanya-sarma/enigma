/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.icrisat.genomicSelection;

import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.nio.file.Files;
import java.nio.file.Paths;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JTextField;
import javax.swing.JButton;
import javax.swing.JOptionPane;
import javax.swing.JScrollPane;
import javax.swing.JFileChooser;
import javax.swing.filechooser.FileNameExtensionFilter;
import javax.swing.GroupLayout.Alignment;
import javax.swing.GroupLayout;
import javax.swing.LayoutStyle.ComponentPlacement;

/**
 *
 * @author Mohan Creating GUI dialog box which genotype,phenotype and result
 * directory are selected Open dialog should be called first when stared project
 * to mainly get result directory and at least one genotype file and phenotype
 * file to start analysis
 */
public class OpenDiall extends JDialog implements ActionListener {

    String lastFileChoosen = "";
    JButton btnOk = new javax.swing.JButton();
    java.awt.Frame frame;
    //stare region
    /*
     * adding the action listener for buttons 
     */
    public OpenDiall(java.awt.Frame parent, boolean modal) {
        super(parent, modal);
        frame = parent;
        initComponents();
        setResizable(false);
        String resultDirPath = Constant.resultdirectory; //getting the result directory path

        //intilizing action lisener of buttons 
        btnOk.addActionListener(this);
    }

    @Override
    @SuppressWarnings("empty-statement")
    public void actionPerformed(ActionEvent evt) {
        Object source = evt.getSource(); //gettig the source (name of the button)

        if (source == btnOk) {
        	System.out.println("Ckiclsked ok");
        }
    }
//end region constructor 

    //start region 
    //initilzing all the GUI components on the dialog box    
    private void initComponents() {
        setDefaultCloseOperation(javax.swing.WindowConstants.DISPOSE_ON_CLOSE);
        setTitle("Load files");
        btnOk.setFont(new java.awt.Font("DejaVu Sans", 0, 15)); // NOI18N
        btnOk.setText("Ok");
        getContentPane().add(btnOk);
     

        pack();
    }
}
