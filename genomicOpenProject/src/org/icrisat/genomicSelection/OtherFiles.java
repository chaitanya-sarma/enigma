/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.icrisat.genomicSelection;

import java.awt.Dimension;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.JFileChooser;
import javax.swing.JOptionPane;
import javax.swing.filechooser.FileNameExtensionFilter;
import org.apache.commons.io.FilenameUtils;

/**
 *
 * @author Mohan
 */
public class OtherFiles extends javax.swing.JDialog {

    /**
     * Dialog box to get the population,pedigree and relationship file in
     * project
     */
    String lastFileChoosen = "";
    //start region
    /*
     * Calling the GUI initiizing method 
     * Adding the action listner for all the buttons
     * ok button action open only the path choosen for the file types
     */

    public OtherFiles(java.awt.Frame parent, boolean modal) {
        super(parent, modal);
        Dimension dim = Toolkit.getDefaultToolkit().getScreenSize();
        this.setLocation(dim.width / 4 - this.getSize().width / 4, dim.height / 3 - this.getSize().height / 3);
        initComponents();
        txtPopulation.setEditable(false);

        btnPopulation.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                JFileChooser chooser;
                lastFileChoosen = Constant.browsepath;
                if (lastFileChoosen.isEmpty()) {
                    chooser = new JFileChooser();
                } else {
                    chooser = new JFileChooser(lastFileChoosen);
                }
                chooser.removeChoosableFileFilter(chooser.getFileFilter());
                chooser.setDialogTitle("Select a population  file");
                FileNameExtensionFilter filter = new FileNameExtensionFilter("CSV", "csv");
                chooser.setFileFilter(filter);
                int r = chooser.showOpenDialog(rootPane);
                if (r == JFileChooser.APPROVE_OPTION) {
                    txtPopulation.setText(chooser.getSelectedFile().getPath());
                    Constant.browsepath = chooser.getSelectedFile().getAbsolutePath();
                    File fileValidation = new File(txtPopulation.getText());
                    boolean poputlationExists = fileValidation.exists();
                    Constant.browsepath = chooser.getSelectedFile().getAbsolutePath();
                    if (poputlationExists == false) {
                        JOptionPane.showMessageDialog(rootPane, "Selected file doesnot exists \nPlease select valid covariate file");
                        txtPopulation.setText("");
                    }
                }
                if (r == JFileChooser.CANCEL_OPTION) {
                    //donothing
                }
            }
        });

        btnHapmap.addActionListener(new ActionListener() {

            @Override
            public void actionPerformed(ActionEvent e) {
                JFileChooser chooser;
                lastFileChoosen = Constant.browsepath;
                if (lastFileChoosen.isEmpty()) {
                    chooser = new JFileChooser();
                } else {
                    chooser = new JFileChooser(lastFileChoosen);
                }
           
                chooser.setDialogTitle("Select a haplotype  file");
                FileNameExtensionFilter filter = new FileNameExtensionFilter("Text", "txt");
                chooser.setFileFilter(filter);
                int r = chooser.showOpenDialog(rootPane);
                if (r == JFileChooser.APPROVE_OPTION) {
                    txtHapmap.setText(chooser.getSelectedFile().getPath());
                    Constant.browsepath = chooser.getSelectedFile().getAbsolutePath();
                    File fileValidation = new File(txtHapmap.getText());
                    boolean hapExists = fileValidation.exists();
                    if (hapExists == false) {
                        JOptionPane.showMessageDialog(rootPane, "Selected file doesnot exists \nPlease select valid hapmap file");
                        txtHapmap.setText("");
                    }
                }
                if (r == JFileChooser.CANCEL_OPTION) {
                    //donothing
                }
            }

        });
        btnOk.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {

                if (!txtPopulation.getText().trim().equals("")) {
                    try {
                        String populatioName = new File(txtPopulation.getText()).getName();
                        String popNameExt = FilenameUtils.removeExtension(populatioName);
                        Date date = new Date();
                        SimpleDateFormat dateFormat = new SimpleDateFormat("mmss");

                        String populationdir = Constant.resultdirectory + popNameExt + dateFormat.format(date) + ".csv";
                        String popluation = txtPopulation.getText(); //getting the result directory path

                        Constant.populationPath = populationdir;
                        CopyFile.copy(popluation, populationdir);
                    } catch (IOException ex) {
                        Logger.getLogger(GsMain.class.getName()).log(Level.SEVERE, null, ex);
                    }
                } else {
                    Constant.populationPath = null;
                }
                if (!txtHapmap.getText().trim().equals("")) {
                    String haploName = new File(txtHapmap.getText()).getName();
                    String removeExt = FilenameUtils.removeExtension(haploName);
                    String haploDir = Constant.resultdirectory + removeExt + ".csv";
                    String haplo = txtHapmap.getText();

                    GsMethods.tabtoCSv(haplo, haploDir);
                    Date date = new Date();
                    SimpleDateFormat dateFormat = new SimpleDateFormat("mmss");
                    String haploColn = Constant.resultdirectory + removeExt + dateFormat.format(date) + ".csv";
                    GsMethods.addColnaTofile(haploDir, haploColn);
                    Constant.haploPath = haploColn;
                    File f = new File(haploDir);
                    boolean delete = f.delete();

                } else {
                    Constant.haploPath = null;
                }
                if (txtPopulation.getText().trim().equals("") && txtHapmap.getText().trim().equals("")) {
                    JOptionPane.showMessageDialog(null, "No file is selected.Please select a file(s) using browse button");
                } else {
                    dispose(); //closing only this
                }

            }
        });
        btnCanel.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                dispose();
                Constant.populationPath = null;
                Constant.haploPath = null;
            }
        });

    }
//end region  constructor 

    /**
     * start region method which generate GUI for the otherFile dialog box
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        lblPopulation = new javax.swing.JLabel();
        txtPopulation = new javax.swing.JTextField();
        btnPopulation = new javax.swing.JButton();
        btnOk = new javax.swing.JButton();
        btnCanel = new javax.swing.JButton();
        txtHapmap = new javax.swing.JTextField();
        btnHapmap = new javax.swing.JButton();
        lblHapmap = new javax.swing.JLabel();

        setDefaultCloseOperation(javax.swing.WindowConstants.DISPOSE_ON_CLOSE);
        setTitle("Import");

        lblPopulation.setFont(new java.awt.Font("DejaVu Sans", 0, 15)); // NOI18N
        lblPopulation.setText("Covariate file :");

        txtPopulation.setFont(new java.awt.Font("DejaVu Sans", 0, 15)); // NOI18N
        txtPopulation.setText(" ");

        btnPopulation.setFont(new java.awt.Font("DejaVu Sans", 0, 15)); // NOI18N
        btnPopulation.setText("Browse");

        btnOk.setFont(new java.awt.Font("DejaVu Sans", 0, 15)); // NOI18N
        btnOk.setText("Ok");

        btnCanel.setFont(new java.awt.Font("DejaVu Sans", 0, 15)); // NOI18N
        btnCanel.setText("Cancel");

        txtHapmap.setFont(new java.awt.Font("DejaVu Sans", 0, 15)); // NOI18N
        txtHapmap.setText(" ");

        btnHapmap.setFont(new java.awt.Font("DejaVu Sans", 0, 15)); // NOI18N
        btnHapmap.setText("Browse");

        lblHapmap.setFont(new java.awt.Font("DejaVu Sans", 0, 15)); // NOI18N
        lblHapmap.setText("Hapmap file (TASSEL) :");

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(getContentPane());
        getContentPane().setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(lblPopulation)
                    .addComponent(lblHapmap))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                    .addComponent(txtPopulation, javax.swing.GroupLayout.DEFAULT_SIZE, 184, Short.MAX_VALUE)
                    .addComponent(txtHapmap))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                    .addComponent(btnPopulation, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(btnHapmap, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                .addContainerGap(159, Short.MAX_VALUE)
                .addComponent(btnOk, javax.swing.GroupLayout.PREFERRED_SIZE, 60, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addGap(39, 39, 39)
                .addComponent(btnCanel)
                .addGap(125, 125, 125))
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(lblPopulation)
                    .addComponent(txtPopulation, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(btnPopulation))
                .addGap(18, 18, 18)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(txtHapmap, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(btnHapmap)
                    .addComponent(lblHapmap))
                .addGap(18, 18, 18)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                    .addComponent(btnCanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(btnOk))
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );

        pack();
    }// </editor-fold>//GEN-END:initComponents
//end region 

    /**
     * main method
     *
     * @param args
     */
    public static void main(String args[]) {
        /*
         * Set the Nimbus look and feel
         */
        //<editor-fold defaultstate="collapsed" desc=" Look and feel setting code (optional) ">
        /*
         * If Nimbus (introduced in Java SE 6) is not available, stay with the
         * default look and feel. For details see
         * http://download.oracle.com/javase/tutorial/uiswing/lookandfeel/plaf.html
         */
        try {
            for (javax.swing.UIManager.LookAndFeelInfo info : javax.swing.UIManager.getInstalledLookAndFeels()) {
                if ("Nimbus".equals(info.getName())) {
                    javax.swing.UIManager.setLookAndFeel(info.getClassName());
                    break;
                }
            }
        } catch (ClassNotFoundException | InstantiationException | IllegalAccessException | javax.swing.UnsupportedLookAndFeelException ex) {
            java.util.logging.Logger.getLogger(OtherFiles.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        }
        //</editor-fold>
        
        //</editor-fold>
    }
    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JButton btnCanel;
    private javax.swing.JButton btnHapmap;
    private javax.swing.JButton btnOk;
    private javax.swing.JButton btnPopulation;
    private javax.swing.JLabel lblHapmap;
    private javax.swing.JLabel lblPopulation;
    private javax.swing.JTextField txtHapmap;
    private javax.swing.JTextField txtPopulation;
    // End of variables declaration//GEN-END:variables
}
