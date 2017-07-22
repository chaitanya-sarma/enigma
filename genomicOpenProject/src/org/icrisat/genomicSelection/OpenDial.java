/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.icrisat.genomicSelection;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.net.URL;
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

import org.icrisat.genomicSelection.PhenoTypeDB.ConnectToPhenotypeDB;
import org.icrisat.genomicSelection.GenoTypeDB.ConnectToGenotypeDB;
import org.icrisat.genomicSelection.Helper.Helper;
import org.icrisat.genomicSelection.Helper.Util;
import javax.swing.GroupLayout.Alignment;
import javax.swing.ImageIcon;
import javax.swing.GroupLayout;
import javax.swing.LayoutStyle.ComponentPlacement;
import java.awt.Font;
import java.awt.Image;

/**
 *
 * @author Mohan Creating GUI dialog box which genotype,phenotype and result
 * directory are selected Open dialog should be called first when stared project
 * to mainly get result directory and at least one genotype file and phenotype
 * file to start analysis
 */
public class OpenDial extends JDialog implements ActionListener {

    public int validate_geno = 0, validate_pheno = 0;
    String lastFileChoosen = "";
    ClosableTabbedPane tabbedPane;
    JLabel lblGeno = new javax.swing.JLabel();
    JTextField txtGenotype;
    JButton btnGenotype = new javax.swing.JButton();
    JLabel lblPheno = new javax.swing.JLabel();
    JTextField txtPhenotype;
    JButton btnPhenotype = new javax.swing.JButton();
    JButton btnBMS = new javax.swing.JButton();
    JButton btnGOBII = new javax.swing.JButton();
    JLabel lblResultDir = new javax.swing.JLabel();
    JTextField txtResultDir;
    JButton btnResultDir = new javax.swing.JButton();
    JButton btnOk = new javax.swing.JButton();
    JButton btnCancle = new javax.swing.JButton();
    java.awt.Frame frame;
    //stare region
    /*
     * adding the action listener for buttons 
     */
    public OpenDial(java.awt.Frame parent, boolean modal) {
        super(parent, modal);
        setSize(new Dimension(900, 900));
        frame = parent;
        tabbedPane = new ClosableTabbedPane();
        txtGenotype = new javax.swing.JTextField();
        txtGenotype.setEditable(false);
        txtPhenotype = new javax.swing.JTextField();
        txtPhenotype.setEditable(false);
        txtResultDir = new javax.swing.JTextField();
        txtResultDir.setEditable(false);
        initComponents();
        setResizable(false);
        String resultDirPath = Constant.resultdirectory; //getting the result directory path
        if (!"empty".equals(resultDirPath)) {
            txtResultDir.setText(resultDirPath);//set the path allready choosen by the user
            txtResultDir.setEditable(false);
            txtResultDir.setOpaque(true); //setting the resultdir text field opaque
            txtResultDir.setMaximumSize(new Dimension(180, 0));
            txtResultDir.setFont(new java.awt.Font("DejaVu Sans", 0, 15));
            btnResultDir.setOpaque(true);//setting result button opaque
            btnResultDir.setEnabled(false); //setting the btnresult button enabled  flase
        } else { //checking result directory is empty or not 
            btnResultDir.addActionListener(this); //if empty enable action listner for btnresult
        }

        //intilizing action lisener of buttons 
        btnGenotype.addActionListener(this);
        btnPhenotype.addActionListener(this);
        btnBMS.addActionListener(this);
        btnGOBII.addActionListener(this);
        btnOk.addActionListener(this);
        btnCancle.addActionListener(this);
        Constant.genotype = null;
        Constant.phenotype = null;
    }

    @Override
    @SuppressWarnings("empty-statement")
    public void actionPerformed(ActionEvent evt) {
        Object source = evt.getSource(); //gettig the source (name of the button)

        if (source == btnOk) {
            String geno = txtGenotype.getText(); //getting the genotype filename 
            String pheno = txtPhenotype.getText(); //getting the phenotype filename 
            String result = txtResultDir.getText(); //getting the result directory path
            String geno_filename = new File(geno).getName();
            String pheno_filename = new File(pheno).getName();
            //checking the genotype,phenotype and result directory is empty
            if (txtGenotype.getText().trim().equals("") || txtPhenotype.getText().trim().equals("") || txtResultDir.getText().trim().equals("")) {
                //displaying a messagedialog box
                JOptionPane.showMessageDialog(rootPane, "Please select the * fields \n* fields are mandatory");
            } else {
                if (Constant.genoList.contains(geno_filename) || Constant.phenoList.contains(pheno_filename)) {
                    if (Constant.genoList.contains(geno_filename) && Constant.phenoList.contains(pheno_filename)) {
                        JOptionPane.showMessageDialog(rootPane, "Genotype & Phenotype  file selected already opened \nPlease select other genotype file");
                    } else {
                        if (Constant.genoList.contains(geno_filename)) {
                            JOptionPane.showMessageDialog(rootPane, "Genotype file selected already opened \nPlease select other genotype file");
                        }
                        if (Constant.phenoList.contains(pheno_filename)) {
                            JOptionPane.showMessageDialog(rootPane, "Phenotype file selected already opened \nPlease select other genotype file");
                        }
                    }
                } else {
                    Constant.resultdirectory = result; //setting the result directory as constant for this session
                    validate_geno = 1;
                    validate_pheno = 1;
                    //checking whether geno and pheno files dn't have any errors
                    if (validate_geno == 1 && validate_pheno == 1) {
                        Constant.genotype = geno;
                        Constant.phenotype = pheno;
                        dispose();
                    } else {
                        setVisible(false);
                        GsMethods gg = new GsMethods();
                        //if geno file have errors then it is 0
                        if (validate_geno == 0) {
                            //reading the geno error file and displaying it
                            JScrollPane textFileReader = gg.textFileReader(Constant.resultdirectory + geno_filename + "_errors.txt");
                            tabbedPane.add(textFileReader);
                        }
                        if (validate_pheno == 0) {
                            //reading the error file and display it
                            JScrollPane textFileReader = gg.textFileReader(Constant.resultdirectory + pheno_filename + "_errors.txt");
                            tabbedPane.add(textFileReader);
                        }
                    }
                }
            }

        }
        if (source == btnCancle) { //if cancle button is selected 
            Constant.genotype = null;
            Constant.phenotype = null;
            dispose(); //close the dialog box and also unallocated the memory also assigned to it 
        }

        if (source == btnGenotype) { //if source is genotype button
            JFileChooser chooser;
            //getting the perviously choosen directory
            lastFileChoosen = Constant.browsepath;
            //checking the last file choosen is empty or not
            if (lastFileChoosen.isEmpty()) {
                chooser = new JFileChooser();
            } else {
                chooser = new JFileChooser(lastFileChoosen);
            }
            chooser.removeChoosableFileFilter(chooser.getFileFilter());
            chooser.setDialogTitle("Select a genotype file");//setting the dialog title
            chooser.setApproveButtonText("Select File");
            //setting the file filter 
            FileNameExtensionFilter filter = new FileNameExtensionFilter("CSV (Comma delimited) (*.csv)", "csv");
            chooser.setFileFilter(filter);
            int r = chooser.showOpenDialog(rootPane);

            if (r == JFileChooser.APPROVE_OPTION) {
                txtGenotype.setText(chooser.getSelectedFile().getPath());
                File fileValidation = new File(txtGenotype.getText());
                boolean genotypeExists = fileValidation.exists();
                Constant.browsepath = chooser.getSelectedFile().getAbsolutePath();
                if (genotypeExists == false) {
                    JOptionPane.showMessageDialog(rootPane, "Selected file doesnot exists \nPlease select valid genotype file");
                    txtGenotype.setText("");
                }
            }

            if (r == JFileChooser.CANCEL_OPTION) {
                //donothing
            }
        }

        if(source == btnBMS){
        	System.out.println("Connect to BMS");
         //   BMSConnect bmsPanel = new BMSConnect(frame);
          // bmsPanel.setLocationRelativeTo(frame);
          //  bmsPanel.setVisible(true);
        	ConnectToPhenotypeDB phenoDB = new ConnectToPhenotypeDB(frame, true);
        	phenoDB.setVisible(true);
        	System.out.println("Connected to BMS");
        }
        
        if(source == btnGOBII){
        	System.out.println("Connect to GOBII");
        //    BMSConnect bmsPanel = new BMSConnect(frame);
        //    bmsPanel.setLocationRelativeTo(frame);
        //    bmsPanel.setVisible(true);
        	ConnectToGenotypeDB genoDB = new ConnectToGenotypeDB(frame, true);
        	genoDB.setVisible(true);
        	System.out.println("Connected to GOBII");
        }
        
        if (source == btnPhenotype) { //if the phenotype button 
            JFileChooser chooser;
            //getting the browsePath if file already choosen  
            lastFileChoosen = Constant.browsepath;
            //check whether it is empty or not 
            if (lastFileChoosen.isEmpty()) {
                chooser = new JFileChooser();
            } else {
                chooser = new JFileChooser(lastFileChoosen);
            }
            chooser.removeChoosableFileFilter(chooser.getFileFilter());
            chooser.setDialogTitle("Select a phenotype  file");//setting the title
            chooser.setApproveButtonText("Select File");
            FileNameExtensionFilter filter = new FileNameExtensionFilter("CSV (Comma delimited) (*.csv)", "csv"); //adding the file filter csv
            chooser.setFileFilter(filter);
            int r = chooser.showOpenDialog(rootPane);

            if (r == JFileChooser.APPROVE_OPTION) {
                txtPhenotype.setText(chooser.getSelectedFile().getPath());
                Constant.browsepath = chooser.getSelectedFile().getAbsolutePath();
                File fileValidation = new File(txtPhenotype.getText());
                boolean phenotypeExists = fileValidation.exists();
                Constant.browsepath = chooser.getSelectedFile().getAbsolutePath();
                if (phenotypeExists == false) {
                    JOptionPane.showMessageDialog(rootPane, "Selected file doesnot exists \nPlease select valid phenotype file");
                    txtPhenotype.setText("");
                }
            }

            if (r == JFileChooser.CANCEL_OPTION) {
                //donothing
            }
        }

        if (source == btnResultDir) { //if the result button selected
            JFileChooser folderchooser = new JFileChooser(); //initilizing the jfilechooser 
            folderchooser.setDialogTitle("Select a folder where result should be stored"); //setting the title 
            folderchooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY); //setting the condiation to select only directorys 
            folderchooser.removeChoosableFileFilter(folderchooser.getFileFilter());
            // disable the "All files" option.
            FileNameExtensionFilter filter = new FileNameExtensionFilter("Any Folder", "FOLDER");
            folderchooser.setFileFilter(filter);
            folderchooser.setApproveButtonText("Select Folder");
            folderchooser.setApproveButtonToolTipText("Select the folder for working and result directoy");
            if (folderchooser.showOpenDialog(this) == JFileChooser.APPROVE_OPTION) {
                OSCheck.OSType os = OSCheck.getOperatingSystemType();
                if ("Window".equals(os.toString())) {
                    txtResultDir.setText(folderchooser.getSelectedFile().toString() + "\\"); //getting the selected directory path and setting to resultdir text field
                    if (Files.isDirectory(Paths.get(txtResultDir.getText()))) {
                        //folder is exists
                        //donothing
                    } else {
                        //folder not exits
                        int option = JOptionPane.showOptionDialog(rootPane, "Folder path doesnot exits \nDo you want to create folder for the path", "New Folder Creation", JOptionPane.YES_NO_OPTION,
                                JOptionPane.WARNING_MESSAGE, null, null, null);
                        if (option == JOptionPane.YES_OPTION) {
                            File fileCreate = new File(txtResultDir.getText());
                            boolean mkdir = fileCreate.mkdir();
                        }
                        if (option == JOptionPane.NO_OPTION) {
                            txtResultDir.setText("");
                        }
                    }
                }
                if ("Linux".equals(os.toString())) {
                    txtResultDir.setText(folderchooser.getSelectedFile().toString() + "/"); //getting the selected directory path and setting to resultdir text field
                    if (Files.isDirectory(Paths.get(txtResultDir.getText()))) {
                        //folder is exists
                    } else {
                        //folder not exits
                        int option = JOptionPane.showOptionDialog(rootPane, "Folder path doesnot exits \nDo you want to create folder for the path", "New Folder Creation", JOptionPane.YES_NO_OPTION,
                                JOptionPane.WARNING_MESSAGE, null, null, null);
                        if (option == JOptionPane.YES_OPTION) {
                            File fileCreate = new File(txtResultDir.getText());
                            boolean mkdir = fileCreate.mkdir();
                        }
                        if (option == JOptionPane.NO_OPTION) {
                            txtResultDir.setText("");
                        }
                    }
                }
                if ("MacOS".equals(os.toString())) {
                    txtResultDir.setText(folderchooser.getSelectedFile().toString() + "/"); //getting the selected directory path and setting to resultdir text field
                }
            }
        }
    }
//end region constructor 

    //start region 
    //initilzing all the GUI components on the dialog box    
    private void initComponents() {
        setDefaultCloseOperation(javax.swing.WindowConstants.DISPOSE_ON_CLOSE);
        setTitle("Load files");
        lblGeno.setFont(new java.awt.Font("DejaVu Sans", 0, 15)); // NOI18N
        lblGeno.setText("Select a genotype file *:");
        txtGenotype.setFont(new java.awt.Font("DejaVu Sans", 0, 15)); // NOI18N
        txtGenotype.setText(" ");

        btnGenotype.setFont(new java.awt.Font("DejaVu Sans", 0, 15)); // NOI18N
        btnGenotype.setText("Browse");

        lblPheno.setFont(new java.awt.Font("DejaVu Sans", 0, 15)); // NOI18N
        lblPheno.setText("Select a phenotype file *:");

        txtPhenotype.setFont(new java.awt.Font("DejaVu Sans", 0, 15)); // NOI18N
        txtPhenotype.setText(" ");

        btnPhenotype.setFont(new java.awt.Font("DejaVu Sans", 0, 15)); // NOI18N
        btnPhenotype.setText("Browse");

        btnBMS.setFont(new java.awt.Font("DejaVu Sans", 0, 15)); // NOI18N
        btnBMS.setText("Connect");
        btnBMS.setBackground(Color.WHITE);
        btnBMS.setIcon(Helper.createIcon("Download_database_Icon.gif"));

        btnGOBII.setFont(new java.awt.Font("DejaVu Sans", 0, 15)); // NOI18N
        btnGOBII = new JButton("Connect");
        btnGOBII.setFont(new Font("Dialog", Font.PLAIN, 15));
        btnGOBII.setBackground(Color.WHITE);
        btnGOBII.setIcon(Helper.createIcon("Download_database_Icon.gif"));
        
        lblResultDir.setFont(new java.awt.Font("DejaVu Sans", 0, 15)); // NOI18N
        lblResultDir.setText("Select a result directory *:");

        txtResultDir.setFont(new java.awt.Font("DejaVu Sans", 0, 15)); // NOI18N
        txtResultDir.setText(" ");

        btnResultDir.setFont(new java.awt.Font("DejaVu Sans", 0, 15)); // NOI18N
        btnResultDir.setText("Browse");

        btnOk.setFont(new java.awt.Font("DejaVu Sans", 0, 15)); // NOI18N
        btnOk.setText("Ok");

        btnCancle.setFont(new java.awt.Font("DejaVu Sans", 0, 15)); // NOI18N
        btnCancle.setText("Cancel");
        
        
        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(getContentPane());
        layout.setHorizontalGroup(
        	layout.createParallelGroup(Alignment.LEADING)
        		.addGroup(layout.createSequentialGroup()
        			.addContainerGap()
        			.addGroup(layout.createParallelGroup(Alignment.TRAILING)
        				.addGroup(layout.createSequentialGroup()
        					.addGroup(layout.createParallelGroup(Alignment.LEADING)
        						.addComponent(lblGeno)
        						.addComponent(lblPheno)
        						.addComponent(lblResultDir))
        					.addGroup(layout.createParallelGroup(Alignment.LEADING, false)
        						.addComponent(txtGenotype)
        						.addComponent(txtPhenotype, GroupLayout.DEFAULT_SIZE, 174, Short.MAX_VALUE)
        						.addComponent(txtResultDir))
        					.addPreferredGap(ComponentPlacement.RELATED)
        					.addGroup(layout.createParallelGroup(Alignment.LEADING)
        						.addComponent(btnGenotype)
        						.addComponent(btnPhenotype)
        						.addComponent(btnResultDir)))
        				.addGroup(layout.createSequentialGroup()
        					.addComponent(btnOk, GroupLayout.PREFERRED_SIZE, 66, GroupLayout.PREFERRED_SIZE)
        					.addGap(18)
        					.addComponent(btnCancle)
        					.addGap(159)))
        			.addPreferredGap(ComponentPlacement.RELATED)
        			.addGroup(layout.createParallelGroup(Alignment.LEADING)
        				.addComponent(btnBMS, GroupLayout.DEFAULT_SIZE, 89, Short.MAX_VALUE)
        				.addComponent(btnGOBII, GroupLayout.DEFAULT_SIZE, 89, Short.MAX_VALUE)))
        );
        layout.setVerticalGroup(
        	layout.createParallelGroup(Alignment.LEADING)
        		.addGroup(layout.createSequentialGroup()
        			.addGap(32)
        			.addGroup(layout.createParallelGroup(Alignment.BASELINE)
        				.addComponent(lblGeno)
        				.addComponent(txtGenotype, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
        				.addComponent(btnGenotype)
        				.addComponent(btnGOBII, GroupLayout.PREFERRED_SIZE, 29, GroupLayout.PREFERRED_SIZE))
        			.addGap(18)
        			.addGroup(layout.createParallelGroup(Alignment.BASELINE)
        				.addComponent(lblPheno)
        				.addComponent(txtPhenotype, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
        				.addComponent(btnPhenotype)
        				.addComponent(btnBMS, GroupLayout.PREFERRED_SIZE, 29, GroupLayout.PREFERRED_SIZE))
        			.addGap(18)
        			.addGroup(layout.createParallelGroup(Alignment.BASELINE)
        				.addComponent(lblResultDir)
        				.addComponent(txtResultDir, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
        				.addComponent(btnResultDir))
        			.addPreferredGap(ComponentPlacement.RELATED, 36, Short.MAX_VALUE)
        			.addGroup(layout.createParallelGroup(Alignment.BASELINE)
        				.addComponent(btnOk)
        				.addComponent(btnCancle))
        			.addGap(30))
        );
        getContentPane().setLayout(layout);

        pack();
    }
//end region 

    //start region 
    //main method 
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
        } catch (ClassNotFoundException ex) {
            java.util.logging.Logger.getLogger(OpenDial.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (InstantiationException ex) {
            java.util.logging.Logger.getLogger(OpenDial.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (IllegalAccessException ex) {
            java.util.logging.Logger.getLogger(OpenDial.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (javax.swing.UnsupportedLookAndFeelException ex) {
            java.util.logging.Logger.getLogger(OpenDial.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        }
        //</editor-fold>
      
    }
    
}
