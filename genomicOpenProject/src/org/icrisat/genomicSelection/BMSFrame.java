/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.icrisat.genomicSelection;

import java.awt.Dimension;
import java.awt.Font;
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
 * @author Chaitanya
 */
public class BMSFrame extends JDialog implements ActionListener {

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
    public BMSFrame(java.awt.Frame parent, boolean modal) {
    	 super(parent, modal);
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
        btnBMS.setText("BMS");

        btnGOBII.setFont(new java.awt.Font("DejaVu Sans", 0, 15)); // NOI18N
        btnGOBII = new JButton("GOBII");
        btnGOBII.setFont(new Font("Dialog", Font.PLAIN, 15));
        
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
        				.addComponent(btnGOBII, GroupLayout.DEFAULT_SIZE, 77, Short.MAX_VALUE)))
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
}
