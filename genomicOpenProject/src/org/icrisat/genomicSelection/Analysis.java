package org.icrisat.genomicSelection;

import com.csvreader.CsvReader;
import java.awt.HeadlessException;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.IOException;
import java.io.File;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.DefaultListModel;
import javax.swing.JOptionPane;

/**
 * @author Mohan get the genotype,phenotype,relationship ,pedigree and
 * population file names selecting the traits , selecting methods names enabling
 * parameters required
 */
public final class Analysis extends javax.swing.JDialog {

    DefaultListModel resultList1 = new DefaultListModel();
    DefaultListModel resultList2 = new DefaultListModel();
    public List genoList_Genotype;
    public List genoList_Phenotype;
    GsMethods gsmethod = new GsMethods();
    public List sortingTraitsSelected = new ArrayList();

    //start region  
    //initilzing the variables and setting action listeners 
    public Analysis(java.awt.Frame parent, boolean modal) {
        super(parent, modal);
        this.setIconImage(Toolkit.getDefaultToolkit().getImage(getClass().getResource("/org/icrisat/genomicSelection/image/log_icrisat.png")));
        this.setLocationRelativeTo(parent);
        //setting the analysis dialog box to the middle application frame based on the present size of the frame
        this.setLocation(parent.getSize().width / 2 - 381, parent.getSize().height / 2 - 296);
        initComponents();
        //setting the analysisCanceled to false before analysis start
        Constant.analysisCancled = false;
        forestsSpin.setEnabled(false);
        //clear the content in then  resultGeratedList 
        Constant.resultGenratedList.clear();

        //clear list selectedTraits 
        Constant.selectedTraits.clear();

        //checking the whether R or Fortran is selected 
        if (Constant.engineSelection == 1) {
            chkFtBayesA.setEnabled(false);
            chkFtRidgeRegression.setEnabled(false);
            tabPaneRFortran.setSelectedIndex(0);
        }
        if (Constant.engineSelection == 2) {
            tabPaneRFortran.setSelectedIndex(1);
            chkBayeCpi.setEnabled(false);
            chkBayesB.setEnabled(false);
            chkBayesLasso.setEnabled(false);
            //   chkCrossValidation.setEnabled(false);
            chkKinshipGauss.setEnabled(false);
            chkRandomForest.setEnabled(false);
            chkRidgeRegression.setEnabled(false);
        }

        boolean check = Constant.traitList.isEmpty();
        if (check == false) {
            Constant.traitList.clear();
        }
        //adding the genotype file name to genocombo
        Iterator genoItr = Constant.genoList.iterator();
        while (genoItr.hasNext()) {
            genoCombo.addItem(genoItr.next());
        }
        genoCombo.setSelectedIndex(1);  //set selected the index 1

        //adding the phenotype file names to phenocombo
        Iterator phenoItr = Constant.phenoList.iterator();
        while (phenoItr.hasNext()) {
            phenoCombo.addItem(phenoItr.next());
        }
        phenoCombo.setSelectedIndex(1);

        String phenoFileName = Constant.phenoList.get(0);
        //setting the first phenotype file traits displayS
        String phenoFileName1 = Constant.phenoList.get(0);
        if ("Select".equals(phenoFileName)) {
            //donothing
        } else {
            phenoTraitReader(phenoFileName1);
        }

        //adding the population file names to populationcombo
        Iterator populationItr = Constant.populationList.iterator();
        while (populationItr.hasNext()) {
            populationCombo.addItem(populationItr.next());
        }
        //setting the values of Bayess,crossvalidation,RandomForest addtional parameters
        Constant.BayessAddParameters = null;
        Constant.RandomForestAddParameteres = 0;
        Constant.crossValidationParameteres = null;
        Constant.noOfCPUs = 1;

        //setting the values to default in analysis
        Constant.perMissing = 0;
        Constant.picValue = 0.0;
        Constant.mafValue = 0.0;

        //getting the no of processors avilable in system
        int cores = Runtime.getRuntime().availableProcessors();
        //setting max value of spiner to aviable cores in the system, default is one 
        coresSpin.setModel(new javax.swing.SpinnerNumberModel(1, 1, cores, 1));

        phenoCombo.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                String phenofile = phenoCombo.getSelectedItem().toString();
                System.out.println("selection file is " + phenofile);
                if ("Select".equals(phenofile)) {
                    resultList1.clear();
                    jList1.setModel(resultList1);
                } else {
                    phenoTraitReader(phenofile);
                }
            }
        });

        chkRandomForest.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if (chkRandomForest.isSelected() == true) {
                    forestsSpin.setEnabled(true);
                }
                if (chkRandomForest.isSelected() == false) {
                    forestsSpin.setEnabled(false);
                }
            }
        });

        btnStart.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
            	System.out.println("Analysys Start button pressed");
                String genofile = genoCombo.getSelectedItem().toString();
                String phenfile = phenoCombo.getSelectedItem().toString();
                String popFile = populationCombo.getSelectedItem().toString(); //coverinet file
                int list2size = jList2.getModel().getSize();

                if ("Select".equals(genofile) || "Select".equals(phenfile)) {
                    JOptionPane.showMessageDialog(null, "Please select the genotype and phenotype files");
                }
                if (!"Select".equals(genofile) && !"Select".equals(phenfile)) {
                    if (list2size == 0) {
                        JOptionPane.showMessageDialog(null, "Please select atleast one triat from trait list for analysis");
                    }
                    if (chkRidgeRegression.isSelected() || chkBayeCpi.isSelected() || chkBayesB.isSelected() || chkBayesLasso.isSelected() || chkRandomForest.isSelected() || chkKinshipGauss.isSelected() || chkFtBayesA.isSelected() || chkFtRidgeRegression.isSelected()) {
                    } else {
                        JOptionPane.showMessageDialog(null, "Please select atleast one method for analysis");
                    }
                }
                if (list2size != 0) {
                    if (chkRidgeRegression.isSelected() || chkBayeCpi.isSelected() || chkBayesB.isSelected() || chkBayesLasso.isSelected() || chkRandomForest.isSelected() || chkKinshipGauss.isSelected() || chkFtBayesA.isSelected() || chkFtRidgeRegression.isSelected()) {

                        if ((!"Select".equals(genofile) && !"Select".equals(phenfile)) && list2size != 0 && (chkRidgeRegression.isSelected() || chkBayeCpi.isSelected() || chkBayesB.isSelected() || chkBayesLasso.isSelected() || chkRandomForest.isSelected() || chkKinshipGauss.isSelected() || chkFtBayesA.isSelected() || chkFtRidgeRegression.isSelected())) {
                            {
                                //assigning whethere which methods are selected 
                                Constant.RidgeReg = chkRidgeRegression.isSelected();
                                Constant.BayesCpi = chkBayeCpi.isSelected();
                                Constant.BayesB = chkBayesB.isSelected();
                                Constant.BayesLasso = chkBayesLasso.isSelected();
                                Constant.RandomeForest = chkRandomForest.isSelected();
                                Constant.KinshipGauss = chkKinshipGauss.isSelected();
                                Constant.BayeAFortran = chkFtBayesA.isSelected();
                                Constant.RidgeRegressionFortran = chkFtRidgeRegression.isSelected();

                                //bayes parameters for fortran script
                                Constant.bayesFortParameters = roundsSpin.getValue() + " " + burninSpin.getValue() + " " + thiningSpin.getValue();

                                //RandomForest validation value
                                Constant.RandomForestAddParameteres = (int) forestsSpin.getValue();

                                //CrossValidation values 
                                Constant.crossValidationParameteres = replicationSpin.getValue() + " " + foldSpin.getValue();
                                Constant.noOfCPUs = (int) coresSpin.getValue();
                                Constant.perMissing = (int) perMissingMarker.getValue();
                                Constant.picValue = (Double) picSpin.getValue();
                                Constant.mafValue = (Double) mafSpin.getValue();

                                //File names selected for analysis geno,pheno and population(Coverient) names
                                Constant.geno_Analysis = genofile;
                                Constant.popu_Analysis = popFile;  //coverient

                                //check if all the triats are selected 
                                int orginalTraitsCount1 = Constant.traitCount;
                                int selecteTraitsCount1 = jList2.getModel().getSize();
                                if (orginalTraitsCount1 - 1 == selecteTraitsCount1) {
                                    try {
                                        File f = File.createTempFile("tmp", ".csv", new File(Constant.resultdirectory));
                                        String tmpFilename = f.getName();
                                        String orgFile1 = Constant.resultdirectory + phenoCombo.getSelectedItem().toString();
                                        String tmpFile1 = Constant.resultdirectory + tmpFilename;
                                        CopyFile.copy(orgFile1, tmpFile1);
                                        Constant.pheno_Analysis = tmpFilename;
                                        for (int i = 0; i < list2size; i++) {
                                            String traitSelected = (String) jList2.getModel().getElementAt(i);
                                            sortingTraitsSelected.add(traitSelected);
                                        }
                                        //sorting the selected traits according to the input file given
                                        sortingTraitNames(orgFile1);
                                    } catch (IOException tmp) {
                                        System.out.println(tmp.getMessage());
                                    }
                                }

                                //if only some traits are selected then creating a temp file with those selected traits colunms 
                                if (orginalTraitsCount1 - 1 != selecteTraitsCount1) {
                                    try {
                                        File f = File.createTempFile("tmp", ".csv", new File(Constant.resultdirectory));
                                        String tmpFilename = f.getName();
                                        String orgFile1 = Constant.resultdirectory + phenoCombo.getSelectedItem().toString();
                                        String tmpFile1 = Constant.resultdirectory + tmpFilename;

                                        for (int i = 0; i < list2size; i++) {
                                            String traitSelected = (String) jList2.getModel().getElementAt(i);
                                            sortingTraitsSelected.add(traitSelected);
                                            Constant.selectedTraits.add(traitSelected);
                                        }
                                        sortingTraitNames(orgFile1); //sorting the selected traits according to the traits

                                        gsmethod.traitRemove(orgFile1, tmpFile1);
                                        Constant.pheno_Analysis = tmpFilename;
                                    } catch (IOException tmp) {
                                        System.out.println(tmp.getMessage());
                                    }
                                }
                                dispose();
                            }
                        }
                    }
                }
            }
        });
        btnCancel.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                closeAnalysis();
            }
        });

        btndoubleiparrow.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                DefaultListModel listModel = (DefaultListModel) jList1.getModel();
                listModel.removeAllElements();
                String phen = Constant.phenotraits;
                DefaultListModel resultList = new DefaultListModel();
                String regexp = "[,]";
                String[] tokens;
                int traitcount = Constant.traitCount;
                for (int i = 1; i < traitcount; i++) {
                    tokens = phen.split(regexp);
                    resultList.addElement(tokens[i]);
                }
                jList2.setModel(resultList);
            }
        });
        btndoubleremovearrow.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                DefaultListModel listModel = (DefaultListModel) jList2.getModel();
                listModel.removeAllElements();
                String phen = Constant.phenotraits;
                DefaultListModel resultList = new DefaultListModel();
                String regexp = "[,]";
                String[] tokens;
                int traitcount = Constant.traitCount;
                for (int i = 1; i < traitcount; i++) {
                    tokens = phen.split(regexp);
                    resultList.addElement(tokens[i]);
                }
                jList1.setModel(resultList);
            }
        });

        btnsingleiparrow.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                DefaultListModel list2 = (DefaultListModel) jList2.getModel();
                DefaultListModel list = (DefaultListModel) jList1.getModel();
                java.util.List select = jList1.getSelectedValuesList();
                for (Object select1 : select) {
                    if (list2.contains(select1)) {
                    } else {
                        list2.addElement(select1);
                        list.removeElement(select1);
                    }
                }
                jList2.setModel(list2);
            }
        });

        btnsingleremovearrow.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                DefaultListModel listModelJlist1 = (DefaultListModel) jList1.getModel();
                DefaultListModel listModel1Jlist2 = (DefaultListModel) jList2.getModel();
                java.util.List selected = jList2.getSelectedValuesList();
                for (Object selected1 : selected) {
                    if (listModelJlist1.contains(selected1)) {
                    } else {
                        listModelJlist1.addElement(selected1);
                        listModel1Jlist2.removeElement(selected1);
                    }
                }
                jList1.setModel(listModelJlist1);
            }
        });
        this.addWindowListener(new WindowAdapter() {
            @Override
            public void windowClosing(WindowEvent e) {

                closeAnalysis();
            }
        });
    }
    public List traits_list;

    public void sortingTraitNames(String phenoPath) {
        try {
            //getting only the headers from it 
            CsvReader reader = new CsvReader(phenoPath, ',');
            reader.readHeaders();
            String[] headers = reader.getHeaders();
            traits_list = new ArrayList();
            traits_list = sortingTraitsSelected;

            FileInputStream stream = new FileInputStream(phenoPath);
            DataInputStream inputStream = new DataInputStream(stream);
            BufferedReader br = new BufferedReader(new InputStreamReader(inputStream));
            String[] str;
            int num = 0;
            String line;
            int headerCount;
            while ((line = br.readLine()) != null) {
                str = line.split(",");
                if (num == 0) {
                    headerCount = str.length;
                    for (int i = 0; i < headerCount; i++) {
                        if (traits_list.contains(headers[i])) {
                            Constant.traitList.add(headers[i]);
                        }
                    }
                }
                num++;
            }
        } catch (FileNotFoundException ex) {
            Logger.getLogger(Analysis.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ex) {
            Logger.getLogger(Analysis.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    //Getting the traits from selected phenotype file from the list 
    public void phenoTraitReader(String filepathSrc) {
        String resultDir = Constant.resultdirectory;
        String filePath = resultDir + filepathSrc;
        int headerCount;
        try {
            CsvReader reader = new CsvReader(filePath, ',');
            reader.readHeaders();
            String[] headers = reader.getHeaders();
            String listnames = ""; //local variable listnames
            String pheno = ""; //local variable pheno
            //Get the count of headers
            headerCount = reader.getHeaderCount();
            Constant.traitCount = headerCount;
            int traitCount = headerCount; // getting the headercount getting no of traits
            //for loop for getting all the headers names and validating the names
            for (int i = 0; i < headerCount; i++) {
                if (headers[i].matches("[a-zA-z0-9_]*")) { //condtion to check the header name is alpha numerical 
                    pheno = pheno + headers[i] + ","; //adding the name to pheno string 
                } else {//variable not alpha numerical 
                    listnames = listnames + headers[i] + ","; // adding to the listname string as it not alpha numerical 
                }
            }
            if (listnames.isEmpty()) { //checking whether any alphanumerical are present 
                Constant.phenotraits = pheno; //adding the header names to the phenotraits constant string
                String phenofiles = Constant.phenotraits;
                String regexp = ",";
                String[] tokens;
                resultList1.clear();
                resultList2.clear();
                for (int i = 1; i < traitCount; i++) {
                    tokens = phenofiles.split(regexp);
                    resultList1.addElement(tokens[i]);
                }
                jList1.setModel(resultList1);
                jList2.setModel(resultList2);
            } else {//non-alpha numerical present 
                //adding all the non-alphanumerical to jdialog boxes and display to the user
                JOptionPane.showMessageDialog(null, "Phenotype file only alpha-numercial characters and '_' is accpeted.Please remove the unwanted characters at" + listnames);
            }
        } catch (IOException | HeadlessException ex) {
            System.out.println(ex);
        }
    }

    private void closeAnalysis() {
        dispose();
        Constant.analysisFileName = null;
        Constant.analysisCancled = true;
    }

    /**
     * Method which build gui by darging and droping the the required
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        fileNameSelectingPanel = new javax.swing.JPanel();
        jLabel1 = new javax.swing.JLabel();
        genoCombo = new javax.swing.JComboBox();
        phenoCombo = new javax.swing.JComboBox();
        jLabel9 = new javax.swing.JLabel();
        populationCombo = new javax.swing.JComboBox();
        lblPheno = new javax.swing.JLabel();
        traitPanel = new javax.swing.JPanel();
        jScrollPane1 = new javax.swing.JScrollPane();
        jList1 = new javax.swing.JList();
        jScrollPane2 = new javax.swing.JScrollPane();
        jList2 = new javax.swing.JList();
        btnsingleiparrow = new javax.swing.JButton();
        btndoubleiparrow = new javax.swing.JButton();
        btndoubleremovearrow = new javax.swing.JButton();
        btnsingleremovearrow = new javax.swing.JButton();
        methodPanel = new javax.swing.JPanel();
        tabPaneRFortran = new javax.swing.JTabbedPane();
        rPanel = new javax.swing.JPanel();
        chkRidgeRegression = new javax.swing.JCheckBox();
        chkBayeCpi = new javax.swing.JCheckBox();
        chkBayesB = new javax.swing.JCheckBox();
        chkBayesLasso = new javax.swing.JCheckBox();
        chkRandomForest = new javax.swing.JCheckBox();
        chkKinshipGauss = new javax.swing.JCheckBox();
        fortranPanel = new javax.swing.JPanel();
        chkFtBayesA = new javax.swing.JCheckBox();
        chkFtRidgeRegression = new javax.swing.JCheckBox();
        additionalParPanel = new javax.swing.JPanel();
        jPanel2 = new javax.swing.JPanel();
        jLabel13 = new javax.swing.JLabel();
        forestsSpin = new javax.swing.JSpinner();
        jPanel3 = new javax.swing.JPanel();
        lblRounds2 = new javax.swing.JLabel();
        roundsSpin = new javax.swing.JSpinner();
        lblRounds3 = new javax.swing.JLabel();
        burninSpin = new javax.swing.JSpinner();
        lblother1 = new javax.swing.JLabel();
        thiningSpin = new javax.swing.JSpinner();
        jPanel4 = new javax.swing.JPanel();
        jLabel7 = new javax.swing.JLabel();
        jLabel8 = new javax.swing.JLabel();
        replicationSpin = new javax.swing.JSpinner();
        foldSpin = new javax.swing.JSpinner();
        jPanel5 = new javax.swing.JPanel();
        jLabel14 = new javax.swing.JLabel();
        coresSpin = new javax.swing.JSpinner();
        btnStart = new javax.swing.JButton();
        btnCancel = new javax.swing.JButton();
        jPanel1 = new javax.swing.JPanel();
        jLabel10 = new javax.swing.JLabel();
        perMissingMarker = new javax.swing.JSpinner();
        jLabel11 = new javax.swing.JLabel();
        picSpin = new javax.swing.JSpinner();
        jLabel12 = new javax.swing.JLabel();
        mafSpin = new javax.swing.JSpinner();

        setDefaultCloseOperation(javax.swing.WindowConstants.DO_NOTHING_ON_CLOSE);
        setTitle("Analysis");
        setFont(new java.awt.Font("Dialog", 0, 15)); // NOI18N

        fileNameSelectingPanel.setBorder(javax.swing.BorderFactory.createTitledBorder(null, "Select file names from combo box", javax.swing.border.TitledBorder.CENTER, javax.swing.border.TitledBorder.DEFAULT_POSITION, new java.awt.Font("SansSerif", 1, 15))); // NOI18N

        jLabel1.setFont(new java.awt.Font("DejaVu Sans", 0, 15)); // NOI18N
        jLabel1.setText("Genotype file name :");

        genoCombo.setFont(new java.awt.Font("DejaVu Sans", 0, 15)); // NOI18N
        genoCombo.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "Select" }));

        phenoCombo.setFont(new java.awt.Font("DejaVu Sans", 0, 15)); // NOI18N
        phenoCombo.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "Select" }));

        jLabel9.setFont(new java.awt.Font("DejaVu Sans", 0, 15)); // NOI18N
        jLabel9.setText("Covariate file name :");

        populationCombo.setFont(new java.awt.Font("DejaVu Sans", 0, 15)); // NOI18N
        populationCombo.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "Select" }));

        lblPheno.setFont(new java.awt.Font("DejaVu Sans", 0, 15)); // NOI18N
        lblPheno.setText("Phenotype file name :");

        javax.swing.GroupLayout fileNameSelectingPanelLayout = new javax.swing.GroupLayout(fileNameSelectingPanel);
        fileNameSelectingPanel.setLayout(fileNameSelectingPanelLayout);
        fileNameSelectingPanelLayout.setHorizontalGroup(
            fileNameSelectingPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(fileNameSelectingPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(fileNameSelectingPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                    .addComponent(jLabel9, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(lblPheno, javax.swing.GroupLayout.Alignment.TRAILING, javax.swing.GroupLayout.DEFAULT_SIZE, 207, Short.MAX_VALUE)
                    .addComponent(jLabel1, javax.swing.GroupLayout.Alignment.TRAILING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                .addGap(18, 18, 18)
                .addGroup(fileNameSelectingPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                    .addComponent(genoCombo, 0, 160, Short.MAX_VALUE)
                    .addComponent(phenoCombo, 0, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(populationCombo, 0, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                .addContainerGap())
        );
        fileNameSelectingPanelLayout.setVerticalGroup(
            fileNameSelectingPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(fileNameSelectingPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(fileNameSelectingPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jLabel1)
                    .addComponent(genoCombo, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(18, 18, 18)
                .addGroup(fileNameSelectingPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(phenoCombo, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(lblPheno))
                .addGap(18, 18, 18)
                .addGroup(fileNameSelectingPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jLabel9)
                    .addComponent(populationCombo, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addContainerGap(15, Short.MAX_VALUE))
        );

        traitPanel.setBorder(javax.swing.BorderFactory.createTitledBorder(null, "Select the trait(s) for analysis", javax.swing.border.TitledBorder.CENTER, javax.swing.border.TitledBorder.DEFAULT_POSITION, new java.awt.Font("SansSerif", 1, 15), java.awt.Color.black)); // NOI18N

        jList1.setFont(new java.awt.Font("DejaVu Sans", 0, 15)); // NOI18N
        jList1.addMouseListener(new java.awt.event.MouseAdapter() {
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                jList1MouseClicked(evt);
            }
        });
        jScrollPane1.setViewportView(jList1);

        jList2.setFont(new java.awt.Font("DejaVu Sans", 0, 15)); // NOI18N
        jList2.addMouseListener(new java.awt.event.MouseAdapter() {
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                jList2MouseClicked(evt);
            }
        });
        jScrollPane2.setViewportView(jList2);

        btnsingleiparrow.setFont(new java.awt.Font("DejaVu Sans", 0, 15)); // NOI18N
        btnsingleiparrow.setText(">");

        btndoubleiparrow.setFont(new java.awt.Font("DejaVu Sans", 0, 15)); // NOI18N
        btndoubleiparrow.setText(">>");

        btndoubleremovearrow.setFont(new java.awt.Font("DejaVu Sans", 0, 15)); // NOI18N
        btndoubleremovearrow.setText("<<");

        btnsingleremovearrow.setFont(new java.awt.Font("DejaVu Sans", 0, 15)); // NOI18N
        btnsingleremovearrow.setText(" <  ");

        javax.swing.GroupLayout traitPanelLayout = new javax.swing.GroupLayout(traitPanel);
        traitPanel.setLayout(traitPanelLayout);
        traitPanelLayout.setHorizontalGroup(
            traitPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(traitPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addComponent(jScrollPane1, javax.swing.GroupLayout.PREFERRED_SIZE, 136, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, 28, Short.MAX_VALUE)
                .addGroup(traitPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                    .addComponent(btndoubleremovearrow, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(btnsingleremovearrow, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(btndoubleiparrow, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(btnsingleiparrow, javax.swing.GroupLayout.PREFERRED_SIZE, 53, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(26, 26, 26)
                .addComponent(jScrollPane2, javax.swing.GroupLayout.PREFERRED_SIZE, 142, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addContainerGap())
        );
        traitPanelLayout.setVerticalGroup(
            traitPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(traitPanelLayout.createSequentialGroup()
                .addGroup(traitPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(traitPanelLayout.createSequentialGroup()
                        .addGap(20, 20, 20)
                        .addComponent(btnsingleiparrow)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addComponent(btndoubleiparrow)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addComponent(btndoubleremovearrow)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addComponent(btnsingleremovearrow)
                        .addGap(0, 0, Short.MAX_VALUE))
                    .addComponent(jScrollPane1, javax.swing.GroupLayout.Alignment.TRAILING)
                    .addComponent(jScrollPane2, javax.swing.GroupLayout.Alignment.TRAILING))
                .addContainerGap())
        );

        methodPanel.setBorder(javax.swing.BorderFactory.createTitledBorder(null, "Select a method(s) to start analysis", javax.swing.border.TitledBorder.CENTER, javax.swing.border.TitledBorder.DEFAULT_POSITION, new java.awt.Font("SansSerif", 1, 15), java.awt.Color.black)); // NOI18N

        tabPaneRFortran.setFont(new java.awt.Font("Tahoma", 1, 15)); // NOI18N

        chkRidgeRegression.setFont(new java.awt.Font("DejaVu Sans", 0, 15)); // NOI18N
        chkRidgeRegression.setText("Ridge Regression BLUP");

        chkBayeCpi.setFont(new java.awt.Font("DejaVu Sans", 0, 15)); // NOI18N
        chkBayeCpi.setText("Bayes Cpi");

        chkBayesB.setFont(new java.awt.Font("DejaVu Sans", 0, 15)); // NOI18N
        chkBayesB.setText("BayesB");

        chkBayesLasso.setFont(new java.awt.Font("DejaVu Sans", 0, 15)); // NOI18N
        chkBayesLasso.setText("Bayes LASSO");

        chkRandomForest.setFont(new java.awt.Font("DejaVu Sans", 0, 15)); // NOI18N
        chkRandomForest.setText("Random Forest");

        chkKinshipGauss.setFont(new java.awt.Font("DejaVu Sans", 0, 15)); // NOI18N
        chkKinshipGauss.setText("Kinship Gauss");

        javax.swing.GroupLayout rPanelLayout = new javax.swing.GroupLayout(rPanel);
        rPanel.setLayout(rPanelLayout);
        rPanelLayout.setHorizontalGroup(
            rPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(rPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(rPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(chkRidgeRegression)
                    .addComponent(chkBayeCpi)
                    .addComponent(chkBayesB)
                    .addComponent(chkBayesLasso)
                    .addComponent(chkRandomForest)
                    .addComponent(chkKinshipGauss))
                .addGap(0, 0, 0))
        );
        rPanelLayout.setVerticalGroup(
            rPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(rPanelLayout.createSequentialGroup()
                .addComponent(chkRidgeRegression)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addComponent(chkBayeCpi)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addComponent(chkBayesB)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addComponent(chkBayesLasso)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addComponent(chkRandomForest)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addComponent(chkKinshipGauss)
                .addGap(0, 0, Short.MAX_VALUE))
        );

        tabPaneRFortran.addTab("R   ", rPanel);

        chkFtBayesA.setFont(new java.awt.Font("DejaVu Sans", 0, 15)); // NOI18N
        chkFtBayesA.setText("BayesA");

        chkFtRidgeRegression.setFont(new java.awt.Font("DejaVu Sans", 0, 15)); // NOI18N
        chkFtRidgeRegression.setText("Ridge Regression");

        javax.swing.GroupLayout fortranPanelLayout = new javax.swing.GroupLayout(fortranPanel);
        fortranPanel.setLayout(fortranPanelLayout);
        fortranPanelLayout.setHorizontalGroup(
            fortranPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(fortranPanelLayout.createSequentialGroup()
                .addGap(18, 18, 18)
                .addComponent(chkFtBayesA)
                .addContainerGap(198, Short.MAX_VALUE))
            .addGroup(fortranPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                .addGroup(fortranPanelLayout.createSequentialGroup()
                    .addGap(16, 16, 16)
                    .addComponent(chkFtRidgeRegression)
                    .addContainerGap(132, Short.MAX_VALUE)))
        );
        fortranPanelLayout.setVerticalGroup(
            fortranPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(fortranPanelLayout.createSequentialGroup()
                .addGap(49, 49, 49)
                .addComponent(chkFtBayesA)
                .addContainerGap(111, Short.MAX_VALUE))
            .addGroup(fortranPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                .addGroup(fortranPanelLayout.createSequentialGroup()
                    .addGap(17, 17, 17)
                    .addComponent(chkFtRidgeRegression)
                    .addContainerGap(143, Short.MAX_VALUE)))
        );

        tabPaneRFortran.addTab("Fortran", fortranPanel);

        javax.swing.GroupLayout methodPanelLayout = new javax.swing.GroupLayout(methodPanel);
        methodPanel.setLayout(methodPanelLayout);
        methodPanelLayout.setHorizontalGroup(
            methodPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(tabPaneRFortran, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
        );
        methodPanelLayout.setVerticalGroup(
            methodPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(methodPanelLayout.createSequentialGroup()
                .addComponent(tabPaneRFortran, javax.swing.GroupLayout.PREFERRED_SIZE, 222, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );

        additionalParPanel.setBorder(javax.swing.BorderFactory.createTitledBorder(null, "Additional Parameters", javax.swing.border.TitledBorder.CENTER, javax.swing.border.TitledBorder.DEFAULT_POSITION, new java.awt.Font("DejaVu Sans", 1, 15), java.awt.Color.black)); // NOI18N

        jPanel2.setBorder(javax.swing.BorderFactory.createTitledBorder(null, "Random Forest", javax.swing.border.TitledBorder.CENTER, javax.swing.border.TitledBorder.DEFAULT_POSITION, new java.awt.Font("DejaVu Sans", 1, 15), java.awt.Color.black)); // NOI18N

        jLabel13.setFont(new java.awt.Font("DejaVu Sans", 0, 15)); // NOI18N
        jLabel13.setText("Forests");

        forestsSpin.setFont(new java.awt.Font("DejaVu Sans", 0, 15)); // NOI18N
        forestsSpin.setModel(new javax.swing.SpinnerNumberModel(Integer.valueOf(10), Integer.valueOf(10), null, Integer.valueOf(1)));
        forestsSpin.setRequestFocusEnabled(false);

        javax.swing.GroupLayout jPanel2Layout = new javax.swing.GroupLayout(jPanel2);
        jPanel2.setLayout(jPanel2Layout);
        jPanel2Layout.setHorizontalGroup(
            jPanel2Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel2Layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(jLabel13)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, 26, Short.MAX_VALUE)
                .addComponent(forestsSpin, javax.swing.GroupLayout.PREFERRED_SIZE, 65, javax.swing.GroupLayout.PREFERRED_SIZE))
        );
        jPanel2Layout.setVerticalGroup(
            jPanel2Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel2Layout.createSequentialGroup()
                .addGap(5, 5, 5)
                .addGroup(jPanel2Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jLabel13, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(forestsSpin, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );

        jPanel3.setBorder(javax.swing.BorderFactory.createTitledBorder(null, "Bayes", javax.swing.border.TitledBorder.CENTER, javax.swing.border.TitledBorder.DEFAULT_POSITION, new java.awt.Font("DejaVu Sans", 1, 15), java.awt.Color.black)); // NOI18N

        lblRounds2.setFont(new java.awt.Font("DejaVu Sans", 0, 15)); // NOI18N
        lblRounds2.setText("Rounds");

        roundsSpin.setFont(new java.awt.Font("DejaVu Sans", 0, 15)); // NOI18N
        roundsSpin.setModel(new javax.swing.SpinnerNumberModel(Integer.valueOf(1000), Integer.valueOf(1), null, Integer.valueOf(1)));

        lblRounds3.setFont(new java.awt.Font("DejaVu Sans", 0, 15)); // NOI18N
        lblRounds3.setText("Burnin");

        burninSpin.setFont(new java.awt.Font("DejaVu Sans", 0, 15)); // NOI18N
        burninSpin.setModel(new javax.swing.SpinnerNumberModel(Integer.valueOf(100), Integer.valueOf(1), null, Integer.valueOf(1)));

        lblother1.setFont(new java.awt.Font("DejaVu Sans", 0, 15)); // NOI18N
        lblother1.setText("Thinning");

        thiningSpin.setFont(new java.awt.Font("DejaVu Sans", 0, 15)); // NOI18N
        thiningSpin.setModel(new javax.swing.SpinnerNumberModel(Integer.valueOf(20), Integer.valueOf(1), null, Integer.valueOf(1)));

        javax.swing.GroupLayout jPanel3Layout = new javax.swing.GroupLayout(jPanel3);
        jPanel3.setLayout(jPanel3Layout);
        jPanel3Layout.setHorizontalGroup(
            jPanel3Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel3Layout.createSequentialGroup()
                .addGroup(jPanel3Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(lblRounds2)
                    .addComponent(lblRounds3)
                    .addComponent(lblother1))
                .addGap(18, 18, 18)
                .addGroup(jPanel3Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(burninSpin)
                    .addComponent(roundsSpin, javax.swing.GroupLayout.DEFAULT_SIZE, 73, Short.MAX_VALUE)
                    .addComponent(thiningSpin))
                .addContainerGap())
        );
        jPanel3Layout.setVerticalGroup(
            jPanel3Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel3Layout.createSequentialGroup()
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addGroup(jPanel3Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(lblRounds2)
                    .addComponent(roundsSpin, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(jPanel3Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(lblRounds3)
                    .addComponent(burninSpin, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(jPanel3Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(lblother1)
                    .addComponent(thiningSpin, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)))
        );

        jPanel4.setBorder(javax.swing.BorderFactory.createTitledBorder(null, "Cross Validation\n", javax.swing.border.TitledBorder.CENTER, javax.swing.border.TitledBorder.DEFAULT_POSITION, new java.awt.Font("DejaVu Sans", 1, 15), java.awt.Color.black)); // NOI18N

        jLabel7.setFont(new java.awt.Font("DejaVu Sans", 0, 15)); // NOI18N
        jLabel7.setText("Replication ");

        jLabel8.setFont(new java.awt.Font("DejaVu Sans", 0, 15)); // NOI18N
        jLabel8.setText("Fold");

        replicationSpin.setFont(new java.awt.Font("DejaVu Sans", 0, 15)); // NOI18N
        replicationSpin.setModel(new javax.swing.SpinnerNumberModel(Integer.valueOf(1), Integer.valueOf(1), null, Integer.valueOf(1)));

        foldSpin.setFont(new java.awt.Font("DejaVu Sans", 0, 15)); // NOI18N
        foldSpin.setModel(new javax.swing.SpinnerNumberModel(Integer.valueOf(2), Integer.valueOf(2), null, Integer.valueOf(1)));

        javax.swing.GroupLayout jPanel4Layout = new javax.swing.GroupLayout(jPanel4);
        jPanel4.setLayout(jPanel4Layout);
        jPanel4Layout.setHorizontalGroup(
            jPanel4Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel4Layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(jPanel4Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(jPanel4Layout.createSequentialGroup()
                        .addComponent(jLabel8, javax.swing.GroupLayout.PREFERRED_SIZE, 44, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addGap(40, 40, 40)
                        .addComponent(foldSpin, javax.swing.GroupLayout.DEFAULT_SIZE, 46, Short.MAX_VALUE))
                    .addGroup(jPanel4Layout.createSequentialGroup()
                        .addComponent(jLabel7)
                        .addGap(8, 8, 8)
                        .addComponent(replicationSpin))))
        );
        jPanel4Layout.setVerticalGroup(
            jPanel4Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel4Layout.createSequentialGroup()
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addGroup(jPanel4Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jLabel7)
                    .addComponent(replicationSpin, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(jPanel4Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                    .addComponent(jLabel8)
                    .addComponent(foldSpin, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(7, 7, 7))
        );

        jPanel5.setBorder(javax.swing.BorderFactory.createTitledBorder(null, "Processor", javax.swing.border.TitledBorder.CENTER, javax.swing.border.TitledBorder.DEFAULT_POSITION, new java.awt.Font("DejaVu Sans", 1, 15), java.awt.Color.black)); // NOI18N

        jLabel14.setFont(new java.awt.Font("DejaVu Sans", 0, 15)); // NOI18N
        jLabel14.setText("Cores");

        coresSpin.setFont(new java.awt.Font("DejaVu Sans", 0, 15)); // NOI18N
        coresSpin.setModel(new javax.swing.SpinnerNumberModel(Integer.valueOf(1), Integer.valueOf(1), null, Integer.valueOf(1)));
        coresSpin.setRequestFocusEnabled(false);
        coresSpin.setValue(1);

        javax.swing.GroupLayout jPanel5Layout = new javax.swing.GroupLayout(jPanel5);
        jPanel5.setLayout(jPanel5Layout);
        jPanel5Layout.setHorizontalGroup(
            jPanel5Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel5Layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(jLabel14)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addComponent(coresSpin, javax.swing.GroupLayout.PREFERRED_SIZE, 65, javax.swing.GroupLayout.PREFERRED_SIZE))
        );
        jPanel5Layout.setVerticalGroup(
            jPanel5Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel5Layout.createSequentialGroup()
                .addGroup(jPanel5Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(coresSpin, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jLabel14))
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );

        javax.swing.GroupLayout additionalParPanelLayout = new javax.swing.GroupLayout(additionalParPanel);
        additionalParPanel.setLayout(additionalParPanelLayout);
        additionalParPanelLayout.setHorizontalGroup(
            additionalParPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(additionalParPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addComponent(jPanel3, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addGap(18, 18, 18)
                .addGroup(additionalParPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                    .addComponent(jPanel2, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(jPanel5, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                .addGap(18, 18, 18)
                .addComponent(jPanel4, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
        );
        additionalParPanelLayout.setVerticalGroup(
            additionalParPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(jPanel4, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
            .addGroup(additionalParPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING, false)
                .addGroup(additionalParPanelLayout.createSequentialGroup()
                    .addComponent(jPanel2, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                    .addComponent(jPanel5, javax.swing.GroupLayout.PREFERRED_SIZE, 61, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addComponent(jPanel3, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
        );

        btnStart.setFont(new java.awt.Font("DejaVu Sans", 1, 15)); // NOI18N
        btnStart.setText("Start");

        btnCancel.setFont(new java.awt.Font("DejaVu Sans", 1, 15)); // NOI18N
        btnCancel.setText("Cancel");

        jPanel1.setBorder(javax.swing.BorderFactory.createTitledBorder(null, "Data Subset\n", javax.swing.border.TitledBorder.CENTER, javax.swing.border.TitledBorder.DEFAULT_POSITION, new java.awt.Font("DejaVu Sans", 1, 15), java.awt.Color.black)); // NOI18N

        jLabel10.setFont(new java.awt.Font("DejaVu Sans", 0, 15)); // NOI18N
        jLabel10.setText("Percentage(%) of missing markers");

        perMissingMarker.setFont(new java.awt.Font("DejaVu Sans", 0, 15)); // NOI18N
        perMissingMarker.setModel(new javax.swing.SpinnerNumberModel(10, 0, 100, 1));

        jLabel11.setFont(new java.awt.Font("DejaVu Sans", 0, 15)); // NOI18N
        jLabel11.setText("PIC value");

        picSpin.setFont(new java.awt.Font("DejaVu Sans", 0, 15)); // NOI18N
        picSpin.setModel(new javax.swing.SpinnerNumberModel(0.0d, 0.0d, 1.0d, 0.1d));

        jLabel12.setFont(new java.awt.Font("DejaVu Sans", 0, 15)); // NOI18N
        jLabel12.setText("Minor allale frequency (MAF)");

        mafSpin.setFont(new java.awt.Font("DejaVu Sans", 0, 15)); // NOI18N
        mafSpin.setModel(new javax.swing.SpinnerNumberModel(0.0d, 0.0d, 1.0d, 0.1d));

        javax.swing.GroupLayout jPanel1Layout = new javax.swing.GroupLayout(jPanel1);
        jPanel1.setLayout(jPanel1Layout);
        jPanel1Layout.setHorizontalGroup(
            jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel1Layout.createSequentialGroup()
                .addGap(0, 0, Short.MAX_VALUE)
                .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(jLabel10)
                    .addComponent(jLabel11)
                    .addComponent(jLabel12))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                    .addComponent(perMissingMarker, javax.swing.GroupLayout.PREFERRED_SIZE, 53, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(picSpin)
                    .addComponent(mafSpin)))
        );
        jPanel1Layout.setVerticalGroup(
            jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel1Layout.createSequentialGroup()
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jLabel10)
                    .addComponent(perMissingMarker, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jLabel11)
                    .addComponent(picSpin, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jLabel12)
                    .addComponent(mafSpin, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)))
        );

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(getContentPane());
        getContentPane().setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(layout.createSequentialGroup()
                        .addContainerGap()
                        .addComponent(additionalParPanel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addGap(18, 18, 18)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                            .addComponent(btnStart, javax.swing.GroupLayout.PREFERRED_SIZE, 165, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(btnCancel, javax.swing.GroupLayout.PREFERRED_SIZE, 165, javax.swing.GroupLayout.PREFERRED_SIZE)))
                    .addGroup(layout.createSequentialGroup()
                        .addGap(12, 12, 12)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                            .addComponent(fileNameSelectingPanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(traitPanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                            .addComponent(methodPanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(jPanel1, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))))
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addGap(1, 1, 1)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(fileNameSelectingPanel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addGap(5, 5, 5)
                        .addComponent(traitPanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(methodPanel, javax.swing.GroupLayout.PREFERRED_SIZE, 256, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(jPanel1, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                    .addComponent(additionalParPanel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(btnStart, javax.swing.GroupLayout.PREFERRED_SIZE, 65, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addComponent(btnCancel, javax.swing.GroupLayout.PREFERRED_SIZE, 63, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addGap(12, 12, 12)))
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );

        pack();
    }// </editor-fold>//GEN-END:initComponents

    /*
     * Mouse click listeners for selecting the traits on double click
     */
    private void jList1MouseClicked(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_jList1MouseClicked
        if (evt.getClickCount() == 2) {
            int index = jList1.locationToIndex(evt.getPoint());
            if (index >= 0) {
                Object o = jList1.getModel().getElementAt(index);
                if (resultList2.contains(o)) {
                    //donothing
                } else {
                    resultList2.addElement(o);
                    jList2.setModel(resultList2);
                    resultList1.removeElement(o);
                    jList1.setModel(resultList1);
                }
            }
        }
    }//GEN-LAST:event_jList1MouseClicked

    private void jList2MouseClicked(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_jList2MouseClicked
        if (evt.getClickCount() == 2) {
            int index = jList2.locationToIndex(evt.getPoint());
            if (index >= 0) {
                Object o = jList2.getModel().getElementAt(index);
                if (resultList1.contains(o)) {
                    //donothing
                } else {
                    resultList1.addElement(o);
                    jList1.setModel(resultList1);
                    resultList2.removeElement(o);
                    jList2.setModel(resultList2);
                }
            }
        }
    }//GEN-LAST:event_jList2MouseClicked

    /**
     * @param args the command line arguments
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
            java.util.logging.Logger.getLogger(Analysis.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        }
        //</editor-fold>
        //</editor-fold>

    }
    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JPanel additionalParPanel;
    private javax.swing.JButton btnCancel;
    private javax.swing.JButton btnStart;
    private javax.swing.JButton btndoubleiparrow;
    private javax.swing.JButton btndoubleremovearrow;
    private javax.swing.JButton btnsingleiparrow;
    private javax.swing.JButton btnsingleremovearrow;
    private javax.swing.JSpinner burninSpin;
    private javax.swing.JCheckBox chkBayeCpi;
    private javax.swing.JCheckBox chkBayesB;
    private javax.swing.JCheckBox chkBayesLasso;
    private javax.swing.JCheckBox chkFtBayesA;
    private javax.swing.JCheckBox chkFtRidgeRegression;
    private javax.swing.JCheckBox chkKinshipGauss;
    private javax.swing.JCheckBox chkRandomForest;
    private javax.swing.JCheckBox chkRidgeRegression;
    private javax.swing.JSpinner coresSpin;
    private javax.swing.JPanel fileNameSelectingPanel;
    private javax.swing.JSpinner foldSpin;
    private javax.swing.JSpinner forestsSpin;
    private javax.swing.JPanel fortranPanel;
    private javax.swing.JComboBox genoCombo;
    private javax.swing.JLabel jLabel1;
    private javax.swing.JLabel jLabel10;
    private javax.swing.JLabel jLabel11;
    private javax.swing.JLabel jLabel12;
    private javax.swing.JLabel jLabel13;
    private javax.swing.JLabel jLabel14;
    private javax.swing.JLabel jLabel7;
    private javax.swing.JLabel jLabel8;
    private javax.swing.JLabel jLabel9;
    private javax.swing.JList jList1;
    private javax.swing.JList jList2;
    private javax.swing.JPanel jPanel1;
    private javax.swing.JPanel jPanel2;
    private javax.swing.JPanel jPanel3;
    private javax.swing.JPanel jPanel4;
    private javax.swing.JPanel jPanel5;
    private javax.swing.JScrollPane jScrollPane1;
    private javax.swing.JScrollPane jScrollPane2;
    private javax.swing.JLabel lblPheno;
    private javax.swing.JLabel lblRounds2;
    private javax.swing.JLabel lblRounds3;
    private javax.swing.JLabel lblother1;
    private javax.swing.JSpinner mafSpin;
    private javax.swing.JPanel methodPanel;
    private javax.swing.JSpinner perMissingMarker;
    private javax.swing.JComboBox phenoCombo;
    private javax.swing.JSpinner picSpin;
    private javax.swing.JComboBox populationCombo;
    private javax.swing.JPanel rPanel;
    private javax.swing.JSpinner replicationSpin;
    private javax.swing.JSpinner roundsSpin;
    private javax.swing.JTabbedPane tabPaneRFortran;
    private javax.swing.JSpinner thiningSpin;
    private javax.swing.JPanel traitPanel;
    // End of variables declaration//GEN-END:variables
}
