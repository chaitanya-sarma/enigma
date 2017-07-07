/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.icrisat.genomicSelection;

import java.awt.Font;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.PrintStream;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JLayer;
import javax.swing.JPanel;
import javax.swing.SwingWorker;
import org.apache.commons.lang3.StringUtils;

/**
 *
 * @author Mohan-60686
 */
public class DataSummary extends JDialog {

    JPanel dataPanel;
    JComboBox genoCombo;
    JComboBox phenoCombo;
    JLabel lblGenotype;
    JLabel lblPhenotype;
    JCheckBox chkMissingPrecent;
    JCheckBox chkbxPICValue;
    JCheckBox chkMaf;
    JCheckBox chkDataSummary;
    JButton btnStart;
    JButton btnCancel;
    private final WaitLayerUI layerUI = new WaitLayerUI();
    Long processStartTime;
    String resultEscp;
    //for running rscript dcn.r
    int perMissing = -1;
    int picValue = -1;
    int maf = -1;
    int engine = 1;  //1 for R and 2 for Fotran
    int calculatesummaryPheno = -1;
    Process processCallR;
    Process processCallFortran;
    boolean processCalled;
    String presentWorkingDir = System.getProperty("user.dir");
    String genoSummaryCsvName, genoHtmlName, phenoHtmlName;
    int runningProcess = 0; //for getting the which process is running ,used when cancel button clicked to kill process
    int exitGenoStatus, exitPhenoStatus;

    public DataSummary(java.awt.Frame parent, boolean modal) {
        super(parent, modal);
        System.out.println("Analysys");
        this.setSize(700, 250);
        dataPanel = new JPanel();
        dataPanel.setSize(700, 250);
        this.setIconImage(Toolkit.getDefaultToolkit().getImage(getClass().getResource("/org/icrisat/genomicSelection/image/log_icrisat.png")));
        init();
        Constant.summaryResultList.clear();
        Constant.summaryResultList.add(0, null);
        Constant.summaryResultList.add(1, null);
        Constant.summaryResultList.add(2, null);
        this.add(new JLayer<>(dataPanel, layerUI));
        //adding genotype and phenotype files to comboBox
        for (String genoList : Constant.genoList) {
            genoCombo.addItem(genoList);
        }
        genoCombo.setSelectedIndex(1);
        for (String phenoList : Constant.phenoList) {
            phenoCombo.addItem(phenoList);
        }
        phenoCombo.setSelectedIndex(1);

        //actionListeners for Genotype and Phenotype combobox 
        genoCombo.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if ("Select".equals(genoCombo.getSelectedItem().toString())) {
                    chkMaf.setEnabled(false);
                    chkMissingPrecent.setEnabled(false);
                    chkbxPICValue.setEnabled(false);
                    chkMaf.setSelected(false);
                    chkMissingPrecent.setSelected(false);
                    chkbxPICValue.setSelected(false);
                } else {
                    chkMaf.setEnabled(true);
                    chkMissingPrecent.setEnabled(true);
                    chkbxPICValue.setEnabled(true);
                    chkMaf.setSelected(true);
                    chkMissingPrecent.setSelected(true);
                    chkbxPICValue.setSelected(true);
                }
            }
        });

        phenoCombo.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if ("Select".equals(phenoCombo.getSelectedItem().toString())) {
                    chkDataSummary.setEnabled(false);
                    chkDataSummary.setSelected(false);
                } else {
                    chkDataSummary.setEnabled(true);
                    chkDataSummary.setSelected(true);
                }
            }
        });

        //setting the all the checkboxs enabled by default 
        chkbxPICValue.setSelected(true);
        chkMissingPrecent.setSelected(true);
        chkMaf.setSelected(true);
        chkDataSummary.setSelected(true);

        OSCheck.OSType os = OSCheck.getOperatingSystemType();
        if ("Window".equals(os.toString())) {
            String resultRemo = StringUtils.removeEnd(Constant.resultdirectory, "\\");
            resultEscp = "\"" + resultRemo + "\"";
        }
        if ("".equals(os.toString())) {
            String resultRemo = StringUtils.removeEnd(Constant.resultdirectory, "/");
            resultEscp = resultRemo;
        }
        btnStart.addActionListener(new ActionListener() {

            @Override
            public void actionPerformed(ActionEvent e) {

                System.out.println("In Start button.");
                System.out.println("Still In Start button.");

                String genoFile = genoCombo.getSelectedItem().toString();
                String phenoFile = phenoCombo.getSelectedItem().toString();
                if (!"Select".equals(genoFile) || !"Select".equals(phenoFile)) {
                    if (chkbxPICValue.isSelected() || chkMissingPrecent.isSelected() || chkMaf.isSelected() || chkDataSummary.isSelected()) {
                        processStartTime = new Date().getTime();
                        layerUI.start();
                        btnStart.setEnabled(false);
                        chkbxPICValue.setEnabled(false);
                        chkMissingPrecent.setEnabled(false);
                        chkMaf.setEnabled(false);
                        chkDataSummary.setEnabled(false);
                        genoCombo.setEnabled(false);
                        phenoCombo.setEnabled(false);
                        try {
                            System.setOut(new PrintStream(new FileOutputStream(Constant.resultdirectory + "log.txt")));
                        } catch (FileNotFoundException ex) {
                            Logger.getLogger(DataSummary.class.getName()).log(Level.SEVERE, null, ex);
                        }

                        SwingWorker<Object, Object[]> swingWorker = new SwingWorker<Object, Object[]>() {

                            @Override
                            protected Object doInBackground() throws Exception {
                                if (chkMissingPrecent.isSelected() == true) {
                                    perMissing = 1;
                                }
                                if (chkbxPICValue.isSelected() == true) {
                                    picValue = 1;
                                }
                                if (chkMaf.isSelected() == true) {
                                    maf = 1;
                                }
                                if (chkDataSummary.isSelected() == true) {
                                    calculatesummaryPheno = 1;
                                }

//*******************************************************************************************************
                                OSCheck.OSType os = OSCheck.getOperatingSystemType();
                                String rscriptPath = "";
                                if ("Window".equals(os.toString())) {
                                    BufferedReader reader = new BufferedReader(new FileReader(presentWorkingDir + "\\os.ini"));
                                    rscriptPath = reader.readLine();
                                }
                                if ("Linux".equals(os.toString())) {
                                    BufferedReader reader = new BufferedReader(new FileReader(presentWorkingDir + "/os.ini"));
                                    rscriptPath = reader.readLine();
                                }

                                if (rscriptPath != null) {
                                    if (Constant.engineSelection == 1) {
                                        engine = 1;
                                    }
                                    if (Constant.engineSelection == 2) {
                                        engine = 2;
                                    }
                                    Date date = new Date();
                                    SimpleDateFormat dateFormat = new SimpleDateFormat("mmss");
                                    genoSummaryCsvName = "Ge_summary" + dateFormat.format(date) + ".csv";

                                    Date date1 = new Date();
                                    SimpleDateFormat dateFormat1 = new SimpleDateFormat("mmss");
                                    genoHtmlName = "Ge_summary" + dateFormat1.format(date1) + ".htm";

                                    Date date2 = new Date();
                                    SimpleDateFormat dateFormat2 = new SimpleDateFormat("mmss");
                                    phenoHtmlName = "Ph_summary" + dateFormat2.format(date2) + ".htm";

                                    if (chkMaf.isSelected() == true || chkMissingPrecent.isSelected() == true || chkbxPICValue.isSelected() == true) {
                                        String cmdGeno = "";
                                        if ("Window".equals(os.toString())) {
                                            cmdGeno = rscriptPath + " \"" + presentWorkingDir + "\\dcN.R\"" + " " + resultEscp + " \"" + genoCombo.getSelectedItem().toString() + "\" " + genoHtmlName
                                                    + " " + perMissing + " " + picValue + " " + maf + " " + "-1" + " " + "2" + " " + "?" + " " + genoSummaryCsvName + " " + Constant.engineSelection + " " + phenoHtmlName + " " + phenoCombo.getSelectedItem().toString();// + " " + calculatesummaryPheno;
                                            System.out.println("dncr call  " + cmdGeno);
                                        }
                                        if ("Linux".equals(os.toString())) {
                                            cmdGeno = rscriptPath + " " + presentWorkingDir + "/dcN.R" + " " + resultEscp + " " + genoCombo.getSelectedItem().toString() + " " + genoHtmlName
                                                    + " " + perMissing + " " + picValue + " " + maf + " " + "-1" + " " + "2" + " " + "?" + " " + genoSummaryCsvName + " " + engine + " " + phenoHtmlName + " " + phenoCombo.getSelectedItem().toString();// + " " + calculatesummaryPheno;

                                        }
                                        runningProcess = 1;
                                        processCalled = true;
                                        processCallR = Runtime.getRuntime().exec(cmdGeno);
                                        exitGenoStatus = processCallR.waitFor();
                                        if (exitGenoStatus == 0) {
                                            Constant.summaryResultList.add(0, genoSummaryCsvName);
                                            Constant.summaryResultList.add(1, genoHtmlName);
                                        }
                                    }
                                    if (chkDataSummary.isSelected() == true) {
                                        String cmdPheno = "";
                                        if ("Window".equals(os.toString())) {
                                            cmdPheno = rscriptPath + " \"" + presentWorkingDir + "\\PhenoS.R\"" + " " + resultEscp + " \"" + genoCombo.getSelectedItem().toString() + "\" " + genoHtmlName
                                                    + " " + perMissing + " " + picValue + " " + maf + " " + "-1" + " " + "2" + " " + "?" + " " + genoSummaryCsvName + " " + Constant.engineSelection + " " + phenoHtmlName + " " + phenoCombo.getSelectedItem().toString();//+ " " + calculatesummaryPheno;
                                        }
                                        if ("Linux".equals(os.toString())) {
                                            cmdPheno = rscriptPath + " " + presentWorkingDir + "/PhenoS.R" + " " + resultEscp + " " + genoCombo.getSelectedItem().toString() + " " + genoHtmlName
                                                    + " " + perMissing + " " + picValue + " " + maf + " " + "-1" + " " + "2" + " " + "?" + " " + genoSummaryCsvName + " " + engine + " " + phenoHtmlName + " " + phenoCombo.getSelectedItem().toString();// + " " + calculatesummaryPheno;
                                        }
                                        processCalled = true;
                                        runningProcess = 2;
                                        processCallFortran = Runtime.getRuntime().exec(cmdPheno);
                                        exitPhenoStatus = processCallFortran.waitFor();
                                        if (exitPhenoStatus == 0) {
                                            Constant.summaryResultList.add(2, phenoHtmlName);
                                        }
                                    }
                                    dispose();
                                }
                                return null;
                            }
                        };
                        swingWorker.execute();
                    }
                }
            }
        });
        //cancel button action lisneter 
        btnCancel.addActionListener(new ActionListener() {

            @Override
            public void actionPerformed(ActionEvent e) {
                try {
                    close();
                } catch (Exception ee) {
                    System.out.println(ee);
                }
            }
        });
        //adding windowListener for closing of dialog box
        this.addWindowListener(
                new WindowAdapter() {
                    @Override
                    public void windowClosing(WindowEvent e) {
                        try {
                            close();
                        } catch (Exception ee) {
                            System.out.println(ee);
                        }
                    }
                });
    }

    private void init() {
        // dataPanel = new JPanel();
        dataPanel.setLayout(null);
        genoCombo = new javax.swing.JComboBox();
        genoCombo.addItem("Select");
        phenoCombo = new javax.swing.JComboBox();
        phenoCombo.addItem("Select");

        lblGenotype = new javax.swing.JLabel("Select a genotype file");
        lblPhenotype = new javax.swing.JLabel("Select a phenotype file");

        chkMissingPrecent = new javax.swing.JCheckBox("Calculate %  missing");
        chkbxPICValue = new javax.swing.JCheckBox("Calculate PIC value");
        chkMaf = new javax.swing.JCheckBox("Calculate MAF");
        chkDataSummary = new javax.swing.JCheckBox("Calculate Data Summary");

        btnStart = new javax.swing.JButton("Start");
        btnCancel = new javax.swing.JButton("Cancel");
        setTitle("Data summary");

        lblGenotype.setFont(new Font("DejaVu Sans", 0, 15));
        lblGenotype.setBounds(20, 10, 200, 20);

        lblPhenotype.setFont(new Font("DejaVu Sans", 0, 15));
        lblPhenotype.setBounds(370, 10, 200, 20);

        genoCombo.setFont(new Font("DejaVu Sans", 0, 15));
        genoCombo.setBounds(180, 10, 150, 20);
        phenoCombo.setFont(new Font("DejaVu Sans", 0, 15));
        phenoCombo.setBounds(530, 10, 150, 20);

        chkMissingPrecent.setFont(new Font("DejaVu Sans", 0, 15));
        chkMissingPrecent.setBounds(20, 50, 200, 20);
        chkbxPICValue.setFont(new Font("DejaVu Sans", 0, 15));
        chkbxPICValue.setBounds(20, 90, 180, 20);
        chkMaf.setFont(new Font("DejaVu Sans", 0, 15));
        chkMaf.setBounds(20, 130, 180, 20);
        chkDataSummary.setFont(new Font("DejaVuSans", 0, 15));
        chkDataSummary.setBounds(370, 50, 200, 20);

        btnStart.setFont(new Font("DejaVu Sans", 0, 15));
        btnStart.setBounds(270, 170, 80, 25);

        btnCancel.setFont(new Font("DejaVu Sans", 0, 15));
        btnCancel.setBounds(360, 170, 90, 25);

        dataPanel.add(lblGenotype);
        dataPanel.add(lblPhenotype);
        dataPanel.add(genoCombo);
        dataPanel.add(phenoCombo);
        dataPanel.add(chkbxPICValue);
        dataPanel.add(chkMissingPrecent);
        dataPanel.add(chkMaf);
        dataPanel.add(chkDataSummary);
        dataPanel.add(btnStart);
        dataPanel.add(btnCancel);
        add(dataPanel);
    }

    //closing the Summary dialog box when clicked on cancel button, window close
    public void close() {
        if (runningProcess == 1 || exitGenoStatus != 0) {
            processCallR.destroy();
            chkDataSummary.setSelected(false); //setting for disabling pheno summary
        }
        if (runningProcess == 2 || exitPhenoStatus != 0) {
            processCallFortran.destroy();
        }

        Constant.summaryCancled = true;
        Constant.summaryResultList.clear();
        Constant.summaryResultList.add(0, null);
        Constant.summaryResultList.add(1, null);
        Constant.summaryResultList.add(2, null);

        File genoSumm = new File(Constant.resultdirectory + genoSummaryCsvName);
        boolean exists = genoSumm.exists();
        if (exists == true) {
            genoSumm.delete();
        }
        File genoSummHtml = new File(Constant.resultdirectory + genoHtmlName);
        boolean existsGenoHtml = genoSummHtml.exists();
        if (existsGenoHtml == true) {
            GsMethods.removePngFromHtml(Constant.resultdirectory + genoHtmlName);
            genoSummHtml.delete();
        }

        File phenoSummHtml = new File(Constant.resultdirectory + phenoHtmlName);
        boolean existsPhenoHtml = phenoSummHtml.exists();
        if (existsPhenoHtml == true) {
            GsMethods.removePngFromHtml(Constant.resultdirectory + phenoHtmlName);
            phenoSummHtml.delete();
        }
        dispose();
    }
}
