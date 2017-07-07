package org.icrisat.genomicSelection;

import com.csvreader.CsvReader;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JLabel;
import javax.swing.JLayer;
import javax.swing.JPanel;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.io.FileWriter;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.PrintStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.SimpleDateFormat;

import java.util.Date;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.JOptionPane;
import javax.swing.SwingWorker;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang3.StringUtils;

/**
 * @author Mohan After selecting the required data in the Analysis class getting
 * the all the data of file names and parameter values Running the script based
 * on the selection of methods, which generates the Html output in result
 * directory
 */
public class ProgressDisplay extends javax.swing.JDialog {

    JPanel mainPanel;
    JLabel lblPerformAnalysis, lblTraitName, lblUsingmethid, lblMethodName;
    JCheckBox chkRidgeRegression, chkBayesCpi, chkBayesB, chkBayesLasso, chkRandomForest, chkKinshipGauss, chkRRFotron, chkBayesAFotron;
    JButton btnCancel;
    private WaitLayerUI1 layerUI1 = new WaitLayerUI1();
    Process p1, p2, p3, p4, p5, p6, runA, abRun, cRun, fotRR, fotBA;
    boolean exit = false;
    int processNum = 0;
    String resultRemo = "";
    String resultDir = "";
    String genotypeFileName = "";
    String phenotypeFileName = "";
    String populationFileName = "";
    String resultFileName = "";
    Double missingPercentage;
    Double minorAllaleFrequnce;
    Double pic;
    int methodName;
    int nofCores;
    String missingCh = "";
    int ridgeReg = -1;
    int kinshipGauss = -1;
    int bayesLasso = -1;
    int bayesCpi = -1;
    int bayesB = -1;
    int randomForest = -1;
    int traitNum;
    int dataType;
    int ridgeFotran = -1;
    int bayesAFortran = -1;
    String alleleSeperator = ":";
    String tmp1 = "";
    String tmp2 = "";
    String tmp3 = "";
    String crossValidation;
    String BayesParameters;
    int RandomForestParameter;
    boolean checkAExe = true;
    int engineSelection = 1;
    int x = 40;
    int y = 60;
    String resultDir1, resultEscp, htmlFileName;
    String inputFolder = "";
    String outputFolder = "";
    String genotypeFile = "";
    String MarkerNames = "Yes";
    String pedigreeFile = "";
    String weightFile = "None";
    String numberOfTraits = "1";
    String individulas = "";
    String inversionOfG = "Ordinary";
    String diagonalFudgeFactor = "0.001";
    String missingPercTreshold = "";
    String MafTreshold = "";
    String picTreshold = "";
    String makeG = "No";
    String makeGinverse = "No";
    String makeA = "No";
    String makeAinverse = "No";
    String makeH = "No";
    String makeHinverse = "No";
    String summaryStatistics = "Yes";
    String dominantMarkers = "No";
    String vanRadenGFullMat = "No";
    String vanRadenGIJA = "No";
    String inverseVanRadenGFullMat = "No";
    String inverseVanRadenGIJA = "No";
    String amatFullMat = "No";
    String amatIJA = "No";
    String inverseAmatFullMat = "No";
    String inverseAmatIJA = "No";
    String hmatFullMat = "No";
    String hmatIJA = "No";
    String inverseHmatFullMat = "No";
    String inverseHmatIJA = "No";
    String presentWorkingDir = System.getProperty("user.dir");
    String convertedFileName;
    long processANStartTime;
    long processANStopTime;
    //start region
    //called after selecting the required option in Analysis dialog box 
    //adding the action listener event for cancel button 

    public ProgressDisplay(java.awt.Frame parent, boolean modal) {

        super(parent, modal);
        this.setSize(400, 350);
        this.setIconImage(Toolkit.getDefaultToolkit().getImage(getClass().getResource("/org/icrisat/genomicSelection/image/log_icrisat.png")));
        this.setLocationRelativeTo(parent);
        mainPanel = new javax.swing.JPanel();
        init();

        //If fotron engine selected
        if (Constant.RidgeRegressionFortran == true) {
            y = y + 30;
            chkRRFotron.setBounds(x, y, 200, 25);
            chkRRFotron.setEnabled(false);
            ridgeFotran = 1;
        }
        if (Constant.BayeAFortran == true) {
            y = y + 30;
            chkBayesAFotron.setBounds(x, y, 200, 25);
            chkBayesAFotron.setEnabled(false);
            bayesAFortran = 1;
        }
        //If R engine selected 
        if (Constant.RidgeReg == true) {
            y = y + 30;
            chkRidgeRegression.setBounds(x, y, 200, 25);
            chkRidgeRegression.setEnabled(false);
            ridgeReg = 1;
        }
        if (Constant.BayesCpi == true) {
            y = y + 30;
            chkBayesCpi.setBounds(x, y, 200, 25);
            chkBayesCpi.setEnabled(false);
            bayesCpi = 1;
        }
        if (Constant.BayesB == true) {
            y = y + 30;
            chkBayesB.setBounds(x, y, 200, 25);
            chkBayesB.setEnabled(false);
            bayesB = 1;
        }
        if (Constant.BayesLasso == true) {
            y = y + 30;
            chkBayesLasso.setBounds(x, y, 200, 25);
            chkBayesLasso.setEnabled(false);
            bayesLasso = 1;
        }

        if (Constant.RandomeForest == true) {
            y = y + 30;
            chkRandomForest.setBounds(x, y, 200, 25);
            chkRandomForest.setEnabled(false);
            randomForest = 1;
        }
        if (Constant.KinshipGauss == true) {
            y = y + 30;
            chkKinshipGauss.setBounds(x, y, 200, 25);
            chkKinshipGauss.setEnabled(false);
            kinshipGauss = 1;
        }

        //emptying before adding the file names
        Constant.resultGenratedList.clear();

        add(new JLayer<>(mainPanel, layerUI1));
        layerUI1.start();
        processANStartTime = System.currentTimeMillis();
        SwingWorker<Object, Object[]> ss = new SwingWorker<Object, Object[]>() {
            @Override
            protected Object doInBackground() throws Exception {
                String presentWorkingDirectory = System.getProperty("user.dir");
                OSCheck.OSType os = OSCheck.getOperatingSystemType();
                BufferedReader reader = null;
                if ("Window".equals(os.toString())) {
                    reader = new BufferedReader(new FileReader(presentWorkingDirectory + "\\os.ini"));
                }
                if ("Linux".equals(os.toString())) {
                    reader = new BufferedReader(new FileReader(presentWorkingDirectory + "/os.ini"));
                }
                String rscriptPath = reader.readLine();
                if (rscriptPath != null) {
                    //adding the rscript log message to log file 
                    System.setOut(new PrintStream(new FileOutputStream(Constant.resultdirectory + "log.txt")));

                    if ("Window".equals(os.toString())) {
                        resultRemo = StringUtils.removeEnd(Constant.resultdirectory, "\\");
                        resultDir = "\"" + resultRemo + "\"";
                    }
                    if ("Linux".equals(os.toString())) {
                        resultRemo = StringUtils.removeEnd(Constant.resultdirectory, "/");
                        resultDir = resultRemo;
                    }

                    genotypeFileName = Constant.geno_Analysis;
                    phenotypeFileName = Constant.pheno_Analysis;
                    populationFileName = Constant.popu_Analysis;
                    int tempMissingValue = Constant.perMissing;
                    int r = 100;
                    missingPercentage = ((double) tempMissingValue) / r;
                    minorAllaleFrequnce = Constant.mafValue;
                    pic = Constant.picValue;
                    methodName = -1;
                    nofCores = Constant.noOfCPUs;
                    missingCh = "?";   //prevoius=-
                    dataType = 2;
                    alleleSeperator = ":";
                    tmp1 = "GsResTmp.txt";
                    tmp2 = "tmp1.txt";
                    tmp3 = "tmp2.txt";
                    crossValidation = Constant.crossValidationParameteres;
                    BayesParameters = Constant.bayesFortParameters;
                    RandomForestParameter = Constant.RandomForestAddParameteres;
                    engineSelection = Constant.engineSelection;
                    //itterating the traitList and performing each selected method      
                    for (int itl = 0; itl < Constant.traitList.size(); itl++) {
                        lblTraitName.setText(Constant.traitList.get(itl));
                        traitNum = itl + 1;//assing the present trait count 

                        //creating a random file name for displaying results
                        Date date = new Date();
                        SimpleDateFormat dateFormat = new SimpleDateFormat("mmss");
                        resultFileName = "res" + Constant.traitList.get(itl) + dateFormat.format(date) + ".htm"; //result file name  //assing the value to the 
                        if (checkAExe == true) {
                            processNum = 1;
                            String aCall = "";
                            if ("Window".equals(os.toString())) {
                                aCall = rscriptPath + " \"" + presentWorkingDirectory + "\\a.R\" " + resultDir + " " + genotypeFileName + " " + phenotypeFileName + " " + resultFileName + " "
                                        + missingPercentage + " " + minorAllaleFrequnce + " " + pic + " " + methodName + " " + nofCores + " " + missingCh + " " + "-1 -1 -1 -1 -1 1 " + traitNum + " " + dataType + " "
                                        + alleleSeperator + " " + tmp1 + " " + tmp2 + " " + tmp3 + " " + engineSelection;
                            }
                            if ("Linux".equals(os.toString())) {
                                aCall = rscriptPath + " " + presentWorkingDirectory + "/a.R " + resultDir + " " + genotypeFileName + " " + phenotypeFileName + " " + resultFileName + " "
                                        + missingPercentage + " " + minorAllaleFrequnce + " " + pic + " " + methodName + " " + nofCores + " " + missingCh + " " + "-1 -1 -1 -1 -1 1 " + traitNum + " " + dataType + " "
                                        + alleleSeperator + " " + tmp1 + " " + tmp2 + " " + tmp3 + " " + engineSelection;
                            }
                            runA = Runtime.getRuntime().exec(aCall);
                            int statA = runA.waitFor();
                            if (statA == 0) {
                                checkAExe = false;
                            }
                        }

                        //running ab
                        String abCall = "";
                        processNum = 2;
                        if ("Window".equals(os.toString())) {
                            abCall = rscriptPath + " \"" + presentWorkingDirectory + "\"" + "\\ab.R " + resultDir + " " + genotypeFileName + " " + phenotypeFileName + " " + resultFileName + " "
                                    + missingPercentage + " " + minorAllaleFrequnce + " " + pic + " " + methodName + " " + nofCores + " " + missingCh + " " + "-1 -1 -1 -1 -1 1 " + traitNum + " " + dataType + " "
                                    + alleleSeperator + " " + tmp1 + " " + tmp2 + " " + tmp3 + " " + engineSelection;// crossValidation + " " + RandomForestParameter + " " + BayesParameters;
                        }
                        if ("Linux".equals(os.toString())) {
                            abCall = rscriptPath + " " + presentWorkingDirectory + "/ab.R" + " " + resultDir + " " + genotypeFileName + " " + phenotypeFileName + " " + resultFileName + " "
                                    + missingPercentage + " " + minorAllaleFrequnce + " " + pic + " " + methodName + " " + nofCores + " " + missingCh + " " + "-1 -1 -1 -1 -1 1 " + traitNum + " " + dataType + " "
                                    + alleleSeperator + " " + tmp1 + " " + tmp2 + " " + tmp3 + " " + engineSelection;// crossValidation + " " + RandomForestParameter + " " + BayesParameters;

                        }
                        Process abRun = Runtime.getRuntime().exec(abCall);
                        int statAB = abRun.waitFor();

                        if (Constant.engineSelection == 2) {
                            if (Constant.RidgeRegressionFortran == true) {
                                processNum = 10;
                                lblMethodName.setText("RidgeRegression...");
                                String bFOrtranCall = "";
                                if ("Window".equals(os.toString())) {
                                    bFOrtranCall = rscriptPath + " \"" + presentWorkingDirectory + "\"" + "\\bFort.R " + resultDir + " " + genotypeFileName + " " + phenotypeFileName + " " + resultFileName
                                            + " " + missingPercentage + " " + minorAllaleFrequnce + " " + pic + " " + methodName + " " + nofCores + " " + missingCh + " " + "1 -1 -1 -1 -1 -1"
                                            + " " + traitNum + " " + dataType + " " + alleleSeperator + " " + tmp1 + " " + tmp2 + " " + tmp3 + " " + Constant.engineSelection + " " + "1" + " " + "-1" + " " + populationFileName + " " + BayesParameters + " " + crossValidation;  //last two values for  RR and BayesA
                                }
                                if ("Linux".equals(os.toString())) {
                                    bFOrtranCall = rscriptPath + " " + presentWorkingDirectory + "/bFort.R" + " " + resultDir + " " + genotypeFileName + " " + phenotypeFileName + " " + resultFileName
                                            + " " + missingPercentage + " " + minorAllaleFrequnce + " " + pic + " " + methodName + " " + nofCores + " " + missingCh + " " + "1 -1 -1 -1 -1 -1"
                                            + " " + traitNum + " " + dataType + " " + alleleSeperator + " " + tmp1 + " " + tmp2 + " " + tmp3 + " " + Constant.engineSelection + " " + "1" + " " + "-1" + " " + populationFileName + " " + BayesParameters + " " + crossValidation;  //last two values for  RR and BayesA
                                }
                                fotRR = Runtime.getRuntime().exec(bFOrtranCall);
                                int status = fotRR.waitFor();
                                if (status == 0) {
                                    chkRRFotron.setSelected(true);
                                }
                            }

                            if (Constant.BayeAFortran == true) {
                                processNum = 11;
                                lblMethodName.setText("BayesA");
                                String bFOrtranCall12 = "";
                                if ("Window".equals(os.toString())) {
                                    bFOrtranCall12 = rscriptPath + " \"" + presentWorkingDirectory + "\"" + "\\bFort.R" + " " + resultDir + " " + genotypeFileName + " " + phenotypeFileName + " " + resultFileName
                                            + " " + missingPercentage + " " + minorAllaleFrequnce + " " + pic + " " + methodName + " " + nofCores + " " + missingCh + " " + "1 -1 -1 -1 -1 -1"
                                            + " " + traitNum + " " + dataType + " " + alleleSeperator + " " + tmp1 + " " + tmp2 + " " + tmp3 + " " + Constant.engineSelection + " " + "-1" + " " + "1" + " " + populationFileName + " " + BayesParameters + " " + crossValidation; //last two RR and BayesA
                                }
                                if ("Linux".equals(os.toString())) {
                                    bFOrtranCall12 = rscriptPath + " " + presentWorkingDirectory + "/bFort.R" + " " + resultDir + " " + genotypeFileName + " " + phenotypeFileName + " " + resultFileName
                                            + " " + missingPercentage + " " + minorAllaleFrequnce + " " + pic + " " + methodName + " " + nofCores + " " + missingCh + " " + "1 -1 -1 -1 -1 -1"
                                            + " " + traitNum + " " + dataType + " " + alleleSeperator + " " + tmp1 + " " + tmp2 + " " + tmp3 + " " + Constant.engineSelection + " " + "-1" + " " + "1" + " " + populationFileName + " " + BayesParameters + " " + crossValidation; //last two RR and BayesA
                                }
                                fotBA = Runtime.getRuntime().exec(bFOrtranCall12);
                                int status = fotBA.waitFor();
                                if (status == 0) {
                                    chkBayesAFotron.setSelected(true);
                                }
                            }
                        }

                        if (Constant.engineSelection == 1) {
                            //itterating for performing method 
                            for (int i = 1; i <= 6; i++) {
                                if (exit == false) {
                                    if (i == 1) {
                                        if (Constant.RidgeReg == true) {
                                            processNum = 3;
                                            lblMethodName.setText("RidgeRegression...");
                                            String bCall = "";
                                            if ("Window".equals(os.toString())) {
                                                bCall = rscriptPath + " \"" + presentWorkingDirectory + "\"" + "\\b.R" + " " + resultDir + " " + genotypeFileName + " " + phenotypeFileName + " " + resultFileName
                                                        + " " + missingPercentage + " " + minorAllaleFrequnce + " " + pic + " " + methodName + " " + nofCores + " " + missingCh + " " + "1 -1 -1 -1 -1 -1"
                                                        + " " + traitNum + " " + dataType + " " + alleleSeperator + " " + tmp1 + " " + tmp2 + " " + tmp3 + " " + engineSelection;//crossValidation + " " + BayesParameters + " " + RandomForestParameter;
                                            }
                                            if ("Linux".equals(os.toString())) {
                                                bCall = rscriptPath + " " + presentWorkingDirectory + "/b.R" + " " + resultDir + " " + genotypeFileName + " " + phenotypeFileName + " " + resultFileName
                                                        + " " + missingPercentage + " " + minorAllaleFrequnce + " " + pic + " " + methodName + " " + nofCores + " " + missingCh + " " + "1 -1 -1 -1 -1 -1"
                                                        + " " + traitNum + " " + dataType + " " + alleleSeperator + " " + tmp1 + " " + tmp2 + " " + tmp3 + " " + engineSelection;//crossValidation + " " + BayesParameters + " " + RandomForestParameter;
                                            }
                                            p1 = Runtime.getRuntime().exec(bCall);

                                            //(***************************************************
                                            //capturing error log when running rscript 
                                            try (BufferedReader reader1 = new BufferedReader(new InputStreamReader(p1.getErrorStream()))) {
                                                String line;
                                                while ((line = reader1.readLine()) != null) {
                                                    System.out.println(line);
                                                }
                                            }
                                            //********************************************************)
                                            int status1 = p1.waitFor();
                                            if (status1 == 0) {
                                                chkRidgeRegression.setSelected(true);
                                            }
                                            if (status1 == 1) {
                                                FileWriter fstream = new FileWriter(Constant.resultdirectory + resultFileName, true);
                                                try (BufferedWriter fbw = new BufferedWriter(fstream)) {
                                                    fbw.write("<hr/><p style=\"color:red\"> <b>Error:</b> <i>Couldn't run the RidgeRegression method.</i> <br>");
                                                    fbw.write("<b>Please check the log file for more details.</b>");
                                                    fbw.write("</p><hr/>");
                                                }
                                            }
                                        }
                                    }

                                    //***********************************************************************************************
                                    //Running bayescpi script 
                                    if (i == 2) {
                                        if (Constant.BayesCpi == true) {
                                            processNum = 4;
                                            lblMethodName.setText("BayesCpi...");
                                            String bCall = "";
                                            if ("Window".equals(os.toString())) {
                                                bCall = rscriptPath + " \"" + presentWorkingDirectory + "\"" + "\\b.R" + " " + resultDir + " " + genotypeFileName + " " + phenotypeFileName + " " + resultFileName
                                                        + " " + missingPercentage + " " + minorAllaleFrequnce + " " + pic + " " + methodName + " " + nofCores + " " + missingCh + " " + "-1 -1 -1 1 -1 -1"
                                                        + " " + traitNum + " " + dataType + " " + alleleSeperator + " " + tmp1 + " " + tmp2 + " " + tmp3 + " " + engineSelection;// crossValidation + " " + BayesParameters + " " + RandomForestParameter;
                                            }
                                            if ("Linux".equals(os.toString())) {
                                                bCall = rscriptPath + " " + presentWorkingDirectory + "/b.R" + " " + resultDir + " " + genotypeFileName + " " + phenotypeFileName + " " + resultFileName
                                                        + " " + missingPercentage + " " + minorAllaleFrequnce + " " + pic + " " + methodName + " " + nofCores + " " + missingCh + " " + "-1 -1 -1 1 -1 -1"
                                                        + " " + traitNum + " " + dataType + " " + alleleSeperator + " " + tmp1 + " " + tmp2 + " " + tmp3 + " " + engineSelection;// crossValidation + " " + BayesParameters + " " + RandomForestParameter;

                                            }
                                            p2 = Runtime.getRuntime().exec(bCall);

                                            //(***********************************************************************
                                            //caputring error stream running rscript
                                            try (BufferedReader reader2 = new BufferedReader(new InputStreamReader(p2.getErrorStream()))) {
                                                String line;
                                                while ((line = reader2.readLine()) != null) {
                                                    System.out.println(line);
                                                }
                                            }
                                            //************************************************************************)
                                            int waitStatus2 = p2.waitFor();
                                            if (waitStatus2 == 0) {
                                                chkBayesCpi.setSelected(true);
                                            }
                                            if (waitStatus2 == 1) {
                                                FileWriter fstream = new FileWriter(Constant.resultdirectory + resultFileName, true);
                                                try (BufferedWriter fbw = new BufferedWriter(fstream)) {
                                                    fbw.write("<hr/><p style=\"color:red\"> <b>Error:</b> <i>Couldn't run the BayesCpi method.</i> <br>");
                                                    fbw.write("<b>Please check the log file for more details.</b>");
                                                    fbw.write("</p><hr/>");
                                                }
                                            }
                                        }
                                    }
                                    if (i == 3) {
                                        if (Constant.BayesB == true) {
                                            processNum = 5;
                                            lblMethodName.setText("BayesB...");
                                            String bCall = "";
                                            if ("Window".equals(os.toString())) {
                                                bCall = rscriptPath + " \"" + presentWorkingDirectory + "\"" + "\\b.R" + " " + resultDir + " " + genotypeFileName + " " + phenotypeFileName + " " + resultFileName
                                                        + " " + missingPercentage + " " + minorAllaleFrequnce + " " + pic + " " + methodName + " " + nofCores + " " + missingCh + " " + "-1 -1 1 -1 -1 -1"
                                                        + " " + traitNum + " " + dataType + " " + alleleSeperator + " " + tmp1 + " " + tmp2 + " " + tmp3 + " " + engineSelection;//crossValidation + " " + BayesParameters + " " + RandomForestParameter;
                                            }
                                            if ("Linux".equals(os.toString())) {
                                                bCall = rscriptPath + " " + presentWorkingDirectory + "/b.R" + " " + resultDir + " " + genotypeFileName + " " + phenotypeFileName + " " + resultFileName
                                                        + " " + missingPercentage + " " + minorAllaleFrequnce + " " + pic + " " + methodName + " " + nofCores + " " + missingCh + " " + "-1 -1 1 -1 -1 -1"
                                                        + " " + traitNum + " " + dataType + " " + alleleSeperator + " " + tmp1 + " " + tmp2 + " " + tmp3 + " " + engineSelection;//crossValidation + " " + BayesParameters + " " + RandomForestParameter;
                                            }
                                            p3 = Runtime.getRuntime().exec(bCall);

                                            //(***************************************************
                                            //caputring error 
                                            try (BufferedReader reader3 = new BufferedReader(new InputStreamReader(p3.getErrorStream()))) {
                                                String line;
                                                while ((line = reader3.readLine()) != null) {
                                                    System.out.println(line);
                                                }
                                            }
                                            //****************************************************)
                                            int waitStatus3 = p3.waitFor();
                                            if (waitStatus3 == 0) {
                                                chkBayesB.setSelected(true);
                                            }
                                            if (waitStatus3 == 1) {
                                                FileWriter fstream = new FileWriter(Constant.resultdirectory + resultFileName, true);
                                                try (BufferedWriter fbw = new BufferedWriter(fstream)) {
                                                    fbw.write("<hr/><p style=\"color:red\"> <b>Error:</b> <i>Couldn't run the BayesB method.</i> <br>");
                                                    fbw.write("<b>Please check the log file for more details.</b>");
                                                    fbw.write("</p><hr/>");
                                                }
                                            }
                                        }
                                    }
                                    if (i == 4) {
                                        if (Constant.BayesLasso == true) {
                                            processNum = 6;
                                            lblMethodName.setText("BayesLasso...");
                                            String bCall = "";
                                            if ("Window".equals(os.toString())) {
                                                bCall = rscriptPath + " \"" + presentWorkingDirectory + "\"" + "\\b.R" + " " + resultDir + " " + genotypeFileName + " " + phenotypeFileName + " " + resultFileName
                                                        + " " + missingPercentage + " " + minorAllaleFrequnce + " " + pic + " " + methodName + " " + nofCores + " " + missingCh + " " + "-1 -1 -1 -1 1 -1"
                                                        + " " + traitNum + " " + dataType + " " + alleleSeperator + " " + tmp1 + " " + tmp2 + " " + tmp3 + " " + engineSelection;//crossValidation + " " + BayesParameters + " " + RandomForestParameter;
                                            }
                                            if ("Linux".equals(os.toString())) {
                                                bCall = rscriptPath + " " + presentWorkingDirectory + "/b.R" + " " + resultDir + " " + genotypeFileName + " " + phenotypeFileName + " " + resultFileName
                                                        + " " + missingPercentage + " " + minorAllaleFrequnce + " " + pic + " " + methodName + " " + nofCores + " " + missingCh + " " + "-1 -1 -1 -1 1 -1"
                                                        + " " + traitNum + " " + dataType + " " + alleleSeperator + " " + tmp1 + " " + tmp2 + " " + tmp3 + " " + engineSelection;//crossValidation + " " + BayesParameters + " " + RandomForestParameter;

                                            }
                                            p4 = Runtime.getRuntime().exec(bCall);

                                            //(************************************************************
                                            //caputring error log
                                            try (BufferedReader reader4 = new BufferedReader(new InputStreamReader(p4.getErrorStream()))) {
                                                String line;
                                                while ((line = reader4.readLine()) != null) {
                                                    System.out.println(line);
                                                }
                                            }
                                            //***************************************************************)
                                            int waitStatus4 = p4.waitFor();
                                            if (waitStatus4 == 0) {
                                                chkBayesLasso.setSelected(true);
                                            }
                                            if (waitStatus4 == 1) {
                                                FileWriter fstream = new FileWriter(Constant.resultdirectory + resultFileName, true);
                                                try (BufferedWriter fbw = new BufferedWriter(fstream)) {
                                                    fbw.write("<hr/><p style=\"color:red\"> <b>Error:</b> <i>Couldn't run the BayesLasso method.</i> <br>");
                                                    fbw.write("<b>Please check the log file for more details.</b>");
                                                    fbw.write("</p><hr/>");
                                                }
                                            }
                                        }
                                    }

                                    if (i == 5) {
                                        if (Constant.RandomeForest == true) {
                                            processNum = 7;
                                            lblMethodName.setText("RandomForest...");
                                            String bCall = "";
                                            if ("Window".equals(os.toString())) {
                                                bCall = rscriptPath + " \"" + presentWorkingDirectory + "\"" + "\\b.R" + " " + resultDir + " " + genotypeFileName + " " + phenotypeFileName + " " + resultFileName
                                                        + " " + missingPercentage + " " + minorAllaleFrequnce + " " + pic + " " + methodName + " " + nofCores + " " + missingCh + " " + "-1 -1 -1 -1 -1 1"
                                                        + " " + traitNum + " " + dataType + " " + alleleSeperator + " " + tmp1 + " " + tmp2 + " " + tmp3 + " " + engineSelection;// crossValidation + " " + BayesParameters + " " + RandomForestParameter;
                                            }
                                            if ("Linux".equals(os.toString())) {
                                                bCall = rscriptPath + " " + presentWorkingDirectory + "/b.R" + " " + resultDir + " " + genotypeFileName + " " + phenotypeFileName + " " + resultFileName
                                                        + " " + missingPercentage + " " + minorAllaleFrequnce + " " + pic + " " + methodName + " " + nofCores + " " + missingCh + " " + "-1 -1 -1 -1 -1 1"
                                                        + " " + traitNum + " " + dataType + " " + alleleSeperator + " " + tmp1 + " " + tmp2 + " " + tmp3 + " " + engineSelection;// crossValidation + " " + BayesParameters + " " + RandomForestParameter;
                                            }
                                            p5 = Runtime.getRuntime().exec(bCall);

                                            //(************************************************************
                                            //caputring inputstream log
                                            try (BufferedReader reader5 = new BufferedReader(new InputStreamReader(p5.getErrorStream()))) {
                                                String line;
                                                while ((line = reader5.readLine()) != null) {
                                                    System.out.println(line);
                                                }
                                            }
                                            //***************************************************************)
                                            int waitStatus5 = p5.waitFor();
                                            if (waitStatus5 == 0) {
                                                chkRandomForest.setSelected(true);
                                            }
                                            if (waitStatus5 == 1) {
                                                FileWriter fstream = new FileWriter(Constant.resultdirectory + resultFileName, true);
                                                try (BufferedWriter fbw = new BufferedWriter(fstream)) {
                                                    fbw.write("<hr/><p style=\"color:red\"> <b>Error:</b> <i>Couldn't run the RandomForest method.</i> <br>");
                                                    fbw.write("<b>Please check the log file for more details.</b>");
                                                    fbw.write("</p><hr/>");
                                                }
                                            }
                                        }
                                    }
                                    if (i == 6) {
                                        if (Constant.KinshipGauss == true) {
                                            processNum = 8;
                                            lblMethodName.setText("KinshipGauss...");
                                            String bCall = "";
                                            if ("Window".equals(os.toString())) {
                                                bCall = rscriptPath + " \"" + presentWorkingDirectory + "\"" + "\\b.R" + " " + resultDir + " " + genotypeFileName + " " + phenotypeFileName + " " + resultFileName
                                                        + " " + missingPercentage + " " + minorAllaleFrequnce + " " + pic + " " + methodName + " " + nofCores + " " + missingCh + " " + "-1 1 -1 -1 -1 -1"
                                                        + " " + traitNum + " " + dataType + " " + alleleSeperator + " " + tmp1 + " " + tmp2 + " " + tmp3 + " " + engineSelection;//crossValidation + " " + BayesParameters + " " + RandomForestParameter;
                                            }
                                            if ("Window".equals(os.toString())) {
                                                bCall = rscriptPath + " " + presentWorkingDirectory + "/b.R" + " " + resultDir + " " + genotypeFileName + " " + phenotypeFileName + " " + resultFileName
                                                        + " " + missingPercentage + " " + minorAllaleFrequnce + " " + pic + " " + methodName + " " + nofCores + " " + missingCh + " " + "-1 1 -1 -1 -1 -1"
                                                        + " " + traitNum + " " + dataType + " " + alleleSeperator + " " + tmp1 + " " + tmp2 + " " + tmp3 + " " + engineSelection;//crossValidation + " " + BayesParameters + " " + RandomForestParameter;

                                            }
                                            p6 = Runtime.getRuntime().exec(bCall);

                                            //(************************************************************
                                            //caputring error  log
                                            try (BufferedReader reader6 = new BufferedReader(new InputStreamReader(p6.getErrorStream()))) {
                                                String line;
                                                while ((line = reader6.readLine()) != null) {
                                                    System.out.println(line);
                                                }
                                            }
                                            //*************************************************************)
                                            int waitStatus6 = p6.waitFor();
                                            if (waitStatus6 == 0) {
                                                chkKinshipGauss.setSelected(true);
                                            }
                                            if (waitStatus6 == 1) {
                                                FileWriter fstream = new FileWriter(Constant.resultdirectory + resultFileName, true);
                                                try (BufferedWriter fbw = new BufferedWriter(fstream)) {
                                                    fbw.write("<hr/><p style=\"color:red\"> <b>Error:</b> <i>Couldn't run the KinshipGauss method.</i> <br>");
                                                    fbw.write("<b>Please check the log file for more details.</b>");
                                                    fbw.write("</p><hr/>");
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        //integrating all the methods outputs
                        processNum = 9;
                        String cCall = "";
                        if ("Window".equals(os.toString())) {
                            cCall = rscriptPath + " \"" + presentWorkingDirectory + "\"" + "\\c.R" + " " + resultDir + " " + genotypeFileName + " " + phenotypeFileName + " " + resultFileName
                                    + " " + missingPercentage + " " + minorAllaleFrequnce + " " + pic + " " + methodName + " " + nofCores + " " + missingCh + " " + "-1 1 -1 -1 -1 -1"
                                    + " " + traitNum + " " + dataType + " " + alleleSeperator + " " + tmp1 + " " + tmp2 + " " + tmp3 + " " + engineSelection;// crossValidation + " " + BayesParameters + " " + RandomForestParameter;

                        }
                        if ("Linux".equals(os.toString())) {
                            cCall = rscriptPath + " " + presentWorkingDirectory + "/c.R" + " " + resultDir + " " + genotypeFileName + " " + phenotypeFileName + " " + resultFileName
                                    + " " + missingPercentage + " " + minorAllaleFrequnce + " " + pic + " " + methodName + " " + nofCores + " " + missingCh + " " + "-1 1 -1 -1 -1 -1"
                                    + " " + traitNum + " " + dataType + " " + alleleSeperator + " " + tmp1 + " " + tmp2 + " " + tmp3 + " " + engineSelection;// crossValidation + " " + BayesParameters + " " + RandomForestParameter;

                        }
                        cRun = Runtime.getRuntime().exec(cCall);
                        int statC = cRun.waitFor();
                        Constant.resultGenratedList.add(resultFileName);

                        if (Constant.RidgeReg == true) {
                            chkRidgeRegression.setSelected(false);
                        }
                        if (Constant.BayesB == true) {
                            chkBayesB.setSelected(false);
                        }
                        if (Constant.BayesCpi == true) {
                            chkBayesCpi.setSelected(false);
                        }
                        if (Constant.BayesLasso == true) {
                            chkBayesLasso.setSelected(false);
                        }
                        if (Constant.RandomeForest == true) {
                            chkRandomForest.setSelected(false);
                        }
                        if (Constant.KinshipGauss == true) {
                            chkKinshipGauss.setSelected(false);
                        }
                        if (Constant.RidgeRegressionFortran == true) {
                            chkRRFotron.setSelected(false);
                        }
                        if (Constant.BayeAFortran == true) {
                            chkBayesAFotron.setSelected(false);
                        }
                        lblMethodName.setText("");
                    }
                    dispose(); //closing the progressDisplay class 
                }
                return "html";
            }

        };
        ss.execute();

        btnCancel.addActionListener(
                new ActionListener() {
                    @Override
                    public void actionPerformed(ActionEvent e) {

                        int option = JOptionPane.showOptionDialog(rootPane, "Do you want to cancel analysis", "Cancel dialog", JOptionPane.YES_NO_OPTION,
                                JOptionPane.WARNING_MESSAGE, null, null, null);
                        if (option == JOptionPane.YES_OPTION) {
                            closeAnalysis();
                        }
                        if (option == JOptionPane.NO_OPTION) {
                            //donothing
                        }
                    }
                });

    }
//end region 

    /**
     * //start region initializing the GUI components
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">                          
    private void init() {
        mainPanel.setLayout(null);

        //    jPanel2 = new javax.swing.JPanel();
        lblPerformAnalysis = new javax.swing.JLabel();
        lblPerformAnalysis.setBounds(30, 10, 200, 20);
        lblTraitName = new javax.swing.JLabel();
        lblTraitName.setBounds(250, 10, 200, 20);
        lblUsingmethid = new javax.swing.JLabel();
        lblUsingmethid.setBounds(30, 40, 200, 20);
        lblMethodName = new javax.swing.JLabel();
        lblMethodName.setBounds(250, 40, 200, 20);
        //  jPanel3 = new javax.swing.JPanel();
        chkRidgeRegression = new javax.swing.JCheckBox();
        //  chkRidgeRegression.setBounds(40, 80, 200, 20);
        chkBayesB = new javax.swing.JCheckBox();
        //  chkBayesB.setBounds(40, 110, 200, 20);
        chkBayesCpi = new javax.swing.JCheckBox();
        //  chkBayesCpi.setBounds(40, 140, 200, 20);
        chkBayesLasso = new javax.swing.JCheckBox();
        //  chkBayesLasso.setBounds(40, 170, 200, 20);
        chkRandomForest = new javax.swing.JCheckBox();
        //  chkRandomForest.setBounds(40, 200, 200, 20);
        chkKinshipGauss = new javax.swing.JCheckBox();
        //  chkKinshipGauss.setBounds(40, 230, 200, 20);
        chkRRFotron = new javax.swing.JCheckBox();

        chkBayesAFotron = new javax.swing.JCheckBox();
        // jButton1 = new javax.swing.JButton();
        //   y = y + 20;
        // jButton1.setBounds(110, 270, 100, 20);
        btnCancel = new javax.swing.JButton();
        btnCancel.setBounds(140, 270, 100, 25);

        setDefaultCloseOperation(javax.swing.WindowConstants.DISPOSE_ON_CLOSE);
        setTitle("Analysis");

        lblPerformAnalysis.setFont(new java.awt.Font("SansSerif", 0, 15)); // NOI18N
        lblPerformAnalysis.setText("Prerforming analysis of trait :");

        lblTraitName.setFont(new java.awt.Font("SansSerif", 0, 15)); // NOI18N
        lblTraitName.setText("trait name");

        lblUsingmethid.setFont(new java.awt.Font("SansSerif", 0, 15)); // NOI18N
        lblUsingmethid.setText("Present running method :");

        lblMethodName.setFont(new java.awt.Font("SansSerif", 0, 15)); // NOI18N
        lblMethodName.setText("Data Processing");

        // javax.swing.GroupLayout jPanel2Layout = new javax.swing.GroupLayout(jPanel2);
        //   jPanel2.setLayout(jPanel2Layout);
        //  jPanel3.setBorder(javax.swing.BorderFactory.createTitledBorder(null, "Selected methods for analysis :", javax.swing.border.TitledBorder.CENTER, javax.swing.border.TitledBorder.DEFAULT_POSITION, new java.awt.Font("SansSerif", 0, 15))); // NOI18N
        chkRidgeRegression.setFont(new java.awt.Font("SansSerif", 0, 15)); // NOI18N
        chkRidgeRegression.setText("RidgeRegression");

        chkBayesB.setFont(new java.awt.Font("SansSerif", 0, 15)); // NOI18N
        chkBayesB.setText("BayesB");

        chkBayesCpi.setFont(new java.awt.Font("SansSerif", 0, 15)); // NOI18N
        chkBayesCpi.setText("BayesCpi");

        chkBayesLasso.setFont(new java.awt.Font("SansSerif", 0, 15)); // NOI18N
        chkBayesLasso.setText("BayesLasso");

        chkRandomForest.setFont(new java.awt.Font("SansSerif", 0, 15)); // NOI18N
        chkRandomForest.setText("RandomForest");

        chkKinshipGauss.setFont(new java.awt.Font("SansSerif", 0, 15)); // NOI18N
        chkKinshipGauss.setText("KinshipGauss");

        chkRRFotron.setFont(new java.awt.Font("SansSerif", 0, 15));
        chkRRFotron.setText("RidgeRegression");

        chkBayesAFotron.setFont(new java.awt.Font("SansSerif", 0, 15));
        chkBayesAFotron.setText("BayesA");

        //   javax.swing.GroupLayout jPanel3Layout = new javax.swing.GroupLayout(jPanel3);
        //   jPanel3.setLayout(jPanel3Layout);
        //jButton1.setFont(new java.awt.Font("Tahoma", 0, 14)); // NOI18N
        //jButton1.setText("Start");
        btnCancel.setFont(new java.awt.Font("SansSerif", 0, 14)); // NOI18N
        btnCancel.setText("Cancel");

        mainPanel.add(lblPerformAnalysis);
        mainPanel.add(lblTraitName);
        mainPanel.add(lblUsingmethid);
        mainPanel.add(lblMethodName);
        mainPanel.add(chkRidgeRegression);
        mainPanel.add(chkBayesCpi);
        mainPanel.add(chkBayesB);
        mainPanel.add(chkBayesLasso);
        mainPanel.add(chkRandomForest);
        mainPanel.add(chkKinshipGauss);
        mainPanel.add(chkBayesAFotron);
        mainPanel.add(chkRRFotron);
        //   mainPanel.add(jButton1);
        mainPanel.add(btnCancel);
        add(mainPanel);

    }// </editor-fold>                        
//end region 

    //main method not used 
    public static void main(String args[]) {
        /* Set the Nimbus look and feel */
        //<editor-fold defaultstate="collapsed" desc=" Look and feel setting code (optional) ">
        /* If Nimbus (introduced in Java SE 6) is not available, stay with the default look and feel.
         * For details see http://download.oracle.com/javase/tutorial/uiswing/lookandfeel/plaf.html 
         */
        try {
            for (javax.swing.UIManager.LookAndFeelInfo info : javax.swing.UIManager.getInstalledLookAndFeels()) {
                if ("Nimbus".equals(info.getName())) {
                    javax.swing.UIManager.setLookAndFeel(info.getClassName());
                    break;
                }
            }
        } catch (ClassNotFoundException ex) {
            java.util.logging.Logger.getLogger(ProgressDisplay.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (InstantiationException ex) {
            java.util.logging.Logger.getLogger(ProgressDisplay.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (IllegalAccessException ex) {
            java.util.logging.Logger.getLogger(ProgressDisplay.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (javax.swing.UnsupportedLookAndFeelException ex) {
            java.util.logging.Logger.getLogger(ProgressDisplay.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        }
        //</editor-fold>
        //</editor-fold>
    }

    //start region
   /*
     * method to kill the running process
     */
    private void closeAnalysis() {
        exit = true;
        try {
            if (processNum == 1) {
                runA.destroy();
            }
            if (processNum == 2) {
                abRun.destroy();
            }
            if (processNum == 3) {
                p1.destroy();
            }
            if (processNum == 4) {
                p2.destroy();
            }
            if (processNum == 5) {
                p3.destroy();
            }
            if (processNum == 6) {
                p4.destroy();
            }
            if (processNum == 7) {
                p5.destroy();
            }
            if (processNum == 8) {
                p6.destroy();
            }
            if (processNum == 9) {
                cRun.destroy();
            }
            if (processNum == 10) {
                fotRR.destroy();
            }
            if (processNum == 11) {
                fotBA.destroy();
            }
        } catch (Exception e) {
            System.out.println(e.getMessage());
        }
        File gsdataFile = new File(Constant.resultdirectory + "gsdata.rdata");
        boolean gsexists = gsdataFile.exists();
        if (gsexists == true) {
            boolean delete = gsdataFile.delete();
        }

        File tempFile = new File(Constant.resultdirectory + Constant.pheno_Analysis);
        boolean existstemp = tempFile.exists();
        if (existstemp == true) {
            boolean tempdelete = tempFile.delete();
        }
        File resAnalysis = new File(Constant.resultdirectory + resultFileName);
        boolean existsGenoHtml = resAnalysis.exists();
        if (existsGenoHtml == true) {
            GsMethods.removePngFromHtml(Constant.resultdirectory + resultFileName);
            boolean delete = resAnalysis.delete();
        }

        File cmsFile = new File(Constant.resultdirectory + "CMImData_For_Analysis.csv");
        boolean cmsexists = cmsFile.exists();
        if (cmsexists == true) {
            boolean cMSdelete = cmsFile.delete();
        }

        File gsResFile = new File(Constant.resultdirectory + "GsResTmp.txt");
        boolean gsResexists = gsResFile.exists();
        if (gsResexists == true) {
            boolean gsResdelete = gsResFile.delete();
        }
        dispose();
        Constant.analysisFileName = null;
        Constant.analysisCancled = true;
    }
    //end region

    public String convertNA(String path) {

        String getname = null;
        PrintWriter pw = null;
        FileWriter fw;

        try {
            FileInputStream stream = new FileInputStream(path);
            DataInputStream inputStream = new DataInputStream(stream);
            BufferedReader br = new BufferedReader(new InputStreamReader(inputStream));
            String[] str;
            int num = 0;
            String line;
            int headerCount = 0;

            String filename = new File(path).getName();
            getname = FilenameUtils.removeExtension(filename) + "-convert.csv";
            fw = new FileWriter(Constant.resultdirectory + getname);
            pw = new PrintWriter(fw);
            //getting the input file name to add in the text file
            while ((line = br.readLine()) != null) {
                str = line.split(",");

                if (num == 0) {
                    headerCount = str.length;
                    for (int i = 0; i < headerCount; i++) {
                        if (i == headerCount - 1) {
                            pw.println(str[i]);
                        } else {
                            pw.print(str[i] + ",");
                        }
                    }
                } else {
                    for (int i = 0; i < headerCount; i++) {
                        if (i == 0) {
                            pw.print(str[i] + ",");
                        }
                        if (i > 0) {
                            if (i == headerCount - 1) {
                                if ("NA".equals(str[i])) {
                                    pw.println("9");
                                } else {
                                    pw.println(str[i]);
                                }
                            } else {
                                if ("NA".equals(str[i])) {
                                    pw.print("9,");
                                } else {
                                    pw.print(str[i] + ",");
                                }
                            }
                        }
                    }
                }
                num++;
            }
        } catch (Exception e) {
            System.out.println(e);
        } finally {
            pw.close();
        }
        return getname;
    }
    //end region

    //start region
    public int headerCount(String path) {
        int countHeader = 0;
        try {
            CsvReader reader = new CsvReader(path, ',');
            reader.readHeaders();
            countHeader = reader.getHeaderCount();
        } catch (FileNotFoundException ex) {
            Logger.getLogger(ProgressDisplay.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ex) {
            Logger.getLogger(ProgressDisplay.class.getName()).log(Level.SEVERE, null, ex);
        }
        return countHeader;
    }
    //end region

    //start region
    //converting the NA in the inoput file to 9 (missing value) for AlphaAGH script
    public String convertTxtFile(String path) {

        String getname = null;
        PrintWriter pw = null;
        FileWriter fw;

        try {
            FileInputStream stream = new FileInputStream(path);
            DataInputStream inputStream = new DataInputStream(stream);
            BufferedReader br = new BufferedReader(new InputStreamReader(inputStream));
            String[] str;
            int num = 0;
            String line;
            int headerCount = 0;

            String filename = new File(path).getName();
            getname = FilenameUtils.removeExtension(filename) + "-subset.csv";
            fw = new FileWriter(Constant.resultdirectory + getname);
            pw = new PrintWriter(fw);
            //getting the input file name to add in the text file
            while ((line = br.readLine()) != null) {
                str = line.split(" +");

                if (num == 0) {
                    headerCount = str.length;
                    for (int i = 0; i < headerCount; i++) {
                        if (i == headerCount - 1) {
                            pw.println(str[i]);
                        } else {
                            pw.print(str[i] + ",");
                        }
                    }
                } else {
                    for (int i = 0; i < headerCount; i++) {
                        if (i == 0) {
                            pw.print(str[i] + ",");
                        }
                        if (i > 0) {
                            if (i == headerCount - 1) {
                                pw.println(str[i]);
                            } else {
                                pw.print(str[i] + ",");
                            }
                        }
                    }
                }
                num++;
            }
        } catch (Exception e) {
            System.out.println(e);
        } finally {
            pw.close();
        }
        return getname;
    }
    //end region      
}
