package org.icrisat.genomicSelection;

import chrriis.dj.nativeswing.NativeSwing;
import chrriis.dj.nativeswing.swtimpl.NativeInterface;
import chrriis.dj.nativeswing.swtimpl.components.JWebBrowser;
import com.csvreader.CsvReader;
import java.awt.BorderLayout;
import java.awt.Dimension;
import java.io.File;
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintWriter;
import java.io.DataInputStream;
import java.io.InputStreamReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.LineNumberReader;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.JTextArea;
import javax.swing.table.DefaultTableModel;
import org.icepdf.ri.common.ComponentKeyBinding;
import org.icepdf.ri.common.SwingController;
import org.icepdf.ri.common.SwingViewBuilder;

/**
 * implementing all the method of GenomicInterface
 */
public class GsMethods implements GenomicInterface {

    JPanel csvpanel;

    //start region 
    //To read the csv file and display 
    @Override
    public JPanel csvreader(String path) {
        try {
            String filename = new File(path).getName();
            CsvReader reader = new CsvReader(path, ',');
            // you have to always call readHeaders first before you do any other operation
            reader.readHeaders();
            //Get CSV file header information
            String[] headers = reader.getHeaders();
            java.util.List<Object[]> list = new ArrayList<>();
            while (reader.readRecord()) {
                //add each row to the list
                list.add(reader.getValues());
            }
            Object[][] data = new String[list.size()][];
            for (int i = 0; i < list.size(); i++) {
                data[i] = list.get(i);
            }

            JTable table = new JTable(data, headers);
            table.setPreferredScrollableViewportSize(new Dimension(5000000, 10000));
            table.setFillsViewportHeight(true);
            table.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
            table.setCellSelectionEnabled(true);
            table.setRowSelectionAllowed(true);
            table.setColumnSelectionAllowed(true);
            //Create the scroll pane and add the table to it.
            JScrollPane scrollPane = new JScrollPane(table);
            scrollPane.setWheelScrollingEnabled(true);

            //Add the scroll pane to this panel.
            csvpanel = new JPanel(new BorderLayout());
            JScrollPane bar = new JScrollPane();

            int horizontalPolicy = JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED;
            bar.setHorizontalScrollBarPolicy(horizontalPolicy);

            csvpanel.setBounds(0, 0, 100, 100);
            csvpanel.setName(filename);
            csvpanel.add(scrollPane);
            csvpanel.repaint();
        } catch (IOException e) {
            e.getMessage();
        }
        return csvpanel;
    }
//end region csvreader

    //start region
    //For reading the text files 
    @Override
    public synchronized JScrollPane textFileReader(String path) {
        JScrollPane scrollpane = new JScrollPane();
        JTextArea txtDisplayArea = new JTextArea();
        txtDisplayArea.setColumns(100);
        txtDisplayArea.setRows(200);
        String name = new File(path).getName();
        scrollpane.setName(name);
        BufferedReader br = null;

        try {
            FileInputStream stream = new FileInputStream(path);
            DataInputStream inputStream = new DataInputStream(stream);
            br = new BufferedReader(new InputStreamReader(inputStream));
            txtDisplayArea.read(br, null);
        } catch (IOException ioe) {
            System.err.println(ioe);
            System.exit(1);
        } finally {
            if (br != null) {
                try {
                    br.close();
                } catch (IOException e) {
                }
            }
            scrollpane.setViewportView(txtDisplayArea);
        }
        return scrollpane;
    }

    //end region textFile reader
    //start region 
    //For opening the pdf file and displaying in tabbed pane
    @Override
    public JPanel pdfviewer(String pdfpath) {
        JPanel viewerComponentPanel = null;
        try {
            // build a controller
            SwingController controller = new SwingController();
            // Build a SwingViewFactory configured with the controller
            SwingViewBuilder factory = new SwingViewBuilder(controller);
            viewerComponentPanel = factory.buildViewerPanel();
            viewerComponentPanel.setName("Gs manual");
            // add copy keyboard command
            ComponentKeyBinding.install(controller, viewerComponentPanel);
            // add interactive mouse link annotation support via callback
            controller.getDocumentViewController().setAnnotationCallback(new org.icepdf.ri.common.MyAnnotationCallback(controller.getDocumentViewController()));
            // Open a PDF document to view
            controller.openDocument(pdfpath);
        } catch (Exception e) {
            e.getMessage();
        }
        return viewerComponentPanel;
    }
//end region pdfviewer

    //start region
    //Method which use to display the html file using DJNative browser 
    @Override
    public JWebBrowser browser(String htmlpath) {
        //System.out.println("calling browser ...");
        JWebBrowser browser = null;
        try {
            NativeInterface.open();
            NativeSwing.initialize();
            // System.out.println("inside the browse");
            browser = new JWebBrowser();
            browser.setPreferredSize(new Dimension(1020, 850));
            browser.navigate("file:///" + htmlpath);
            // browser.navigate("http://www.google.com");
            browser.setBarsVisible(false);
            browser.setButtonBarVisible(false);
            browser.setDefaultPopupMenuRegistered(true);
            browser.setStatusBarVisible(false);
            String nam = new File(htmlpath).getName();
            browser.setName(nam);
        } catch (Exception browseException) {
            browseException.getMessage();
        }
        return browser;
    }
    //end region browser

    //start region
    //adding print feature to jwebbrowser not working 
    public void print(boolean showDialog) {
        JWebBrowser nativeWebBrowser = new JWebBrowser();
        nativeWebBrowser.print(showDialog);
    }
    //end region 

    //start region 
    //Creating the text file i.e writting the data from outputs to the file
    @Override
    public void writeToFile(java.util.List<String> porSavelist, String path) {
        BufferedWriter out;
        try {
            File file = new File(path);
            out = new BufferedWriter(new FileWriter(file, true));
            for (String s : porSavelist) {
                out.write(s);
                out.newLine();
            }
            out.close();
        } catch (IOException e) {
            System.err.println(e);
        }
    }
    //end region  writeToFile
    //star region
    //reading the csv file only but method can used to read SNP file also
    JTable tablePheno;
    JPanel csvpanelPheno;

    public JPanel csvReaderSNP(String path1) {
        int count2 = 0;
        try {
            try {
                DefaultTableModel modelTest;
                String filename = new File(path1).getName();
                //initialized a CsvReader object with file path and delimiter
                CsvReader reader = new CsvReader(path1, ',');
                // you have to always call readHeaders first before you do any other operation
                reader.readHeaders();
                //Get CSV file header information
                String[] headers = reader.getHeaders();
                tablePheno = new JTable();
                modelTest = new DefaultTableModel(headers, 0) {
                };
                File f = new File(path1);
                InputStream is = new FileInputStream(f);
                Scanner scan = new Scanner(is);
                String[] arr = null;
                while (scan.hasNextLine()) {
                    String line = scan.nextLine();
                    if (line.contains(",")) {
                        arr = line.split(",");
                    }
                    if (count2 == 0) {
                    } else {
                        Object[] dat1 = new Object[arr.length];
                        System.arraycopy(arr, 0, dat1, 0, arr.length);
                        modelTest.addRow(dat1);
                    }
                    count2++;
                }
                tablePheno.setModel(modelTest);
                tablePheno.setAutoCreateRowSorter(true);
                tablePheno.setPreferredScrollableViewportSize(tablePheno.getPreferredSize());
                tablePheno.setFillsViewportHeight(true);
                tablePheno.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);

                tablePheno.setAutoCreateRowSorter(true);
                tablePheno.setRowSelectionAllowed(true);
                tablePheno.setColumnSelectionAllowed(true);
                tablePheno.setCellSelectionEnabled(true);
                // table.setRowSorter(sorter);
                //Create the scroll pane and add the table to it.
                JScrollPane scrollPane = new JScrollPane(tablePheno);

                scrollPane.setWheelScrollingEnabled(true);
                //Add the scroll pane to this panel.
                csvpanelPheno = new JPanel(new BorderLayout());
                csvpanelPheno.setBounds(0, 0, 100, 100);
                csvpanelPheno.setName(filename);
                csvpanelPheno.add(scrollPane);
                csvpanelPheno.repaint();
            } catch (IOException | NumberFormatException e) {
                e.getMessage();
            }
        } catch (Exception ee) {
            ee.getMessage();
        }
        return csvpanelPheno;
    }
    //end region csvreaderSNP
    //start region
    //creating a temp file for selected trait of phenotype file selected in analysis dialog box  
    public List traits_list;
    PrintWriter pw = null;
    FileWriter fw;

    public void traitRemove(String phenoPath, String tmpFilename) {
        try {
            CsvReader reader = new CsvReader(phenoPath, ',');
            reader.readHeaders();
            String[] headers = reader.getHeaders();
            traits_list = new ArrayList();
            traits_list = Constant.selectedTraits;
            //adding the genotype the first column data which will not be selected
            traits_list.add("Genotype");
            FileInputStream stream = new FileInputStream(phenoPath);
            DataInputStream inputStream = new DataInputStream(stream);
            BufferedReader br = new BufferedReader(new InputStreamReader(inputStream));
            String[] str;
            int row = 0;
            int num = 0;
            String line;
            int headerCount = 0;
            fw = new FileWriter(tmpFilename);
            pw = new PrintWriter(fw);
            while ((line = br.readLine()) != null) {
                row++;
                str = line.split(",");
                if (num == 0) {
                    headerCount = str.length;
                    for (int i = 0; i < headerCount; i++) {
                        if (traits_list.contains(headers[i])) {
                            if (i == headerCount - 1) {  //checking if it is the last trait
                                pw.print(str[i]);
                                pw.println();
                            } else {
                                pw.print(str[i] + ",");
                            }
                        } else {
                            if (i == headerCount - 1) {
                                pw.println();
                            }
                        }
                    }
                } else {
                    for (int i = 0; i < headerCount; i++) {
                        if (traits_list.contains(headers[i])) {
                            if (i == headerCount - 1) {
                                pw.print(str[i]);
                                pw.println();
                            } else {
                                pw.print(str[i] + ",");
                            }
                        } else {
                            if (i == headerCount - 1) {
                                pw.println();
                            }
                        }
                    }
                }
                num++;
            }
        } catch (FileNotFoundException ex) {
            Logger.getLogger(GsMethods.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ex) {
            Logger.getLogger(GsMethods.class.getName()).log(Level.SEVERE, null, ex);
        } finally {
            pw.flush();
            pw.close();
        }
    }
    //end region trait remove

    
    //start region
    //removing the png file when html file deleted
    public  static void removePngFromHtml(String path) {
        //   try {
        try (Scanner scaner = new Scanner(new File(path))) {
            while (scaner.hasNextLine()) {
                String next = scaner.nextLine();
                // System.out.println(next);
                if (next.contains("png")) { //getting the line of images 
                    //    System.out.println(next);
                    String[] split = next.split("="); //spliting the string obtained 
                    for (int j = 0; j < split.length; j++) { //itterating the string 
                        if (j == 2) {       //we get the file path in 2 line of split
                            // System.out.println(split[j]);
                            String inputString = split[j];  //converting to string the data
                            inputString = inputString.substring(0, inputString.lastIndexOf(" ")) + ""; // removing the extra data in the path at end
                            String replaced = inputString.replaceAll("'", ""); //replacing the ' with empty space
                            //System.out.println(replaceAll); 
                            File deleteFile = new File(replaced);  //delete images (png) of particular related html file
                            boolean delRes = deleteFile.delete();
                            System.out.println("deleted png files " + delRes);
                        }
                    }
                }
            }
            scaner.close();
        } catch (FileNotFoundException ex) {
            Logger.getLogger(GsMethods.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
//end region

    //start region
    //convert tab delimited to csv file
    public static void tabtoCSv(String path, String outPath) {
        PrintWriter printWriter1 = null;
        FileWriter fileWrite;

        try {
            FileInputStream stream = new FileInputStream(path);
            DataInputStream inputStream = new DataInputStream(stream);
            BufferedReader br = new BufferedReader(new InputStreamReader(inputStream));
            String[] str;
            String line;
            int headerCount;

            // System.out.println("getname  " + getname);
            fileWrite = new FileWriter(outPath);
            printWriter1 = new PrintWriter(fileWrite);
            //getting the input file name to add in the text file
            while ((line = br.readLine()) != null) {
                // row++;
                str = line.split("\t");

                headerCount = str.length;
                for (int i = 0; i < headerCount; i++) {
                    if (i >= 1 & i < 11) {

                    } else {
                        if (i == headerCount - 1) {
                            printWriter1.println(str[i]);
                            // System.out.println(str[i]);
                        } else {
                            //System.out.print(str[i] + ",");
                            printWriter1.print(str[i] + ",");
                        }
                    }
                }
            }
        } catch (IOException e) {
            System.out.println(e);
        } finally {
            printWriter1.close();
        }

    }
     //end region

    //start region
    //convert haplotype file to genotype format
    public static void transpose(String path, String outpath) throws IOException {
        int countRows = countRows(path);
        int headerCount = headerCount(path);
        System.out.println("count " + countRows + "  header   " + headerCount);
        Scanner scain;
        int Rowc = 0;
        String inputline;
        String[][] arrayValues = new String[countRows][headerCount];
        try {
            scain = new Scanner(new BufferedReader(new FileReader(path)));
            while (scain.hasNextLine()) {
                inputline = scain.nextLine();
                String[] inarray = inputline.split(",");
                for (int i = 0; i < inarray.length; i++) {
                    arrayValues[Rowc][i] = (inarray[i]);
                }
                Rowc++;
            }
        } catch (FileNotFoundException e) {
            e.getMessage();
        }
        //transpose matrix 
        String[][] temp = new String[arrayValues[0].length][arrayValues.length];
        for (int i = 0; i < arrayValues.length; i++) {
            for (int j = 0; j < arrayValues[0].length; j++) {
                temp[j][i] = arrayValues[i][j];
            }
        }
        try ( //creating a csv file of transpose
                PrintWriter printwriter1 = new PrintWriter(new FileWriter(outpath))) {
            for (String[] temp1 : temp) {
                for (int j = 0; j < temp1.length; j++) {
                    if (j == temp1.length - 1) {
                        printwriter1.print(temp1[j]);
                    } else {
                        printwriter1.print(temp1[j] + ",");
                    }
                }
                printwriter1.println();
            }
        }
    }
    //end region
    
    //start region
    //count number of headers
    public static int countRows(String path) {
        FileInputStream stream = null;
        try {
            stream = new FileInputStream(path);
            DataInputStream inputStream = new DataInputStream(stream);
            BufferedReader br = new BufferedReader(new InputStreamReader(inputStream));
            int count = 0;
            while ((br.readLine()) != null) {
                count++;
            }
            return count;
        } catch (FileNotFoundException ex) {
            Logger.getLogger(GsMethods.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ex) {
            Logger.getLogger(GsMethods.class.getName()).log(Level.SEVERE, null, ex);
        } finally {
            try {
                stream.close();
            } catch (IOException ex) {
                Logger.getLogger(GsMethods.class.getName()).log(Level.SEVERE, null, ex);
            }

        }
        return 0;
    }
       //end region

    //start region
    //count noof columns in  the file
    public static int headerCount(String path) {
        // System.out.println("get header count");
        int countHeader = 0;
        try {
            CsvReader reader = new CsvReader(path, ',');
            reader.readHeaders();
            //  String[] headers = reader.getHeaders();
            countHeader = reader.getHeaderCount();
            //  System.out.println(countHeader);
            //    numberOfSnp = countHeader + "";

        } catch (FileNotFoundException ex) {
            Logger.getLogger(GsMethods.class
                    .getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ex) {
            Logger.getLogger(GsMethods.class
                    .getName()).log(Level.SEVERE, null, ex);
        }
        return countHeader;
    }

    //start region
    //convert the haplo files with : in the fields
    public static void addColnaTofile(String path, String outPath) {
        PrintWriter printWriter1 = null;
        FileWriter fileWrite;

        try {
            FileInputStream stream = new FileInputStream(path);
            DataInputStream inputStream = new DataInputStream(stream);
            BufferedReader br = new BufferedReader(new InputStreamReader(inputStream));
            String[] str;
            String line;
            int headerCount;
            int rowCount = 0;
            fileWrite = new FileWriter(outPath);
            printWriter1 = new PrintWriter(fileWrite);
            //getting the input file name to add in the text file
            while ((line = br.readLine()) != null) {
                if (rowCount == 0) {
                    //leaving the first line as it is trait names
                    str = line.split(",");

                    headerCount = str.length;
                    for (int i = 0; i < headerCount; i++) {
                        if (i == headerCount - 1) {
                            printWriter1.println(str[i]);
                        } else {
                            printWriter1.print(str[i] + ",");
                        }
                    }
                } else {
                    str = line.split(",");
                    headerCount = str.length;
                    for (int i = 0; i < headerCount; i++) {
                        if (i == 0) {
                            //leaving the first column as it is genotype names
                            printWriter1.print(str[i] + ",");
                        }
                        if (i == headerCount - 1) {
                            // printing last column data 
                            if ("NA".equals(str[i])) {
                                printWriter1.println(str[i]);
                            } else {
                                if ("NN".equals(str[i])) {
                                    printWriter1.println("NA");
                                } else {
                                    String addCol = str[i];
                                    addCol = new StringBuilder(addCol).insert(addCol.length() - 1, ":").toString();
                                    printWriter1.println(addCol);
                                }
                            }
                        } else if (i > 0) {
                            if ("NA".equals(str[i])) {
                                printWriter1.print(str[i] + ",");
                            } else {
                                if ("NN".equals(str[i])) {
                                    printWriter1.print("NA,");
                                } else {
                                    String addCol = str[i];
                                    addCol = new StringBuilder(addCol).insert(addCol.length() - 1, ":").toString();
                                    printWriter1.print(addCol + ",");
                                }
                            }
                        }
                    }
                }
                rowCount++;
            }
        } catch (IOException e) {
            System.out.println(e);
        } finally {
            printWriter1.close();
        }
    }

    //getting the time required for completing analysis and summary
    public static void addingProcessTime2HTMlFile(int nooflines, String htmlPath, String timeTaken) {
        try {
            List<String> lines = Files.readAllLines(Paths.get(htmlPath), StandardCharsets.UTF_8);
            String charCountOfLine = lines.get(nooflines - 3);
            int length = charCountOfLine.length() - 4;
            StringBuilder sb = new StringBuilder(lines.get(nooflines - 3));
            sb.insert(length, timeTaken);
            lines.remove(nooflines - 3);
            lines.add(nooflines - 3, sb.toString()); 
            Files.write(Paths.get(htmlPath), lines, StandardCharsets.UTF_8);
        } catch (IOException ex) {
            Logger.getLogger(GsMethods.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    //getting number of lines in text file (Used for html now)
    public static int noofLinesInHTML(String htmmlPath) {
        int linenumber = 0;
        File file = new File(htmmlPath);

        if (file.exists()) {

            FileReader fr = null;
            try {
                fr = new FileReader(file);
                try (LineNumberReader lnr = new LineNumberReader(fr)) {
                    while (lnr.readLine() != null) {
                        linenumber++;
                    }
                }
            } catch (FileNotFoundException ex) {
                Logger.getLogger(GsMethods.class.getName()).log(Level.SEVERE, null, ex);
            } catch (IOException ex) {
                Logger.getLogger(GsMethods.class.getName()).log(Level.SEVERE, null, ex);
            } finally {
                try {
                    fr.close();
                } catch (IOException ex) {
                    Logger.getLogger(GsMethods.class.getName()).log(Level.SEVERE, null, ex);
                }
            }
        } else {
            System.out.println("File does not exists!");
        }
        return linenumber;
    }

    //for summary html file
    //getting the time required for completing analysis and summary
    public static void addingProcessTime2HTMlFileSummary(int nooflines, String htmlPath, String timeTaken) {
        try {
            List<String> lines = Files.readAllLines(Paths.get(htmlPath), StandardCharsets.UTF_8);
            //   StringBuilder bbb = new StringBuilder(lines.get(nooflines - 2));
            String charCountOfLine = lines.get(nooflines - 2);
            StringBuilder sb = new StringBuilder(lines.get(nooflines - 2));
            int length = charCountOfLine.length() - 4;
            sb.insert(length, timeTaken);

            lines.remove(nooflines - 2);
            lines.add(nooflines - 2, sb.toString()); // for example
            Files.write(Paths.get(htmlPath), lines, StandardCharsets.UTF_8);
        } catch (IOException ex) {
            Logger.getLogger(GsMethods.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    
    public static String stripFileExtension(String fileName) {
        int dotInd = fileName.lastIndexOf('.');
        // if dot is in the first position,
        // we are dealing with a hidden file rather than an extension
        return (dotInd > 0) ? fileName.substring(0, dotInd) : fileName;
    }
}
