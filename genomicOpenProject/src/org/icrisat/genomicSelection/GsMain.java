package org.icrisat.genomicSelection;

import chrriis.dj.nativeswing.swtimpl.components.JWebBrowser;
import com.jtattoo.plaf.aero.AeroLookAndFeel;
import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.HeadlessException;
import java.awt.Image;
import java.awt.Insets;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.InputEvent;
import java.awt.event.KeyEvent;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.Iterator;
import java.util.List;
import java.util.Properties;
import java.util.concurrent.TimeUnit;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.imageio.ImageIO;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JPopupMenu;
import javax.swing.JScrollPane;
import javax.swing.JSplitPane;
import javax.swing.JTree;
import javax.swing.KeyStroke;
import javax.swing.UIManager;
import javax.swing.UnsupportedLookAndFeelException;
import javax.swing.event.TreeSelectionEvent;
import javax.swing.event.TreeSelectionListener;
import javax.swing.filechooser.FileNameExtensionFilter;
import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.DefaultTreeModel;
import javax.swing.tree.MutableTreeNode;
import javax.swing.tree.TreePath;
import org.apache.commons.io.FilenameUtils;
import org.jdesktop.application.Application;
import org.jdesktop.application.SingleFrameApplication;

/**
 *
 * @author mohan
 */
/**
 * GsMain class main class for the application. Which is extended with
 * SingleFrameApplication which has start and stop methods
 */
public class GsMain extends SingleFrameApplication implements ActionListener {

    JFrame frame; //frame for whole application
    JSplitPane splitPane; //split pane which used to add tabbedpane and tree
    DynamicTree dyntree; //Class for all operation related to tree structure
    Analysis analysisobj; //class for the analysis to perform
    OtherFiles otherFile; //class for loading the other files(relationship matrix ,pedigree,population structure)
    ClosableTabbedPane tabbedPane;  //class for closing the selected tab in window
    private BrowserToolBar toolbar; //class for add the toolbar menu which can be pluged any where 
    GsMethods gsMethods = new GsMethods(); //initilizing the Gsmehtod class for calling the methods csvreader and other methods
    DefaultMutableTreeNode genotypeNode, phenotypeNode, realtionshipNode, pedigreeNode, populationNode, resultsNode, logNode; //root nodes for the dynamic tree class 
    JMenu fileMenu = new JMenu("File"); //declaring the filemenu item for jmenu
    JMenu anlysisMenu = new JMenu("Analysis"); //declaring the analysis menu item for jmenu
    JMenu helpMenu = new JMenu("Help"); //declaring the help menu item
    //declaring menuitems for menu bar
    private final JMenuItem OpenItem = new JMenuItem("Open");
    private final JMenuItem OtherFiles = new JMenuItem("Import");
    private final JMenuItem genotypeFile = new JMenuItem("Genotype file");
    private final JMenuItem phenotypeFile = new JMenuItem("Phenotype file");
    private final JMenuItem openProject = new JMenuItem("Open project"); //menuitmes in menu bar
    private final JMenuItem exit = new JMenuItem("Exit"); //menuitmes in menubar
    private final JMenuItem genomicItem = new JMenuItem("Genomic selection"); //menuitmes in menu bar
    private final JMenuItem dataSummaryItem = new JMenuItem("Data summary");
    private final JMenuItem engine = new JMenuItem("Engine"); //for selecting R or Fortran 
    private final JMenuItem saveFile = new JMenuItem("Save As"); //menuitem in menu bar
    private final JMenuItem saveProject = new JMenuItem("Save Project"); //menuitem in menu bar
    private final JMenuItem about = new JMenuItem("About");
    private final JMenuItem manual = new JMenuItem("User manual");
    public java.util.List<String> porSaveList = new ArrayList<>();//arraylist for saving whole project files
    private JButton btnOpen, btnOtherfiles, btnAnalysis, btnDataSummary; //buttons added
    //removed the btnOpenProject,btnRemove,btnSave button
    int analysisCount = 0; //for getting the count of the no of time analysis is runned 
    String lastFileChoosen = ""; //Maintaning the last selected file path for genotype,phenotype,pedigree,population and relationship files
    int StatusForSaving = 0; //0-if yes and save clicked ,1 -if yes and cancel the save  //used in save dialog
    int genoc = 0, phenoc = 0, populationc = 0, resultc = 0;  //using when saveporject done and exiting with and without changes
    boolean statusProSave = false;

    /**
     * startup method to launch the application
     */
    @Override
    protected void startup() {
        System.out.println("Yellow");

        //start region 
        //look and feel of appilication
        Properties props = new Properties();
        props.put("logoString", "ISMU");
        props.put("dynamicLayout", "on");
        props.put("linuxStyleScrollBar", "on");
        props.put("windowDecoration", "on");
        props.put("tooltipBackgroundColor", "");
        props.put("menuTextFont", "Arial BOLD 15");
        props.put("textAntiAliasing", "on");
        AeroLookAndFeel.setCurrentTheme(props);
        try {
            UIManager.setLookAndFeel("com.jtattoo.plaf.aero.AeroLookAndFeel");
        } catch (ClassNotFoundException | InstantiationException | IllegalAccessException | UnsupportedLookAndFeelException ex) {
            Logger.getLogger(GsMain.class.getName()).log(Level.SEVERE, null, ex);
        }
       //end region of look and feel 

        //start region 
        //setting the content of menubar 
        JMenuBar menuBar = new JMenuBar(); //initilizing the jmenubar
        menuBar.add(fileMenu);
        fileMenu.setFont(new Font("DejaVu Sans", Font.PLAIN, 14));
        menuBar.add(anlysisMenu);
        anlysisMenu.setFont(new Font("DejaVu Sans", Font.PLAIN, 14));
        menuBar.add(helpMenu);
        helpMenu.setFont(new Font("DejaVu Sans", Font.PLAIN, 14));

        //adding the menuitems to the file menu
        fileMenu.add(OpenItem); //adding open menu item for file menu 
        OpenItem.setFont(new Font("DejaVu Sans", Font.PLAIN, 14));
        OpenItem.setMnemonic(KeyEvent.VK_O);
        OpenItem.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_O, InputEvent.CTRL_DOWN_MASK));

        fileMenu.add(OtherFiles); //adding open menu item for file menu 
        OtherFiles.setFont(new Font("DejaVu Sans", Font.PLAIN, 14));
        OtherFiles.setMnemonic(KeyEvent.VK_T);
        OtherFiles.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_T, InputEvent.CTRL_DOWN_MASK));

        fileMenu.add(genotypeFile);
        genotypeFile.setFont(new Font("DejaVu Sans", Font.PLAIN, 14));

        fileMenu.add(phenotypeFile);
        phenotypeFile.setFont(new Font("DejaVu Sans", Font.PLAIN, 14));

        fileMenu.add(openProject); //adding open menu item for file menu 
        openProject.setFont(new Font("DejaVu Sans", Font.PLAIN, 14));

        fileMenu.add(saveFile);
        saveFile.setFont(new Font("DejaVu Sans", Font.PLAIN, 14));
        saveFile.setMnemonic(KeyEvent.VK_S);
        saveFile.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_S, InputEvent.CTRL_DOWN_MASK));

        fileMenu.add(saveProject);
        saveProject.setFont(new Font("DejaVu Sans", Font.PLAIN, 14));
        fileMenu.setFont(new Font("DejaVu Sans", Font.PLAIN, 14));

        fileMenu.add(exit);  //adding exit menu item for file menu
        exit.setFont(new Font("DejaVu Sans", Font.PLAIN, 14));
        exit.setMnemonic(KeyEvent.VK_F4);
        exit.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_F4, InputEvent.CTRL_DOWN_MASK));

        anlysisMenu.add(dataSummaryItem);
        dataSummaryItem.setFont(new Font("DejaVu Sans", Font.PLAIN, 14));
        dataSummaryItem.setMnemonic(KeyEvent.VK_D);
        dataSummaryItem.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_D, InputEvent.CTRL_DOWN_MASK));

        //adding the menuitmes to the anlaysis menu 
        anlysisMenu.add(genomicItem);  //adding the gs menu  analysis menu
        genomicItem.setFont(new Font("DejaVu Sans", Font.PLAIN, 14));
        genomicItem.setMnemonic(KeyEvent.VK_G);
        genomicItem.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_G, InputEvent.CTRL_DOWN_MASK));

        anlysisMenu.add(engine);
        engine.setFont(new Font("DejaVu Sans", Font.PLAIN, 14));
        engine.setMnemonic(KeyEvent.VK_E);
        engine.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_E, InputEvent.CTRL_DOWN_MASK));

        helpMenu.add(manual);
        manual.setFont(new Font("DejaVu Sans", Font.PLAIN, 14));
        manual.setMnemonic(KeyEvent.VK_M);
        manual.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_M, InputEvent.CTRL_DOWN_MASK));

        //adding the menuitems to the help menu
        helpMenu.add(about);
        about.setFont(new Font("DejaVu Sans", Font.PLAIN, 14));
        about.setMnemonic(KeyEvent.VK_A);
        about.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_A, InputEvent.CTRL_DOWN_MASK));

        //adding action listener to the all menuitems of  Menu bar
        OpenItem.addActionListener(this);
        openProject.addActionListener(this);
        OtherFiles.addActionListener(this);
        genotypeFile.addActionListener(this);
        phenotypeFile.addActionListener(this);
        saveFile.addActionListener(this);
        saveProject.addActionListener(this);
        exit.addActionListener(this); //adding action listener for exit menu item
        genomicItem.addActionListener(this);
        about.addActionListener(this);
        manual.addActionListener(this);
        dataSummaryItem.addActionListener(this);
        engine.addActionListener(this);
        //end region menubar

        //start region 
        //setting the content of the toolbar in appicaltion
        toolbar = new BrowserToolBar(); //toolbar 
        btnOpen = new JButton("Open");
        btnOpen.setToolTipText("Selecting the Genotype,Phenotype files and result directory");
        //adding the image to the button
        try {
            Image img = ImageIO.read(getClass().getResource("/org/icrisat/genomicSelection/image/open-file-icon.png"));
            btnOpen.setIcon(new ImageIcon(img));
        } catch (IOException ex) {
            Logger.getLogger(GsMain.class.getName()).log(Level.SEVERE, null, ex);
        }
        btnOtherfiles = new JButton("Import");
        btnOtherfiles.setToolTipText("Selecting the population file");
        //adding the image to the button
        try {
            Image img = ImageIO.read(getClass().getResource("/org/icrisat/genomicSelection/image/other.gif"));
            btnOtherfiles.setIcon(new ImageIcon(img));
        } catch (IOException ex) {
            Logger.getLogger(GsMain.class.getName()).log(Level.SEVERE, null, ex);
        }

        btnDataSummary = new JButton("Data summary");
        btnDataSummary.setToolTipText("Cleaning of data");
        //adding the image to the button
        try {
            Image img = ImageIO.read(getClass().getResource("/org/icrisat/genomicSelection/image/filter-icon1.png"));
            btnDataSummary.setIcon(new ImageIcon(img));
        } catch (IOException ex) {
            Logger.getLogger(GsMain.class.getName()).log(Level.SEVERE, null, ex);
        }

        btnAnalysis = new JButton("Analysis");
        btnAnalysis.setToolTipText("Analysis");
        //adding the image to the button
        try {
            Image img1 = ImageIO.read(getClass().getResource("/org/icrisat/genomicSelection/image/run_cropped.jpg"));
            btnAnalysis.setIcon(new ImageIcon(img1));
        } catch (IOException ex) {
            Logger.getLogger(GsMain.class.getName()).log(Level.SEVERE, null, ex);
        }
        btnOpen.addActionListener(this);
        btnOtherfiles.addActionListener(this);
        btnDataSummary.addActionListener(this);
        btnAnalysis.addActionListener(this);
        toolbar.add(btnOpen);
        toolbar.add(btnOtherfiles);
        toolbar.add(btnDataSummary);
        toolbar.add(btnAnalysis);
        //end region

        //start region
        //initialization of jframe 
        //adding the menubar,toolbar,home panel in it 
        //setting the size based on the screen size
        frame = new JFrame();
        frame.setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
        frame.setTitle("ISMU 2.0 pipeline");
        frame.setIconImage(Toolkit.getDefaultToolkit().getImage(getClass().getResource("/org/icrisat/genomicSelection/image/log_icrisat.png")));
        //action listener for the jframe
        frame.setBackground(Color.WHITE);

        //listener to close the application when click on close button 
        frame.addWindowListener(new WindowAdapter() {
            @Override
            public void windowClosing(WindowEvent e) {
                shutdown();
            }
        });

        //setting the size of the appliaction based on the screen resolution
        frame.pack();
        Insets insets = frame.getInsets();
        //System.out.println(new Dimension(insets.left + insets.right, insets.top + insets.bottom));
        frame.setSize(new Dimension(insets.left + insets.right, insets.top + insets.bottom));
        splitPane = new JSplitPane(); //initlizing the splitPanel 
        splitPane.setBackground(Color.WHITE);
        splitPane.setDividerLocation(200);//setting the divider location at 200 pixcel location
        dyntree = new DynamicTree(); //initlizing the Dynamic tree class 

        //method adding the default nodes to the tree 
        populateTree(dyntree);
        try {
            dyntree.addObject(logNode, "log.txt"); //node to see the log text file
        } catch (Exception e) {
            e.getMessage();
        }
        //initilizing closableTabbedPane class  
        tabbedPane = new ClosableTabbedPane();
        tabbedPane.setBackground(Color.WHITE);
        //adding the home panel to the tabbedpane
        Home homePanel = new Home();
        tabbedPane.add(homePanel);

        frame.setJMenuBar(menuBar); //adding the menu bar to the frame
        frame.add(toolbar, BorderLayout.NORTH); //adding the toolbar to frame on north of the panel
        //splitPane.setLeftComponent(dyntree); //is added at the dynamic tree initilization //added at the tree initilization
        splitPane.setRightComponent(tabbedPane); //adding the tabbedpane to the right side split panel
        frame.add(splitPane); //adding the split pane to the frame
        //resize the frame which adjust the components 
        frame.pack();
        //setting the frame to visiable
        frame.setVisible(true);
        System.out.println(Constant.cmdOpenProject);
        if ("noargsismu".equals(Constant.cmdOpenProject)) {
            //donothing
        } else {
            commanndLineProjectOpening(Constant.cmdOpenProject);
        }
    }
//end region of frame

    //start region 
    //method to save the present working state of application
    @Override
    protected void shutdown() {
        Constant.genoCount = Constant.genoList.size();
        Constant.phenCount = Constant.phenoList.size();
        Constant.populationCount = Constant.populationList.size();
        Constant.resultFileCount = Constant.resultList.size();
        if (statusProSave == true) {
            if (genoc == Constant.genoCount && phenoc == Constant.phenCount && resultc == Constant.resultFileCount) {
                System.exit(0);
            } else {
                appenedSavedChangesToIsmu();
                System.exit(0);
            }
        } else {
            if ("empty".equals(Constant.resultdirectory)) {
                System.exit(0); //to close the applicaiton-
            } else {
                if (!"".equals(Constant.resultdirectory) && Constant.genoList.isEmpty() && Constant.phenoList.isEmpty() && Constant.populationList.isEmpty() && Constant.outputResultsList.isEmpty()) {
                    //as all files are deleted nothing to save
                    System.exit(0);
                } else {
                    if (Constant.projectName == null) {
                        int option = JOptionPane.showOptionDialog(frame, "Do you want to save the project", "Save Dialog", JOptionPane.YES_NO_CANCEL_OPTION,
                                JOptionPane.WARNING_MESSAGE, null, null, null);
                        if (option == JOptionPane.YES_OPTION) {
                            saveProject();
                            if (StatusForSaving == 0) {  ///check whethere the save or cancel is clicked when exit and save
                                System.exit(0);
                            }
                        }
                        if (option == JOptionPane.NO_OPTION) {
                            System.exit(0); //close application without saving state
                        }
                        if (option == JOptionPane.CANCEL_OPTION) {
                            //donothing
                        }
                    } else {
                        if (genoc == Constant.genoCount && phenoc == Constant.phenCount && resultc == Constant.resultFileCount) {
                            System.exit(0);
                        } else {
                            int option = JOptionPane.showOptionDialog(frame, "Do you want to save changes made in project", "Save Dialog", JOptionPane.YES_NO_CANCEL_OPTION,
                                    JOptionPane.WARNING_MESSAGE, null, null, null);
                            PrintWriter printwriter1 = null;
                            try {
                                printwriter1 = new PrintWriter(Constant.projectName);
                            } catch (FileNotFoundException ex) {
                                Logger.getLogger(GsMain.class.getName()).log(Level.SEVERE, null, ex);
                            }
                            if (option == JOptionPane.YES_OPTION) {
                                String pojPath = Constant.projectName;
                                if (!pojPath.isEmpty()) {
                                    String resultdir = Constant.resultdirectory;
                                    if (!"".equals(resultdir)) {
                                        porSaveList.add(resultdir);
                                        printwriter1.println(resultdir);
                                    }
                                    String browsepath = Constant.browsepath;
                                    if (!"".equals(browsepath)) {
                                        porSaveList.add(browsepath);
                                        printwriter1.println(browsepath);
                                    }
                                    String genofilenames = "";
                                    Iterator genoItr = Constant.genoList.iterator();
                                    Constant.genoCount = Constant.genoList.size();
                                    while (genoItr.hasNext()) {
                                        genofilenames = genofilenames + genoItr.next() + ",";
                                    }
                                    porSaveList.add(genofilenames);
                                    printwriter1.println(genofilenames);

                                    String phenofilenames = "";
                                    Iterator phenoItr = Constant.phenoList.iterator();
                                    Constant.phenCount = Constant.phenoList.size();
                                    while (phenoItr.hasNext()) {
                                        phenofilenames = phenofilenames + phenoItr.next() + ",";
                                    }
                                    porSaveList.add(phenofilenames);
                                    printwriter1.println(phenofilenames);

                                    String populationfilenames = "";
                                    Iterator populationItr = Constant.populationList.iterator();
                                    Constant.populationCount = Constant.populationList.size();
                                    while (populationItr.hasNext()) {
                                        populationfilenames = populationfilenames + populationItr.next() + ",";
                                    }
                                    porSaveList.add(populationfilenames);
                                    printwriter1.println(populationfilenames);

                                    String output_files = "";
                                    Iterator outItr = Constant.outputResultsList.iterator();
                                    Constant.resultFileCount = Constant.outputResultsList.size();
                                    while (outItr.hasNext()) {
                                        output_files = output_files + outItr.next() + ",";
                                    }
                                    porSaveList.add(output_files);
                                    printwriter1.println(output_files);

                                    int enigneSelected = Constant.engineSelection;
                                    printwriter1.println(enigneSelected);
                                }
                                printwriter1.close();
                                StatusForSaving = 0;
                                System.exit(0);
                            }
                            if (option == JOptionPane.NO_OPTION) {
                                System.exit(0);
                            }
                            if (option == JOptionPane.CANCEL_OPTION) {
                                statusProSave = false;
                            }
                        }
                    }
                }
            }
        }
    }
//end region shutdown()

    //start region 
    //save current state of the application 
    public void saveProject() {
        PrintWriter printwriter = null;
        try {
            JFileChooser chooser = new JFileChooser();
            chooser.removeChoosableFileFilter(chooser.getFileFilter());
            FileNameExtensionFilter filter = new FileNameExtensionFilter("ISMU (Integrated SNP Mining & Utilization Pipeline)(*.ismu)", "ismu");
            chooser.setFileFilter(filter);

            chooser.setCurrentDirectory(new File(Constant.resultdirectory));
            int retrival = chooser.showSaveDialog(null);
            if (retrival == JFileChooser.APPROVE_OPTION) {
                String pojPath = chooser.getSelectedFile() + ".ismu";
                Constant.projectName = pojPath;
                printwriter = new PrintWriter(new FileWriter(chooser.getSelectedFile() + ".ismu"));
                if (!pojPath.isEmpty()) {
                    String resultdir = Constant.resultdirectory;
                    if (!"".equals(resultdir)) {
                        porSaveList.add(resultdir);
                        printwriter.println(resultdir);
                    }
                    String browsepath = Constant.browsepath;
                    if (!"".equals(browsepath)) {
                        porSaveList.add(browsepath);
                        printwriter.println(browsepath);
                    }
                    String genofilenames = "";
                    Iterator genoItr = Constant.genoList.iterator();
                    Constant.genoCount = Constant.genoList.size();
                    while (genoItr.hasNext()) {
                        genofilenames = genofilenames + genoItr.next() + ",";
                    }
                    porSaveList.add(genofilenames);
                    printwriter.println(genofilenames);

                    String phenofilenames = "";
                    Iterator phenoItr = Constant.phenoList.iterator();
                    Constant.phenCount = Constant.phenoList.size();
                    while (phenoItr.hasNext()) {
                        phenofilenames = phenofilenames + phenoItr.next() + ",";
                    }
                    porSaveList.add(phenofilenames);
                    printwriter.println(phenofilenames);

                    String populationfilenames = "";
                    Iterator populationItr = Constant.populationList.iterator();
                    Constant.populationCount = Constant.populationList.size();
                    while (populationItr.hasNext()) {
                        populationfilenames = populationfilenames + populationItr.next() + ",";
                    }
                    porSaveList.add(populationfilenames);
                    printwriter.println(populationfilenames);

                    String output_files = "";
                    Iterator outItr = Constant.outputResultsList.iterator();
                    Constant.resultFileCount = Constant.outputResultsList.size();
                    while (outItr.hasNext()) {
                        output_files = output_files + outItr.next() + ",";
                    }
                    porSaveList.add(output_files);
                    printwriter.println(output_files);
                    int enigneSelected = Constant.engineSelection;
                    printwriter.println(enigneSelected);
                    printwriter.println(Constant.projectName);
                    printwriter.println(Constant.genoCount);
                    printwriter.println(Constant.phenCount);
                    printwriter.println(Constant.populationCount);
                    printwriter.println(Constant.resultFileCount);
                }
                printwriter.close();
                StatusForSaving = 0;
                statusProSave = true;
            }
            if (retrival == JFileChooser.CANCEL_OPTION) {
                StatusForSaving = 1;
                statusProSave = false;
            }
        } catch (IOException ex) {
            Logger.getLogger(GsMain.class.getName()).log(Level.SEVERE, null, ex);
        }

    }
//end region saveproject

    //appeding changes for the saved project 
    public void appenedSavedChangesToIsmu() {
        try {
            PrintWriter printwriter = null;
            printwriter = new PrintWriter(new FileWriter(Constant.projectName));
            String resultdir = Constant.resultdirectory;
            if (!"".equals(resultdir)) {
                porSaveList.add(resultdir);
                printwriter.println(resultdir);
            }
            String browsepath = Constant.browsepath;
            if (!"".equals(browsepath)) {
                porSaveList.add(browsepath);
                printwriter.println(browsepath);
            }
            String genofilenames = "";
            Iterator genoItr = Constant.genoList.iterator();
            Constant.genoCount = Constant.genoList.size();
            while (genoItr.hasNext()) {
                genofilenames = genofilenames + genoItr.next() + ",";
            }
            porSaveList.add(genofilenames);
            printwriter.println(genofilenames);

            String phenofilenames = "";
            Iterator phenoItr = Constant.phenoList.iterator();
            Constant.phenCount = Constant.phenoList.size();
            while (phenoItr.hasNext()) {
                phenofilenames = phenofilenames + phenoItr.next() + ",";
            }
            porSaveList.add(phenofilenames);
            printwriter.println(phenofilenames);

            String populationfilenames = "";
            Iterator populationItr = Constant.populationList.iterator();
            Constant.populationCount = Constant.populationList.size();
            while (populationItr.hasNext()) {
                populationfilenames = populationfilenames + populationItr.next() + ",";
            }
            porSaveList.add(populationfilenames);
            printwriter.println(populationfilenames);

            String output_files = "";
            Iterator outItr = Constant.outputResultsList.iterator();
            Constant.resultFileCount = Constant.outputResultsList.size();
            while (outItr.hasNext()) {
                output_files = output_files + outItr.next() + ",";
            }
            porSaveList.add(output_files);
            printwriter.println(output_files);
            //  System.out.println("appneded file list out " + output_files+"   , ");

            int enigneSelected = Constant.engineSelection;
            printwriter.println(enigneSelected);

            printwriter.println(Constant.projectName);
            printwriter.println(Constant.genoCount);
            printwriter.println(Constant.phenCount);
            printwriter.println(Constant.populationCount);
            printwriter.println(Constant.resultFileCount);
            printwriter.close();
            StatusForSaving = 0;
            statusProSave = true;
        } catch (IOException ex) {
            Logger.getLogger(GsMain.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    //start region 
    //setting the actionlistners for the all the toolbar buttons and menuitems
    @Override
    public void actionPerformed(ActionEvent e) {
        //get which button or menuitem clicked 
        Object source = e.getSource();

        /*
         *This opens a Load file dilog box to get the path of genotype,phenotype and result directory
         * Then checking the path adding the genotype and phenotye file to the tree,tabbedpane
         * copy genotype and phenotype file to the result directory
         * also add the file names to the lists 
         */
        if (source == OpenItem || source == btnOpen) {
            try {
                OpenDial dialog1 = new OpenDial(frame, true);
                dialog1.setLocationRelativeTo(frame);
                dialog1.setVisible(true);

                if (!"null".equals(Constant.genotype)) {
                    String genopath = Constant.genotype;
                    String phenopath = Constant.phenotype;

                    String genoname = new File(genopath).getName();
                    String phenoname = new File(phenopath).getName();

                    Date date = new Date();
                    SimpleDateFormat dateFormat = new SimpleDateFormat("mmss");
                    SimpleDateFormat dateFormat1 = new SimpleDateFormat("mmss");
                    String genotemp = GsMethods.stripFileExtension(genoname) + dateFormat.format(date) + ".csv";
                    String phenotemp = GsMethods.stripFileExtension(phenoname) + dateFormat1.format(date) + ".csv";

                    String genoresultdir = Constant.resultdirectory + genotemp;
                    String phenoresultdir = Constant.resultdirectory + phenotemp;

                    //copying data to result folder selected
                    if (!"".equals(genopath)) {
                        try {
                            CopyFile.copy(genopath, genoresultdir);
                        } catch (IOException ex) {
                            Logger.getLogger(GsMain.class.getName()).log(Level.SEVERE, null, ex);
                        }
                    }
                    if (!"".equals(phenopath)) {
                        try {
                            CopyFile.copy(phenopath, phenoresultdir);
                        } catch (IOException ex) {
                            Logger.getLogger(GsMain.class.getName()).log(Level.SEVERE, null, ex);
                        }
                    }

                    JPanel csvreaderGeno = gsMethods.csvReaderSNP(genoresultdir);
                    tabbedPane.add(csvreaderGeno);
                    Constant.iconNode = genotemp;
                    dyntree.addObject(genotypeNode, genotemp, Boolean.TRUE);
                    int tabCountGeno = tabbedPane.getTabCount();
                    tabbedPane.setSelectedIndex(tabCountGeno - 1);

                    JPanel csvreaderPheno = gsMethods.csvReaderSNP(phenoresultdir);
                    tabbedPane.add(csvreaderPheno);
                    Constant.iconNode = phenotemp;
                    dyntree.addObject(phenotypeNode, phenotemp, Boolean.TRUE);
                    int tabCountPheno = tabbedPane.getTabCount();
                    tabbedPane.setSelectedIndex(tabCountPheno - 1);

                    Constant.genoList.add(genotemp);
                    Constant.phenoList.add(phenotemp);

                    //create a log file for pipeline
                    File logFile = new File(Constant.resultdirectory + "log.txt");
                    try {
                        logFile.createNewFile();
                    } catch (IOException ex) {
                        Logger.getLogger(GsMain.class.getName()).log(Level.SEVERE, null, ex);
                    }
                }
            } catch (Exception openException) {
                openException.getMessage();
            }
        }
        /*
         * This opens OtherFiles dialog box.
         * Display the selected files in tabbedpane and adds the names of the file to the coresspending node 
         * copys the selected file data to the result directory 
         */
        if (source == OtherFiles || source == btnOtherfiles) {
            try {
                if (Constant.genoList.size() > 0 && Constant.phenoList.size() > 0) {
                    otherFile = new OtherFiles(frame, true);
                    otherFile.setVisible(true);
                    if (Constant.populationPath != null) {
                        String populatioName = new File(Constant.populationPath).getName();
                        String populationdir = Constant.resultdirectory + populatioName;
                        dyntree.addObject(populationNode, populatioName, Boolean.TRUE);
                        Constant.populationList.add(populatioName);
                        JPanel csvreaderPopulation = gsMethods.csvReaderSNP(populationdir);
                        tabbedPane.add(csvreaderPopulation);
                        int count = tabbedPane.getTabCount();
                        tabbedPane.setSelectedIndex(count - 1);
                    }
                    if (Constant.haploPath != null) {
                        String haploName = new File(Constant.haploPath).getName();
                        String haploDir = Constant.resultdirectory + haploName;
                        dyntree.addObject(genotypeNode, haploName, Boolean.TRUE);
                        Constant.genoList.add(haploName);
                        JPanel csvreaderHaplo = gsMethods.csvReaderSNP(haploDir);
                        tabbedPane.add(csvreaderHaplo);
                        int count = tabbedPane.getTabCount();
                        tabbedPane.setSelectedIndex(count - 1);
                    }
                } else {
                    JOptionPane.showMessageDialog(frame, "Import atleast a Genotype and Phenotype files");
                }
            } catch (HeadlessException otherException) {
                otherException.getMessage();
            }
        }

        /*
         * Opens a Analysis dialog box 
         * Then checks for the if analysis cancled if not opens the another dailog box Analysis progress display. 
         * Then display the genarated html files in tabbedpanes and also the names of the files to the tree under Results node
         */
        if (source == btnAnalysis || source == genomicItem) {
        	System.out.println("blue");
            try {
                if (Constant.genoList.size() > 0 && Constant.phenoList.size() > 0) {
                    analysisobj = new Analysis(frame, true);
                    analysisobj.setVisible(true);
                    if (Constant.analysisCancled == false) {
                        long startTime = System.currentTimeMillis();
                        ProgressDisplay progressDisplay = new ProgressDisplay(frame, true);
                        progressDisplay.setVisible(true);
                        long stopTime = System.currentTimeMillis();
                        long difference = stopTime - startTime;
                        String res = "";
                        long days = TimeUnit.MILLISECONDS.toDays(difference);
                        long hours1 = TimeUnit.MILLISECONDS.toHours(difference) - TimeUnit.DAYS.toHours(TimeUnit.MILLISECONDS.toDays(difference));
                        long minutes1 = TimeUnit.MILLISECONDS.toMinutes(difference) - TimeUnit.HOURS.toMinutes(TimeUnit.MILLISECONDS.toHours(difference));
                        long seconds1 = TimeUnit.MILLISECONDS.toSeconds(difference) - TimeUnit.MINUTES.toSeconds(TimeUnit.MILLISECONDS.toMinutes(difference));
                        if (days == 0) {
                            res = String.format(" (Total Time: %02d Hr : %02d Min : %02d Sec)", hours1, minutes1, seconds1);
                        } else {
                            res = String.format(" (Total Time: %dd Days : %02d Hr : %02d Min : %02d Sec)", days, hours1, minutes1, seconds1);
                        }

                        //itterating the loop of results generated by running the analysis on traits 
                        String presentWorkingDirectory = System.getProperty("user.dir");
                        for (String storeName : Constant.resultGenratedList) {
                            if ("error.htm".equals(storeName)) {
                                JWebBrowser errorHtm = gsMethods.browser(presentWorkingDirectory + "/error.htm");
                                tabbedPane.add(errorHtm);
                                tabbedPane.setSelectedIndex(tabbedPane.getTabCount() - 1);
                            } else {
                                File file = new File(Constant.resultdirectory + storeName);
                                if (file.exists()) {
                                    int noofLinesInHTML = GsMethods.noofLinesInHTML(Constant.resultdirectory + storeName);
                                    GsMethods.addingProcessTime2HTMlFile(noofLinesInHTML, Constant.resultdirectory + storeName, res);
                                    JWebBrowser htmlPage = gsMethods.browser(Constant.resultdirectory + storeName);
                                    tabbedPane.add(htmlPage);
                                    dyntree.addObject(resultsNode, storeName, Boolean.TRUE);
                                    Constant.outputResultsList.add(storeName);
                                    Constant.resultList.add(storeName);
                                    tabbedPane.setSelectedIndex(tabbedPane.getTabCount() - 1);
                                }
                            }
                        }
                    }
                } else {
                    JOptionPane.showMessageDialog(frame, "Import atleast a Genotype and Phenotype files");
                }
            } catch (HeadlessException anlaysisException) {
                anlaysisException.getMessage();
            }
        }

        /*
         *  open save dialog box in which we can specify required file name and directory 
         *  if file extension is  html then it convert to pdf and save in directory specified
         *  if file file extension is .csv then it copys the file to the selected directory with given name 
         */
        if (source == saveFile) {
            try {
                if (Constant.genoList.size() > 0 && Constant.phenoList.size() > 0) {
                    int tabIndexNumber = tabbedPane.getSelectedIndex();
                    String resultdir = Constant.resultdirectory;
                    String name = tabbedPane.getTabTitleAt(tabIndexNumber);
                    String ext = FilenameUtils.getExtension(name);
                    if (!"Home".equals(name)) {
                        if ("htm".equals(ext) || "html".equals(ext)) {
                            String path = resultdir + name;
                            JFileChooser fileChooser = new JFileChooser();
                            fileChooser.setDialogTitle("Save as");
                            fileChooser.setSelectedFile(new File(FilenameUtils.removeExtension(path)));
                            FileNameExtensionFilter filter = new FileNameExtensionFilter("PDF", "pdf");
                            fileChooser.setFileFilter(filter);
                            int userSelection = fileChooser.showSaveDialog(frame);
                            if (userSelection == JFileChooser.APPROVE_OPTION) {
                                String fileNameGiven = "";
                                String presentWorkingDirectory = System.getProperty("user.dir");
                                OSCheck.OSType os = OSCheck.getOperatingSystemType();
                                String cmd = null;
                                if ("Window".equals(os.toString())) {
                                    fileNameGiven = "\"" + fileChooser.getSelectedFile().toString() + ".pdf\"";
                                    cmd = "\"" + presentWorkingDirectory + "\\wkhtmltopdf\" " + "\"file:///" + Constant.resultdirectory + name + "\" " + fileNameGiven;
                                }
                                if ("Linux".equals(os.toString())) {
                                    fileNameGiven = fileChooser.getSelectedFile().toString() + ".pdf";
                                    cmd = presentWorkingDirectory + "/wkhtmltopdf " + "file:///" + Constant.resultdirectory + name + " " + fileNameGiven;
                                }
                                try {
                                    Process processCall = Runtime.getRuntime().exec(cmd);
                                    try {
                                        processCall.waitFor();
                                        JOptionPane.showMessageDialog(frame, "  Pdf conversion done  ");

                                    } catch (InterruptedException ex) {
                                        Logger.getLogger(GsMain.class.getName()).log(Level.SEVERE, null, ex);
                                    }
                                } catch (IOException ex) {
                                    Logger.getLogger(GsMain.class.getName()).log(Level.SEVERE, null, ex);
                                }
                            }
                        } else {
                            String path = resultdir + name;
                            Constant.tabbselected = path;
                            JFileChooser fileChooser = new JFileChooser();
                            fileChooser.setDialogTitle("Save as");
                            fileChooser.setSelectedFile(new File(path));
                            int userSelection = fileChooser.showSaveDialog(frame);
                            if (userSelection == JFileChooser.APPROVE_OPTION) {
                                File fileToSave = fileChooser.getSelectedFile();
                                CopyFile.copy(path, fileToSave.getAbsolutePath());
                            }
                            if (userSelection == JFileChooser.CANCEL_OPTION) {
                                //donothing in this
                            }
                        }
                    }
                } else {
                    JOptionPane.showMessageDialog(frame, "Open any file to save");
                }
            } catch (HeadlessException | IOException saveException) {
                saveException.getMessage();
            }
        }

        /*
         * for the menu item exit 
         * calls the shutdown method  
         */
        if (source == exit) {
            if ("empty".equals(Constant.resultdirectory)) {
                System.exit(0);
            } else {
                try {
                    shutdown();
                } catch (Exception exitException) {
                    exitException.getMessage();
                }
            }
        }

        /*
         * for menu item saveProject
         * calls the saveproject method 
         */
        if (source == saveProject) {
            if (Constant.genoList.isEmpty() == false && Constant.phenoList.isEmpty() == false) {
                saveProject();
            }
            genoc = Constant.genoCount;
            phenoc = Constant.phenCount;
            populationc = Constant.populationCount;
            resultc = Constant.resultFileCount;
        }

        /*
         * opens dilaog box to select the project folder 
         * Reads the pro.ismu files 
         * opens the prevoius state of the project 
         */
        if (source == openProject) {
            if ("empty".equals(Constant.resultdirectory)) {
                try {
                    statusProSave = true; //setting the project as saved beacuse it is loading existing project
                    JFileChooser chooser;
                    lastFileChoosen = Constant.browsepath;
                    chooser = new JFileChooser();
                    chooser.setDialogTitle("Select the project file from the ISMU 2.0 directory");
                    chooser.removeChoosableFileFilter(chooser.getFileFilter()); //remove all file selection
                    FileNameExtensionFilter filter = new FileNameExtensionFilter("ISMU (Integrated SNP Mining & Utilization Pipeline)(*.ismu)", "ismu");
                    chooser.setFileFilter(filter);  //selecting the .ismu project
                    chooser.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
                    int r = chooser.showOpenDialog(frame);
                    if (r == JFileChooser.APPROVE_OPTION) {
                        String fileSelected = chooser.getSelectedFile().getPath();
                        int counter = 0;
                        File f = new File(fileSelected);
                        boolean exists = f.exists();
                        if (exists == true) {
                        }
                        BufferedReader in = null;
                        try {
                            OSCheck.OSType os = OSCheck.getOperatingSystemType();
                            if ("Window".equals(os.toString())) {
                                in = new BufferedReader(new FileReader(fileSelected));
                            }
                            if ("Linux".equals(os.toString())) {
                                in = new BufferedReader(new FileReader(fileSelected));
                            }
                            String[] str;
                            try {
                                while (in.ready()) {
                                    String dataLine = in.readLine();
                                    if (!dataLine.isEmpty()) {
                                        if (counter == 0) {
                                            Constant.resultdirectory = dataLine;
                                            counter++;
                                        } else if (counter >= 1) {
                                            if (counter == 1) {
                                                Constant.browsepath = dataLine;
                                            }
                                            str = dataLine.split(",");
                                            if (counter == 2) {
                                                for (String str1 : str) {
                                                    Constant.genotype = str1;
                                                    String genopath = Constant.resultdirectory + str1;
                                                    JPanel csvreaderGeno = gsMethods.csvreader(genopath);
                                                    tabbedPane.add(csvreaderGeno);
                                                    dyntree.addObject(genotypeNode, str1, Boolean.TRUE);
                                                    int count1 = tabbedPane.getTabCount();
                                                    tabbedPane.setSelectedIndex(count1 - 1);
                                                    Constant.genoList.add(str1);
                                                }
                                            }
                                            if (counter == 3) {
                                                for (String str1 : str) {
                                                    Constant.phenotype = str1;
                                                    String phenoPath = Constant.resultdirectory + str1;
                                                    JPanel csvreaderPheno = gsMethods.csvreader(phenoPath);
                                                    tabbedPane.add(csvreaderPheno);
                                                    dyntree.addObject(phenotypeNode, str1, Boolean.TRUE);
                                                    int count2 = tabbedPane.getTabCount();
                                                    tabbedPane.setSelectedIndex(count2 - 1);
                                                    Constant.phenoList.add(str1);
                                                }
                                            }
                                            if (counter == 4) {
                                                for (String str1 : str) {
                                                    Constant.population = str1;
                                                    String populationPath = Constant.resultdirectory + str1;
                                                    JPanel csvreaderPopulation = gsMethods.csvreader(populationPath);
                                                    tabbedPane.add(csvreaderPopulation);
                                                    dyntree.addObject(populationNode, str1, Boolean.TRUE);
                                                    int count5 = tabbedPane.getTabCount();
                                                    tabbedPane.setSelectedIndex(count5 - 1);
                                                    Constant.populationList.add(str1);
                                                }
                                            }
                                            if (counter == 5) {
                                                for (String str1 : str) {
                                                    String path = Constant.resultdirectory + str1;
                                                    String extension = FilenameUtils.getExtension(path);
                                                    if ("htm".equals(extension) || "html".equals(extension)) {
                                                        JWebBrowser browserOpen = gsMethods.browser(path);
                                                        tabbedPane.add(browserOpen);
                                                        dyntree.addObject(resultsNode, str1, Boolean.TRUE);
                                                        int count6 = tabbedPane.getTabCount();
                                                        tabbedPane.setSelectedIndex(count6 - 1);
                                                        Constant.outputResultsList.add(str1);
                                                        Constant.resultList.add(str1);
                                                    }
                                                    if ("csv".equals(extension)) {
                                                        JPanel csvreader = gsMethods.csvreader(path);
                                                        tabbedPane.add(csvreader);
                                                        dyntree.addObject(resultsNode, str1, Boolean.TRUE);
                                                        int count6 = tabbedPane.getTabCount();
                                                        tabbedPane.setSelectedIndex(count6 - 1);
                                                        Constant.outputResultsList.add(str1);
                                                        Constant.resultList.add(str1);
                                                    }
                                                }
                                            }
                                            if (counter == 6) {
                                                Constant.engineSelection = Integer.parseInt(dataLine);
                                            }
                                            if (counter == 7) {
                                                Constant.projectName = dataLine;
                                            }
                                            if (counter == 8) {
                                                genoc = Integer.parseInt(dataLine);
                                            }
                                            if (counter == 9) {
                                                phenoc = Integer.parseInt(dataLine);
                                            }
                                            if (counter == 10) {
                                                //  Constant.RandomForestCount = Integer.parseInt(dataLine);
                                            }
                                            if (counter == 11) {
                                                resultc = Integer.parseInt(dataLine);
                                            }
                                            counter++;
                                        }
                                    } else {
                                        counter++;
                                    }
                                }
                            } catch (IOException ex) {
                                Logger.getLogger(GsMain.class.getName()).log(Level.SEVERE, null, ex);
                            }
                            try {
                                in.close();
                            } catch (IOException ex) {
                                Logger.getLogger(GsMain.class.getName()).log(Level.SEVERE, null, ex);
                            }
                        } catch (FileNotFoundException ex) {
                            Logger.getLogger(GsMain.class.getName()).log(Level.SEVERE, null, ex);
                        } finally {
                            try {
                                in.close();
                            } catch (IOException ex) {
                                Logger.getLogger(GsMain.class.getName()).log(Level.SEVERE, null, ex);
                            }
                        }
                    }
                    if (r == JFileChooser.CANCEL_OPTION) {
                        //do nothing when clicked 
                    }
                } catch (HeadlessException | NumberFormatException openExecption) {
                    openExecption.getMessage();
                }
            } else {
                JOptionPane.showMessageDialog(frame, "Already a project is running.\nPlease close the existing and open the saved project");
            }
        }
        /*
         * Clicked on datasummary button or datasummary item 
         * opens Datasummary dialog box 
         * checks if summary is canceled or not
         * if exits,then display genarated output file in tabbed pane and adds the file name to the result node 
         */
        if (source == btnDataSummary || source == dataSummaryItem) {
            System.out.println("Analyzing");            System.out.println("Analyzing");
            try {
                if (Constant.genoList.size() > 0 && Constant.phenoList.size() > 0) {
                    long startTime = System.currentTimeMillis();
                    DataSummary dataSummary = new DataSummary(frame, true);
                    dataSummary.setLocationRelativeTo(frame);
                    dataSummary.setVisible(true);
                    long stopTime = System.currentTimeMillis();
                    long difference = stopTime - startTime;
                    long seconds = difference / 1000;
                    long minutes = seconds / 60;
                    long hours = minutes / 60;
                    String resSummary = "";
                    long days = TimeUnit.MILLISECONDS.toDays(difference);
                    long hours1 = TimeUnit.MILLISECONDS.toHours(difference) - TimeUnit.DAYS.toHours(TimeUnit.MILLISECONDS.toDays(difference));
                    long minutes1 = TimeUnit.MILLISECONDS.toMinutes(difference) - TimeUnit.HOURS.toMinutes(TimeUnit.MILLISECONDS.toHours(difference));
                    long seconds1 = TimeUnit.MILLISECONDS.toSeconds(difference) - TimeUnit.MINUTES.toSeconds(TimeUnit.MILLISECONDS.toMinutes(difference));
                    if (days == 0) {
                        resSummary = String.format(" (Total Time: %02d Hr : %02d Min : %02d Sec)", hours1, minutes1, seconds1);
                    } else {
                        resSummary = String.format(" (Total Time: %dd Days : %02d Hr : %02d Min : %02d Sec)", days, hours, minutes, seconds);
                    }
                    System.out.println(resSummary);
                    if (null != Constant.summaryResultList.get(0)) {
                        String genosummarCsv = Constant.resultdirectory + Constant.summaryResultList.get(0);
                        System.out.println("output path csv " + genosummarCsv);
                        JPanel browserDataSummary = gsMethods.csvReaderSNP(genosummarCsv);
                        tabbedPane.add(browserDataSummary);
                        tabbedPane.setSelectedIndex(tabbedPane.getTabCount() - 1);
                        dyntree.addObject(resultsNode, Constant.summaryResultList.get(0), true);
                        Constant.outputResultsList.add(Constant.summaryResultList.get(0).toString());
                        Constant.resultList.add(Constant.summaryResultList.get(0).toString());
                    }

                    if (null != Constant.summaryResultList.get(1)) {
                        String genosummaryHtml = Constant.resultdirectory + Constant.summaryResultList.get(1);
                        int noofLinesInHTML = GsMethods.noofLinesInHTML(genosummaryHtml);
                        GsMethods.addingProcessTime2HTMlFileSummary(noofLinesInHTML, genosummaryHtml, resSummary);
                        JWebBrowser browser = gsMethods.browser(genosummaryHtml);
                        tabbedPane.add(browser);
                        tabbedPane.setSelectedIndex(tabbedPane.getTabCount() - 1);
                        dyntree.addObject(resultsNode, Constant.summaryResultList.get(1), true);
                        Constant.outputResultsList.add(Constant.summaryResultList.get(1).toString());
                        Constant.resultList.add(Constant.summaryResultList.get(1).toString());
                    }

                    if (null != Constant.summaryResultList.get(2)) {
                        String phenoSummaryHtml = Constant.resultdirectory + Constant.summaryResultList.get(2);
                        int noofLinesInHTML = GsMethods.noofLinesInHTML(phenoSummaryHtml);
                        GsMethods.addingProcessTime2HTMlFileSummary(noofLinesInHTML, phenoSummaryHtml, resSummary);
                        JWebBrowser browe = gsMethods.browser(phenoSummaryHtml);
                        tabbedPane.add(browe);
                        tabbedPane.setSelectedIndex(tabbedPane.getTabCount() - 1);
                        dyntree.addObject(resultsNode, Constant.summaryResultList.get(2));
                        Constant.outputResultsList.add(Constant.summaryResultList.get(2).toString());
                        Constant.resultList.add(Constant.summaryResultList.get(2).toString());
                    }
                } else {
                    JOptionPane.showMessageDialog(frame, "Import atleast a Genotype and Phenotype files");
                }
            } catch (HeadlessException we) {
                we.getMessage();
            }
        }

        /*
         selecting the R or Fortran for script running
         */
        if (source == engine) {
            EngineSelection engSelection = new EngineSelection(frame, true);
            engSelection.setLocationRelativeTo(frame);
            engSelection.setVisible(true);
        }

        /*
         * Opens About dialog box 
         */
        if (source == about) {
            System.out.println("about called");
            try {
                About abut = new About(frame, true);
                abut.setVisible(true);
            } catch (Exception abt) {
                abt.getMessage();
            }
        }
        /*
         * Opens user manual in tabbedpane 
         */
        if (source == manual) {
            try {
                String presentWorkingDir = System.getProperty("user.dir");
                OSCheck.OSType os = OSCheck.getOperatingSystemType();
                if ("Window".equals(os.toString())) {
                    JPanel pdfviewer = gsMethods.pdfviewer(presentWorkingDir + "\\ISMU 2 0 User Manual.pdf");
                    tabbedPane.add(pdfviewer);
                    int count = tabbedPane.getTabCount();
                    tabbedPane.setSelectedIndex(count - 1);
                }
                if ("Linux".equals(os.toString())) {
                    JPanel pdfviewer = gsMethods.pdfviewer(presentWorkingDir + "/ISMU 2 0 User Manual.pdf");
                    tabbedPane.add(pdfviewer);
                    int count = tabbedPane.getTabCount();
                    tabbedPane.setSelectedIndex(count - 1);
                }
            } catch (Exception e1) {
                e1.getMessage();
            }
        }

        /*
         * Opens file chooser box 
         * the selected file is displayed in tabbedpane and name added in genotype node
         * copy the selected file to the result directory 
         */
        if (source == genotypeFile) {
            if (Constant.genoList.size() > 0 && Constant.phenoList.size() > 0) {
                lastFileChoosen = Constant.browsepath;
                JFileChooser genoChooser = new JFileChooser(lastFileChoosen);
                genoChooser.setDialogTitle("Select a genotype file");
                FileNameExtensionFilter filter = new FileNameExtensionFilter("CSV", "csv");
                genoChooser.setFileFilter(filter);
                int r = genoChooser.showOpenDialog(frame);
                if (r == JFileChooser.APPROVE_OPTION) {
                    try {
                        Constant.browsepath = genoChooser.getSelectedFile().getAbsolutePath();
                        String selectedFile = genoChooser.getSelectedFile().getPath();
                        if (Constant.genoList.contains(selectedFile)) {
                            JOptionPane.showMessageDialog(frame, "Genotype file selected already opened \n Please select other genotype file");
                        } else {
                            Constant.genoFileIsValid = true;
                            String genoName = new File(selectedFile).getName();
                            Date date = new Date();
                            SimpleDateFormat dateformat = new SimpleDateFormat("mmss");
                            String genotmp = GsMethods.stripFileExtension(genoName) + dateformat.format(date) + ".csv";
                            String genoResultDir = Constant.resultdirectory + genotmp;
                            CopyFile.copy(selectedFile, genoResultDir);
                            JPanel csvTable = gsMethods.csvReaderSNP(genoResultDir);
                            tabbedPane.add(csvTable);
                            dyntree.addObject(genotypeNode, genotmp, Boolean.TRUE);
                            int tabCount = tabbedPane.getTabCount();
                            tabbedPane.setSelectedIndex(tabCount - 1);
                            Constant.genoList.add(genotmp);
                        }
                    } catch (IOException ex) {
                        Logger.getLogger(GsMain.class.getName()).log(Level.SEVERE, null, ex);
                    }
                }
                if (r == JFileChooser.CANCEL_OPTION) {
                    //donothing 
                }
            } else {
                JOptionPane.showMessageDialog(frame, "For the first time use the open button or \nopen in menu to choose files and result directory");
            }

        }

        /*
         *opens file chooser dialog box 
         * display the selected file in tabbedpane and display the name under phenotype node
         * copy the file to result directory 
         */
        if (source == phenotypeFile) {
            if (Constant.genoList.size() > 0 && Constant.phenoList.size() > 0) {
                lastFileChoosen = Constant.browsepath;
                JFileChooser phenoChooser = new JFileChooser(lastFileChoosen);
                phenoChooser.setDialogTitle("Select a phenotype file");
                FileNameExtensionFilter filter = new FileNameExtensionFilter("CSV", "csv");
                phenoChooser.setFileFilter(filter);
                int r = phenoChooser.showOpenDialog(frame);
                if (r == JFileChooser.APPROVE_OPTION) {
                    Constant.browsepath = phenoChooser.getSelectedFile().getAbsolutePath();
                    try {
                        String selectedFile = phenoChooser.getSelectedFile().getPath();
                        if (Constant.phenoList.contains(selectedFile)) {
                            JOptionPane.showMessageDialog(null, "Phenotype file selected already opened \n Please select other genotype file");
                        } else {
                            Constant.phenoFileIsValid = true;
                            String genoName = new File(selectedFile).getName();
                            Date date = new Date();
                            SimpleDateFormat dateformat = new SimpleDateFormat("mmss");
                            String phenotmp = GsMethods.stripFileExtension(genoName) + dateformat.format(date) + ".csv";
                            String phenoResultDir = Constant.resultdirectory + phenotmp;
                            CopyFile.copy(selectedFile, phenoResultDir);
                            JPanel csvTable = gsMethods.csvReaderSNP(phenoResultDir);
                            tabbedPane.add(csvTable);
                            dyntree.addObject(phenotypeNode, phenotmp, Boolean.TRUE);
                            int tabCount = tabbedPane.getTabCount();
                            tabbedPane.setSelectedIndex(tabCount - 1);
                            Constant.phenoList.add(phenotmp);
                        }
                    } catch (IOException ex) {
                        Logger.getLogger(GsMain.class.getName()).log(Level.SEVERE, null, ex);
                    }
                }
                if (r == JFileChooser.CANCEL_OPTION) {
                    //donothing 
                }
            } else {
                JOptionPane.showMessageDialog(frame, "For the first time use the open button or \nopen in menu to choose files and result directory");
            }
        }
    }
//end region of action listener 

//open project 
    /**
     *
     * @param pathProject
     */
    public void commanndLineProjectOpening(String pathProject) {
        statusProSave = true; //setting the project as saved beacuse it is loading existing project

        String fileSelected = pathProject;
        int counter = 0;
        BufferedReader in = null;
        try {
            OSCheck.OSType os = OSCheck.getOperatingSystemType();
            if ("Window".equals(os.toString())) {
                in = new BufferedReader(new FileReader(fileSelected));
            }
            if ("Linux".equals(os.toString())) {
                in = new BufferedReader(new FileReader(fileSelected));
            }
            String[] str;
            try {
                while (in.ready()) {
                    String dataLine = in.readLine();
                    if (!dataLine.isEmpty()) {
                        if (counter == 0) {
                            Constant.resultdirectory = dataLine;
                            counter++;
                        } else if (counter >= 1) {
                            if (counter == 1) {
                                Constant.browsepath = dataLine;
                            }
                            str = dataLine.split(",");
                            if (counter == 2) {
                                for (String str1 : str) {
                                    Constant.genotype = str1;
                                    String genopath = Constant.resultdirectory + str1;
                                    JPanel csvreaderGeno = gsMethods.csvreader(genopath);
                                    tabbedPane.add(csvreaderGeno);
                                    dyntree.addObject(genotypeNode, str1, Boolean.TRUE);
                                    int count1 = tabbedPane.getTabCount();
                                    tabbedPane.setSelectedIndex(count1 - 1);
                                    Constant.genoList.add(str1);
                                }
                            }
                            if (counter == 3) {
                                for (String str1 : str) {
                                    Constant.phenotype = str1;
                                    String phenoPath = Constant.resultdirectory + str1;
                                    JPanel csvreaderPheno = gsMethods.csvreader(phenoPath);
                                    tabbedPane.add(csvreaderPheno);
                                    dyntree.addObject(phenotypeNode, str1, Boolean.TRUE);
                                    int count2 = tabbedPane.getTabCount();
                                    tabbedPane.setSelectedIndex(count2 - 1);
                                    Constant.phenoList.add(str1);
                                }
                            }
                            if (counter == 4) {
                                for (String str1 : str) {
                                    Constant.population = str1;
                                    String populationPath = Constant.resultdirectory + str1;
                                    JPanel csvreaderPopulation = gsMethods.csvreader(populationPath);
                                    tabbedPane.add(csvreaderPopulation);
                                    dyntree.addObject(populationNode, str1, Boolean.TRUE);
                                    int count5 = tabbedPane.getTabCount();
                                    tabbedPane.setSelectedIndex(count5 - 1);
                                    Constant.populationList.add(str1);
                                }
                            }
                            if (counter == 5) {
                                for (String str1 : str) {
                                    String path = Constant.resultdirectory + str1;
                                    String extension = FilenameUtils.getExtension(path);
                                    if ("htm".equals(extension) || "html".equals(extension)) {
                                        JWebBrowser browserOpen = gsMethods.browser(path);
                                        tabbedPane.add(browserOpen);
                                        dyntree.addObject(resultsNode, str1, Boolean.TRUE);
                                        int count6 = tabbedPane.getTabCount();
                                        tabbedPane.setSelectedIndex(count6 - 1);
                                        Constant.outputResultsList.add(str1);
                                        Constant.resultList.add(str1);
                                    }
                                    if ("csv".equals(extension)) {
                                        JPanel csvreader = gsMethods.csvreader(path);
                                        tabbedPane.add(csvreader);
                                        dyntree.addObject(resultsNode, str1, Boolean.TRUE);
                                        int count6 = tabbedPane.getTabCount();
                                        tabbedPane.setSelectedIndex(count6 - 1);
                                        Constant.outputResultsList.add(str1);
                                        Constant.resultList.add(str1);
                                    }
                                }
                            }
                            if (counter == 6) {
                                Constant.engineSelection = Integer.parseInt(dataLine);
                            }
                            if (counter == 7) {
                                Constant.projectName = dataLine;
                            }
                            if (counter == 8) {
                                genoc = Integer.parseInt(dataLine);
                            }
                            if (counter == 9) {
                                phenoc = Integer.parseInt(dataLine);
                            }
                            if (counter == 10) {
                                //  Constant.RandomForestCount = Integer.parseInt(dataLine);
                            }
                            if (counter == 11) {
                                resultc = Integer.parseInt(dataLine);
                            }
                            counter++;
                        }
                    } else {
                        counter++;
                    }
                }
            } catch (IOException ex) {
                Logger.getLogger(GsMain.class.getName()).log(Level.SEVERE, null, ex);
            }
            try {
                in.close();
            } catch (IOException ex) {
                Logger.getLogger(GsMain.class.getName()).log(Level.SEVERE, null, ex);
            }
        } catch (FileNotFoundException ex) {
            Logger.getLogger(GsMain.class.getName()).log(Level.SEVERE, null, ex);
        } finally {
            try {
                in.close();
            } catch (IOException ex) {
                Logger.getLogger(GsMain.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
    }
//star region
//main method of the class
//call the start() 

    public static void main(String args[]) {

        String path = "";//args[0];
        System.out.println("path saved " + path);
        File f1 = new File(path);

        if (f1.exists()) {
            System.out.println("File exist");
            Constant.cmdOpenProject = path;
        } else {
            System.out.println("File not found!");
            Constant.cmdOpenProject = "noargsismu";
        }
        Application.launch(GsMain.class, args);
    }
//end region main

    //start region 
    //class for handling the tree in application
    //adding and removing nodes in tree
    //opening the closed files in tabbedpane when double clicked on node
    final class DynamicTree implements TreeSelectionListener, MouseListener {

        private final DefaultMutableTreeNode rootNode;
        private final DefaultTreeModel treeModel;
        private final JTree tree;

        /*
         * constructor of dynmictree class
         * initilizing rootNode,treeModel,tree 
         * adding the tree to leftComponent of splitPane
         */
        public DynamicTree() {
            rootNode = new DefaultMutableTreeNode("Data");
            treeModel = new DefaultTreeModel(rootNode);
            tree = new JTree(treeModel);
            JScrollPane scrollPane = new JScrollPane(tree);
            splitPane.setLeftComponent(scrollPane);
            tree.addMouseListener(this);
        }

        //remove all the node of the tree
        public void clear() {
            rootNode.removeAllChildren();
            treeModel.reload();
        }

        //method to remove the current node selected
        public void removeCurrentNode() {
            TreePath currentSelection = tree.getSelectionPath();
            if (currentSelection != null) {
                DefaultMutableTreeNode currentNode = (DefaultMutableTreeNode) (currentSelection.getLastPathComponent());
                MutableTreeNode parent = (MutableTreeNode) (currentNode.getParent());
                if (parent != null) {
                    treeModel.removeNodeFromParent(currentNode);
                }
            } // Either there was no selection, or the root was selected.
        }

        //add the child node to the tree
        public DefaultMutableTreeNode addObject(DefaultMutableTreeNode parent, Object child) {
            return addObject(parent, child, false);
        }

        /*
         * method called when adding nodes to tree GsMain class
         * parameters passed are parent node,child name and visibility
         */
        public DefaultMutableTreeNode addObject(DefaultMutableTreeNode parent, Object child, boolean shouldBeVisible) {
            DefaultMutableTreeNode childNode = new DefaultMutableTreeNode(child);
            if (parent == null) {
                parent = rootNode;
            }
            // It is key to invoke this on the TreeModel, and NOT DefaultMutableTreeNode
            treeModel.insertNodeInto(childNode, parent, parent.getChildCount());

            tree.setCellRenderer(null);
            // Make sure the user can see the lovely new node.
            if (shouldBeVisible) {
                tree.scrollPathToVisible(new TreePath(childNode.getPath()));
            }
            return childNode;
        }

        /*
         * actionListeners of tree 
         */
        @Override
        //not used
        public void valueChanged(TreeSelectionEvent tse) {
        }
        //not used 

        @Override
        public void mouseClicked(MouseEvent e) {
        }

        /*
         * when double clicked on child node on tree 
         * it highlights the tab if already file is opened in tabbedpane
         * or 
         * if file name selected is closed, then it opens the file in tabbedpane 
         * when right clicked on child node shows a popup menu (remove)
         */
        @Override
        public void mousePressed(MouseEvent e) {
            int nodeSelection = tree.getRowForLocation(e.getX(), e.getY());
            TreePath selPath = tree.getPathForLocation(e.getX(), e.getY());
            if (nodeSelection != -1) {
                if (e.getClickCount() == 1) {
                    //donothing
                } else if (e.getClickCount() == 2) {
                    DefaultMutableTreeNode selectedNode = ((DefaultMutableTreeNode) selPath.getLastPathComponent());
                    String selectedNodeName = selectedNode.toString();
                    if (selectedNode.isLeaf()) {
                        if (!"Genotype".equals(selectedNodeName) || !"Phenotype".equals(selectedNodeName)
                                || !"Relationship matrix".equals(selectedNodeName) || !"Pedigree data".equals(selectedNodeName)
                                || !"Covariate".equals(selectedNodeName) || !"Result Directory".equals(selectedNodeName) || !"Log".equals(selectedNodeName)) {
                            try {
                                int count = tabbedPane.getTabCount(); //storing the tab count to the variable
                                String resultdir = Constant.resultdirectory; //getting the result directory path
                                boolean checkCondition = true; //boolean variable for checking the condiation 
                                for (int i = 0; i < count; i++) { //for loop
                                    if (selectedNodeName.equals(tabbedPane.getTabTitleAt(i))) { //checking the condition whether selectednode name and tabbedpane title is same
                                        tabbedPane.setSelectedIndex(i); //setting the tabbedpane i
                                        i = count + 1; //adding the 1 to the i value
                                        checkCondition = false; //setting the boolean variable to false for not entring into the else part
                                    }
                                }
                                if (checkCondition == true) {
                                    //checking local variable is true or false 
                                    //if the selected node file is not opened in closedtabbedpane then her it opens
                                    String resultStoringdir = resultdir + selectedNodeName; //adding the selectednode name to result directory
                                    String extension = selectedNodeName.substring(selectedNodeName.lastIndexOf(".") + 1, selectedNodeName.length());
                                    if ("html".equals(extension) || "htm".equals(extension)) {
                                        JWebBrowser browser = gsMethods.browser(resultStoringdir);
                                        tabbedPane.add(browser);
                                    }
                                    if ("txt".equals(extension)) {
                                        if (Constant.resultdirectory != null) {
                                            JScrollPane textFileReaderScroll = gsMethods.textFileReader(resultStoringdir);
                                            tabbedPane.add(textFileReaderScroll);
                                        }
                                    }
                                    if ("csv".equals(extension)) {
                                        JPanel csvreaderResult = gsMethods.csvreader(resultStoringdir);
                                        tabbedPane.add(csvreaderResult);
                                    }
                                    tabbedPane.setSelectedIndex(count);
                                    tabbedPane.tabAboutToClose(count);
                                }
                            } catch (Exception ie) {
                                System.out.println(ie.getMessage());
                            }
                        }
                    }
                }
            }
            if (((e.getModifiers() & InputEvent.BUTTON3_MASK) != 0) && (tree.getSelectionCount() > 0)) {
                DefaultMutableTreeNode selectedNode = ((DefaultMutableTreeNode) selPath.getLastPathComponent());
                if ("Genotype".equals(selectedNode.toString()) || "Phenotype".equals(selectedNode.toString())
                        || "Relationship matrix".equals(selectedNode.toString()) || "Pedigree data".equals(selectedNode.toString())
                        || "Covariate".equals(selectedNode.toString()) || "Result Directory".equals(selectedNode.toString())
                        || "Log".equals(selectedNode.toString()) || "log.txt".equals(selectedNode.toString()) || "Data".equals(selectedNode.toString())) {
                    //donothing
                } else {
                    showMenu(e.getX(), e.getY());
                }
            }
        }

        /*
         * shows popup menu remove
         * clicked on remove button calls the deleteSeletctedItems()
         */
        protected void showMenu(int x, int y) {
            JPopupMenu removePopup = new JPopupMenu();
            JMenuItem menuItem = new JMenuItem("Remove");
            TreePath path = tree.getSelectionPath();
            Object node1 = path.getLastPathComponent();
            if (node1 == tree.getModel().getRoot()) {
                menuItem.setEnabled(false);
            }
            removePopup.add(menuItem);
            menuItem.addActionListener(new ActionListener() {
                @Override
                public void actionPerformed(ActionEvent event) {
                    del();
                    deleteSelectedItems();
                }
            });
            removePopup.show(tree, x, y);
        }

        //delteing files 
        //WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
        List dele = new ArrayList();
        DefaultMutableTreeNode n;

        public void del() {
            dele.clear();
            try {
                DefaultTreeModel model = (DefaultTreeModel) (tree.getModel());
                TreePath[] paths = tree.getSelectionPaths();

                for (TreePath path : paths) {
                    n = (DefaultMutableTreeNode) (path.getLastPathComponent());
                    model.removeNodeFromParent(n);   //removing the node from the tree 
                    dele.add(n.toString()); //adding the list of node names for deleting the cressponding files in the folder
                    System.out.println("nodes list " + n.toString());
                }
            } catch (Exception e) {
                System.out.println(e);
            }
        }
        //*********************************************************
        //??????????????????????????????????????????????????????????????
        /*
         * removes the node ,delete the file from result directory and close the tab
         */
        //   List deleteFileList = new ArrayList(); //get the list of nodes selected for deleting 
        //  DefaultMutableTreeNode node;

        protected void deleteSelectedItems() {
            /*   deleteFileList.clear();
             try {
             DefaultTreeModel model = (DefaultTreeModel) (tree.getModel());
             TreePath[] paths = tree.getSelectionPaths();

             for (TreePath path : paths) {
             node = (DefaultMutableTreeNode) (path.getLastPathComponent());
             model.removeNodeFromParent(node);   //removing the node from the tree 
             deleteFileList.add(node.toString()); //adding the list of node names for deleting the cressponding files in the folder
             System.out.println("nodes list " + node.toString());
             }
             */ for (int count = 0; count <= dele.size() - 1; count++) { //itterating the array list 
                String listitt = dele.get(count).toString();
                //removing the file name lists 
                if (Constant.genoList.contains(listitt)) {
                    System.out.println("geno exists deleted");
                    Constant.genoList.remove(listitt);
                }
                if (Constant.phenoList.contains(listitt)) {
                    System.out.println("pheno exists deleted");
                    Constant.phenoList.remove(listitt);
                }
                if (Constant.populationList.contains(listitt)) {
                    System.out.println("population exists deleted");
                    Constant.populationList.remove(listitt);
                }
                if (Constant.outputResultsList.contains(listitt)) {
                    System.out.println("results exists deleted");
                    Constant.outputResultsList.remove(listitt);
                }

                //if html file delted cresspondig png files deletion method 
                if ("htm".equals(FilenameUtils.getExtension(listitt)) || "html".equals(FilenameUtils.getExtension(listitt))) {
                    GsMethods.removePngFromHtml(Constant.resultdirectory + listitt);
                }

                //deleting the file in the directory 
                File file = new File(Constant.resultdirectory + listitt);
                boolean delete = file.delete();
                
                System.out.println("delete status " + listitt + " " + delete);
                if (delete == false) {
                    boolean delete1 = file.delete();
                    System.out.println("delete try 1" + delete1);
                }
                
                //removing the tabs from the tabbedpane 
                for (int i = 0; i <= tabbedPane.getTabCount(); i++) {
                    if (listitt.equals(tabbedPane.getTabTitleAt(i))) { //checking the condition whether selectednode name and tabbedpane title is same
                        tabbedPane.remove(i);
                        i = tabbedPane.getTabCount();
                    } else {
                        if (i == tabbedPane.getTabCount() - 1) {
                            i = tabbedPane.getTabCount();
                        }
                    }
                }
            }
            // } catch (Exception e) {
            //    System.out.println(e);
            // }
        }
//not used listener

        @Override
        public void mouseReleased(MouseEvent e) {
            //    throw new UnsupportedOperationException("Not supported yet.");
        }

        //not used listener
        @Override
        public void mouseEntered(MouseEvent e) {
            //   throw new UnsupportedOperationException("Not supported yet.");
        }

        //not used listener
        @Override
        public void mouseExited(MouseEvent e) {
            //   throw new UnsupportedOperationException("Not supported yet.");
        }

        /**
         * Returns an ImageIcon, or null if the path was invalid.
         */
        ImageIcon createImageIcon(String path) {
            java.net.URL imgURL = DynamicTree.class.getResource(path);
            if (null != imgURL) {
                return new ImageIcon(imgURL);
            } else {
                System.err.println("Couldn't find file: " + path);
                return null;
            }
        }
    }
//end region DynamicTree class 
//start region
//to add the child nodes to the root node
//which is initilized when starting of appilcation

    void populateTree(DynamicTree treePanel) {
        String genoNode = "Genotype";
        String phenoNode = "Phenotype";
        String popNode = "Covariate";
        String resNode = "Result Directory";
        String logNode1 = "Log";

        genotypeNode = treePanel.addObject(null, genoNode);
        phenotypeNode = treePanel.addObject(null, phenoNode);
        populationNode = treePanel.addObject(null, popNode);
        resultsNode = treePanel.addObject(null, resNode);
        logNode = treePanel.addObject(null, logNode1);
    }
}
//end region GsMain class 
