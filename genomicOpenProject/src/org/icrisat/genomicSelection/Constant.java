 /*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
/*
 * 
 *In the class all the constant variables are decalred where it can be used in several places
 * It helps in easy way of asinging values
 * 
 */
package org.icrisat.genomicSelection;

import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author Mohan
 */
//class which maintances the present selected value for the field 
public class Constant {

    static String genotype = null; //stores the present genotype file choosen
    static String phenotype = null;//stores the present phenotype file choosen
    static String population = "";//stores the present pooulation file choosen
    static String resultdirectory = "empty"; //stores the working directory where the results and input files in these directory
    static String phenotraits = ""; //stores the count of the phenotraits for choosen file
    static int traitCount; //getting the trait count of the phenotype file choosen 
    static String browsepath = ""; //to store the present selected file directory of geno,pheno,...
    static String tabbselected = "";// to save the particular file selected (i.e the selected)
    static String projectdir = ""; //where the project should be saved 
    static String dataSummaryFileName = "";//name of the html file genarated when data summary is runned based on the input file name
    static String dataSummaryFileNameR = "";
    static boolean summaryCancled = false;
    static String analysisFileName = ""; //file name of the html file output genarated when 
    static String populationPath = "";
    static String haploPath = ""; //dir path of haplotype
    static public java.util.List<String> genoList = new ArrayList<>();
    static public java.util.List<String> phenoList = new ArrayList<>();
    static public java.util.List<String> populationList = new ArrayList<>();
    static public java.util.List<String> resultList = new ArrayList<>();
    static public java.util.List<String> porSaveList = new ArrayList<>();//for saving the state of the application and writing into text file
    static public java.util.List<String> outputResultsList = new ArrayList<>();
    static boolean genoFileIsValid;
    static String genoFilePathOnly = "";
    static boolean genoFileSelectionCanceled;
    static boolean phenoFileIsValid;
    static String phenoFilePathOnly = "";
    static boolean phenoFileSelectionCanceled = false;
    //In analysis which methods are selected
    static boolean RidgeReg = false;
    static boolean BayesCpi = false;
    static boolean BayesB = false;
    static boolean BayesLasso = false;
    static boolean PedigreeBlup = false;
    static boolean RandomeForest = false;
    static boolean KinshipGauss = false;
    static boolean RidgeRegressionFortran = false;
    static boolean BayeAFortran = false;
    //selected genotype and phenotype file in  analysis
    static String geno_Analysis = "";
    static String pheno_Analysis = "";
    static String popu_Analysis = "";
    static String BayessAddParameters = null;
    static int RandomForestAddParameteres = 0;
    static String crossValidationParameteres = null;
    static int noOfCPUs = 1;
    static boolean analysisCancled = false;//for is analysis is cancled
    //adding counter for the all the methods
    static int RidgeRegressionCount = 0;
    static int BayesCpiCount = 0;
    static int BayesBCount = 0;
    static int BayesLassoCount = 0;
    static int PedigreeBlupCount = 0;
    static int RandomForestCount = 0;
    static int KinshipGaussCount = 0;
    static public java.util.List<String> traitList = new ArrayList<>(); //getting all the trait selected names to progress display dialog
    //adding all the result genrated to arraylist
    static public List<String> resultGenratedList = new ArrayList<>();
    //adding additional parameters
    static int perMissing;
    static Double picValue;
    static Double mafValue;
    static List selectedTraits = new ArrayList(); //only for if some traits are only selected
    //  static List sortingForSelectedTraits=new ArrayList(); //for sorting the selected traits in analysis
    static public String iconNode = ""; //for adding the icon for particular node 
    static public int engineSelection = 1;  //1=R and 2=Fortran    //Engine Selection 
    static public String phnotypeSummary = "";
    //saving project 
    static public String projectName = null;
    static public int genoCount = 0;
    static public int phenCount = 0;
    static public int populationCount = 0;
    static public int resultFileCount = 0;
    static public String bayesFortParameters = null;
    static public List summaryResultList = new ArrayList();
    
    static public String cmdOpenProject=null;
}
