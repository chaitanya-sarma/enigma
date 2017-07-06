args<-commandArgs(TRUE)
#source("bayesBparallel.R")
#source("bayesCpi.R")
#source("bayesBparallel.R")
#source("bayesLasso.R")
#source("BLRnonVerbose.R")
#source("pedigreeBLUP.R")
TFile<-tempfile()
#TFile
sink(TFile)


StartTime<-Sys.time()
#path = "c:\Program Files\r\R-3.0.1\bin\x64"
#cd D:\ICRISAT\Work\Analysis\Work\ISMU\ISMU  

GetRFileName <- function (start = "RES_")
{
  paste(start, substring(date(), 5, 7), substring(date(), 9, 10), "_", substring(date(), 12, 13), substring(date(), 15, 16), substring(date(), 18, 19), abs(round(rnorm(1),3)), sep = "")  
}


require(R2HTML)

#args <- c("C:\\Users\\Arathore\\Desktop\\Abayes", "Genotype.csv", "Phenotype.csv", "Result.htm" ,"0.3", "-1", "-1", -1, 2, "-", 1, 1, -1, -1, -1, 1, 2,2,":",  "GSResTMP.txt", "TMP1.txt", "TMP2.txt","dummy", "-1",  "1")
#rscript bFort.R C:\ABC Genotype.csv Phenotype.csv Result.htm 0.3 -1 -1 -1 2 - 1 1 -1 -1 -1 1 1 2 1 GSResTMP.txt TMP1.txt TMP2.txt 2 1 -1 a.txt 100 110 22

AWD       <- args[1]
AGDATA    <- args[2]
APDATA    <- args[3]
ARFILE    <- args[4]
AMISP      <- as.numeric(args[5])
APIC  	  <- as.numeric(args[6])
AMAF		  <- as.numeric(args[7])
AGSMETHOD	<- args[8]
ANCORE		<- as.numeric(args[9])
AMISSCH		<- args[10]
ARR	<- args[11]
AKG	<- args[12]
ABB	<- args[13]
ABCPI	<- args[14]
ABL	<- args[15]
ARF	<- args[16]
ATraitNo <-as.numeric( args[17])
AMTYPE  <- as.numeric( args[18])
ASEP    <- args[19]

ATEMP_FILE_CORR    <- args[20]   # Will Store: Method, Coorel and  P value
ATEMP_FILE_OTHERS1  <- args[21]  # 
ATEMP_FILE_OTHERS2  <- args[22]  #
Engine  <-  as.numeric(args[23])   #
FRR	    <-  as.numeric(args[24])   # Fortran Method : RR
FBA     <-  as.numeric(args[25])   # Fortran Method : BayesA
CovFileName    <-  args[26]   # Fortran Method : BayesA

NIter <- args[27]  # This is Rounds
NBurnIn<-args[28]  # This is Burn iN
NThin <- args[29]  # Thining is not in Fortran, But we kept for later use

Replications = as.numeric(args[30])
Folds= as.numeric(args[31])


if(is.na(NIter)) NIter <-1000
if(is.na(NBurnIn)) NBurnIn <-100
if(is.na(CovFileName)) CovFileName <- "a.txt"

CovTrainTxt <-"None"
CovTestTxt <- "None"
nCov <-0

# if (Replications > 1) CrossValidation =1
CrossValidation <- 0
if (Folds > 1) CrossValidation <- 1   # Cross Validation is only possible if we have Folds more than 1


oldwd<-setwd(AWD)
THRes<-source(paste(oldwd, "/FortGS.R", sep=""))

load("gsdata.rdata")


########## Prepare Filef for AlphaBayes ################################
# Extracting phenotypic data for genotype in training population 
# phenoTra<-FinalBLUP[complete.cases(FinalBLUP),]
i<-ATraitNo
FinalPheno<-FinalBLUP
pheno<-FinalPheno[,c(1,1+i)]
phenoTraining<-pheno[!is.na(pheno[,2]),]
# phenoTraining<-setNames(phenoTra, rep(" ", length(phenoTra)))

# Extracting genotype names in testing population 
phenoTes<-FinalBLUP[!(FinalBLUP[,1]%in%phenoTraining[,1]),][1]

# Extracting genotypic data in training population  
genoTra<-geno[, c("Marker",phenoTraining[,1])]
#genoTraining<-data.frame(colname=rownames(genoTra),genoTra)

# Extracting genotypic data in testing population  
genoTes<-geno
#genoTesting<-data.frame(colname=rownames(genoTes),genoTes)

if(! (CovFileName =="Select")) # i.e. Covariate Selected, as Defaul Value is "Select"
{
  # Read Covariate File
  Covariates <- read.csv(CovFileName, header=T, row.names=1)
  
  # Covariate file for training population
  CovSubSetTrain <- Covariates[rownames(Covariates) %in% colnames(genoTra)[-1],] # -1 is to remove "marker" from col names
  CovSubSetTrain <-CovSubSetTrain [order(rownames(CovSubSetTrain )),]
  
  # Covariate file for testing population
  CovSubSetTest <- Covariates[rownames(Covariates) %in% colnames(genoTes)[-1],] # -1 is to remove "marker" from col names
  CovSubSetTest<-CovSubSetTest[order(rownames(CovSubSetTest)),]
  nCov<-ncol(CovSubSetTest)
  
  # Writing Covariate file for geno training
  write.table(CovSubSetTrain, "CovTrain.txt",                 quote=FALSE, na="9", col.names=F)
  write.table(CovSubSetTest, "CovTest.txt",                   quote=FALSE, na="9", col.names=F)
  
  CovTrainTxt <-"CovTrain.txt"
  CovTestTxt <- "CovTest.txt"
  
}


write.table(phenoTraining, "phenoTraining.txt", row.names=FALSE,col.names=FALSE,quote=FALSE, na="9")
write.table(genoTra, "genoTraining.txt",        row.names=FALSE,                quote=FALSE, na="9")
write.table(genoTes, "genoTesting.txt",         row.names=FALSE,                quote=FALSE, na="9")

### Get Number of SNps and Genotypes
nSNPs <- nrow(genoTra)
nGenoTr <- ncol(genoTra)-1  # Same will be number of genotype for Geno Testing
# nSNPs
# nGenoTr
# genoTra
################################################################################
# Data Conversion Done #
################################################################################

# As we will do it trait wise no Loop is required
i<-ATraitNo

#for (i in 1:NTraits )
#{
pheno<-FinalPheno[,c(1,1+i)]
ppi<-300  # Resolution
IR <- .5  # Image Ratio
SG <- 8   #Significant Decimal Digits


cat("\nValue of FRR", FRR, "\n")
# ridge regression

if( ! (FRR == -1))
{  
  SelectedMethod<- "RidgeRegression"
  # Get Temporary and Unique Graph File Name
  PNGName <- paste(getwd(), "/", GetRFileName(), ".png", sep="");
  THRes<-png(PNGName, width=IR*10*ppi, height=IR*10*ppi, res=ppi)
  
  Seed = rnorm(1)
  # Start RR BLUP 
  FortGS (oldwd=oldwd,   	FilePATH		= getwd(), 
          GenotypeTraining	= "genoTraining.txt",
          GenotypeTesting		= "genoTesting.txt",
          MarkerNames			= "Yes",
          MarkersInRows      	= "Yes",
          PhenotypeTraining	= "phenoTraining.txt",
          PhenotypeTesting			= "None",
          CovariateTraining = CovTrainTxt,
          CovariateTestin  =  CovTestTxt,             
          FixedSnpFile		= "None",
          nSnp				= nSNPs,
          nCov        = nCov ,
          nAnisTr				= nGenoTr,
          nAnisTe				= nGenoTr,
          nRound			 	= NIter,
          nBurn				= NBurnIn,
          VareA				= 0.5,
          VarE				= 0.5,
          NumberOfProcessors	= 3,
          ScalingOption		= 2,
          MissingGenoCode		= 9,
          MarkerSolver		= "Ridge"
  )
  
  
  #########################################################
  #   Now Process Results and put in a format as of R
  #########################################################
  
  gsOut.rr<-read.table("Ebv.txt", header=F)
  names(gsOut.rr)<-c("geno","predVals")
  
  correl<-cor.test(gsOut.rr$predVals,pheno[,2])
  THRes<-plot(gsOut.rr$predVals, pheno[,2], col="red", type="p", pch=19, main="Scatter Diagram", xlab="GEBV", ylab="Observed Phenotype")
  
  GEBV<-data.frame(Genotype=colnames(geno)[-1], Observed= pheno[,2], GEBV=round(gsOut.rr$predVals,SG), row.names=NULL)
  hasPheno <- !is.na(pheno[,2])
  
  #   if(sum(hasPheno)<length(pheno[,2]))
  #   {
  #     TGEBV <- GEBV[hasPheno,]
  #     BGEBV <-GEBV[!hasPheno,]
  #   } else {
  #     TGEBV <- GEBV[hasPheno,]
  #     BGEBV<-NULL
  #   }
  
  # Close Device PNG  
  dev.off() # Graph Ready
  
  # Prepare and write Result to File  
  
  ResultSTR <- paste("\nRidgeRegression\t",  round(correl$estimate, SG), "\t",round(correl$p.value, SG), sep="", collapse="\n" )     
  GSRes <- file(ATEMP_FILE_CORR, "at")   #Open in APPEND Mode to write results of RR
  cat(ResultSTR, file=GSRes)
  
  
  if(CrossValidation==1)
  {
  THRes<-source(paste(oldwd, "/crossVFort.R", sep=""))
  CVRes <- crossVFort( Replications =Replications, Folds =Folds)
  ResultSTR <- paste("\nCrossValidation(R=",Replications, ",F=", Folds, ")\t",  round(CVRes[1], SG), "\t",round(CVRes [2], SG), sep="", collapse="\n" )     
  cat(ResultSTR, file=GSRes)
  }
  
    
  close(GSRes)

    
}



# BayesB
# NOTE: because BayesB uses multicore, it has to be done from the command line, not the GUI
#   unless you specify n.core=1
if( ! (FBA == -1))
{	
  SelectedMethod<-"BayesA"
  
  # Get Temporary and Unique Graph File Name
  PNGName <- paste(getwd(), "/", GetRFileName(), ".png", sep="");
  THRes<-png(PNGName, width=IR*10*ppi, height=IR*10*ppi, res=ppi)
  
  Seed = rnorm(1)
  # Start Bayes A 
  FortGS (oldwd=oldwd,     FilePATH		= getwd(), 
          GenotypeTraining	= "genoTraining.txt",
          GenotypeTesting		= "genoTesting.txt",
          MarkerNames			  = "Yes",
          MarkersInRows      	= "Yes",
          PhenotypeTraining	= "phenoTraining.txt",
          PhenotypeTesting			= "None",
          CovariateTraining = CovTrainTxt,
          CovariateTestin  =  CovTestTxt,           
          FixedSnpFile		= "None",
          nSnp				= nSNPs,
          nCov        = nCov ,
          nAnisTr				= nGenoTr,
          nAnisTe				= nGenoTr,
          nRound  		 	= NIter,
          nBurn				= NBurnIn,
          VareA				= 0.5,
          VarE				= 0.5,
          NumberOfProcessors	= 3,
          ScalingOption		= 2,
          MissingGenoCode		= 9,
          MarkerSolver  	= "BayesA"
  )
  # Start BA 
  
  #########################################################
  #   Now Process Results and put in a format as of R
  #########################################################
  gsOut.ba<-read.table("Ebv.txt", header=F)
  
  names(gsOut.ba)<-c("geno","predVals")
  
  # If in worst case, AlphaBayes do not converge and all values are NaN, handle here
  LoopI <- gsOut.ba [1,2]
  if(is.nan(LoopI)) 
  { # i.e Model not solved and Ebv.txt has NaN
    correl <- list(estimate=NaN, p.value=1)    
  }else
  {
    correl<- cor.test(gsOut.ba$predVals,pheno[,2])
    plot(gsOut.ba$predVals, pheno[,2], col="maroon", type="p", pch=19, main="Scatter Diagram", xlab="GEBV", ylab="Observed Phenotype")
  }
  
  GEBV<-data.frame(Genotype=colnames(geno)[-1],Observed= pheno[,2],GEBV=round(gsOut.ba$predVals,SG),row.names=NULL)
  
  hasPheno <- !is.na(pheno[,2])
  
  #   if(sum(hasPheno)<length(pheno[,2]))
  #   {
  #     TGEBV <- GEBV[hasPheno,]
  #     BGEBV <-GEBV[!hasPheno,]
  #   } else {
  #     TGEBV <- GEBV[hasPheno,]
  #     BGEBV<-NULL
  #   }
  # Close Device PNG  
  dev.off() # Graph Ready
  
  # Prepare and write Result to File  
  
  ResultSTR <- paste("\nBayesA\t",  round(correl$estimate, SG), "\t",round(correl$p.value, SG), sep="", collapse="\n" )  
  GSRes <- file(ATEMP_FILE_CORR, "at")   #Open in APPEND Mode to write results of RR
  cat(ResultSTR, file=GSRes)

  if(CrossValidation==1)
  {
    THRes<-source(paste(oldwd, "/crossVFort.R", sep=""))
    CVRes <- crossVFort( Replications =Replications, Folds =Folds)
    ResultSTR <- paste("\nCrossValidation(R=",Replications, ",F=", Folds, ")\t",  round(CVRes[1], SG), "\t",round(CVRes [2], SG), sep="", collapse="\n" )     
    cat(ResultSTR, file=GSRes)
  }
  
  
  
  close(GSRes)
  
} 

# Now Write Results to HTML Page ;-)
# Calculate Time

EndTime <- Sys.time()
DiffTime <- EndTime-StartTime 
DiffTimeInString <- paste("\nTime take for ", SelectedMethod, " : ", round(c(DiffTime), 4), "Seconds")

# 
# # 1 : Start Phenotype Summary 
# PSummary <- summary(FinalPheno[,i])
#   TRes<-HTML.title("Phenotypic Summary", HR=3, file=ARFILE, Append=T)
#   TRes<-HTML(PSummary, HR=1, file=ARFILE,  Append=T)
# Make PNG of Historgam
# 
# PNGNameTrait <- paste(getwd(), "/", GetRFileName(), ".png", sep="");
# THRes<-png(PNGNameTrait, width=IR*10*ppi, height=IR*10*ppi, res=ppi)
# TRes<-hist(FinalPheno[,i], xlab=NameTrait[i])
# dev.off()
#   TRes<- HTMLInsertGraph(PNGNameTrait, file=ARFILE)

# Now put GS Results
# First Heading with Selected method
TRes<-HTML.title(paste("Genomic Selection Method: ", SelectedMethod, "(Fortran)"), HR=3, file=ARFILE, Append=T)
TRes<- HTMLInsertGraph(PNGName, file=ARFILE)

# Write Graph Title

TRes<-HTML("<div align=\"center\">", file=ARFILE, Append=T)
TRes<-HTML(paste("<b align=\"center\">Observed Phenotype Vs GEBV (", SelectedMethod,")</b>"), HR=4, file=ARFILE, Append=T, Align="center")
TRes<-HTML("<br/>", file=ARFILE, Append=T)

# Write Correlations

TRes<-data.frame(Method= SelectedMethod, Correlation = round(correl$estimate, SG), "p.value" = round(correl$p.value, SG), row.names=NULL)
Res<-HTML(TRes, file=ARFILE, Append=T)

 if(CrossValidation==1)
  {
  TRes<-data.frame("Cross Validation"= paste("Repl=",Replications, ", Folds=", Folds, sep="" ), "Correlation" = round(CVRes[1,1], SG), "p.value" = round(CVRes [2], SG), row.names=NULL)
  Res<-HTML(TRes, file=ARFILE, Append=T)
  }
  

# Write GEBVs
TRes<-HTML("<br/>", file=ARFILE, Append=T)
TRes<-HTML(paste("<b align=\"center\">GEBV for ", colnames(pheno[2]),"</b>"), HR=4, file=ARFILE, Append=T, Align="center")
Res<-HTML(GEBV, file=ARFILE, Append=T,  digits=6, nsmall = 6)

# Close Center Alignement
TRes<-HTML("</div>", file=ARFILE, Append=T)


# This and below may go to c.r 

# Write All Correlation and P Value

# Res<-read.table(ATEMP_FILE_CORR, header=T)
# TRes<-HTML(Res, file=ARFILE, Append=T)

# Write Time and Logo
# TRes<-HTML("<hr/><p align=center> ", file=ARFILE, Append=T)
# TRes<-HTML(paste("Results Generated by ISMU 2.0 : ", format(Sys.time(), "%a %b %d %X %Y")), file=ARFILE, Append=T)
# TRes<-HTMLEndFile(file , file=ARFILE)

# TRes<-HTML("<table border=1><tr class= firstline><td>Method</td><td>Correlation</td><td>Prob>t</td></tr>", file=ARFILE, Append=T)
# TRes<-HTML(paste("<tr><td>", SelectedMethod, "</td><td>", round(correl$estimate, SG), "</td><td>", round(correl$p.value, SG),"</td></tr>"), file=ARFILE, Append=T)

sink()

setwd(oldwd)
file.remove(TFile)
