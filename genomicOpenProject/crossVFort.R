# We have All files name with us and we need to do following:
# 1. INPUT : FOld and Reps
# 2. Loop REP
# 3. Read All files in R, Geno, Pheno and Covariate also if required
# 4. Sample Size = nrow(PhenoData)*(1-1/Folds)
# 5. Take Sample and subset all files
# 6. Recreate files and call alphabayes
# 

crossVFort <- function( Replications =2, Folds =4)
{
                     
sink("ErrorCVFort.txt")
# Relplications <-2
# Folds <-4

phenoTraining <- read.table("phenoTraining.txt", na.strings = "9", header=F)#, row.names=1)
genoTra <-read.table("genoTraining.txt",        na.strings = "9", header=T)#, row.names=1)
#genoTes <-read.table ("genoTesting.txt",         na.strings = "9", header=T)#, row.names=1)
CVVector <- NULL

for( i in 1:(Folds*Replications))
{
# Take Sample of Phenotypic Data
SubsetCV <- sample(nrow(phenoTraining), nrow(phenoTraining)*(1-1/Folds))   # Make a subset Index  needed to fit GS model
SubsetCV
# Make Subset Training
SubsetCV_phenoTraining  <- phenoTraining [SubsetCV,]   # Take subset from Pheno
SubsetCV_genoTra        <- genoTra[,c(1,SubsetCV+1)]        # Take subset from Training 



# Make Subset Testing
SubsetCV_phenoTesting   <- phenoTraining [-SubsetCV,]   # Take subset from Pheno
SubsetCV_genoTes        <- genoTra[,- (SubsetCV+1)]     # Take subset from Testing 

# Initialize other Parameters

nGenoTr    <-nrow(SubsetCV_phenoTraining)
nGenoTest  <-nrow(phenoTraining) - nGenoTr
nSnp  			= nrow(genoTra)

# Work on Covariate

if(! (CovFileName =="Select")) # i.e. Covariate Selected, as Defaul Value is "Select"
{
  # Read Covariate File
  
  CovSubSetTrain <-  read.table("CovTrain.txt",                 na.strings ="9") #, col.names=F)
  CovSubSetTest <- read.table("CovTest.txt",                   na.strings ="9") #, col.names=F)
  
  SubsetCV_CovSubSetTrain <- CovSubSetTrain[SubsetCV,]
  SubsetCV_CovSubSetTest  <- CovSubSetTrain[-SubsetCV,] 
    
  nCov<-ncol(CovSubSetTrain)-1
  
  # Writing Covariate file for geno training

  CovTrainTxt <-  paste("CovTrainCV", i, ".txt", sep="")
  CovTestTxt <-   paste("CovTestCV", i, ".txt", sep="")
  write.table(SubsetCV_CovSubSetTrain, CovTrainTxt,                 quote=FALSE, na="9", col.names=F, row.names=F)
  write.table(SubsetCV_CovSubSetTest, CovTestTxt,                   quote=FALSE, na="9", col.names=F, row.names=F)

}

# Now write data sets
write.table(SubsetCV_phenoTraining, "phenoTrainingCV.txt", row.names=FALSE, col.names=FALSE,quote=FALSE, na="9")
write.table(SubsetCV_genoTra, "genoTrainingCV.txt",        row.names=FALSE,                quote=FALSE, na="9")
write.table(SubsetCV_genoTes, "genoTestingCV.txt",         row.names=FALSE,                quote=FALSE, na="9")

# Now Call AlphaBayes
if( ! (FRR == -1))  MarkerSolverM <-"Ridge"
if( ! (FBA == -1))  MarkerSolverM <-"BayesA"

# oldwd <- setwd(AWD)
# THRes <- source(paste(oldwd, "/FortGS.R", sep=""))

FortGS (oldwd=oldwd,     FilePATH		= getwd(), 
        GenotypeTraining	= "genoTrainingCV.txt",
        GenotypeTesting		= "genoTestingCV.txt",
        MarkerNames			= "Yes",
        MarkersInRows      	= "Yes",
        PhenotypeTraining	= "phenoTrainingCV.txt",
        PhenotypeTesting			= "None",
        CovariateTraining = CovTrainTxt,
        CovariateTestin  =  CovTestTxt,             
        FixedSnpFile		= "None",
        nSnp				= nSNPs,
        nCov        = nCov ,
        nAnisTr				= nGenoTr,
        nAnisTe				= nGenoTest,
        nRound			 	= NIter,
        nBurn				= NBurnIn,
        VareA				= 0.5,
        VarE				= 0.5,
        NumberOfProcessors	= 3,
        ScalingOption		= 2,
        MissingGenoCode		= 9,
        MarkerSolver		= MarkerSolverM 
)


gsOut.rr<-read.table("Ebv.txt", header=F)
names(gsOut.rr)<-c("geno","predVals")

correl<-cor.test(gsOut.rr$predVals,SubsetCV_phenoTesting[,2] )

CVVector <- c(CVVector,   correl$estimate)

print(correl)
#cat (correl)
}
print(CVVector)
#cat(CVVector)
res<-t.test(CVVector)


sink()
CVP = data.frame(CV=mean(CVVector), p.value=res$p.value)
return(CVP)
}
