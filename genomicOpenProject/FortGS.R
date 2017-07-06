# GenotypeTraining	= "GenoTrainingTr.txt"
# GenotypeTesting		= "GenoTestingTr.txt"
# MarkerNames			= "Yes"
# MarkersInRows      	= "Yes"
# PhenotypeTraining	= "PhenoTraining.csv"
# TbvTesting			= "None"
# FixedSnpFile		= "None"
# nSnp				= 2920
# nAnisTr				= 315
# nAnisTe				= 10
# nRound			 	= 1000
# nBurn				= 100
# VareA				= 0.5
# VarE				= 0.5
# NumberOfProcessors	= 3
# ScalingOption		= 2
# MissingGenoCode		= 9
# MarkerSolver		= "BayesA"

FortGS <- function(	oldwd=oldwd,   	FilePATH		= FilePATH, 
					GenotypeTraining	= "GenoTrainingTr.txt",
					GenotypeTesting		= "GenoTestingTr.txt",
					MarkerNames			= "Yes",
					MarkersInRows      	= "Yes",
					PhenotypeTraining	= "PhenoTraining.csv",
					PhenotypeTesting  = "None",
					CovariateTraining = "CovTrain.txt",
					CovariateTesting  = "CovTest.txt",

          
          
          FixedSnpFile		= "None",
					nSnp				= 2920,
					nCov        = 3 ,
					
					nAnisTr				= 315,
					nAnisTe				= 10,
					nRound			 	= 1000,
					nBurn				= 100,
					VareA				= 0.5,
					VarE				= 0.5,
					NumberOfProcessors	= 3,
					ScalingOption		= 2,
					MissingGenoCode		= 9,
					MarkerSolver		= "Ridge",
					Seed = rnorm(1)
                  )
{
cat("\nMethod used:\t",MarkerSolver, "\n")  
FilePATH <- FilePATH  # Same as AWD, however add a \ for Fortarn Program to work

IPFileGS <- "AlphaBayesSpec.txt"

# Make input file
FRes <- file(IPFileGS, "wt")

#cat("##############################Input parameters############################", fill=T, file=FRes)
cat("InputFolder\t,\"", 	FilePATH, "/\""	, fill=T, file=FRes, sep="")
cat("OutputFolder\t,\"",   	FilePATH, "/\""	, fill=T, file=FRes, sep="")
cat("GenotypeTraining\t,", 	GenotypeTraining 	, fill=T, file=FRes, sep="")
cat("GenotypeTesting\t,", 	GenotypeTesting 	, fill=T, file=FRes, sep="")
cat("MarkerNames\t,", 		MarkerNames			, fill=T, file=FRes, sep="")
cat("MarkersInRows\t,", 	MarkersInRows 		, fill=T, file=FRes, sep="")
cat("PhenotypeTraining\t,", PhenotypeTraining 	, fill=T, file=FRes, sep="")
cat("PhenotypeTesting       \t,", 		PhenotypeTesting 			, fill=T, file=FRes, sep="")     # This is TbvTesting earlier
cat("CovariateTraining\t,",   	CovariateTraining 			, fill=T, file=FRes, sep="")   # new
cat("CovariateTestin\t,",   	CovariateTesting 			, fill=T, file=FRes, sep="")       # new
cat("FixedSnpFile\t,", 		FixedSnpFile 		, fill=T, file=FRes, sep="")
cat("nSnp\t,", 				nSnp 				, fill=T, file=FRes, sep="")
cat("nCov\t,",   			nCov 				, fill=T, file=FRes, sep="")                         # new
cat("nAnisTr\t,", 			nAnisTr 			, fill=T, file=FRes, sep="")
cat("nAnisTe\t,", 			nAnisTe 			, fill=T, file=FRes, sep="")
cat("nRound\t,", 			nRound 				, fill=T, file=FRes, sep="")
cat("nBurn\t,", 			nBurn 				, fill=T, file=FRes, sep="")
cat("VareA\t,", 			VareA 				, fill=T, file=FRes, sep="")
cat("VarE\t,", 				VarE 				, fill=T, file=FRes, sep="")
cat("NumberOfProcessors\t,",NumberOfProcessors 	, fill=T, file=FRes, sep="")
cat("ScalingOption\t,", 	ScalingOption 		, fill=T, file=FRes, sep="")
cat("MissingGenoCode\t,", 	MissingGenoCode 	, fill=T, file=FRes, sep="")
cat("MarkerSolver\t,", 		MarkerSolver 		, fill=T, file=FRes, sep="")

close(FRes)

INFName <- paste( FilePATH,"/", IPFileGS, sep="")
##########################
#  This line it till jenaz remove bug from his code. This will ensure Specification file is at 
#  Both Places. I will aso make Seed.txt at both places
#
file.copy(from=INFName, to=paste( oldwd,"/", IPFileGS, sep=""), overwrite=T)

FRes <- file( paste(oldwd, "/Seed.txt", sep=""),  "wt")
cat(Seed, fill=T, file=FRes)
close(FRes)

FRes <- file( paste(FilePATH, "/Seed.txt", sep=""),  "wt")
cat(Seed, fill=T, file=FRes)
close(FRes)

cat("\nSeed = ", Seed, "\n\n")
####################################################

AGHExePath <- paste(oldwd, "/", "AlphaBayes", sep="")

file.remove("Ebv.txt")

a<-system2(AGHExePath, args=c("--par", shQuote (INFName) ), minimized = FALSE, invisible = T, stdout=F)
#QUIT after a maximum number of tries say 10?
  


}

#  oldwd <- "C:\\Users\\Arathore\\Desktop\\Abayes"
#  setwd(oldwd)
#  FortGS(oldwd=oldwd, FilePATH=oldwd)
# setwd("C:\\Users\\Arathore\\Desktop\\Abayes")
# FortGS(oldwd="D:/ICRISAT/Work/Analysis/Work/ISMU/ISMU", FilePATH="C:\\Users\\Arathore\\Desktop\\Abayes")

# 
# 
# 
# 
# 
# while(is.nan(LoopI))
# {
#   if(Times>5)
#   {
#     cat("\nError: AlphaBayes could not solve model!!!\n Quiting GS.\n")
#     #     EBV_FP <- file("Ebv.txt", "wt")
#     #     cat("0 0", file=EBV_FP)
#     
#     break
#   }
#   
#   a<-system2(AGHExePath, args=c("--par", shQuote (INFName) ), minimized = FALSE, invisible = T, stdout="")
#   THRes<-read.table("Ebv.txt", header=F)
#   Times <- Times+1
#   cat("\nSystem2 Return Value:\t", a, "\n")
#   cat("\nTime: ", Times, "  AlphaBayes in search of  correct seed.\n")
#   LoopI <- THRes [1,2]
#   
#   #QUIT after a maximum number of tries say 10?
#   
#   
# }
