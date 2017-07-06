args<-commandArgs(TRUE)
args<-c("D:/ICRISAT/Work/Analysis/Work/ISMU/ISMU","Genotype.csv","Result.htm")

setwd(args[1])

FileName1 <- paste(rnorm(1), "_1.png", sep="")
FileName2 <- paste(rnorm(1), "_2.png", sep="")

needPackages <- c("genetics","R2HTML")
needPackages <- needPackages[!needPackages %in% installed.packages()]
if (length(needPackages) > 0) install.packages(needPackages, repos="http://lib.stat.cmu.edu/R/CRAN")
require(R2HTML)
require(genetics)



# Reading Genotypic Data from directory
Marker<-read.csv(args[2],stringsAsFactors=F)

MisRatio<-sapply (as.data.frame(Marker[,-1]), function(X) sum(is.na(X))*100 / length(X))
nbreaks<-seq(0,100,by=5)
a=hist(MisRatio,breaks=nbreaks,plot=F)
per_miss<-round((a$count/length(MisRatio)*100),2)


lab<-paste(a$count, "\n(", per_miss, "%)", sep="")
ppi<-300
png(FileName1, width=10*ppi, height=6*ppi, res=ppi)

par(bg="lightyellow",cex=0.6 )

hist(MisRatio,breaks=nbreaks,col="green",axes=F,labels=lab,
     main="Frequency of Missing %", 
     xlab="%Missing", ylab="Number of Markers", cex.lab=1.2, cex.main=1, xpd=T)

rug(MisRatio)
axis(2,cex.axis=1)
axis(1, at=nbreaks,cex.axis=1,las=2)
dev.off()



MRes<-sapply( as.data.frame(Marker[,-1]), function (X) summary(genotype(X, X)  ))
PIC <- as.numeric(MRes[5,])
names(PIC)<-names(MisRatio)
nbreaks<-seq(0,0.5,by=0.025)
b=hist(PIC,breaks=nbreaks,plot=F)
per_PIC<-round((b$count/length(PIC)*100),3)

lab<-paste(b$count,"\n(", per_PIC, "%)", sep="")
png(FileName2, width=10*ppi, height=6*ppi, res=ppi)
par(bg="lightyellow", cex=0.6)
hist(PIC,breaks=nbreaks,col="sky blue",bg="gray",axes=F,labels=lab,main="Polymorphics Information Content", xlab="PIC", ylab="Number of Markers", xpd=T)
rug(PIC)
axis(2,cex.axis=0.8)
axis(1, at=nbreaks,cex.axis=0.8,las=2)
dev.off()

HTML.title("Histogram", file=args[3],append=F) 
HTML("Summary of Missing Genotypes",file=args[3])
# HTML(a,file=args[3])


HTMLInsertGraph(FileName1,file=args[3],Align="left", WidthHTML="100%",append=T)
HTML("Summary of Polymorphic Information Content",file=args[3])
HTMLInsertGraph(FileName2,file=args[3],Align="left", WidthHTML="100%", append=T)  
