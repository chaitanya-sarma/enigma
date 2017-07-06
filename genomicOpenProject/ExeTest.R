a<-system2("AlphaBayes.exe" , args=c("--par", shQuote ("AlphaBayesSpec.txt") ), minimized = T, invisible =T, stdout="")
a
a<-system(paste("AlphaBayes.exe " , " --par AlphaBayesSpec.txt"))
a