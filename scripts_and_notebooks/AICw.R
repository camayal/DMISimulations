##Setup 
# Load libraries
library(AICcmodavg); library(minpack.lm)


# Name of files list
files <- list("orr/orr_timeComparison.csv",
              "orr-inf/orr-inf_timeComparison.csv",
              "orr-inf+rev/orr-inf+rev_timeComparison.csv",
              "dem/dem_timeComparison.csv",
              "dem+rev/dem+rev_timeComparison.csv",
              "dem+sel/dem+sel_timeComparison.csv",
              "dem+rev+sel/dem+rev+sel_timeComparison.csv"
)


# Define output for svg file
svg(filename="Plot_fitCurves_GeneralCurve.svg", 
    width=8, 
    height=10)
par(mfrow=c(4,2))


mainTableResult <- data.frame ("Model"= character(), 
                               "lin_AICWt" = double(), 
                               "quad_AICWt" = double(), 
                               "log_AICWt" = double()
                               #"bert_AICWt" = double()
                               )


for (i in seq_along(files)) {
  # Get name
  prename = strsplit(files[[i]], "/")
  prename = strsplit(prename[[1]][length(table(prename))], "_")
  if (any(prename[[1]]=="pureDataFrame.csv")){
    rawText = "_raw" 
    rawData = TRUE
    }
  else {
    rawText = ""
    rawData = FALSE
  }
  name = paste( unlist(prename[[1]][1]), rawText, collapse="")
  name = gsub(" ", "", name, fixed = TRUE)
  #print (name)

  # Load file
  dat <- read.csv(file=files[[i]], header=TRUE, sep=",")

  if (rawData){
    x = dat$generations
    y = dat$avgDMIsperHyb
  }else{
    x = dat[[1]]
    y = dat[[2]]
  }
    

  nlc <- nls.lm.control(maxiter = 1024)

  # Create list of models
  print (paste0("Models AICs for ", name, ":"))
  
  models <- list(
    nlsLM(y ~ I((a * x) + b),  data = dat, start = list(a = -1, b = 0), control=nlc), #linear
    nlsLM(y ~ I((a * x ^ 2) + (b * x) + c), data = dat, start = list(a = 0, b = 0, c = 0), control=nlc),  #quadratic
    #nlsLM(y ~ I(a * exp(b * x) + c), data = dat, start = list(a = 0,b = 0, c = 0.1), control=nlc), #exponential
    nlsLM(y ~ a/(1 + (b * exp(-c * x))), data = dat, start = list(a = max(y),b = 0.5, c = 0), control=nlc) #logistic
    #nlsLM(y ~ a * (1 -  (b * exp(-c * x))), data = dat, start = list(a = max(y),b = 0, c = 0), control=nlc) #Bertalanffy 
  )
  
  #linears <- lm(y ~ poly(x,3))

  op <- par(cex = 0.5)
  
  plot(x, y, main=name, xlab="Num. of generations", ylab="Avg. number of DMIs", pch=19)
  segments(x,y-dat[[3]],x,y+dat[[3]])
  lines(x, fitted(models[[1]]), col = 2, lty=2)
  lines(x, fitted(models[[2]]), col = 3, lty=3)
  lines(x, fitted(models[[3]]), col = 4, lty=4)
  #lines(x, fitted(models[[4]]), col = 5, lty=5)
  
  #lines(x, fitted(linears), col = 5, lty=5)
  
  #Print AIC table
  Modnames <-c("linear", "quadratic", "logistic") #, "bertalanffy") 
  aicTable = aictab(cand.set = models, modnames = Modnames, sort = FALSE, second.ord = FALSE) ## with AIC
  aiccTable = aictab(cand.set = models, modnames = Modnames, sort = FALSE, second.ord = TRUE) ## with AICc
  
  #Show result in console
  print (aicTable)
  #cat ("polynomic",AIC(linears), "\n\n")
 
  
  mainTableResult <- rbind(mainTableResult, data.frame("Model"= name, 
                                                       "lin_AICWt" = round(aicTable$AICWt[1], digits = 2), 
                                                       "quad_AICWt" = round(aicTable$AICWt[2], digits = 2), 
                                                       "log_AICWt" = round(aicTable$AICWt[3], digits = 2)
                                                       # "bert_AICWt" = round(aicTable$AICWt[4], digits = 2)
                                                       ))
  
  aicslinear <- paste0("Linear (AIC: ",round(aicTable$AIC[1], digits = 2)," - AICwt: ", round(aicTable$AICWt[1], digits = 2), ")")
  aicsquad <- paste0("Quadratic (AIC: ",round(aicTable$AIC[2], digits = 2)," - AICwt: ", round(aicTable$AICWt[2], digits = 2), ")") 
  aicslog <- paste0("Logistic (AIC: ",round(aicTable$AIC[3], digits = 2)," - AICwt: ", round(aicTable$AICWt[3], digits = 2), ")")
  # aicsbert <- paste0("Bertalanffy (AIC: ",round(aicTable$AIC[4], digits = 2)," - AICwt: ", round(aicTable$AICWt[4], digits = 2), ")")
  
  legend("bottomright", legend=c(aicslinear,aicsquad,aicslog), lty = c(2:5), col = c(2:5),  bty = 'n')

  
  
  cat("Parameters for " , name , " \n" ,
      "Linear\n",
      "a=" , coef(models[[1]])[1], "\n" ,
      "b=" , coef(models[[1]])[2], "\n" ,
      "Quadratic\n",
      "a=" , coef(models[[2]])[1], "\n" ,
      "b=" , coef(models[[2]])[2], "\n" ,
      "c=" , coef(models[[2]])[3], "\n" ,
      "Logistic\n",
      "a=" , coef(models[[3]])[1], "\n" ,
      "b=" , coef(models[[3]])[2], "\n" ,
      "c=" , coef(models[[3]])[3], "\n" ,
      # "Bertalanffy\n",         
      # "a=" , coef(models[[4]])[1], "\n" ,
      # "b=" , coef(models[[4]])[2], "\n" ,
      # "c=" , coef(models[[4]])[3], "\n" ,
      sep="")
  
  cat("\n=======================================\n\n")
  
}



# Stop svg 
dev.off()




print (mainTableResult)






