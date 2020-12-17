##Setup 
# Load libraries
library(AICcmodavg); library(minpack.lm)


# Name of files list
files <- list("orr/orrinfoHybs_AVG_DnMis_intime_withErrors.csv",
              "orr-inf/orr-infinfoHybs_AVG_DnMis_intime_withErrors.csv",
              "orr-inf+rev/orr-inf+revinfoHybs_AVG_DnMis_intime_withErrors.csv",
              "ind/indinfoHybs_AVG_DnMis_intime_withErrors.csv",
              "ind+rev/ind+revinfoHybs_AVG_DnMis_intime_withErrors.csv",
              "ind+sel/ind+selinfoHybs_AVG_DnMis_intime_withErrors.csv",
              "ind+rev+sel/ind+rev+selinfoHybs_AVG_DnMis_intime_withErrors.csv"
)


# Define output for svg file
svg(filename="Plot_fitCurves_DD_DACurves.svg", 
    width=8, 
    height=20)
# Set a 7x2 grid for plots
par(mfrow=c(7,2))

# Setup the main result table
mainTableResult <- data.frame ("Model"= character(), 
                               "lin_dd_AICWt" = double(), 
                               "quad_dd_AICWt" = double(), 
                               #"exp_dd_AICWt" = double(), 
                               "log_dd_AICWt" = double(),
                               "lin_da_AICWt" = double(), 
                               "quad_da_AICWt" = double(), 
                               #"exp_da_AICWt" = double(), 
                               "log_da_AICWt" = double()
                               )

# Fit every set of data using three functions (linear, quadratic, and logistic)
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
    y2 = dat[[4]]
  }
    

  nlc <- nls.lm.control(maxiter = 1024)

  # Create list of models
  print (paste0("Models AICs for ", name, " - ",names(dat)[2],":"))
  


  
  
  modelsdd <- list(
    nlsLM(y ~ I((a * x) + b),  data = dat, start = list(a = -1, b = 0), control=nlc), #linear
    nlsLM(y ~ I((a * x ^ 2) + (b * x) + c), data = dat, start = list(a = 0, b = 0, c = 0), control=nlc),  #quadratic
    #nlsLM(y ~ I(a * exp(b * x) + c), data = dat, start = list(a = 0,b = 0, c = 0.1), control=nlc), #exponential
    nlsLM(y ~ a/(1 + (b * exp(-c * x))), data = dat, start = list(a = max(y),b = 0.1, c = 0), control=nlc) #logistic
    #nlsLM(y ~ a * exp(-((x - b)^2)/(2 * c^2)), data = dat, start = list(a = max(y),b = 0, c = max(y)), control=nlc, trace = TRUE) #gauss from past3
    #nlsLM(y ~ a*(b^((x-c)^2)), data = dat, start = list(a = max(y),b = 1, c = min(y)), control=nlc, trace =  TRUE) #gauss from tutorial
  )

  print (paste0("Models AICs for ", name," - ",names(dat)[4], ":"))
  
  modelsda <- list(
    nlsLM(y2 ~ I((a * x) + b),  data = dat, start = list(a = -1, b = 0), control=nlc), #linear
    nlsLM(y2 ~ I((a * x ^ 2) + (b * x) + c), data = dat, start = list(a = 0, b = 0, c = 0), control=nlc),  #quadratic
    #nlsLM(y ~ I(a * exp(b * x) + c), data = dat, start = list(a = 0,b = 0, c = 0.1), control=nlc), #exponential
    nlsLM(y2 ~ a/(1 + (b * exp(-c * x))), data = dat, start = list(a = max(y),b = 0.1, c = 0), control=nlc) #logistic
    #nlsLM(y2 ~ a * exp(-((x - b)^2)/(2 * c^2)), data = dat, start = list(a = 0,b = max(y), c = 0), control=nlc) #gauss from past3
  )
  
  
  #fits <- nls_multstart(y ~ a * exp(-((x - b)^2)/(2 * c^2)), data = dat, control=nlc, iter = 500)
  
  
  #linears <- lm(y ~ poly(x,3))

  op <- par(cex = 0.5)
  
  
  
  
  #Print AIC table
  Modnames <-c("linear", "quadratic", "logistic") #, "gaussian") #, "quad2", "quad1", "exp no C")
  aicTabledd = aictab(cand.set = modelsdd, modnames = Modnames, sort = FALSE, second.ord = FALSE) ## with AIC
  aiccTabledd = aictab(cand.set = modelsdd, modnames = Modnames, sort = FALSE, second.ord = TRUE) ## with AICc
  
  aicTableda = aictab(cand.set = modelsda, modnames = Modnames, sort = FALSE, second.ord = FALSE) ## with AIC
  aiccTableda = aictab(cand.set = modelsda, modnames = Modnames, sort = FALSE, second.ord = TRUE) ## with AICc
  
  #Show result in console
  print (aicTabledd)
  print (aicTableda)
  #cat ("polynomic",AIC(linears), "\n\n")
 
  
  mainTableResult <- rbind(mainTableResult, data.frame("Model"= name, 
                                                       "lin_dd_AICWt" = round(aicTabledd$AICWt[1], digits = 2), 
                                                       "quad_dd_AICWt" = round(aicTabledd$AICWt[2], digits = 2), 
                                                       #"exp_dd_AICWt" = round(aicTable$AICWt[3], digits = 2), 
                                                       "log_dd_AICWt" = round(aicTabledd$AICWt[3], digits = 2),
                                                       "lin_da_AICWt" = round(aicTableda$AICWt[1], digits = 2), 
                                                       "quad_da_AICWt" = round(aicTableda$AICWt[2], digits = 2), 
                                                       #"exp_da_AICWt" = round(aicTable$AICWt[3], digits = 2), 
                                                       "log_da_AICWt" = round(aicTableda$AICWt[3], digits = 2)
                                                       ))
  
  
  plot(x, y, main=paste0(name," - ",names(dat)[2]), xlab="Num. of generations", ylab="Avg. number of DMIs", pch=19)
  #segments(x,y-dat[[3]],x,y+dat[[3]])
  lines(x, fitted(modelsdd[[1]]), col = 2, lty=2)
  lines(x, fitted(modelsdd[[2]]), col = 3, lty=3)
  lines(x, fitted(modelsdd[[3]]), col = 4, lty=4)
  #lines(x, fitted(linears), col = 5, lty=5)

  aicslinear <- paste("Linear (AIC: ",round(aicTabledd$AIC[1], digits = 2)," - AICwt: ", round(aicTabledd$AICWt[1], digits = 2), ")")
  aicsquad <- paste("Quadratic (AIC: ",round(aicTabledd$AIC[2], digits = 2)," - AICwt: ", round(aicTabledd$AICWt[2], digits = 2), ")") 
  aicslog <- paste("Logistic (AIC: ",round(aicTabledd$AIC[3], digits = 2)," - AICwt: ", round(aicTabledd$AICWt[3], digits = 2), ")")
  
  legend("bottomright", legend=c(aicslinear,aicsquad,aicslog), lty = c(2:4), col = c(2:4),  bty = 'n')
  
  
  plot(x, y2, main=paste0(name," - ",names(dat)[4]), xlab="Num. of generations", ylab="Avg. number of DMIs", pch=19)
  #segments(x,y-dat[[3]],x,y+dat[[3]])
  lines(x, fitted(modelsda[[1]]), col = 2, lty=2)
  lines(x, fitted(modelsda[[2]]), col = 3, lty=3)
  lines(x, fitted(modelsda[[3]]), col = 4, lty=4)
  #lines(x, fitted(linears), col = 5, lty=5)

  aicslinear <- paste("Linear (AIC: ",round(aicTableda$AIC[1], digits = 2)," - AICwt: ", round(aicTableda$AICWt[1], digits = 2), ")")
  aicsquad <- paste("Quadratic (AIC: ",round(aicTableda$AIC[2], digits = 2)," - AICwt: ", round(aicTableda$AICWt[2], digits = 2), ")") 
  aicslog <- paste("Logistic (AIC: ",round(aicTableda$AIC[3], digits = 2)," - AICwt: ", round(aicTableda$AICWt[3], digits = 2), ")")
  
  legend("bottomright", legend=c(aicslinear,aicsquad,aicslog), lty = c(2:4), col = c(2:4),  bty = 'n')
  
  
  cat("Parameters for " , paste0(name," - ",names(dat)[2]) , " \n" ,
      "Linear\n",
      "a=" , coef(modelsdd[[1]])[1], "\n" ,
      "b=" , coef(modelsdd[[1]])[2], "\n" ,
      "Quadratic\n",
      "a=" , coef(modelsdd[[2]])[1], "\n" ,
      "b=" , coef(modelsdd[[2]])[2], "\n" ,
      "c=" , coef(modelsdd[[2]])[3], "\n" ,
      # "Exponential\n",
      # "a=" , coef(models[[3]])[1], "\n" ,
      # "b=" , coef(models[[3]])[2], "\n" ,
      "Logistic\n",         
      "a=" , coef(modelsdd[[3]])[1], "\n" ,
      "b=" , coef(modelsdd[[3]])[2], "\n" ,
      "c=" , coef(modelsdd[[3]])[3], "\n\n" ,
      sep="")
 
   
  
  cat("Parameters for " , paste0(name," - ",names(dat)[3]) , " \n" ,
      "Linear\n",
      "a=" , coef(modelsda[[1]])[1], "\n" ,
      "b=" , coef(modelsda[[1]])[2], "\n" ,
      "Quadratic\n",
      "a=" , coef(modelsda[[2]])[1], "\n" ,
      "b=" , coef(modelsda[[2]])[2], "\n" ,
      "c=" , coef(modelsda[[2]])[3], "\n" ,
      # "Exponential\n",
      # "a=" , coef(models[[3]])[1], "\n" ,
      # "b=" , coef(models[[3]])[2], "\n" ,
      "Logistic\n",         
      "a=" , coef(modelsda[[3]])[1], "\n" ,
      "b=" , coef(modelsda[[3]])[2], "\n" ,
      "c=" , coef(modelsda[[3]])[3], "\n" ,
      sep="")
  
  cat("\n=======================================\n\n")
  
}



# Stop svg 
dev.off()

# Print final table
print (mainTableResult)
#print(xtable(mainTableResult, type = "latex"))





