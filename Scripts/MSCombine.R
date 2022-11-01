library(MScombine)
library(tidyverse)

#Load in the data - They have been filtered down to only feautres found thru the GLMM and RMANOVA to save computational time
negtrim <- read_csv('~/Documents/Projects/poopsoup/MScombine/neg_trimmed.csv')
postrim <- read_csv('~/Documents/Projects/poopsoup/MScombine/pos_trimmed.csv')
adducts <- read_csv('~/Documents/Projects/poopsoup/MScombine/data/adducts.csv')

#Format properly to work with MScombine
neg <- negtrim %>%
  #Use either the given neutral mass or calculate by multiplying the m/z by the charge (z)
  mutate(mass = ifelse(is.na(neutral_mass), mz*charge, neutral_mass)) %>%
  dplyr::select(compound, mass, rt, starts_with('neg')) %>%
  #Fix the names to work with MScombine
  dplyr::rename('Compound Name' = compound, 'Mass' = mass, 'RT' = rt)

pos <- postrim %>%
  #Use either the given neutral mass or calculate by multiplying the m/z by the charge (z)
  mutate(mass = ifelse(is.na(neutral_mass), mz*charge, neutral_mass)) %>%
  dplyr::select(compound, mass, rt, starts_with('pos')) %>%
  dplyr::rename('Compound Name' = compound, 'Mass' = mass, 'RT' = rt)

#Find common features between the two ionization modes:
commons <- FindCommon(POSITIVE = pos, NEGATIVE = neg, ADDUCTS = adducts, Masstolerance = 0.002, RTtolerance = 0.2)

#Clean up our data to remove mismatched compounds:
commonsbetter <- RemoveMismatch(CommonEntities = commons)

#The StudyRTdiff() function from MSCombine would not work so I pulled out the source code as a solution
r <- commonsbetter
colnames(r)[colnames(r) == 'RT+'] <- 'x'
colnames(r)[colnames(r) == 'RT-'] <- 'y'

plot(r$x, r$y)

linearfit <- lm(`RT-` ~ `RT+`, data = commonsbetter)
plot(linearfit)

res <- residuals(linearfit)
colnames(commonsbetter)[colnames(commonsbetter) == 'RT+'] <- 'x'
colnames(commonsbetter)[colnames(commonsbetter) == 'RT-'] <- 'y'
commonsbetter$Residuals <- res
#MScombine::FilterbyRT(commonsbetter, MaxResidual = 0.2, MinResidual = -0.2)

#Once again not working so pulling source code:
FilterbyRT<-function(CommonEntitiesImproved,MaxResidual,MinResidual) {
  CommonEntitiesPreFiltered<-CommonEntitiesImproved[(CommonEntitiesImproved$Residuals)<MaxResidual,]
  CommonEntitiesFiltered<-CommonEntitiesPreFiltered[(CommonEntitiesPreFiltered$Residuals)>MinResidual,]
  pdf("RTpositive_vs_RTnegative_Filtered.pdf")
  plot(CommonEntitiesFiltered$x,CommonEntitiesFiltered$y,main="RT+ vs RT-")
  dev.off()
  LinearFitFiltered <- lm(y ~ x, data = CommonEntitiesFiltered)
  par(mfrow=c(1,4))
  pdf("LinearRegressionCharacteristics_Filtered.pdf")
  graph3<-plot(LinearFitFiltered)
  dev.off()
  pdf("Histogram_of_residuals_Filtered.pdf")
  hist(residuals(LinearFitFiltered),breaks="FD")
  dev.off()
  colnames(CommonEntitiesFiltered)[colnames(CommonEntitiesFiltered)=="x"] <- "RT+"
  colnames(CommonEntitiesFiltered)[colnames(CommonEntitiesFiltered)=="y"] <- "RT-"
  write.table(CommonEntitiesFiltered,file="CommonEntitiesFiltered.csv",sep=",",row.names=FALSE,na="")
  CommonEntitiesFiltered
}

rtfilt <- FilterbyRT(commonsbetter, MaxResidual = 0.2, MinResidual = -0.2)

#For the last step, I want to combine the two polarities - When MSCombine does this polarity is lost
#I want to be able to map compound back to polarity for annotation so I will pull and use their source code again
#From Source code
e <- 2
f <- 2
k <- dim(rtfilt)
nentities <- k[1]
d <- 1

IDtoDel=c("CpdIDtoDelete");
IDtoDelNeg = rbind(IDtoDel,c(0));
IDtoDelPos = rbind(IDtoDel,c(0));

while(d<(nentities+1)){
  if(rtfilt[d,9]>rtfilt[d,10]){
    IDtoDelNeg[e] = rtfilt[d,2];
    e=e+1;
    d=d+1;
  } 
  else {
    IDtoDelPos[f] = rtfilt[d,1];
    f=f+1;
    d=d+1;
  }
}

IDtoDelNegDef<-unique(IDtoDelNeg)
IDtoDelPosDef<-unique(IDtoDelPos)

POSITIVE <- as.data.frame(pos)
NEGATIVE <- as.data.frame(neg) 

POSITIVEDef<-POSITIVE[-which(POSITIVE[,1] %in% IDtoDelPosDef),]
NEGATIVEDef<-NEGATIVE[-which(NEGATIVE[,1] %in% IDtoDelNegDef),]
###End source code section

#Save out the final 
#write.csv(POSITIVEDef, file = 'mscPos.csv')
#write.csv(NEGATIVEDef, file = 'mscNeg.csv')






