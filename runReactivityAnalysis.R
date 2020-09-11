require(readxl)
require(ggpubr)
require(stringr)
require(gplots)
require(brglm)
require(ggrepel)

### To be edited

setwd("D:/Wechseldatentraeger/08072018/Reactivity/ForRevision/")
myPath <- setwd(getwd())
myFileName_UniformMass <- "uniformMassData.xlsx"
mySheetName_UniformMass <- 1
myFileName_UniformSurface <- "uniformSurfaceData.xlsx"
mySheetName_UniformSurface <- 1
myFileName_UnitNames <- "unitNames.xlsx"
mySheetName_UnitNames <- 1

drawBoxplots <- function(myMeanColname, myAssay, myUnitNames){
  
  myColnameForPlot <- str_replace_all(myMeanColname, c("\\.x" = "_mass","\\.y" = "_surface", "AUC" = "AUC_mass", "FRAS_Mean" = "FRAS_Mean_surface"))
  myUnitName <- myUnitNames$Unit[which(myUnitNames$Assay == myColnameForPlot)]
  
  if (grepl(".x", myMeanColname, fixed = TRUE)){
    myMeanColnameInOtherMetric <- str_replace_all(myMeanColname, c("\\.x" = "\\.y"))
  } else{
    myMeanColnameInOtherMetric <- str_replace_all(myMeanColname, c("\\.y" = "\\.x"))
  } 
  
  
  ### pairwise complete
  
  myPosInds_pairwisecomplete <- which(eval(parse(text=paste("myData_merged$", myAssay, ".x", sep = ""))) == "active" & !is.na(eval(parse(text=paste("myData_merged$", myMeanColname, sep = "")))) & !is.na(eval(parse(text=paste("myData_merged$", myMeanColnameInOtherMetric, sep = "")))))
  myNegInds_pairwisecomplete <- which(eval(parse(text=paste("myData_merged$", myAssay, ".x", sep = ""))) == "passive" & !is.na(eval(parse(text=paste("myData_merged$", myMeanColname, sep = "")))) & !is.na(eval(parse(text=paste("myData_merged$", myMeanColnameInOtherMetric, sep = "")))))
  
  myDataForTest <- data.frame(myValues=c(eval(parse(text=paste("myData_merged$", myMeanColname, sep = "")))[myPosInds_pairwisecomplete], eval(parse(text=paste("myData_merged$", myMeanColname, sep = "")))[myNegInds_pairwisecomplete]),myClass=c(rep("active",length(eval(parse(text=paste("myData_merged$", myMeanColname, sep = "")))[myPosInds_pairwisecomplete])), rep("passive",length(eval(parse(text=paste("myData_merged$", myMeanColname, sep = "")))[myNegInds_pairwisecomplete]))))
  mySigTest <- compare_means(myValues ~ myClass, myDataForTest, method = "wilcox.test", paired = FALSE)
  
  
  if (length(unlist(strsplit(x = myUnitName, "(?<=blank |control |l ])", perl=T)))==1){
    p <- ggboxplot(myDataForTest, x = "myClass", y = "myValues",add = "jitter", xlab=paste(myAssay, "categorization", sep=" "), ylab=myUnitName, yscale="log10",color = "myClass", palette =c("black", "purple"),legend = "none",bxp.errorbar=T,font.label = list(size = 24))
  } else{
    if (!grepl("ug/ml",myUnitName)){
      p <- ggboxplot(myDataForTest, x = "myClass", y = "myValues",add = "jitter", xlab=paste(myAssay, "categorization", sep=" "), ylab=paste(unlist(strsplit(x = myUnitName, "(?<=blank |control |l ])", perl=T))[1],unlist(strsplit(x = myUnitName, "(?<=blank |control |l ])", perl=T))[2], sep="\n") ,yscale="log10",color = "myClass", palette =c("black", "purple"),legend = "none",bxp.errorbar=T,font.label = list(size = 24))
    } else{
      if (!grepl("deposited", myUnitName)){
        p <- ggboxplot(myDataForTest, x = "myClass", y = "myValues",add = "jitter", xlab=paste(myAssay, "categorization", sep=" "), ylab=eval(bquote(expression(atop(.(unlist(strsplit(x = myUnitName, "(?<=blank |control |l ])", perl=T))[1]),paste(.(str_remove(unlist(strsplit(x = myUnitName, "(?<=blank |control |l ])", perl=T))[2],"ug/ml")),mu,g/ml))))) ,yscale="log10",color = "myClass", palette =c("black", "purple"),legend = "none",bxp.errorbar=T,font.label = list(size = 24))
      } else{
        p <- ggboxplot(myDataForTest, x = "myClass", y = "myValues",add = "jitter", xlab=paste(myAssay, "categorization", sep=" "), ylab=eval(bquote(expression(atop(.(unlist(strsplit(x = myUnitName, "(?<=blank |control |l ])", perl=T))[1]),paste(.(str_remove(unlist(strsplit(x = myUnitName, "(?<=blank |control |l ])", perl=T))[2],"ug/ml - deposited")),mu,g/ml, " - deposited"))))) ,yscale="log10",color = "myClass", palette =c("black", "purple"),legend = "none",bxp.errorbar=T,font.label = list(size = 24))
      }
    }
  }
  #  Add p-value
  p2 <- p + stat_compare_means(method = "wilcox.test",label.y=log10(max(myDataForTest$myValues))+log10(10),label.x = 1.2, size=6)+
    stat_compare_means(label = "..p.signif..", method = "wilcox.test",label.y=log10(max(myDataForTest$myValues))+log10(5),label.x = 1.5, size=6)  
  p3 <- ggpar(p2, font.y=list(size=18),font.x=list(size=18),font.tickslab=list(size=18))
  ggsave(filename=paste0("Boxplot_",myAssay,myColnameForPlot,"_pairwiseComplete.png"), plot=p3, dpi="retina")
  
  
  
  ### complete
  
  myPosInds_complete <- which(eval(parse(text=paste("myData_merged$", myAssay, ".x", sep = ""))) == "active" & !is.na(myData_merged$CPH_Mean.x) & !is.na(myData_merged$DMPO_Mean.x) & !is.na(myData_merged$FRAS_AUC) & !is.na(myData_merged$Carbonyls_Mean.x) & !is.na(myData_merged$Carbonyls_deposited_Mean.x) & !is.na(myData_merged$CPH_Mean.y) & !is.na(myData_merged$DMPO_Mean.y) & !is.na(myData_merged$FRAS_Mean) & !is.na(myData_merged$Carbonyls_Mean.y) & !is.na(myData_merged$Carbonyls_deposited_Mean.y))
  myNegInds_complete <- which(eval(parse(text=paste("myData_merged$", myAssay, ".x", sep = ""))) == "passive" & !is.na(myData_merged$CPH_Mean.x) & !is.na(myData_merged$DMPO_Mean.x) & !is.na(myData_merged$FRAS_AUC) & !is.na(myData_merged$Carbonyls_Mean.x) & !is.na(myData_merged$Carbonyls_deposited_Mean.x) & !is.na(myData_merged$CPH_Mean.y) & !is.na(myData_merged$DMPO_Mean.y) & !is.na(myData_merged$FRAS_Mean) & !is.na(myData_merged$Carbonyls_Mean.y) & !is.na(myData_merged$Carbonyls_deposited_Mean.y))
  
  myDataForTest <- data.frame(myValues=c(eval(parse(text=paste("myData_merged$", myMeanColname, sep = "")))[myPosInds_complete], eval(parse(text=paste("myData_merged$", myMeanColname, sep = "")))[myNegInds_complete]),myClass=c(rep("active",length(eval(parse(text=paste("myData_merged$", myMeanColname, sep = "")))[myPosInds_complete])), rep("passive",length(eval(parse(text=paste("myData_merged$", myMeanColname, sep = "")))[myNegInds_complete]))))
  mySigTest <- compare_means(myValues ~ myClass, myDataForTest, method = "wilcox.test", paired = FALSE)
  
  
  if (length(unlist(strsplit(x = myUnitName, "(?<=blank |control |l ])", perl=T)))==1){
    p <- ggboxplot(myDataForTest, x = "myClass", y = "myValues",add = "jitter", xlab=paste(myAssay, "categorization", sep=" "), ylab=myUnitName, yscale="log10",color = "myClass", palette =c("black", "purple"),legend = "none",bxp.errorbar=T,font.label = list(size = 24))
  } else{
    if (!grepl("ug/ml",myUnitName)){
      p <- ggboxplot(myDataForTest, x = "myClass", y = "myValues",add = "jitter", xlab=paste(myAssay, "categorization", sep=" "), ylab=paste(unlist(strsplit(x = myUnitName, "(?<=blank |control |l ])", perl=T))[1],unlist(strsplit(x = myUnitName, "(?<=blank |control |l ])", perl=T))[2], sep="\n") ,yscale="log10",color = "myClass", palette =c("black", "purple"),legend = "none",bxp.errorbar=T,font.label = list(size = 24))
    } else{
      if (!grepl("deposited", myUnitName)){
        p <- ggboxplot(myDataForTest, x = "myClass", y = "myValues",add = "jitter", xlab=paste(myAssay, "categorization", sep=" "), ylab=eval(bquote(expression(atop(.(unlist(strsplit(x = myUnitName, "(?<=blank |control |l ])", perl=T))[1]),paste(.(str_remove(unlist(strsplit(x = myUnitName, "(?<=blank |control |l ])", perl=T))[2],"ug/ml")),mu,g/ml))))) ,yscale="log10",color = "myClass", palette =c("black", "purple"),legend = "none",bxp.errorbar=T,font.label = list(size = 24))
      } else{
        p <- ggboxplot(myDataForTest, x = "myClass", y = "myValues",add = "jitter", xlab=paste(myAssay, "categorization", sep=" "), ylab=eval(bquote(expression(atop(.(unlist(strsplit(x = myUnitName, "(?<=blank |control |l ])", perl=T))[1]),paste(.(str_remove(unlist(strsplit(x = myUnitName, "(?<=blank |control |l ])", perl=T))[2],"ug/ml - deposited")),mu,g/ml, " - deposited"))))) ,yscale="log10",color = "myClass", palette =c("black", "purple"),legend = "none",bxp.errorbar=T,font.label = list(size = 24))
      }
    }
  }
  #  Add p-value
  p2 <- p + stat_compare_means(method = "wilcox.test",label.y=log10(max(myDataForTest$myValues))+log10(10),label.x = 1.2, size=6)+
    stat_compare_means(label = "..p.signif..", method = "wilcox.test",label.y=log10(max(myDataForTest$myValues))+log10(5),label.x = 1.5, size=6)  
  p3 <- ggpar(p2, font.y=list(size=18),font.x=list(size=18),font.tickslab=list(size=18))
  ggsave(filename=paste0("Boxplot_",myAssay,myColnameForPlot,"_complete.png"), plot=p3, dpi="retina")
  
}



drawHeatmap <- function(myData_merged){
  
  myNumericScaledData <- data.frame(CPH=scale(myData_merged$CPH_Mean.y),DMPO=scale(myData_merged$DMPO_Mean.y),FRAS=scale(myData_merged$FRAS_Mean),Carbonyls=scale(myData_merged$Carbonyls_deposited_Mean.y))
  rownames(myNumericScaledData) <- myData_merged$NM
  
  myNAcounts <- apply(myNumericScaledData,1,function(x) length(which(is.na(x))))
  myNumericScaledData_All_complete <- myNumericScaledData[which(myNAcounts == 0),]
  colnames(myNumericScaledData_All_complete) <- c("ESR - CPH","ESR - DMPO","FRAS","Carbonyls")
  
  png("Heatmap_surface_complete.png", res=500, width=3000, height=2000)
  par(cex.main=0.6)
  heatmap.2(as.matrix(myNumericScaledData_All_complete), col=rev(heat.colors(20)), breaks = c(seq(-5,5,by=0.5)), Colv = F, dendrogram="row", na.rm=T, na.color="lightgrey", margins=c(7,11),trace="none",cexRow=0.7, cexCol=0.7, keysize=1.5, key.par=list(mgp=c(1, 0.4, 0), mar=c(2.5, 2, 1, 1), cex.lab=0.8, cex.axis=0.5))
  dev.off()
  
}


computeCorrelations <- function(myColnames){
  
  myCor <- cor.test(eval(parse(text=paste("myData_merged$", myColnames[1], sep = ""))),eval(parse(text=paste("myData_merged$", myColnames[2], sep = ""))), alternative = "two.sided", method = "spearman", exact = T, conf.level = 0.95)
  write(str_replace(myColnames[1],"\\.y",""), file = "Correlations.txt", append = T)
  write(str_replace(myColnames[2],"\\.y",""), file = "Correlations.txt", append = T)
  write(paste("rho:",myCor$estimate,sep=" "), file = "Correlations.txt", append = T)
  write(paste("p:", myCor$p.value, sep=" "), file = "Correlations.txt", append = T)
  write("\n", file = "Correlations.txt", append = T)
  
  # add complete correlations
 
  myData_uniformSurface_allComplete <- myData_merged[which(!is.na(myData_merged$CPH_Mean.y) & !is.na(myData_merged$DMPO_Mean.y) & !is.na(myData_merged$FRAS_Mean) & !is.na(myData_merged$Carbonyls_Mean.y) & !is.na(myData_merged$Carbonyls_deposited_Mean.y)),]
  
  myCor <- cor.test(eval(parse(text=paste("myData_uniformSurface_allComplete$", myColnames[1], sep = ""))),eval(parse(text=paste("myData_uniformSurface_allComplete$", myColnames[2], sep = ""))), alternative = "two.sided", method = "spearman", exact = T, conf.level = 0.95)
  write(str_replace(myColnames[1],"\\.y",""), file = "Correlations_completeCases.txt", append = T)
  write(str_replace(myColnames[2],"\\.y",""), file = "Correlations_completeCases.txt", append = T)
  write(paste("rho:",myCor$estimate,sep=" "), file = "Correlations_completeCases.txt", append = T)
  write(paste("p:",myCor$p.value,sep=" "), file = "Correlations_completeCases.txt", append = T)
  write("\n", file = "Correlations_completeCases.txt", append = T)
   
}


computeRegression <- function(myColnames, myAssay, i){


  myPosInds <- which(eval(parse(text=paste("myData_merged$", myAssay, ".x", sep = ""))) == "active" & !is.na(myData_merged$CPH_Mean.y) & !is.na(myData_merged$DMPO_Mean.y) & !is.na(myData_merged$FRAS_Mean) & !is.na(myData_merged$Carbonyls_deposited_Mean.y))
  myNegInds <- which(eval(parse(text=paste("myData_merged$", myAssay, ".x", sep = ""))) == "passive" & !is.na(myData_merged$CPH_Mean.y) & !is.na(myData_merged$DMPO_Mean.y) & !is.na(myData_merged$FRAS_Mean) & !is.na(myData_merged$Carbonyls_deposited_Mean.y))
  
  myPosData <- data.frame(CPH_Mean=myData_merged$CPH_Mean.y[myPosInds], DMPO_Mean=myData_merged$DMPO_Mean.y[myPosInds], FRAS_Mean=myData_merged$FRAS_Mean[myPosInds], Carbonyls_deposited_Mean=myData_merged$Carbonyls_deposited_Mean.y[myPosInds] ,Class = rep("active",length(myPosInds)))
  rownames(myPosData) <- myData_merged$NM[myPosInds]
  myNegData <- data.frame(CPH_Mean=myData_merged$CPH_Mean.y[myNegInds], DMPO_Mean=myData_merged$DMPO_Mean.y[myNegInds], FRAS_Mean=myData_merged$FRAS_Mean[myNegInds], Carbonyls_deposited_Mean=myData_merged$Carbonyls_deposited_Mean.y[myNegInds] ,Class = rep("passive",length(myNegInds)))
  rownames(myNegData) <- myData_merged$NM[myNegInds]
  
  myReactivity_All <- rbind(myPosData, myNegData)
  
  
  myClass_numeric <- as.numeric(myReactivity_All[,5])
  myClass_numeric[which(as.numeric(myReactivity_All[,5]) == 2)] <- 0
  myReactivity_All_numeric <- data.frame(myReactivity_All,Class_numeric=myClass_numeric)
  myReactivity_All_numeric_withoutNARows <- myReactivity_All_numeric[apply(myReactivity_All_numeric, 1,function(x) !any(is.na(x))),]
  
  
  myReactivity_All_numeric_withoutNARows <- myReactivity_All_numeric_withoutNARows[-which(rownames(myReactivity_All_numeric_withoutNARows)=="Fe2O3 larger"),] ### remove Fe2O3 larger
  
  
  doLOOCV <- function(myRow){
    output <- brglm(formula = as.formula(paste("Class_numeric", paste(str_replace(myColnames,"\\.y",""), collapse =" + "), sep=" ~ ")), data=myReactivity_All_numeric_withoutNARows[-myRow,], family=binomial)
    yweight <- predict(output, myReactivity_All_numeric_withoutNARows[myRow,][c(1:4)],type="response")
    return(c(output$aic,yweight))
  }
  
  myPredictions <- sapply(1:dim(myReactivity_All_numeric_withoutNARows)[1],doLOOCV)
  
  myMeanAIC <- mean(myPredictions[1,])
  mySDAIC <- sd(myPredictions[1,])
  mySensitivity <- sum(myPredictions[2,] >= 0.5 & myReactivity_All_numeric_withoutNARows$Class_numeric == 1) / sum(myReactivity_All_numeric_withoutNARows$Class_numeric == 1)
  mySpecificity <- sum(myPredictions[2,] <= 0.5 & myReactivity_All_numeric_withoutNARows$Class_numeric == 0) / sum(myReactivity_All_numeric_withoutNARows$Class_numeric == 0)
  myBalancedAccuracy <- (mySensitivity + mySpecificity)/2 
  
  
  write(str_replace(myColnames,"\\.y",""), file = paste0("Regression_", myAssay, ".txt"), append = T)
  write(paste0("AIC mean: ", myMeanAIC), file = paste0("Regression_", myAssay, ".txt"), append = T)
  write(paste0("AIC sd: ", mySDAIC), file = paste0("Regression_", myAssay, ".txt"), append = T)
  write(paste0("Balanced accuracy: ", myBalancedAccuracy), file = paste0("Regression_", myAssay, ".txt"), append = T)
  write(paste0("Sensitivity: ", mySensitivity), file = paste0("Regression_", myAssay, ".txt"), append = T)
  write(paste0("Specificity: ", mySpecificity), file = paste0("Regression_", myAssay, ".txt"), append = T)
  write("\n", file = paste0("Regression_", myAssay, ".txt"), append = T)
    
}


createFamilyBarplots <- function(myMaterialClass){
  
  if (myMaterialClass == "Silica"){
    myData <- data.frame(myMaterials=c("BaSO4 NM-220 (neg. control)", "Bentonite", "Kaolin", "Mn2O3 (pos. control)", "SiO2 Aerosil 200", "SiO2 Levasil 100", "SiO2 Levasil 300", "SiO2 Levasil 50", "SiO2 NM-203", "SiO2_15_Amino", "SiO2_15_Phospho", "SiO2_15_unmod"), myColors = c("grey","cadetblue2","deepskyblue","black","mediumpurple3","mediumorchid3","mediumorchid4","mediumorchid1","palevioletred","maroon2","deeppink3","hotpink1"))
    myData_complete <- merge(myData_merged,myData,by.x="NM",by.y="myMaterials")
  } else if (myMaterialClass == "IronOxide"){
    myData <- data.frame(myMaterials=c("BaSO4 NM-220 (neg. control)", "Fe2O3 larger", "Fe2O3 nanoform A", "Fe2O3 nanoform B", "Mn2O3 (pos. control)"), myColors = c("grey","darkolivegreen1", "darkolivegreen3","darkolivegreen","black"))
    myData_complete <- merge(myData_merged,myData,by.x="NM",by.y="myMaterials")
  } else if (myMaterialClass == "ZincOxide"){
    myData <- data.frame(myMaterials=c("BaSO4 NM-220 (neg. control)", "Mn2O3 (pos. control)", "ZnO NM-110", "ZnO NM-111"), myColors = c("grey","black","chocolate1","chocolate3"))
    myData_complete <- merge(myData_merged,myData,by.x="NM",by.y="myMaterials")
  }
           
  myMeanColnames_surface <- myMeanCols[which(grepl(paste(c("\\.y","FRAS_Mean"),collapse="|"),myMeanCols))]
  myMeanColnames_surface_withoutCarbonyls <- myMeanColnames_surface[-which(myMeanColnames_surface == "Carbonyls_Mean.y")]  
  
     
  createFamilyBarplot <- function(myColname){
    
    myNAInds <- which(is.na(eval(parse(text=paste("myData_complete$", myColname, sep = "")))))
    if(length(myNAInds) != 0){
      myOrder <- order(eval(parse(text=paste("myData_complete$", myColname, sep = "")))[-myNAInds])
      myMeans <- eval(parse(text=paste("myData_complete$", myColname, sep = "")))[-myNAInds][myOrder]
      mySDs <- eval(parse(text=paste("myData_complete$", gsub("Mean","SD",myColname), sep = "")))[-myNAInds][myOrder]
      myNames <- myData_complete$NM[-myNAInds][myOrder]
      myColors <- as.character(myData_complete$myColors[-myNAInds])[myOrder]
    } else{
      myOrder <- order(eval(parse(text=paste("myData_complete$", myColname, sep = ""))))
      myMeans <- eval(parse(text=paste("myData_complete$", myColname, sep = "")))[myOrder]
      mySDs <- eval(parse(text=paste("myData_complete$", gsub("Mean","SD",myColname), sep = "")))[myOrder]
      myNames <- myData_complete$NM[myOrder]
      myColors <- as.character(myData_complete$myColors)[myOrder]
    }
    
    
    myColnameForPlot <- str_replace_all(myColname, c("\\.x" = "_mass","\\.y" = "_surface", "AUC" = "AUC_mass", "FRAS_Mean" = "FRAS_Mean_surface"))
    myUnitName <- myUnitNames$Unit[which(myUnitNames$Assay == myColnameForPlot)]
   
    png(paste0("Barplot_", myMaterialClass, "_",str_replace(myColname,"\\.y",""),".png"), res=500, width=2000, height=2000)
    par(mar=c(9,6,4.1,2.1))
    par(cex=1)
    barCenters <- barplot(myMeans, beside = TRUE,log="y",cex.names =1,cex.axis=0.8,las=2,col=myColors,names=myNames,ylab=paste(unlist(strsplit(x = myUnitName, "(?<=blank |control |l ])", perl=T))[1],unlist(strsplit(x = myUnitName, "(?<=blank |control |l ])", perl=T))[2], sep="\n"))
    arrows(barCenters, myMeans+mySDs, barCenters, myMeans-mySDs, lwd = 1.5, angle = 90, code = 3, length = 0.01)
    dev.off()
    
  }
                 
  sapply(myMeanColnames_surface_withoutCarbonyls, createFamilyBarplot)

}


drawLumoPlots <- function(myData_merged){
  
  myDataToPlot <- myData_merged[which(!is.na(myData_merged$LUMO_Ec_mixed) & !is.na(myData_merged$CPH_Mean.y)),]
  
  p <- ggplot(myDataToPlot, aes(x=LUMO_Ec_mixed, y=log2(CPH_Mean.y)),col="darkgrey") +
    geom_point() + 
    xlim(-5, 0.5) + 
    xlab("LUMO/Ec") + 
    ylab("log2(ESR CPH)") +
    geom_vline(xintercept = -2.4, col="blue") +
    geom_vline(xintercept = -0.2, col="lightblue") +
    geom_vline(xintercept = -2.39, col="lightblue") + 
    geom_errorbar(aes(ymin = log2(CPH_Mean.y - CPH_SD.y), ymax = log2(CPH_Mean.y + CPH_SD.y),width=0.05,col="darkgrey")) +
    geom_rect(aes(xmin = -4.84, xmax = -4.12, ymin = -Inf, ymax = Inf),fill = "green", alpha = 0.03) +
    geom_text_repel(label=myDataToPlot$NM, size=3)
  ggsave(filename="Scatterplot_LUMOAndEcvsCPH.png", plot=p, dpi="retina")
  
  
  myDataToPlot <- myData_merged[which(!is.na(myData_merged$LUMO_Ec_mixed) & !is.na(myData_merged$DMPO_Mean.y)),]
  
  p <- ggplot(myDataToPlot, aes(x=LUMO_Ec_mixed, y=log2(DMPO_Mean.y)),col="darkgrey") +
    geom_point() + 
    xlim(-5, 0.5) + 
    xlab("LUMO/Ec") + 
    ylab("log2(ESR DMPO)") +
    geom_vline(xintercept = -2.4, col="blue") +
    geom_vline(xintercept = -0.2, col="lightblue") +
    geom_vline(xintercept = -2.39, col="lightblue") + 
    geom_errorbar(aes(ymin = log2(DMPO_Mean.y - DMPO_SD.y), ymax = log2(DMPO_Mean.y + DMPO_SD.y),width=0.05,col="darkgrey")) +
    geom_rect(aes(xmin = -4.84, xmax = -4.12, ymin = -Inf, ymax = Inf),fill = "green", alpha = 0.03) +
    geom_text_repel(label=myDataToPlot$NM, size=3)
  ggsave(filename="Scatterplot_LUMOAndEcvsDMPO.png", plot=p, dpi="retina")
  
}



runReactivityAnalysis <- function(myPath, myFileName_UniformMass, mySheetName_UniformMass, myFileName_UniformSurface, mySheetName_UniformSurface, myFileName_UnitNames, mySheetName_UnitNames){

  
  myOutDir <- "outputDirectory"

  if (file.exists(myOutDir)){
    setwd(file.path(myPath, myOutDir))
  } else {
    dir.create(file.path(myPath, myOutDir))
    setwd(file.path(myPath, myOutDir))
  }

  
  ### read data
  
  myData_uniformMass <- read_excel(paste0(myPath, "/", myFileName_UniformMass),mySheetName_UniformMass)
  myData_uniformSurface <- read_excel(paste0(myPath, "/", myFileName_UniformSurface),mySheetName_UniformSurface)
  myUnitNames <- read_excel(paste0(myPath, "/", myFileName_UnitNames),mySheetName_UnitNames)
  
  myData_merged <- merge(myData_uniformMass,myData_uniformSurface,by = "NM",all=T)
  
  
  ### define macrophage assay inds 
  
  myMacrophageAssayPosInds_Mass <- which(myData_merged$MacrophageAssay.x == "active")
  myMacrophageAssayNegInds_Mass <- which(myData_merged$MacrophageAssay.x == "passive")
  
  myMacrophageAssayPosInds_Surface <- which(myData_merged$MacrophageAssay.y == "active")
  myMacrophageAssayNegInds_Surface <- which(myData_merged$MacrophageAssay.y == "passive")
  
  
  ### define STIS inds 
  
  mySTISPosInds_Mass <- which(myData_merged$STIS.x == "active")
  mySTISNegInds_Mass <- which(myData_merged$STIS.x == "passive")
  
  mySTISPosInds_Surface <- which(myData_merged$STIS.y == "active")
  mySTISNegInds_Surface <- which(myData_merged$STIS.y == "passive")
  
  
  
  ### draw boxplots
  
  myMeanCols <- colnames(myData_merged)[which(grepl(paste(c("_Mean","AUC"),collapse="|"),colnames(myData_merged)))]
  
  myBoxplots_MacrophageAssay <- sapply(myMeanCols, drawBoxplots, "MacrophageAssay", myUnitNames)
  myBoxplots_STIS <- sapply(myMeanCols, drawBoxplots, "STIS", myUnitNames)
  
  
  ### draw heatmap
  
  myHeatmap <- drawHeatmap(myData_merged)
  
  
  ### compute correlations
  
  myMeanColnames_surface <- myMeanCols[which(grepl(paste(c("\\.y","FRAS_Mean"),collapse="|"),myMeanCols))]
  myMeanColnames_surface_withoutCarbonyls <- myMeanColnames_surface[-which(myMeanColnames_surface == "Carbonyls_Mean.y")]  
  myMeanColnames_surface_combinations <- combn(myMeanColnames_surface_withoutCarbonyls, 2)
  
  myCorrelations <- apply(myMeanColnames_surface_combinations,2,computeCorrelations)
  
  
  ### compute regression
  
  myMeanColnames_surface <- myMeanCols[which(grepl(paste(c("\\.y","FRAS_Mean"),collapse="|"),myMeanCols))]
  myMeanColnames_surface_withoutCarbonyls <- myMeanColnames_surface[-which(myMeanColnames_surface == "Carbonyls_Mean.y")]  
  
  for (i in 1:4){
    myRegression_MacrophageAssay <- apply(combn(myMeanColnames_surface_withoutCarbonyls, i), 2, computeRegression, "MacrophageAssay", i)
    myRegression_STIS <- apply(combn(myMeanColnames_surface_withoutCarbonyls, i), 2, computeRegression, "STIS", i)
  }

  
  ### Create barplots for NM families
  
  myFamilyBarplots <- sapply(c("Silica","IronOxide","ZincOxide"), createFamilyBarplots)

  
  ### Create LUMO plots
  
  myLumoPlots <- drawLumoPlots(myData_merged)
    
}


runReactivityAnalysis(myPath, myFileName_UniformMass, mySheetName_UniformMass, myFileName_UniformSurface, mySheetName_UniformSurface, myFileName_UnitNames, mySheetName_UnitNames)
