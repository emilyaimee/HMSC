# Bibliotecas
library(ape)
library(MASS)
library(Hmsc)
library(usethis)
library(corrplot)
library(readxl)
library(Rcpp)
library(dplyr)
library(labdsv)

################################
#(2) Examining MCMC convergence
################################
#####EXAMINAR TODOS OS PAR?METROS
#converte o posterior hmsc em uma lista de objetos

x = 

modelo = readRDS("D:/Dissertacao/modelo/modeloOK/normal_todos.rds")
mpost = convertToCodaObject(modelo)
mpost$Omega
#Tamanho defetivo da amostra (Principais: Beta, Gamma e Omega)
## Interpretation: Species niches
ess.beta = effectiveSize(mpost$Beta)
psrf.beta = gelman.diag(mpost$Beta,multivariate=FALSE)$psrf
par(mfrow=c(1,2))
hist(ess.beta)
hist(psrf.beta)
par(mfrow=c(1,1))
par(mar = rep(4,4))
plot(mpost$Beta)
plot(ess.beta)
lines(ess.beta)
summary(mpost$Beta)
x = summary(mpost$Beta)
mean(ess.beta)

## Interpretation: Influence of traits on niches
ess.gamma = effectiveSize(mpost$Gamma)
psrf.gamma = gelman.diag(mpost$Gamma,multivariate=FALSE)$psrf
par(mfrow=c(1,2))
hist(ess.gamma)
hist(psrf.gamma)
par(mfrow=c(1,1))
plot(ess.gamma)
lines(ess.gamma)
mean(ess.gamma)

## Interpretation: Species associations
ess.omega = effectiveSize(mpost$Omega[[1]])
psrf.omega = gelman.diag(mpost$Omega[[1]],multivariate=FALSE)$psrf
par(mfrow=c(1,2))
hist(ess.omega)
hist(psrf.omega)
mean(ess.omega)
#For all effects ESS are high (near 4000 which is the theoretical maximum (we 
#sampled 4 chains x1000 samples)) and PSRF are near 1, so diagnostics are good.

## Interpretation: Phylogenetic signal in species niches
ess.rho = effectiveSize(mpost$Rho)
psrf.rho = gelman.diag(mpost$Rho,multivariate=FALSE)$psrf
ess.rho
summary(mpost$Beta)

## Interpretation: Residual covariance of species niches
ess.V = effectiveSize(mpost$V)
psrf.V = gelman.diag(mpost$V,multivariate=FALSE)$psrf
hist(ess.V)
hist(psrf.V)

##########################
#(3) Evaluating ALL_SUMMER fit
##########################
#probit avalia com o RMSE, AUC e Tjur R?
#outros com RMSE e R?

#Para analisar o poder explicat?rio do ALL_SUMMERo (evaluateALL_SUMMERFit), ser? 
#computado os valores da distribui??o posterior simulada
predY = computePredictedValues(modelo)
MF = evaluateModelFit(hM=modelo, predY=predY)
MF
df = data.frame(MF)
write.csv(df, file = "D:/teste/df4.csv", row.names = FALSE)
print(df)
hist(MF$TjurR2, xlim = c(0,1), main=paste0("Mean = ", round(mean(MF$TjurR2),2)))
MF$R2
#LER: https://raw.githubusercontent.com/hmsc-r/HMSC/master/vignettes/vignette_1_univariate.pdf
#Avalia o ajuste dos dados
#Avalia o ajuste dos dados
partition = createPartition(modelo, nfolds = 4, column = "SITES")
partition
#Y = (ALL_SUMMER$Y)
#partition.sp = Y[,c(4,5,15,20,27,28,36,38,42,44,49,52)] #SUMMER
#partition.sp
preds = computePredictedValues(modelo, partition = partition)
MFCV = evaluateModelFit(hM=modelo, predY=pred)
MFCV
df = data.frame(MFCV)
write.csv(df, file = "D:/teste/df8.csv", row.names = FALSE)

saveRDS(preds, "./pred_summer.rds" )
pred = readRDS("D:/Dissertacao/modelo/predicao/pred_probit_exo.rds")
x = data.frame(pred)
hist(MFCV$TjurR2, xlim = c(0,1), main=paste0("Mean = ", round(mean(MFCV$TjurR2),2)))


#O MF > MFCV
WAIC = computeWAIC(modelo)
par(mfrow=c(1,1))
boxplot(MF)
boxplot(MFCV)

par(mar=c(10,2,1.5,1))

VP = computeVariancePartitioning(modelo)
plotVariancePartitioning(modelo, VP=VP, args.legend=list(x="bottom", bty="n", cex=0.1, y.intersp = 0.8, inset=-0.42, ncol = 2))


plotVariancePartitioning(modelo, VP=VP, #xlab = "Species", 
                         names = c(1:26), 
                         cex.names = 1,cex.lab = 1,
                         args.legend=list(x="bottom", bty="n", cex=0.7, y.intersp = 0.8, inset=-0.40, ncol = 2))

par(mar=c(8,4,2,1))
install.packages("RColorBrewer")
library("RColorBrewer")
plotVariancePartitioning(modelo, VP=VP,names = c(1:72), col=brewer.pal(n=7, name='RdBu'), 
                         cex.names = 0.7,
                         cex.lab=0.1, cex.axis = 0.7,cex.main = 0.9,
                         args.legend=list(x="bottom", 
                                          bty="n",
                                          cex=0.8, 
                                          y.intersp = 0.8,
                                          x.intersp = 0.3,
                                          inset=-0.40, ncol = 2))

VP$R2T$Beta
VP$R2T$Y

#-----------------------------------------------------
OmegaCor = computeAssociations(modelo)
supportLevel = 0.90
for (r in 1:modelo$nr){
  toPlot = ((OmegaCor[[r]]$support>supportLevel) +
              (OmegaCor[[r]]$support<(1-supportLevel))>0)*OmegaCor[[r]]$mean
  corrplot(toPlot, method = 'color',
         col = colorRampPalette(c('blue','white','red'))(200),
         tl.cex = .6,tl.col = 'black',
         title = paste('random effect level'),mar = c(0,0,1,0))}

###################################
#(4) Exploring parameter estimates
###################################
# Phylogenetic patterns in environmental response
postBeta = getPostEstimate(modelo, parName="Beta")
plotBeta(modelo, post=postBeta, param = "Support",
         plotTree = T, supportLevel = 0.90, split = 0.4, spNamesNumbers = c(F,T),
         marTree = c(4,0,0,1),mar = c(4,0,0,1))

plotBeta(modelo, post=postBeta, param="Support", 
         supportLevel = 0.90, split = 0.4,
         plotTree=TRUE, spNamesNumbers=c(FALSE,FALSE),
         cex = c(0.7,0.7,0.8),marTree = c(4, 3, 2, 0),
         mar = c(9,2,0,1))

plotBeta(modelo, post=postBeta, param="Support", 
         supportLevel = 0.90, plotTree=TRUE, split = 0.4,
         marTree = c(8, 0, 0, 1),
         spNamesNumbers=c(F,T), cex = c(0.5,0.5,0.7),
         mar = c(8,2,0,0), newplot = TRUE)


dev.off()


# Visualise effects of traits
postGamma = getPostEstimate(modelo, parName="Gamma")
plotGamma(modelo, post=postGamma, param="Support",  
          supportLevel = 0.90, trNamesNumbers=c(TRUE,TRUE),
          cex = c(0.8,0.6,0.8), 
          mar = c(14,14,3,0), 
          main = ("Trait effects on species environmental responses"))

par(mfrow=c(1,2))
# Species associations
OmegaCor = computeAssociations(modelo)
supportLevel = 0.90
for (r in 1:modelo$nr){
  plotOrder = corrMatOrder(OmegaCor[[r]]$mean,order="AOE")
  toPlot = ((OmegaCor[[r]]$support>supportLevel) +
              (OmegaCor[[r]]$support<(1-supportLevel))>0)*OmegaCor[[r]]$mean
  par(xpd=T)
  colnames(toPlot)=rownames(toPlot)=gsub("_"," ",x=colnames(toPlot))
  corrplot(toPlot[plotOrder,plotOrder], method = "color",
           col=colorRampPalette(c("blue","white","red"))(200),
           title="",type="lower",tl.col="black",tl.cex=.4, mar=c(0,0,0,0))
}

for (r in 1:modelo$nr){
  plotOrder = corrMatOrder(OmegaCor[[r]]$mean,order="AOE")
  toPlot = ((OmegaCor[[r]]$support>supportLevel) +
              (OmegaCor[[r]]$support<(1-supportLevel))>0)*OmegaCor[[r]]$mean
  par(xpd=T)
  colnames(toPlot)=rownames(toPlot)=gsub("_"," ",x=colnames(toPlot))
  corrplot(toPlot[plotOrder,plotOrder], method = "color",
           col=colorRampPalette(c("blue","white","red"))(200),
           title="",type="lower",tl.col="black",tl.cex=.5, mar=c(0,0,0,0))
}

## (5) Making predictions
Gradient = constructGradient(modelo, focalVariable="TSM", ngrid = 50)
Gradient$XDataNew

predY = predict(modelo, XData=Gradient$XDataNew, studyDesign=Gradient$studyDesignNew,
                ranLevels=Gradient$rLNew, expected=TRUE)

par(mfrow=c(1,1))
plotGradient(modelo, Gradient, pred=predY, measure="S", las=1,
             showData = TRUE, main='Species richness (measure="S")')

