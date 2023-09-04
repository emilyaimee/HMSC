#Diretorio
setwd("D:/Dissertacao/Matrixes")
dir()

# Bibliotecas
library(ape)
library(MASS)
library(Hmsc)
library(corrplot)
library(readxl)
library(Rcpp)
library(dplyr)
library(labdsv)


set.seed(1)

data = read_xlsx("m_Y.xlsx")

Y = as.matrix(data[,7:78])

# codigo para matriz de presença e ausencia #
#Y = replace(x = Y, list = is.na(Y), values = 0)
#Y = as.matrix(1*(Y>0))
# ------------------------------------------ #

Y = log(Y+1)
#Y[Y==0]=NA
#Y = 1*(Y>0)
# Pegando apenas os dados referentes ?s vari?veis ambientais
XData = data.frame(data[,79:81])
names(XData)

XData = XData[,c(1,2,3)]

colnames(XData) = c("CHL", "TSM", "CMV")


# Abrindo a matriz com os dados das caracter?sticas das esp?cies
TrData = read_xlsx("Matrix_T_2.xlsx")
rownames(TrData)
colnames(Y)
rownames(TrData) = colnames(Y)
rownames(TrData) == colnames(Y)

#Dados Filogen?ticos
phyloTree <- read.tree("Mean2.tree")
plot(phyloTree)
phyloTree$tip.label = colnames(Y)
phyloTree$tip.label == colnames(Y)
plot(phyloTree)

# STUDY DESIGN
# Ano > Ilha > Sitio (Temporalmente expl?cito)
studyDesign = matrix(NA,nrow(Y),4)
studyDesign[,1] = data$Locality
studyDesign[,2] = data$Sites
studyDesign[,3] = data$Month
studyDesign[,4] = data$Year
studyDesign = as.data.frame(studyDesign)
colnames(studyDesign) = c("LOCALITY", "SITES", "MONTH","YEAR")
studyDesign[,1]=as.factor(studyDesign[,1])
studyDesign[,2]=as.factor(studyDesign[,2])
studyDesign[,3]=as.factor(studyDesign[,3])
studyDesign[,4]=as.factor(studyDesign[,4])
head(studyDesign)

# RANDOM EFFECT STRUCTURE

data$Lat = as.numeric(as.character(data$Lat))
data$Lon = as.numeric(as.character(data$Lon))
#xy-coordenadas
sites = levels(studyDesign[,2])
sites
nsites = length(sites)
nsites
xy = matrix(0, nrow = nsites, ncol = 2)
for (i in 1:nsites){
  rows=studyDesign[,2]==sites[[i]]
  xy[i,1] = mean(data[rows,]$Lon)
  xy[i,2] = mean(data[rows,]$Lat)
}
colnames(xy) = c("x","y")
sRL = xy
sRL
rownames(sRL) = sites
plot(xy, asp=1)

rL = HmscRandomLevel(sData=sRL)
rL2 = HmscRandomLevel(units = levels(studyDesign$LOCALITY))
rL3 = HmscRandomLevel(units = levels(studyDesign$MONTH))
rL4 = HmscRandomLevel(units = levels(studyDesign$YEAR))

names(TrData)

TrFormula = ~ body_size_max + diet + depth_range + Geographic_range_index 


XFormula = ~ poly(CHL, degree=1, raw = TRUE) + 
  poly(TSM, degree=1, raw = TRUE) + 
  poly(CMV, degree=1, raw = TRUE) 


# check that community data are numeric and have finite numbers. If the script
# writes "Y looks OK", you are ok.
if (is.numeric(as.matrix(Y)) || is.logical(as.matrix(Y)) && is.finite(sum(Y, na.rm=TRUE))) {
  print("Y looks OK")
} else {
  print("Y should be numeric and have finite values")	}
# Check that the stydy design data do not have missing values (they are allowed for Y but not S, X, P or Tr)
if (any(is.na(studyDesign))) {
  print("S has NA values - not allowed for")
} else {
  print("S looks ok")	}
# Check that the covariate data do not have missing values (they are allowed for Y but not S, X, P or Tr)
if (any(is.na(XData))) {
  print("X has NA values - not allowed for")
} else {
  print("X looks ok")	}

# Check that the Tr data do not have missing values (they are allowed for Y but not S, X, P or Tr)
if (any(is.na(TrData))) {
  print("Tr has NA values - not allowed for")
} else {
  print("Tr looks ok")	}
# Check that the phylogenetic/taxonomic data do not have missing values (they are allowed for Y but not S, X, P or Tr)
if (any(is.na(phyloTree))) {
  print("P has NA values - not allowed for")
} else {
  print("P looks ok")	}

model = Hmsc(Y=Y,
              XData = XData, XFormula = XFormula,
              TrData = TrData, TrFormula = TrFormula, 
              phyloTree = phyloTree,
              distr = "normal",
              studyDesign = studyDesign, 
              ranLevels = list("SITES" =rL, "LOCALITY"=rL2, 'MONTH'=rL3, "YEAR" =rL4))


nChains = 4
#nParallel = 1 #(pode dar bug - e-mail do ovaskainen)
samples = 1000
thin = 100
transient = 50*thin

ptm = proc.time()

model = sampleMcmc(model, thin = thin, samples = samples, transient = transient,
                   nChains = nChains)

computational.time = proc.time() - ptm
computational.time

saveRDS(model, "normal_exo.rds")

