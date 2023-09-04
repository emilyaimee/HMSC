library(sp)
library(raster)
library(ncdf4)
library(rgdal)
library(oceanmap)
library(wind_direction)

# carregando os dados .shp para as boxes e os dados .nc para TSM, CHL e Componente do vento
p_dentro = readOGR('D:/Dissertacao/Boxes/chl/dentro_chl.shp')
p_fora = readOGR("D:/Dissertacao/Boxes/chl/fora_chl.shp")
box = readOGR('D:/vento/box/box.shp')
setwd("D:/vento")
files = list.files(pattern = "*.nc")
data.sat = stack(files[1], varname = "v10")

#-------------------------------------------------#
box= c(-45, -40, -25, -21)
data.sat = crop(data.sat, box)
data.sat[data.sat < 0] = NA

v1 = crop(data.sat, box)
plot(data.sat[[1:5]])

data.sat = data.sat - 273.15

plot(data.sat[[308:312]])
data.sat[[1105:1128]]

#----recorte da area dentro da box e calculo da media------------------#
box_d = crop(data.sat[[308:312]], p_dentro)
min(box_d)

box_f = crop(data.sat[[308:312]], p_fora)
mean(box_f)

#-------------------------testando media ponderada--------------------#
test = as.data.frame(data.sat[[1:5]])
test = na.omit(test)
summary(test)

x = subset(test, test['analysed.sea.surface.temperature.5'] < 20.00100)
y = subset(test, test['analysed.sea.surface.temperature.5'] > 19.99999)
a = sum(x)/17
b = sum(y)/5033
c = a*17
d = b*5033
z = c+d
s = a+b/5049

med.pon = (((sum(x))*17)+((sum(y))*5033))/5049

#--------Delimitando periodo de 5 dias para os dados----------------#
idx = seq(as.Date('2017-11-04'), as.Date('2017-11-08'), 'day')
names(box_f) = idx

indices = format(as.Date(names(box_f), format = "X%Y.%m.%d"), format = "%y")
indices
indices = as.numeric(indices) 

# extraindo a media #
media = stackApply(box_f, indices, fun=mean)
plot(media)

valores = as.data.frame(media)
valores[1:5]
valores = na.omit(valores) 
summary(valores)

# extraindo a media de 25% dos valores maximos #
x = subset(valores, valores['index_19'] < 20.00000)
y = subset(valores, valores['index_17'] > 19.99940)

max25 = subset(valores, valores['index_17'] > 0.6535)
media25 = mean(max25$index_17)
media25

# media circular apenas para os dados de vento #
x = circ.mean(v1$X2017.03.12)
