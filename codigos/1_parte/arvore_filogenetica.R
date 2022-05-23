# ---- 1° etapa Pré-processamento dos dados de censo visual ---- #
# ---- 1.2. Cosntrução da árvore filogenética ---- #

# a) Primeira parte
# Alterando o diretório de trabalho
setwd("C:/Users/emily/OneDrive/?rea de Trabalho/Mestrado/filogenia")
dir()

# chamando a biblioteca filogenética: http://www.phytools.org/Cordoba2017/ex/2/Intro-to-phylogenies.html
library(phytools)

# Antes é preciso baixar a árvore filogenética do Rabosky: https://fishtreeoflife.org/downloads/
# Lnedo o dado
fishurl2 = "actinopt_full.trees.xz"
fish_tree_full = ape::read.tree(fishurl2) 
str(fish_tree_full)
a <- fish_tree_full[[1]]

#"Numeros de nos" - os indices da borda do objeto "phylo", uma matriz que contem os indices de inicio e fim dos nos 
fish_tree_full[[1]]$edge

#Um objeto de classe "phylo" tambem pode ter outros componentes, o mais comum e o comprimento da borda (um vetor de classe "numerico" contendo todo o comprimento da borda)
fish_tree_full[[1]]$edge.length

#Label - especies
fish_tree_full[[1]]$tip.label

#Numero de nos internos da arvore
fish_tree_full[[1]]$Nnode

#Numero total de especies
Ntip(fish_tree_full[[1]])

#Remover os 'phylos' iguais (nao removeu nenhum)
uni<-unique(obj, incomparables = FALSE,
          use.edge.length = TRUE,
          use.tip.label = TRUE)

# Extrair especies que eu quero
species = c('Abudefduf_saxatilis',
             'Acanthurus_bahianus',
             'Acanthurus_chirurgus',
             'Acanthurus_coeruleus',
             'Acanthostracion_polygonius',
             'Aluterus_monoceros',
             'Amblycirrhitus_pinos',
             'Apogon_americanus',
             'Anisotremus_virginicus',
             'Balistes_capriscus',
             'Bodianus_pulchellus',
             'Bodianus_rufus',
             'Bothus_ocellatus',
             'Canthigaster_figueiredoi',
             'Cantherhines_pullus',
             'Caranx_latus',
             'Centropyge_aurantonotus',
             'Chaetodon_sedentarius',
             'Chaetodon_striatus',
             'Chilomycterus_spinosus',
             'Chromis_multilineata',
             'Coryphopterus_glaucofraenum',
             'Cryptotomus_roseus',
             'Dactylopterus_volitans',
             'Decapterus_macarellus',
             'Decapterus_punctatus',
             'Diplodus_argenteus',
             'Dules_auriga',
             'Elacatinus_figaro',
             'Epinephelus_margitus',
             'Gtholepis_thompsoni',
             'Gymnothorax_moringa',
             'Haemulon_aurolineatum',
             'Haemulon_plumierii',
             'Halichoeres_brasiliensis',
             'Halichoeres_poeyi',
             'Holocentrus_adscensionis',
             'Kyphosus_sectatrix',
             'Kyphosus_vaigiensis',
             'Labrisomus_nuchipinnis',
             'Malacoctenus_delalandii',
             'Mycteroperca_acutirostris',
             'Pagrus_pagrus',
             'Pempheris_schomburgkii',
             'Pinguipes_brasilianus',
             'Plectrypops_retrospinis',
             'Pronotogrammus_martinicensis',
             'Pareques_acumitus',
             'Parablennius_marmoreus',
             'Parablennius_pilicornis',
             'Pomacanthus_paru',
             'Priacanthus_aretus',
             'Pseudupeneus_maculatus',
             'Sargocentron_bullisi',
             'Stephanolepis_hispidus',
             'Scarus_zelindae',
             'Scorpae_plumieri',
             'Serranus_baldwini',
             'Sparisoma_axillare',
             'Sparisoma_frondosum',
             'Sparisoma_radians',
             'Sparisoma_tuiupiranga',
             'Sphoeroides_spengleri',
             'Stegastes_fuscus',
             'Stegastes_pictus',
             'Synodus_intermedius',
             'Synodus_synodus')
species 
class(species)

#sapply = pega uma lista, vetor ou dataframe e transforma em vetor ou matriz
#grep - busca de correspond?ncias de uma express?o/padr?o regular em um vetor 
#de caracteres. Ele retorna os ?ndices no vetor de caracteres que cont?m uma 
#correspond?ncia. grepl() retorna um vetor VERDADEIRO/FALSO indicando quais 
#elementos do vetor de caracteres cont?m uma correspond?ncia

# Lendo as especies do Rabosky
ii<-sapply(species,grep,fish_tree_full[[1]]$tip.label)
ii = unlist(ii)
class(ii)

#match entre as especies (minhas especies devem ser as mesma do Rabosky)
species<-fish_tree_full[[1]]$tip.label[ii]
species

# For para cortar para todas as especies
obj = list()
for (i in 1:100) {
  obj[[i]] <- keep.tip(fish_tree_full[[i]],species)
   
}

# Média de nós 
mean <- averageTree(obj, start = NULL, method="branch.score.difference")
mean

plotTree(mean,ftype="i")
plotTree(a,offset=1)
tiplabels()
nodelabels()

#Visualizando a árvore filogenética das minhas espécies
a <- drop.clade(mean, species)
plotTree(a)

# Salvando a árvore
write.tree(mean, file = "Mean.tree")

#-------------------------------------------------------------------------#
mean1 <- averageTree(obj, start = NULL, method="symmetric.difference")

#1. Symmetric.difference:**********
#2. Branch.score.difference (!!!)
#3. Path.difference: *********
#4. Quadratic.path.difference:

#TRUE or FALSE - Are species equal?
all.equal(anolis.noPR,anolis.noPR1)

comparePhylo(anolis.noPR, anolis.noPR1, plot = FALSE, force.rooted = FALSE,
             use.edge.length = FALSE)