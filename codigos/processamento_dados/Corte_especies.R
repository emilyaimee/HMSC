# ---- Etapa de pre-processamento dos dados de censo visual ---- #
# ---- Criterio de corte das especies ---- #

# a) Primeira parte
setwd("~/Documentos/Emily/corte_especie/Data") 
dir()

# Carregando o dado
d <- read.csv("data_arraial_esp_2017.csv", header = TRUE, sep = ',')
str(d)

# Separando por regiao
region <- unique(d$sites)
region

data = unique(d$data2)
data

# Inspecionando as regioes (dentro e fora)
data.frame(dplyr::summarise(dplyr::group_by(d[d$sites == "anequim",], code), 
                            abun = sum(abun, na.rm = TRUE), size = max(size_cm)))

# Salvando como RDS na pasta R_Objects
saveRDS(d, "D:/Dissertacao/corte_especie/R_Objects/Planilha_Nova_allspecies4.rds")

# b) Segunda parte 

#Alterando o diretório de trabalho para R_Objects
setwd("D:/Dissertacao/corte_especie/R_Objects") 
dir()
d <- readRDS("Planilha_Nova_allspecies4.rds")

# Obter apenas abundancia

# flat_to_matrix #
#Esta funsao remodela os dados registrados em censos visuais na comunidade
#formato matricial. Informa tambem o local, ano, profundidade e observador(es) de um censo.
#Se fornecida com uma lista de constantes a e b, a funsao tambem pode construir
#uma matriz comunitária baseada em biomassa (deve informar biom = TRUE)
#Esta funsao foi projetada para o banco de dados PELD ILOC e requer as colunas
#site, transect_id, ano e codigo (como codigo para nomes de especies). Para aplica-lo em
#outros conjuntos de dados, basta alterar os nomes das colunas.

flat_to_matrix <- function(Dt, tr_method = "tot", 
                           transform = TRUE) {
  
  library(magrittr)
# Check data
  if (!is.data.frame(Dt)) '<-'(Dt, as.data.frame(Dt))
  
# Soma das abundancias por sites e por ano.
  abun <- Dt %>% dplyr::group_by(sites, data2, code) %>%
    dplyr::summarise(abun = mean(abun)) %>%
    #Media das abundâncias
    reshape2::dcast(sites + data2 ~ code, 
                    value.var = "abun", fun.aggregate = mean, fill = 0)
  
  # Adicionando nomes
  rownames(abun) <- apply(abun[,c("sites", "data2")], 1,
                          paste, collapse = "_")
  abun <- abun[ , !colnames(abun) %in% c("sites", "data2")]
  
  
  # Se a biomassa estiver disponivel, devolva-a, caso contrario, devolva apenas abundancia
  retrn <- tryCatch(list("abun" = abun,),
                    error = function(e) list("abun" = abun))
  
  # Adicionando informacoes extras
  if (any(colnames(Dt) == "size_cm")) {
    # Adicioanando o tamanho a matrix 
    retrn$size <- dplyr::group_by(Dt, sites, data2, code) %>%
      dplyr::summarise(size = weighted.mean(size_cm, abun, na.rm = TRUE)) %>%
      reshape2::dcast(sites + data2 ~ code, 
                      value.var = "size", fun.aggregate = mean, fill = 9999)
    retrn$size <- retrn$size[,!colnames(retrn$size) %in% c("sites", "data2")]
  }
  
  
  if (any(colnames(Dt) == "depth_m")) {
    
    # Calculando a profundidade media 
    retrn$pred$depth <- dplyr::summarise(dplyr::group_by(Dt, sites, data2),
                                         depth = mean(depth_m, na.rm = TRUE))$depth
  }
  
  if (any(colnames(Dt) == "observer")) {
    # Use a funsao de *summarise* para obter cada observador do censo
    retrn$pred$obsvr <- dplyr::summarise(dplyr::group_by(Dt, sites, data2), 
                                         observer = paste(unique(observer), 
                                                          collapse = "_"))$observer
  }
  
  if (any(colnames(Dt) == "lat")) {
    
    # Summarise dados calculando a profundidade media.
    retrn$pred$lat <- dplyr::summarise(dplyr::group_by(Dt, sites, data2),
                                       lat = mean(lat, na.rm = TRUE))$lat
  }
  
  if (any(colnames(Dt) == "lon")) {
    
    # Summarise dados calculando a profundidade media.
    retrn$pred$lon <- dplyr::summarise(dplyr::group_by(Dt, sites, data2),
                                       lon = mean(lon, na.rm = TRUE))$lon
  }
  attr(retrn, "class") <- "censusftm"
  
  return(retrn)
}

# Dividido por regiao
d <- split(d, d$sites)

Data <- lapply(d, flat_to_matrix)

# Remove censos sem peixes
Data <- lapply(Data, function(x) {
  x[['mass']] <- x[['mass']][!apply(x[['abun']], 1, function(y) all(y == 0)),]
  x[['abun']] <- x[['abun']][!apply(x[['abun']], 1, function(y) all(y == 0)),]
  x[['size']] <- x[['size']][!apply(x[['abun']], 1, function(y) all(y == 0)),]
  x
})

# c) CORTE DAS ESPÉCIES (3 critérios)

require(magrittr)
# Manter apenas as especies principais por meio de criterios de *detecao* e *abundancia*.
Data <- Map(function(x, Dt, i) {
  
  # Criterio de presenca
  crt_abs <- data.frame(Dt %>% 
                          dplyr::group_by(code, year) %>% 
                          dplyr::summarise(.groups = "keep")) %>% 
    reshape2::dcast(code ~ year, fun.aggregate = function(x) {1}, fill = 0,
                    value.var = "year")
  
  # Crie os nomes das linhas dos codigos
  rownames(crt_abs) <- crt_abs$code
  crt_abs <- crt_abs[,-1]
  #A especie deve estar presente em pelo menos 70% da serie e ausente em apenas 30
  crt_abs <- rowSums(crt_abs) == ncol(crt_abs)
  
  # Criterio de abundancia
  crt_det <- data.frame(Dt %>% 
                          dplyr::group_by(code, year) %>% 
                          dplyr::summarise(n = rep(1, length(year)),
                                           .groups = "keep")) %>% 
    reshape2::dcast(code ~ year, value.var = "n", fun.aggregate = sum, fill = 0)
  # Crie os nomes das linhas dos codigos
  rownames(crt_det) <- crt_det$code
  crt_det <- crt_det[,-1]
  
  # Criterios minimos de registros
  # Quantos anos registraram uma determinada especie em pelo menos 5 transectos com 0.2 individuos por 40 m2?
  crt_det <- rowSums(crt_det > 5) > ((50/100)*ncol(crt_det))
  
  crt_abun <- data.frame(Dt %>% 
                           dplyr::group_by(transect_id, year, code) %>% 
                           dplyr::summarise(abun = sum(abun), .groups = "keep") %>%
                           reshape2::dcast(transect_id + year ~ code, fill = 0, value.var = "abun") %>%
                           reshape2::melt(id.vars = c("transect_id", "year"), value.var = "abun") %>%
                           dplyr::group_by(year, variable) %>%
                           dplyr::summarise(abun = mean(value), .groups = "keep") %>%
                           dplyr::group_by(variable) %>%
                           dplyr::summarise(abun = sd(abun)), .groups = "keep")
  crt_abun <- `names<-`(crt_abun$abun, as.character(crt_abun$variable))
  # 0.2 (4) individuos por 40 m2
  crt_abun <- crt_abun > 0.2 
  
  #Print  criterios
  cat(i,"\n")
  print(paste("Total of species:", length(crt_abs)))
  print(paste("Presence:", sum(crt_abs)))
  print(paste("Detectability:", sum(crt_det)))
  print(paste("Abundance:", sum(crt_abun)))
  print(paste("Combined:", sum(crt_abs | crt_abun | crt_det)))
  
  # Selecione as especies que atenderam a pelo menos 1 criterio
  RM <- crt_abs | crt_abun | crt_det
  x[['mass']] <- x[['mass']][,RM[names(x[["mass"]])]]
  x[['abun']] <- x[['abun']][,RM[names(x[["abun"]])]]
  x[['size']] <- x[['size']][,RM[names(x[["size"]])]]
  x
}, Data, d, names(Data))

#Transforme tamanhos 9999 para NAs
Data <- lapply(Data, function(x) {
  x$size[x$size == 9999] <- NA
  x
})

# Salvando as especies aprovadas como RDS na pasta R_Objects
saveRDS(Data, 'R_Objects/Data.rds')

# Lendo os valores de abunância para as especies filtradas por regiao 
a <- Data[['dentro']][['abun']]

# Salvando como csv.
write.table(a, file='ALL_dentro_02.csv', sep=';', dec=',', row.names=TRUE)

