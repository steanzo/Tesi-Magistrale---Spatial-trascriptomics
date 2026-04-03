library(dplyr)
library(INLA)
library(ggplot2)
library(ggregplot)
library(patchwork)
load("dati/dati_finali.RData")


m <- c(21, 41, 61, 81, 101)
n <- c(40, 60, 80, 100, 118)

load("modelli/GFP e WGA/NB-SPDE comp geni 1-20.RData")

ris <- lapply(prova.multi.BN.SPDE, function(x){
  x$summary.fixed
})

for (el in 1:length(m)){
  cat(sprintf("\r Esecuzione... %d - %d", m[el], n[el]))
  load(paste0("modelli/GFP e WGA/NB-SPDE comp geni ",m[el],"-",n[el],".RData"))
  tmp.ris <- lapply(prova.multi.BN.SPDE, function(x){
    x$summary.fixed
  })
  ris <- c(ris, tmp.ris)
}


# Ossero significatività logGFP ------------------------------------------------
mtr.logGFP <- lapply(ris, function(x){
  tmp <- x["logGFP",]
  names(tmp) <- colnames(x)
  tmp
}) 
mtr.logGFP <- as.data.frame(do.call(rbind, mtr.logGFP))

mtr.logGFP$signif <- 0
mtr.logGFP$signif[sign(mtr.logGFP$`0.025quant`) == sign(mtr.logGFP$`0.975quant`)] <- 1
mtr.logGFP[,c(1,3,5,8)]

# osservo grafici e guardo la significatività se reale o meno
g1 <- data %>% 
  filter(img == "c26t16") %>%
  ggplot() +
  geom_point(aes(x = X1, y= -Y1, col= log(GFP_ch1_mean+1))) +
  scale_color_gradientn('intensità',       # titolo della legenda
                        colors = c("black", "red", "yellow"))
g2 <- data %>% 
  filter(img == "c26t17") %>%
  ggplot() +
  geom_point(aes(x = X1, y= -Y1, col= log(GFP_ch1_mean+1))) +
  scale_color_gradientn('intensità',       # titolo della legenda
                        colors = c("black", "red", "yellow"))

g.logGFP <- g1 + g2 

list.grafici = list()
plot.img.logGFP.vs.Gene <- function(gene){
  segno <- sign(mtr.logGFP[gene, "mean"])
  if(mtr.logGFP[gene, "signif"] == 0){
    titolo <- paste0(gene, ": non signfificativo (", segno, ")")
  } else{
    titolo <- paste0(gene, ": signfificativo (", segno, ")")
  }
  titolo
  g3 <- data %>%
    filter(img == "c26t16") %>%
    ggplot() +
    geom_point(aes(x = X1, y = -Y1, col = data.gene[data$img == "c26t16", gene])) +
    scale_color_gradientn('intensità',       # titolo della legenda
                          colors = c("black", "red", "yellow"))  # scala di colori: nero -> rosso -> giallo
  
  g4 <- data %>%
    filter(img == "c26t17") %>%
    ggplot() +
    geom_point(aes(x = X1, y= -Y1, col = (data.gene[data$img == "c26t17", gene]))) +
    scale_color_gradientn('intensità',       # titolo della legenda
                          colors = c("black", "red", "yellow"))
 g.gene <- g3 + g4 
 
 print(g.logGFP/g.gene + plot_annotation(titolo))
 #readline(prompt = "Premi Invio per vedere il prossimo grafico...")
}





tmp <- cbind(data.gene[,rownames(mtr.logGFP)], log(data$GFP_ch1_mean+1))
relazioni <- cor(tmp)[,119][-119]

confronto.cor.signif <- cbind(mtr.logGFP$signif, relazioni, sign(mtr.logGFP$mean)) %>% as.data.frame()
colnames(confronto.cor.signif) <- c("signif", "corr", "segno")
str(confronto.cor.signif)


# osservo significativi positivi 
sign.pos <- confronto.cor.signif %>% filter(signif == 1, segno == 1) 
range(sign.pos$corr)
hist(sign.pos$corr)


# questi sono i significativi positivi con maggiore correlazione
rownames(sign.pos)[which(sign.pos$corr>=0.25)]
plot.img.logGFP.vs.Gene("Cryab")
plot.img.logGFP.vs.Gene("Eno3")
plot.img.logGFP.vs.Gene("Hspb8")

#Geni significativi con minore corr
rownames(sign.pos)[which(sign.pos$corr<0.15)]
plot.img.logGFP.vs.Gene("Car3") # non molto chiaro
plot.img.logGFP.vs.Gene("Tpm2") 
plot.img.logGFP.vs.Gene("Myl1")


# osservo non significativi positivi 
no.sign.pos <- confronto.cor.signif %>% filter(signif == 0, segno == 1) 
range(no.sign.pos$corr)
hist(no.sign.pos$corr)

# quali geni hanno corr positiva ma non sono considerati segnali
rownames(no.sign.pos)[which(no.sign.pos$corr>=0.2)]

plot.img.logGFP.vs.Gene("Des") # sembra esseci relazione positiva
plot.img.logGFP.vs.Gene("Slc25a4")

# osservo significativi negativi 
sign.neg <- confronto.cor.signif %>% filter(signif == 1, segno == -1) 

plot.img.logGFP.vs.Gene("Foxo1")


# osservo non significativi negitivi 
no.sign.neg <- confronto.cor.signif %>% filter(signif == 0, segno == -1) 
range(no.sign.neg$corr)
hist(no.sign.neg$corr)

# quali geni hanno corr negitiva ma non sono considerati segnali
rownames(no.sign.neg)[which(no.sign.neg$corr<=-0.1)]
no.sign.neg$corr[no.sign.neg$corr<=-0.1]
# Provaiamo a giustificare alcuni casi:

plot.img.logGFP.vs.Gene("Trim63")
  #' ha correlazione -0.24.
  #' ho correlazione negativa, ma bassa. possiamo giustificare la 
  #' non sigificatività poichè solo per 3 tessuti si edivenzia
  #' una relazione negativa (anche se non su tutta l'are del tessuto)


plot.img.logGFP.vs.Gene("Fbxo32")
plot.img.logGFP.vs.Gene("Serpina3n")
plot.img.logGFP.vs.Gene("Foxo3") #non veod nulla



#+ Commentare ben:
#+ Osservando i grafici non si individuano relazioni logGFP e espressione genica perfette
#+ Nella maggior parte dei casi se utilizzassimo solo alcuni tessuti si osserverebbero 
#+ delle plausibili relazioni 


# Osservo quante variabili mantenere -------------------------------------------
mtr.Var <- lapply(ris, function(y){
  apply(y[-1,], 1, function(x){
    sign(x[3]) == sign(x[5])
  })
}) 
mtr.Var <- as.data.frame(do.call(rbind, mtr.Var))

colSums(mtr.Var)

# elimino dapi (1) e count.nucleus (7) sono significative in pochi geni
