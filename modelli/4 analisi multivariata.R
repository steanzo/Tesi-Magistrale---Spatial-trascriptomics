#+ In questo file applico il mdoello spaziale solo per 50 geni. Come faccio:
#+ 1. applico modello bn (no SPDE) e provo ad individuare quali geni sono
#+    differenzialmente esperessi
#+ 2. applico modello BN con SPDE su tali geni


library(ggplot2)
library(dplyr)
library(Voyager)
library(SingleCellExperiment)
library(SpatialFeatureExperiment)
library(scater)
library(spatialreg)
library(spdep)
library(edgeR)
library(INLA)
library(INLAutils)
library(MASS)
load("dati/dati_finali.RData")

rownames(data) = data$X
rownames(data.gene) = data.gene$X
data$X2 <- NULL
data$Y2 <- NULL

# Creo oggetto SingleCellExperiment e trasformo variabili
count = t(data.gene[,-1])
max((data$GFP_ch1_mean))
data$logGFP<- log(data$GFP_ch1_mean+1)
data$logGFP[data$logGFP>0] <- 1

data$dapi <- sqrt( 0.2126*data$dapi_ch0_mean + 0.7152*data$dapi_ch1_mean + 
                     0.0722*data$dapi_ch2_mean)

sce <- SingleCellExperiment(assays = list(counts = count),
                            colData = data[,-1])

sce <- addPerCellQC(sce)

# faccio filtraggio geni e normalizzazione

filter <- rowMeans(assay(sce))>=10
geni.interesse <- c(rownames(assay(sce))[filter], "Foxo1", "Foxo3")
sce <- sce[geni.interesse,]

dge <- calcNormFactors(sce, method = "TMM")

design <- model.matrix(~ logGFP + WGA_ch0_mean, 
                       data = colData(sce))
dge <- estimateDisp(dge, design)


# adatto modello BN multivariato
fit <- glmFit(dge, design)

# testiamo nullità coefficiente logGFP
fit$coefficients
res <- glmLRT(fit, coef = 2)

top <- topTags(res, n=Inf)$table
head(top)
table(top$FDR<=0.05)



# Modelli BN-SPDE per diversi geni ---------------------------------------------
# Partiamo considerando i primi 20 per tempi computazionali.

# parte che è in comune per ogni modello

coords <- data[,c("X1", "Y1", "id_tissue")]

for (el in unique(coords$id_tissue)){
  x_1 <- mean(coords[coords$id_tissue == el, "X1"])
  y_1 <- mean(coords[coords$id_tissue == el, "Y1"])
  coords[coords$id_tissue == el, "X1"] <- coords[coords$id_tissue == el, "X1"] - x_1
  coords[coords$id_tissue == el, "Y1"] <- coords[coords$id_tissue == el, "Y1"] - y_1
}

coords$id_tissue <- NA
coords <- as.matrix(coords[,c("X1", "Y1")])
mesh0 <- inla.mesh.2d(loc.domain = coords, max.edge = c(20,30))

Groups = "id_tissue"

NGroups <- length(unique(data[,Groups])) 


A <- inla.spde.make.A(mesh0, # Leave
                       loc = coords, # Leave
                       group = as.numeric(as.factor(data[,Groups])),
                       n.group = NGroups) 

spde <- inla.spde2.matern(mesh = mesh0, alpha = 2)

w <- inla.spde.make.index(
  name    = 'w', 
  n.spde  = spde$n.spde,
  n.group = NGroups) 

tmm_factors <- calcNormFactors(assay(sce), method = "TMM")

data.nb <-  data.frame(logGFP = data$logGFP, 
                       WGA = data$WGA_ch0_mean, dapi = data$dapi,
                       count.nucleus = data$count.nucleus, tmm_factors = tmm_factors
) 

X0 <- model.matrix(~ -1 + logGFP + count.nucleus + WGA + dapi +
                     tmm_factors, 
                   data = data.nb) 
X <- as.data.frame(X0) 
N <- nrow(data)


f.comp = y ~ -1 + Intercept + logGFP + WGA  +
  f(w, model = spde, group = w.group, control.group = list(model = 'iid')) +
  f(id_tissue, model = "iid")
#+ offset(tmm_factors)

f.null = y ~ -1 + Intercept + WGA +
  f(w, model = spde, group = w.group, control.group = list(model = 'iid')) +
  f(id_tissue, model = "iid") + offset(tmm_factors)



#+ Creo un funzione che inserirò nell'apply. Richiede:
#+ - vettore risposta (geni interesse)
#+ -
modello.SPDE <- function(risp, data.gene, data, X, N, w, spde, A, f, tmm_factors){
  Stack.est <- inla.stack( 
    data = list(y = data.gene[, risp]), # Leave
    A = list(1,1, 1,A), # Change the A matrix to the new one
    effects = list(
      Intercept = rep(1, N), # Leave
      X = X, # Leave
      id_tissue = data$id_tissue,
      w = w),
    tag = "est") # CHANGE
  Stack.pred <- inla.stack( 
    data = list(y = rep(NA, length(data.gene[, risp]))), # Leave
    A = list(1,1,1, A), # Change the A matrix to the new one
    effects = list(
      Intercept = rep(1, N), 
      X = X, 
      id_tissue = data$id_tissue,
      w = w),
    tag = "pred")
  stack <- inla.stack(Stack.est, Stack.pred)
  
  IM6 <- inla(f,
              family = "nbinomial",
              data = inla.stack.data(stack), # Don't forget to change the stack!
              control.compute = list(dic = TRUE, cpo = TRUE,
                                     return.marginals.predictor=TRUE),
              control.predictor = list(A = inla.stack.A(stack), compute=TRUE, 
                                       link = 1),
              control.inla = list(strategy = "laplace")
  )
  return(IM6)
}
length(rownames(top))
m <- 21
n <- 40
geni.prova <- rownames(top)[m:n]

#geni.prova <- c("Foxo1", "Foxo3")
counter <- 0
execution_time <- system.time({
  prova.multi.BN.SPDE <- lapply(geni.prova,
                                function(risp){
                                  counter <<- counter + 1  # Incrementa il contatore
                                  #cat(risp," --> Elemento numero", counter, "di", length(geni.prova), "\n")
                                  cat(sprintf("\r %s ... %d di %d", risp, counter, length(geni.prova)))
                                  modello.SPDE(risp, data.gene, data, X, N, w, spde, A, f.comp, tmm_factors)
                                }
                               )
  })

nomi.geni = geni.prova
names(prova.multi.BN.SPDE) <- geni.prova
save(prova.multi.BN.SPDE, nomi.geni, file = paste0("modelli/no_offset/NB-SPDE comp geni ",m,"-",n,".RData"))

prova.multi.BN.SPDE$Cryab$summary.fixed

#save(prova.multi.BN.SPDE, nomi.geni, file = paste0("NB-SPDE.comp geni Foxo1-3 .RData"))


# faccio con i modelli nullo senza logGFP
#load("NB-SPDE-null geni_21-41.RData")
# execution_time <- system.time({
#   prova.multi.BN.SPDE.null <- lapply(rownames(top)[m:n],
#                                 function(risp){
#                                   counter <<- counter + 1  # Incrementa il contatore
#                                   cat(risp," --> Elemento numero", counter, "di", length(rownames(top)[m:n]), "\n")
#                                   modello.SPDE(risp, data.gene, data, X, N, w, spde, A, f.null)
#                                 }
#   )
# })
# 
# names(prova.multi.BN.SPDE.null) <- rownames(top)[m:n]
# save(prova.multi.BN.SPDE.null, nomi.geni, file = paste0("NB-SPDE-null geni_",m,"-",n,".RData"))


#+ estraggo p-value (come se fosse frequentista) e estremi superiore e inferiore
#+ con livell osignfiicatività 0.05
prova.risultati <- lapply(prova.multi.BN.SPDE, 
                          function(x){
                            coef.logGFP <- x$summary.fixed["logGFP",]
                            t.logGFP <- coef.logGFP$mean/coef.logGFP$sd
                            p.value <- 2 * (1 - pnorm(abs(t.logGFP)))
                            estremo.inf <- x$summary.fixed["logGFP",3]
                            estremo.sup <- x$summary.fixed["logGFP",5]
                            ris <- c(p.value, estremo.inf, estremo.sup)
                            names(ris) <- c('p.value', 'estremo.inf', 'estremo.sup')
                            return (ris)
                          })

prova.multi.BN.SPDE$Pdk4$summary.fixed


matrix.result <- as.data.frame(do.call(rbind, prova.risultati))
matrix.result <- matrix.result[order(matrix.result$p.value),]

matrix.result$p.adj <- p.adjust(matrix.result$p.value, "BH")
matrix.result
rownames(matrix.result)[which(matrix.result$p.adj<0.05)]


# altro metodo Bayesian factor

bayes.factor <- sapply(names(prova.multi.BN.SPDE), 
       function(x){
         log_marginal_likelihood_full <- prova.multi.BN.SPDE[[x]]$mlik[1,]
         log_marginal_likelihood_null <- prova.multi.BN.SPDE.null[[x]]$mlik[1,]
         
         bayes_factor <- exp(log_marginal_likelihood_full - log_marginal_likelihood_null) 
       })
names(bayes.factor) <- names(prova.multi.BN.SPDE)

matrix.result$bayes.factors <- bayes.factor[rownames(matrix.result)]
matrix.result
sort(bayes.factor)


# probabilistic direction 
library(bayestestR)

pd.value <- sapply(names(prova.multi.BN.SPDE), 
                       function(x){
                         posterior_logGFP <- prova.multi.BN.SPDE[[x]]$marginals.fixed$logGFP
                         # Calculate the Bayesian p-value equivalent
                         if (mean(prova.multi.BN.SPDE[[x]]$summary.fixed["logGFP", "mean"]) > 0) {
                           # If the posterior mean is positive, calculate P(beta_X1 <= 0)
                           p <- 1- inla.pmarginal(0, posterior_logGFP)
                         } else {
                           # If the posterior mean is negative, calculate P(beta_X1 >= 0)
                           p <- inla.pmarginal(0, posterior_logGFP)
                         }
                         pd_to_p(p)
                       })
matrix.result$pd.value <- pd.value[rownames(matrix.result)]
matrix.result$pd.value.adj <- p.adjust(matrix.result$pd.value, "BH")
matrix.result

save(matrix.result, file =  paste0("risultati geni ",m,"-",n,".RData"))

View(matrix.result)


# eseguendo su server di risso escono i seguenti risultati:
#' significativi con correzione pd.value (approssimato in p.value):
#' "Cryab"  "Hspb8"  "Mb"     "Fhl1"   "Eef1b2"