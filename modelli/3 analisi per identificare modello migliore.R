# Provo a seguire processo normale per analisi di dati omici

load("dati/dati_finali.RData")

library(dplyr)
library(SingleCellExperiment)
library(SpatialFeatureExperiment)
library(scater)
library(Voyager)
library(spatialreg)
library(spdep)
library(edgeR)
library(INLA)
library(INLAutils)

risp = "Trim63"
#risp ="Cryab"


# Osserviamo se necassitiamo di normalizzare i dati ----------------------------

rownames(data) = data$X
rownames(data.gene) = data.gene$X
data$X2 <- NULL
data$Y2 <- NULL
count = t(data.gene[,-1])

data$logGFP<- log(data$GFP_ch1_mean+1)
plot(density(data$logGFP))

plot(density(data$WGA_ch0_mean))
data$WGA_transfrom <- log(data$WGA_ch0_mean+100)

plot(density(data$WGA_transfrom))

data$dapi <- sqrt( 0.2126*data$dapi_ch0_mean + 0.7152*data$dapi_ch1_mean + 
                     0.0722*data$dapi_ch2_mean)

# creo oggetto spatial data
sce <- SingleCellExperiment(assays = list(counts = count),
                            colData = data[,-1])


# Osserviamo il conteggi otoale di geni come si distribuzisce 
sce <- addPerCellQC(sce)
barplot(sce$total, col=sce$growth_cond, las=2)
  # il numero di geni rilevati non è uniforme tra i campioni

# Filtriamo per i geni poco espressi e conserviamo Foxo1 e Foxo3
filter <- rowMeans(assay(sce))>=10
sum(filter)
rownames(assay(sce))[grepl("Foxo", rownames(assay(sce)))]
geni.interesse <- c(rownames(assay(sce))[filter], "Foxo1", "Foxo3")

sce <- sce[geni.interesse,]

# Normalizzazione TMM
assay(sce, "logcounts") <- log1p(assay(sce))
plotRLE(sce[,1:50], style="full") + geom_hline(yintercept = 0)
tmm_factors <- calcNormFactors(assay(sce), method = "TMM")

sizeFactors(sce) <- tmm_factors * sce$total
sce <- logNormCounts(sce, name="tmm", transform="none")

boxplot(log1p(assay(sce, "tmm")[,1:20]), las=2)
plotRLE(sce[,1:20], style="full", exprs_value = "tmm", exprs_logged = FALSE) + 
  ggtitle("TMM")  + geom_hline(yintercept = 0)



# ho dei fattori di normalizzazione


# Sovradispersione nei dati? ---------------------------------------------------
design <- model.matrix(~ logGFP+ count.nucleus + WGA_ch0_mean + dapi, 
                       data = colData(sce))

dge <- calcNormFactors(sce, method = "TMM")
dge <- estimateDisp(dge, design)
plotMeanVar(dge, show.raw.vars = TRUE, show.tagwise.vars = TRUE,
            show.ave.raw.vars = FALSE)

 # Vi è una lieve sovradispersione



# Modello binomiale negativo ---------------------------------------------------

#+ Poichè siamo in presenza di sovradispersione andremo ad utilizzare un modello
#+ binomiale negativo ed per consoiderare la normalizzazioen aggiungeremo un 
#+ offset 

# -Usiamo modello binomiale negativo con offset
tmm_factors*colData(sce)$total

library(MASS)


data.nb <-  data.frame(risp = data.gene[,risp], logGFP = data$logGFP, 
                       WGA = data$WGA_ch0_mean, dapi = data$dapi,
                       count.nucleus = data$count.nucleus, tmm_factors = tmm_factors
                       ) 
mod1 <- glm.nb(risp ~ logGFP + WGA+ dapi + count.nucleus + offset(tmm_factors),
               data = data.nb)
summary(mod1)

mod2 <- glm.nb(risp ~ logGFP + WGA+ dapi + count.nucleus,
               data = data.nb)
summary(mod2)


par(mfrow=c(2,2))
plot(mod1)
par(mfrow=c(1,1))


# Considero effetto spaziale ---------------------------------------------------

## Tratto le coordinate --------------------------------------------------------
coords <- data[,c("X1", "Y1", "id_tissue")]
coords %>%
  ggplot() +
  geom_point(aes(x = X1, y = Y1))


for (el in unique(coords$id_tissue)){
  x_1 <- mean(coords[coords$id_tissue == el, "X1"])
  y_1 <- mean(coords[coords$id_tissue == el, "Y1"])
  coords[coords$id_tissue == el, "X1"] <- coords[coords$id_tissue == el, "X1"] - x_1
  coords[coords$id_tissue == el, "Y1"] <- coords[coords$id_tissue == el, "Y1"] - y_1
}


coords %>%
  ggplot() +
  geom_point(aes(x = X1, y = Y1))


plot(density(dist(coords)))

hist(dist(data[data$id_tissue == "c26t16_1",c("X1", "Y1")]), breaks = 500)
abline(v = 50, col = "green")

library(reshape2)
x <- dist(data[data$id_tissue == "c26t16_1",c("X1", "Y1")])
df <- melt(as.matrix(x), varnames = c("row", "col"))

ggplot(df, aes( x = value)) +
  geom_histogram(bins =400) +
  geom_vline(xintercept = 100, col = "lightgreen")


coords$id_tissue <- NA
coords <- as.matrix(coords[,c("X1", "Y1")])

## Modello SPDE ----------------------------------------------------------------
data$tmm_factors = tmm_factors
mesh0 <- inla.mesh.2d(loc.domain = coords, max.edge = c(20,30))
plot(mesh0)

Groups = "id_tissue"

NGroups <- length(unique(data[,Groups])) 

A2 <- inla.spde.make.A(mesh0, # Leave
                       loc = coords, # Leave
                       group = as.numeric(as.factor(data[,Groups])),
                       n.group = NGroups) 

spde <- inla.spde2.matern(mesh = mesh0, alpha = 2)

w <- inla.spde.make.index(
  name    = 'w', 
  n.spde  = spde$n.spde,
  n.group = NGroups)  


X0 <- model.matrix(~ -1 + logGFP + count.nucleus + WGA + dapi +
                     tmm_factors, 
                   data = data.nb) 
X <- as.data.frame(X0) 
N <- nrow(data)

Stack.est <- inla.stack( 
  data = list(y = data.gene[, risp]), # Leave
  A = list(1,1, 1,A2), # Change the A matrix to the new one
  effects = list(
    Intercept = rep(1, N), # Leave
    X = X, # Leave
    id_tissue = data$id_tissue,
    w = w),
  tag = "est") # CHANGE

Stack.pred <- inla.stack( 
  data = list(y = rep(NA, length(data.gene[,risp]))), # Leave
  A = list(1,1,1, A2), # Change the A matrix to the new one
  effects = list(
    Intercept = rep(1, N), # Leave
    X = X, # Leaveù
    id_tissue = data$id_tissue,
    w = w),
  tag = "pred")

stack <- inla.stack(Stack.est, Stack.pred)


# Modello BN 
f1 = y ~ -1 + Intercept + logGFP + count.nucleus + WGA + dapi +
         offset(tmm_factors)

IM1 <- inla(f1,
            family = "nbinomial",
            data = inla.stack.data(stack), # Don't forget to change the stack!
            control.compute = list(dic = TRUE, cpo = TRUE, 
                                   return.marginals.predictor=TRUE),
            control.predictor = list(A = inla.stack.A(stack), compute=TRUE, 
                                     link = 1) 
)

IM1$summary.fixed
summary(IM1)
library(INLAOutputs)
PredPValue(IM1)



IM1$marginals.fixed$logGFP %>%
  ggplot(aes(x,y)) +
  geom_line()+xlab("logit")+ylab("density")+ggtitle("Intercept")

IM1$marginals.fixed$logGFP %>% 
  as.data.frame()%>% 
  reframe(x = -x, y) %>%
  ggplot(aes(x,y)) +
  geom_line()+xlab("logit")+ylab("density")+ggtitle("Intercept")


# p-value H0: coef= 0 vs H1 coef != 0ù
    # 1 - (P(coef > 0) + P(coef < 0)) =
    # 1 - (1- P(coef <= 0) + 1- P(coef >= 0))

1-(1-inla.pmarginal(0, IM1$marginals.fixed$logGFP) +
     1-inla.pmarginal(0, IM1$marginals.fixed$logGFP %>%
                        as.data.frame() %>% 
                        reframe(x = -x, y))
   )




# Modello BN con effetto spaziale
f2 = y ~ -1 + Intercept + logGFP+ count.nucleus + WGA + dapi +
         f(w, model = spde, group = w.group, control.group = list(model = 'iid')) +
         f(id_tissue, model = "iid") + offset(tmm_factors)

IM6 <- inla(f2,
            family = "nbinomial",
            data = inla.stack.data(stack), # Don't forget to change the stack!
            control.compute = list(dic = TRUE, cpo = TRUE,
                                   return.marginals.predictor=TRUE),
            control.predictor = list(A = inla.stack.A(stack), compute=TRUE, 
                                     link = 1)
)
IM6$summary.fixed


IM6$marginals.fixed$logGFP %>%
  ggplot(aes(x,y)) +
  geom_line()+xlab("logit")+ylab("density")+ggtitle("Intercept")

IM6$marginals.fixed$logGFP %>% 
  as.data.frame()%>% 
  reframe(x = -x, y) %>%
  ggplot(aes(x,y)) +
  geom_line()+xlab("logit")+ylab("density")+ggtitle("Intercept")


# p-value H0: coef= 0 vs H1 coef != 0

1-(1-inla.pmarginal(0, IM6$marginals.fixed$logGFP) +
  1-inla.pmarginal(0, IM6$marginals.fixed$logGFP %>%
                   as.data.frame() %>% 
                   reframe(x = -x, y)
  ))


IM7 <- inla(f2,
            family = "zeroinflatednbinomial0",
            data = inla.stack.data(stack), # Don't forget to change the stack!
            control.compute = list(dic = TRUE, cpo = TRUE,
                                   return.marginals.predictor=TRUE),
            control.predictor = list(A = inla.stack.A(stack), compute=TRUE, 
                                     link = 1)
)

#names(inla.models()$likelihood)[grepl("zeroinflatednbin", names(inla.models()$likelihood))]
#inla.models()$likelihood$zeroinflatednbinomial0
 
hist(data.gene[, risp])

# residui modello
index.pred <- inla.stack.index(stack, "pred")$data
fitted.values <- IM1$summary.fitted.values[index.pred, "mean"]
residuals <- data.gene[, risp] - fitted.values

ggplot(data, aes(x = fitted.values, y = residuals)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", col = "red") +
  labs(x = "Fitted Values", y = "Residuals", title = "Residuals vs Fitted Values")

  # residui non hanno pattern particolari




  # non è tanto normale





# Guardiamo il grafico DIC
SpatialList <- list(IM1, IM6, IM7)
INLADICFig(SpatialList, ModelNames = c("Base", "SPDE_bn ", 
                                       "SPDE_bn infl zeri"))

range(data.gene[,risp])
range(fitted.values)
range(mod1$fitted.values)



summary(mod1)
mod1$deviance
IM1$summary
summary(IM1)



top
