#############################
#     ANALYSES PROSTATE     #
#############################

source("src/prostate_libraries.R")


roc <- readRDS("data/roc_sextant.rds")
roc_lobe <- readRDS("data/roc_lobe.rds")
roc_pat <- readRDS("data/roc_pat.rds")
#cst <- readRDS("data/cst.rds")
#pstat <- readRDS("data/pstat.rds")


#Description
head(pstat)
head(roc)
head(cst)

#tableau roc
#en histo, 3 patients n'ont pas de cancer, 34 patients ont un cancer, dont 21 patients sur un seul lobe et 13 patients sur les 2 lobes 
tmp <- roc %>% group_by(patient,lobeD) %>% summarise(K=max(ADK_histo)) %>% group_by(patient) %>% summarise(tot=sum(K))
table(tmp$tot)

#85 sextants avec un K en histo
table(roc$ADK_histo)

# moyennes et écart-type sur l'ensemble des sextants, interet?
mean(roc$AL_DCE0)
mean(roc$AL_DCE1)
mean(roc$RP_DCE0)
mean(roc$RP_DCE1)




#Sensibilité_Specificité

#fonctions utilisees:
# get_threshold : 
#   pour chaque variable dit si likert dépasse le seuil de 3 et seuil de 4
# compute_se_sp :
#   pour chaque variable dépassement du seuil 0/1, calcul Se et Sp par rapport à ADK_histo

#1-choix de la taille de la tumeur
#toute tailles
roc1 <- roc
roc2 <- roc_lobe
roc3 <- roc_pat
#taille >=5
roc1 <- roc[roc$taille_IRM>=5 & !is.na(roc$taille_IRM), ]
roc2 <- roc_lobe[roc_lobe$taille_IRM>=5 & !is.na(roc_lobe$taille_IRM), ]
roc3 <- roc_pat[roc_pat$taille_IRM>=5 & !is.na(roc_pat$taille_IRM), ]
#taille <10
roc1 <- roc[roc$taille_IRM<10 & !is.na(roc$taille_IRM), ]
roc2 <- roc_lobe[roc_lobe$taille_IRM<10 & !is.na(roc_lobe$taille_IRM), ]
roc3 <- roc_pat[roc_pat$taille_IRM<10 & !is.na(roc_pat$taille_IRM), ]
#taille >=10
roc1 <- roc[roc$taille_IRM>=10 & !is.na(roc$taille_IRM), ]
roc2 <- roc_lobe[roc_lobe$taille_IRM>=10 & !is.na(roc_lobe$taille_IRM), ]
roc3 <- roc_pat[roc_pat$taille_IRM>=10 & !is.na(roc_pat$taille_IRM), ]

#2-Calcul Se_sp et comparaison Se DCE0/DCE1
#sextant
sespbis <- lapply(3:4, function(x) compute_se_sp(roc1, seuil=x, unit="sextant", R=1000, type="bca")) #bca marche avec R=1000 mais pas 100
sespsext <- do.call(rbind,sespbis)
#lobe
sespbis <- lapply(3:4, function(x) compute_se_sp (roc2, seuil=x, unit="lobe", R=1000, type="bca"))
sesplobe <- do.call(rbind,sespbis)
#pat
sespbis <- lapply(3:4, function(x) compute_se_sp (roc3, seuil=x, unit="patient", R=1000, type="bca"))
sesppat <- do.call(rbind,sespbis)
write.table(print(rbind(sespsext,sesplobe,sesppat)),file="clipboard",sep="\t", row.names=F)








#kappa
#sextant

data(expsy)
ckappa(expsy[,c(12,14)])          # Cohen's kappa for binary diagnosis
#to obtain a 95%confidence interval:
#library(boot)
#ckappa.boot <- function(data,x) {ckappa(data[x,])[[2]]}
#res <- boot(expsy[,c(12,14)],ckappa.boot,1000)
#quantile(res$t,c(0.025,0.975))    # two-sided bootstrapped confidence interval of kappa
#boot.ci(res,type="bca")         # adjusted bootstrap percentile (BCa) confidence interval (better)
#ckappa(expsy[,c(11,13)])          #Cohen's kappa for non binary diagnosis


#pstat

head(pstat)
table(table(pstat$patient)) #37 patients, tous ont info pour les 2 lobes 
table(table(pstat$Lobe))
pstat$LG_pos <- ifelse(pstat$Lobe=="G" & pstat$BP_positives==1,1,0)
pstat$LD_pos <- ifelse(pstat$Lobe=="D" & pstat$BP_positives==1,1,0)
tmp <- pstat[,c("patient","LG_pos","LD_pos")]
tmp2 <- tmp %>% group_by(patient) %>% summarise(LG_pos = max(LG_pos), LD_pos = max(LD_pos))
table(tmp2$LD_pos) #12 patients ont biopsie droitpositive
table(tmp2$LG_pos) #20 patients ont biopsie gauche positive
table(tmp2$LG_pos==1 & tmp2$LD_pos==1) #5patients sont positif pour lobe gauche et lobe droit



