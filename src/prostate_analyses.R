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
df <- data.frame(df = c("roc1", "roc2", "roc3"),unit=c("sextant","lobe","patient"), stringsAsFactors=FALSE)

for (i in 1:3){
  
  .l <- lapply(c(3:4), function(seuil){
    res <- compute_se_sp(.roc=get(df[i,"df"]), seuil, unit=df[i,"unit"], R=1000, type="bca") 
  })
  .l <- do.call(rbind, .l)
  
  sesp <- if (i==1) .l else rbind(sesp, .l)
}

sesp <- sesp %>% select(-N)
write.table(print(sesp),file="clipboard",sep="\t", row.names=F)

# #sextant
# sespbis <- lapply(3:4, function(x) compute_se_sp(roc1, seuil=x, unit="sextant", R=1000, type="bca")) #bca marche avec R=1000 mais pas 100
# sespsext <- do.call(rbind,sespbis)
# #lobe
# sespbis <- lapply(3:4, function(x) compute_se_sp (roc2, seuil=x, unit="lobe", R=1000, type="bca"))
# sesplobe <- do.call(rbind,sespbis)
# #pat
# sespbis <- lapply(3:4, function(x) compute_se_sp (roc3, seuil=x, unit="patient", R=1000, type="bca"))
# sesppat <- do.call(rbind,sespbis)
# write.table(print(rbind(sespsext,sesplobe,sesppat)),file="clipboard",sep="\t", row.names=F)






#3-Cohen's Kappa

df <- data.frame(df = c("roc1", "roc2", "roc3"),unit=c("sextant","lobe","patient"), stringsAsFactors=FALSE)

for (i in 1:3){
  .l <- lapply(c(3:4), function(seuil){
    res <- get_kappa_CI_conc (.roc=get(df[i,"df"]), seuil, unit=df[i,"unit"], R=1000, type="bca") 
  })
  .l <- do.call(rbind, .l)
  kappa_tot <- if (i==1) .l else rbind(kappa_tot, .l)
}
write.table(print(kappa_tot),file="clipboard",sep="\t", row.names=F)


# #pstat
# 
# head(pstat)
# table(table(pstat$patient)) #37 patients, tous ont info pour les 2 lobes 
# table(table(pstat$Lobe))
# pstat$LG_pos <- ifelse(pstat$Lobe=="G" & pstat$BP_positives==1,1,0)
# pstat$LD_pos <- ifelse(pstat$Lobe=="D" & pstat$BP_positives==1,1,0)
# tmp <- pstat[,c("patient","LG_pos","LD_pos")]
# tmp2 <- tmp %>% group_by(patient) %>% summarise(LG_pos = max(LG_pos), LD_pos = max(LD_pos))
# table(tmp2$LD_pos) #12 patients ont biopsie droitpositive
# table(tmp2$LG_pos) #20 patients ont biopsie gauche positive
# table(tmp2$LG_pos==1 & tmp2$LD_pos==1) #5patients sont positif pour lobe gauche et lobe droit
# 

#4- AUC estimation
var <- colnames(roc1[-c(1:4,9)])

get_rocobj <- function (data, var1, unit){
  #browser()
  tryCatch({
  print(var1)
  data$vartmp <- data[,var1]
  gm1 <- glmer(ADK_histo ~ vartmp + (1 | patient), data = data,
               family = binomial)
  # gm1 <- glmer(data$ADK_histo ~ data[ ,var1] + (1 | data$patient),
  #              family = binomial)
  p <- as.numeric(predict(gm1, type="response"))
  if (var1 %in% c("AL_DCE1", "AL_DCE0", "RP_DCE1", "RP_DCE0") & unit!="patient" ) rocobj <- roc(data$ADK_histo, p, smooth=TRUE)
  else rocobj <- roc(data$ADK_histo, p, smooth=FALSE)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}


df <- data.frame(df = c("roc1", "roc2", "roc3"),unit=c("sextant","lobe","patient"), stringsAsFactors=FALSE)

for (i in 1:3) {
  data_tmp <- get(df[i,"df"])
  unit <- df[i,"unit"]
  
  data_tmp <- get_threshold(data_tmp)
  data_tmp$patient <- as.character(data_tmp$patient)
  
  .l <- lapply(var, function(x) get_rocobj(data=data_tmp, var1=x, unit=unit))
  
  if (unit=="patient"){
    dfAUC <- data.frame(variable=var, 
                        do.call(rbind,lapply(1:12, function(x) round(as.numeric(ci.auc(.l[[x]], boot.n=2000, boot.stratified=FALSE)),3))))
  } else {
    dfAUC <- data.frame(variable=var, 
                        do.call(rbind,lapply(1:12, function(x) round(as.numeric(ci.auc(.l[[x]], method="bootstrap", boot.n=2000, boot.stratified=FALSE)),3))))
  }
  dfAUC$auc.95CI <- paste0(dfAUC$X2, "[", dfAUC$X1, "-", dfAUC$X3, "]")
  dfAUC <- dfAUC[,c(1,5)]
  dfAUC$unit <- unit
  if (i ==1 ) tab <- dfAUC
  else tab <- rbind(tab, dfAUC)
  
  pdf(paste0("data/AUC_curve",unit,".pdf"))
  par(mfrow=c(2,2))
  for (i in 1:length(var)){
    plot(.l[[i]], main = paste0(var[i],unit))
  }
  dev.off()
}

write.table(print(tab), file="clipboard", sep ="\t")


#----------------------
#Analyse complémentaire

d <- read.csv2("data/BPC_vs_BPST.csv")

d$BPST_lgmax <- as.numeric(as.character(d$BPST_lgmax))
d$BPC_lgmax <- as.numeric(as.character(d$BPC_lgmax))

table(d$BPST_positif)
table(d$BPC_positif)

round(prop.table(table(d$BPST_positif)),2)
round(prop.table(table(d$BPC_positif)),2)
fisher.test(d$BPST_positif, d$BPC_positif)

wilcox.test(d$BPST_lgmax, d$BPC_lgmax)

xBPST<-d$BPST_lgmax
xBPC<-d$BPC_lgmax

keep <- c(xBPST,xBPC)
obs0<- mean(xBPST)
obs1<- mean(xBPC)
#diff.obs <- abs(obs0 - obs1)
diff.obs <- obs0 - obs1

.gpe <- rep(c("gBPST","gBPC"),c(length(xBPST),length(xBPC)))

perm.test<- function(keepIT=keep,.gpe=.gpe){
  mixgpe <- sample(.gpe,replace = FALSE)
  g0 <- keepIT[mixgpe=="gBPST"]
  g1 <- keepIT[mixgpe=="gBPC"]
  m0<-mean(g0)
  m1<-mean(g1)
  #diff<-abs(m0-m1)
  diff<-m0-m1
  return(diff)
}

many.samp<- replicate (100000, perm.test(keep,.gpe))

p.val <-length(many.samp[abs(many.samp)>= abs(diff.obs)]) / length(many.samp)

# #pour tracer courbe
hist(many.samp,main=paste0("Difference de moyenne ",x))
abline(v=diff.obs,lwd=2,col=2)

cat(paste0("moyenne longueur carotte BPST = ", round(obs0, 2), "\nmoyenne longueur carotte BPC = ", round(obs1, 2), "\npval test de permutation = ", p.val ))
cat(sd(xBPST), sd(xBPC))  
