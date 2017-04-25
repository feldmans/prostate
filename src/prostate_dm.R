source("src/prostate_libraries.R")


#-----------
#ancienne version

# # chargement
#roc <- read.csv2("data/ROC_sextant.csv") #premier tableau roc envoyé
#roc_new <- read.csv2("data/tableau_pour_courbr_ROC2_sextant.csv")#tableau envoyé le 16/01/2017 par AL
# microk <- read.csv2("data/microcancer.csv")
# size <- read.csv2("data/taille.csv")

# # dm
# colnames(roc_new) <- colnames(roc)
# roc <- roc_new
# rm(roc_new)
# 
# #mise en forme du tableau roc à partir de roc sextant, des microK et de la taille
# roc$lobeD <- ifelse(roc$sextant %in% c("AD","MD","BD"),1,0)
# roc <- merge(roc, microk[ ,c("patient","sextants","microcancer")], by = c("patient","sextants"), all.x=T, all.y=F)
# roc <- merge(roc, size[ ,c("patient","sextants","taille_IRM")], by = c("patient","sextants"), all.x=T, all.y=F)
# roc$ADK_histo <- ifelse(roc$microcancer==1 & !is.na(roc$microcancer), 1, roc$ADK_histo)
# roc$taille_IRM <- ifelse(roc$microcancer==1 & !is.na(roc$microcancer), 4, roc$taille_IRM)
# 
# roc <- roc[,c(1,8,2:7,10)]
# 
# saveRDS(roc,"data/roc_sextant.rds")
#----------
#version 2017 04 07

roc <- read.csv2("data/sextant_20170407.csv") 

#dm fichier "data/sextant_20170407.csv"
roc[roc$taille_IRM=="?", "taille_IRM"] <- NA
roc$sextants <- as.character(roc$sextants)
roc$taille_IRM <- as.numeric(as.character(roc$taille_IRM))
roc$microcancer[is.na(roc$microcancer)] <- 0
roc$lobeD <- ifelse(roc$sextant %in% c("AD","MD","BD"),1,0)

#imputation taille IRM de "data/sextant_20170407.csv"
roc %<>% left_join(roc %>% group_by(patient) %>% summarise(tmax = max(taille_IRM,na.rm=T))) %>%
  mutate(taille_IRM=tmax) %>% select(-tmax)
saveRDS(roc,"data/roc_sextant.rds")


#TABLEAUX lobe, sextant, patient (à partir de roc_sextant)

#tableau sextant
#deja cree

#tableau lobe
roc_lobe <- data.frame(roc %>% group_by(patient,lobeD) %>% 
                         select(ADK_histo, grep("_DCE",colnames(roc)), taille_IRM) %>%
                         summarise_each(funs(max(.,na.rm=T))))
saveRDS(roc_lobe,"data/roc_lobe.rds")

#tableau patient
roc_pat <- data.frame(roc %>% group_by(patient) %>% 
                        select(ADK_histo, grep("_DCE",colnames(roc)), taille_IRM) %>%
                               summarise_each(funs(max(.,na.rm=T))))
saveRDS(roc_pat,"data/roc_pat.rds")


#-------
#verif 
# roc_new <- read.csv2("data/tableau_pour_courbr_ROC2_sextant.csv")#tableau envoyé le 16/01/2017 par AL
# roc <- read.csv2("data/ROC_sextant.csv") #premier tableau roc envoyé
# colnames(roc_new) <- colnames(roc)
# all.equal(roc_new,roc)
# vT <- apply(roc_new==roc,1,sum)<7
# roc_new[vT,]
# roc[vT,]
# #il y a des différences, je relance avec roc_new
#----------------------------------

#autres bases, non utilisees
pstat <- read.csv2("data/prostate_stat.csv") 
#cst <- read.csv2("data/BPC_vs_BPST.csv") #utilise mais loadé directement dans l'analyse car pas de dm
roc_l <- read.csv2("data/ROC_lobe.csv")

#data_magment sur tableau inutilise finalement
roc_l[c(FALSE,TRUE), "patient"] <- roc_l[c(TRUE,FALSE), "patient"]
roc_l$lobeD <- ifelse(roc_l$Lobe=="D",1,0)
roc_l <- roc_l[ , c(1,8,3:7)]

pstat[c(FALSE,TRUE), "patient"] <- pstat[c(TRUE,FALSE), "patient"]
saveRDS(pstat, "data/pstat.rds")

cst[,c("BPST_lgmax","BPC_lgmax")] <- apply(cst[,c("BPST_lgmax","BPC_lgmax")],2,function(x)as.numeric(x))
saveRDS(cst,"data/cst.rds")


#-------------------------------

#difference entre tableau lobe Anne luzurier (roc_l) et tableau lobe refait sous r a partir de sextant
print_diff <- function(x) roc_lobe[roc_l[ , x]!= roc_lobe[ , x], c("patient","lobeD")]
names_diff <- unique(do.call(rbind, lapply(names(roc_l),print_diff)))
names_diff <- names_diff[order(names_diff$patient), ]

tab_diff <- apply(names_diff,1,function(i){
  cond <- roc_l$patient==i[1] & roc_l$lobeD==i[2]
  cbind(original="",roc_l[cond, ],recalc="",roc_lobe[cond, ])
}) 
tab_diff <- do.call(rbind,tab_diff)
write.table(print(tab_diff), file="clipboard",sep="\t",row.names=FALSE)

tab_diff_sext <- apply(names_diff,1,function(i){
  cond <- roc$patient==i[1] & roc$lobeD==i[2]
  cbind(roc[cond, ])
}) 
tab_diff_sext <- do.call(rbind,tab_diff_sext)
write.table(print(tab_diff_sext), file="clipboard",sep="\t",row.names=FALSE)




