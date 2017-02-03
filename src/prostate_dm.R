source("src/prostate_libraries.R")

roc <- read.csv2("data/ROC_sextant.csv")
microk <- read.csv2("data/microcancer.csv")
size <- read.csv2("data/taille.csv")

#mise en forme du tableau roc à partir de roc sextant, des microK et de la taille
roc$lobeD <- ifelse(roc$sextant %in% c("AD","MD","BD"),1,0)

roc <- merge(roc, microk[ ,c("patient","sextants","microcancer")], by = c("patient","sextants"), all.x=T, all.y=F)
roc <- merge(roc, size[ ,c("patient","sextants","taille_IRM")], by = c("patient","sextants"), all.x=T, all.y=F)
roc$ADK_histo <- ifelse(roc$microcancer==1 & !is.na(roc$microcancer), 1, roc$ADK_histo)
roc$taille_IRM <- ifelse(roc$microcancer==1 & !is.na(roc$microcancer), 4, roc$taille_IRM)

roc <- roc[,c(1,8,2:7,10)]

saveRDS(roc,"data/roc_sextant.rds")

#recreation des tableaux lobe, sextant, patient à partir de roc_sextant

#tableau lobe
roc_lobe <- data.frame(roc %>% group_by(patient,lobeD) %>% select(4:9) %>% summarise_each(funs(max(.,na.rm=T))))
saveRDS(roc_lobe,"data/roc_lobe.rds")

#tableau patient
roc_pat <- data.frame(roc %>% group_by(patient) %>% select(4:9) %>% summarise_each(funs(max(.,na.rm=T))))
saveRDS(roc_pat,"data/roc_pat.rds")



#----------------------------------

#autres bases, non utilisees
pstat <- read.csv2("data/prostate_stat.csv") 
cst <- read.csv2("data/BPC_vs_BPST.csv")
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




