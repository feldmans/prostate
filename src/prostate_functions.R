get_threshold <- function(.roc){
  roc2 <- .roc 
  #ajout variable 0/1 pour seuil 3
  for (i in colnames(roc2)[str_sub(colnames(roc2),-5,-5)=="_"]){
    roc2[ ,paste0(i,"_3")] <- ifelse(roc2[ ,i] >= 3, 1, 0)
  }
  #ajout variable 0/1 pour seuil 4
  for (i in colnames(roc2)[str_sub(colnames(roc2),-5,-5)=="_"]){
    roc2[ ,paste0(i,"_4")] <- ifelse(roc2[ ,i] >= 4, 1, 0)
  }
  return(roc2)
}

compute_se_sp <- function(.roc, seuil, unit, R=1000, type="bca") {
  roc2 <- data.frame(.roc)
  #selection des colonnes "dépassement du seuil oui ou non"
  col_to_test <- colnames(roc2)[grep(seuil, colnames(roc2))] #"RP_DCE0_4" a tiret en 2e position en partant de la fin
  #pour chacune de ces colonnes, calcul se sp en comparant avec ADK_histo
  se_sp <- lapply(col_to_test, function(x){
    print(x)
    tmp <- roc2[,c("ADK_histo",x)]
    tab <- table(tmp[ ,x],tmp$ADK_histo,deparse.level = 2)
    #print(tab)
    #se : proba (signe/maladie)= proba (AL_DCE1_3=1 sachant ADK=1)
    #se <- paste0(round(tab[2,2]/sum(tab[,2]),2)*100,"%")
    #browser()
    se <- paste0(round(tab[rownames(tab)==1,colnames(tab)[colnames(tab)==1]]/sum(tab[,colnames(tab)[colnames(tab)==1]]),2)*100,"%")
    #sp : proba(pas de signe sachant pas malade = proba (AL_DCE1_3=0 sachant ADK=0))
    
    se_CI <- BootseCi(tmp, "se", R, type)
    #if(x=="RP_DCE0_3")browser()

    if(any(rownames(tab)==0) == FALSE) sp_CI <- NA 
    #if(any(colnames(tab)==0) == FALSE) sp <- "No M- with histo" 
    #if(any(colnames(tab)==0) == FALSE) browser()
    if((any(colnames(tab)==0)== FALSE & any(rownames(tab)==0)== FALSE) == TRUE) sp_CI <- NA 
    else {
      #sp <- paste0(round(tab[1,1]/sum(tab[,1]),2)*100,"%")
      sp <- paste0(round(tab[rownames(tab)==0,colnames(tab)[colnames(tab)==0]]/sum(tab[,colnames(tab)[colnames(tab)==0]]),2)*100,"%")
      sp_CI <- BootseCi(tmp, mes="sp" , R, type)
      #browser()
    }
    N <- nrow(tmp)
    #return(c(se,sp))
    #browser()
    nM0S0 <- tab[rownames(tab)==0 , colnames(tab)==0]
    nM0S1 <- tab[rownames(tab)==1 , colnames(tab)==0]
    nM1S0 <- tab[rownames(tab)==0 , colnames(tab)==1]
    nM1S1 <- tab[rownames(tab)==1 , colnames(tab)==1]
    for (i in c("nM0S0",  "nM0S1", "nM1S0", "nM1S1")){
    #lapply(c(nM0S0,  nM0S1, nM1S0, nM1S1), function(i){
      obj <- get(i)
      if (length(obj)==0) obj <- 0
      assign(i, obj)
    }
    
    #if(x=="AL_DCE0_3")browser()
    #.df <- cbind(N, se_CI, sp_CI, nM0S0, nM0S1, nM1S0, nM1S1) ; colnames(.df) <- c("N","se_CI","sp_CI", "M-S-", "M-S+", "M+S-", "M+S+") #obligatoire pour que les NA ait un colnames aussi, sinon rbind impossible
    .df <- cbind(N, se_CI, sp_CI, nM0S0, nM0S1, nM1S0, nM1S1) ; colnames(.df) <- c("N","se_CI","sp_CI", "nM0S0", "nM0S1", "nM1S0", "nM1S1") #obligatoire pour que les NA ait un colnames aussi, sinon rbind impossible
    return(.df)
  })
  #browser()
  se_sp <- do.call(rbind, se_sp)
  #chaque liste/ligne est une colonne "dépassement du seuil oui ou non" qui nous donne dans le nom de la colonne DCE, seuil et juge.
  #sesp <- data.frame(mesure = col_to_test, N = se_sp[ , 1], Se = se_sp[ , 2], Sp = se_sp[ , 3])
  sesp <- data.frame(mesure = col_to_test, se_sp)
  sesp$unit <- unit
  sesp$DCE <- str_sub(sesp$mesure, 7, 7)
  sesp$seuil <- str_sub(sesp$mesure, -1, -1)
  sesp$juge <- str_sub(sesp$mesure, 1, 2)
  
  #sesp <- sesp[ ,c(4,6,7,5,2:3)]
  sesp <- sesp[ ,c(9:12,2:8)]
  return(sesp)
}

# rocbis <- get_threshold(roc_lobe)
# sesp_lobe <- compute_se_sp (data.frame(rocbis), seuil=3, unit="lobe")
# sesp_lobe <- lapply(3:4, function(x) compute_se_sp (roc_lobe, seuil=x, unit="lobe"))
# sesp_lobe <- do.call(rbind,sesp_lobe)

#sextant
# #seuil
# roc2 <- roc 
# #ajout variable 0/1 pour seuil 3
# for (i in colnames(roc2)[5:8]){
#   roc2[ ,paste0(i,"_3")] <- ifelse(roc2[ ,i] >= 3, 1, 0)
# }
# #ajout variable 0/1 pour seuil 4
# for (i in colnames(roc2)[5:8]){
#   roc2[ ,paste0(i,"_4")] <- ifelse(roc2[ ,i] >= 4, 1, 0)
# }
# 
# #calcul se sp
# #selection des colonnes "dépassement du seuil oui ou non"
# col_to_test <- colnames(roc2)[c(grep("DCE0_",colnames(roc2)),grep("DCE1_",colnames(roc2)))]
# #pour chacune de ces colonnes, calcul se sp en comparant avec ADK_histo
# se_sp <- lapply(col_to_test, function(x){
#          tmp <- roc2[,c("ADK_histo",x)]
#          tab <- table(tmp[ ,x],tmp$ADK_histo,deparse.level = 2)
#          #se : proba (signe/maladie)= proba (AL_DCE1_3=1 sachant ADK=1)
#          se <- paste0(round(tab[2,2]/sum(tab[,2]),2)*100,"%")
#          #sp : proba(pas de signe sachant pas malade = proba (AL_DCE1_3=0 sachant ADK=0))
#          sp <- paste0(round(tab[1,1]/sum(tab[,1]),2)*100,"%")
#          return(c(se,sp))
#        })
# se_sp <- do.call(rbind, se_sp)
# #chaque liste/ligne est une colonne "dépassement du seuil oui ou non" qui nous donne dans le nom de la colonne DCE, seuil et juge.
# sesp_sextant <- data.frame(mesure = col_to_test, Se = se_sp[,1], Sp = se_sp[,2])
# sesp_sextant$DCE <- str_sub(sesp_sextant$mesure, 7, 7)
# sesp_sextant$seuil <- str_sub(sesp_sextant$mesure, -1, -1)
# sesp_sextant$juge <- str_sub(sesp_sextant$mesure, 1, 2)
# 
# sesp_sextant <- sesp_sextant[ ,c(1,4:6,2:3)]

# #lobe
# roc2 <- roc_lobe 
# #ajout variable 0/1 pour seuil 3
# for (i in colnames(roc2)[str_sub(colnames(roc2),-5,-5)=="_"]){
#   roc2[ ,paste0(i,"_3")] <- ifelse(roc2[ ,i] >= 3, 1, 0)
# }
# #ajout variable 0/1 pour seuil 4
# for (i in colnames(roc2)[str_sub(colnames(roc2),-5,-5)=="_"]){
#   roc2[ ,paste0(i,"_4")] <- ifelse(roc2[ ,i] >= 4, 1, 0)
# }
# 
# #calcul se sp
# #selection des colonnes "dépassement du seuil oui ou non"
# col_to_test <- colnames(roc2)[str_sub(colnames(roc2),-2,-2)=="_"] #"RP_DCE0_4" a tiret en 2e position en partant de la fin
# #pour chacune de ces colonnes, calcul se sp en comparant avec ADK_histo
# se_sp <- lapply(col_to_test, function(x){
#   tmp <- roc2[,c("ADK_histo",x)]
#   tab <- table(tmp[ ,x],tmp$ADK_histo,deparse.level = 2)
#   #se : proba (signe/maladie)= proba (AL_DCE1_3=1 sachant ADK=1)
#   se <- paste0(round(tab[2,2]/sum(tab[,2]),2)*100,"%")
#   #sp : proba(pas de signe sachant pas malade = proba (AL_DCE1_3=0 sachant ADK=0))
#   sp <- paste0(round(tab[1,1]/sum(tab[,1]),2)*100,"%")
#   return(c(se,sp))
# })
# se_sp <- do.call(rbind, se_sp)
# #chaque liste/ligne est une colonne "dépassement du seuil oui ou non" qui nous donne dans le nom de la colonne DCE, seuil et juge.
# sesp_lobe <- data.frame(mesure = col_to_test, Se = se_sp[,1], Sp = se_sp[,2])
# sesp_lobe$DCE <- str_sub(sesp_lobe$mesure, 7, 7)
# sesp_lobe$seuil <- str_sub(sesp_lobe$mesure, -1, -1)
# sesp_lobe$juge <- str_sub(sesp_lobe$mesure, 1, 2)
# 
# sesp_lobe <- sesp_lobe[ ,c(1,4:6,2:3)]



se.boot<- function(data, .mes, indices){ 
  .dat<- data[indices,]
  tab <- table(.dat[ ,2],.dat$ADK_histo,deparse.level = 2)
  if(.mes=="se") {
    if (any(colnames(tab)==1) == FALSE | any(rownames(tab)==1) == FALSE) mesure <- NA
    else mesure <- tab[rownames(tab)==1,colnames(tab)[colnames(tab)==1]]/sum(tab[,colnames(tab)[colnames(tab)==1]])
  }
  if(.mes=="sp") {
    if (any(colnames(tab)==0) == FALSE | any(rownames(tab)==0) == FALSE) mesure <- NA
    else mesure <- tab[rownames(tab)==0,colnames(tab)[colnames(tab)==0]]/sum(tab[,colnames(tab)[colnames(tab)==0]])
  }
  return(mesure)
}

bootse<- function (tmp, mes, R)  {    #x = "1":6 ou x="tot"
  .df <- tmp
  .res<- boot(data=.df, statistic = se.boot , .mes=mes, R=R)
  return(.res)
}



BootseCi <- function(tmp, mes, R, type)  { #type="bca" ou "perc"
  .bootres <- bootse (tmp=tmp, mes, R=R)
  if (all(na.omit(.bootres$t)==na.omit(.bootres$t)[1])) return(NA)
  .n <- length (.bootres$t0) #donne le nombre de resultat boot realise : 1 pour internet, 1 pour telephone
  .list.ci <- lapply(1:.n, function(x) boot.ci(.bootres,index=x,type=type)) #fct boot.ci : intervalle de confiance pour chaque boot
  if (type=="perc") type2 <- "percent" 
  else type2 <- type
  .res <- data.frame (t (sapply (.list.ci, function (x) x[[type2]][4:5]))) #selectionne les valeur de IC
  rownames (.res) <- names (.bootres$t0)
  #print(.res)
  #print(str(.res))
  colnames (.res) <- c ("CI_L", "CI_U")
  .res$est <- as.numeric (.bootres$t0)
  .res$n <- nrow(tmp)
  .res <- .res[, c (4,3, 1, 2)]
  .ans <- round (.res, 2) #fait un arrondi sur chaque valeur
  #if(mes=="se").ans <- data.frame (N=.res$n, Se_CI=paste0 (.ans$est*100, "% [", .ans$CI_L*100, "%-", .ans$CI_U*100, "%]")) #met en forme les valeurs
  if(mes=="se").ans <- data.frame (Se_CI=paste0 (.ans$est*100, "% [", .ans$CI_L*100, "%-", .ans$CI_U*100, "%]")) #met en forme les valeurs
  if(mes=="sp").ans <- data.frame (Sp_CI=paste0 (.ans$est*100, "% [", .ans$CI_L*100, "%-", .ans$CI_U*100, "%]")) #met en forme les valeurs
  .ans <- as.vector(.ans)
  return (.ans)
}
