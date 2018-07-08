

get.div.gl <- function(data, qvalue, div.measure){
if(div.measure=="Alpha"|div.measure=="Beta"|div.measure == "Gamma"){
if(qvalue==0|qvalue==1|qvalue==2){
   require(tidyverse)
  
   Pop <- nPop(data)
   Loc <- nLoc(data)

   #prop of reference allele 
   al.prop <- function(x) (sum(x, na.rm = T))/ (2*(length(x>=0)))

   #groups by population, creates allele proportion per locus in each group, then filters to unique populations/loci  
   AlProp <- as.data.frame(cbind(Site = data$pop, as.data.frame(data))) %>%
            gather(2:(Loc +1), key = "Locus",  value = "Prop") %>%
            group_by(Site, Locus) %>% 
            mutate_all(funs(al.prop)) %>%
            distinct()
   
#Calculate Alpha Diversity (H - Entropy (Diversity), D - Numbers Equivalent Diversity)
if(div.measure == "Alpha"){
  if(qvalue == 0){
    
  AlphaValues <- AlProp %>%
                 group_by(Site) %>%
                 summarise(
                 Alpha.q0.H = ((sum((Prop>0&Prop<1), na.rm = T))*2 + #sum in this case is logical (e.g. T =1)
                              sum((Prop==0|Prop==1), na.rm = T))/(sum((Prop>-1), na.rm = T)), #removes NAs from counts
                 Alpha.q0.D = Alpha.q0.H)
  }
  
  if(qvalue == 1){  
      AlProp <- AlProp %>%
                group_by(Site) %>%
                mutate(Sh = ifelse(Prop==1|Prop==0, 0, -((Prop)*log(Prop))-((1-Prop)*log(1-Prop)))) %>%
                ungroup()
    
    
      AlphaValues <- AlProp %>%
                     group_by(Site) %>%
                     summarise(
              Alpha.q1.H = mean(Sh, na.rm =T), 
              Alpha.q1.Var = var(Sh, na.rm=T),
              Alpha.q1.D = exp(Alpha.q1.H))   #exponetial of averaged Shannon values give 1Da  
  }
  
  if(qvalue == 2){
      AlProp <- AlProp %>%
            group_by(Site) %>%
            mutate(Sh = ifelse(Prop==1|Prop==0, 0, -((Prop)*log(Prop))-((1-Prop)*log(1-Prop))),
            He = 2*Prop*(1-Prop)) %>%
            ungroup()
      
      AlphaValues <- AlProp %>%
                     group_by(Site) %>%
                     summarise(
              Alpha.q2.H = mean(He, na.rm =T),
              Alpha.q2.Var = var(He, na.rm=T),
              Alpha.q2.D = (1/(1-Alpha.q2.H)))#Calculates 2Da from averaged He values
  
  }
  
  div_data <- AlphaValues  
}

#Calculate Beta Diversity (H - Entropy (Diversity), D - Numbers Equivalent Diversity) 
if(div.measure == "Beta"){
    
  AlProp <- AlProp %>%
            group_by(Site) %>%
            mutate(Sh = ifelse(Prop==1|Prop==0, 0, -((Prop)*log(Prop))-((1-Prop)*log(1-Prop))),
            He = 2*Prop*(1-Prop)) %>%
            ungroup()
  
  #Create a list of all pairwise site combinations
  AlProp2 <- tidyr::expand(AlProp, Site, Site2 = Site)
  
  #Function for calculating beta values for each pairwise site
  #Maybe could be faster
      if(qvalue == 0){
          calcbeta <- function(Site1, Site2, data, qvalue) {
    
    data1 <-  filter(data, Site == Site1) #all loci and proportion data for Site 1, as data frame
    data2 <-  filter(data, Site == Site2) #all loci and proportion data for Site 2, as data frame
    
    Alpool <-  0.5*(data1$Prop + data2$Prop) #make vector of average proportions per loci
    

    sharedtypes <- left_join(data1, data2, by = "Locus") %>%
      mutate(sharedtypes = ifelse((is.na(Prop.x)|is.na(Prop.y)), NA, #if there are any NAs, then NA
                                  ifelse(                            #else
                                    ((Prop.x==0&Prop.y==1)|(Prop.x==1&Prop.y==0)), 0, #if opposite fixed, then 0 
                                    ifelse(((Prop.x>0&Prop.x<1)&(Prop.y>0&Prop.y<1)), 2, 
                                           1))),
             notypes = ifelse((is.na(Prop.x)|is.na(Prop.y)), NA,
                              ifelse((Prop.x==1&Prop.y==1)|(Prop.x==0&Prop.y==0), 1, 2)),
             types1 = ifelse((is.na(Prop.x)|is.na(Prop.y)), NA,
                             ifelse((Prop.x==1)|(Prop.x==0), 1, 2)),
             types2 = ifelse(is.na(Prop.y), NA,
                             ifelse((Prop.y==1)|(Prop.y==0), 1, 2))
      ) %>%
      mutate(Jac = 1-((sharedtypes)/(notypes)),
             Sor = 1-((2*sharedtypes)/(types1 + types2))) %>%
      dplyr::select(Jac, Sor)

    Sor = mean(sharedtypes$Sor, na.rm = T)
    Sorvar = var(sharedtypes$Sor, na.rm = T)
    
    Jac = mean(sharedtypes$Jac, na.rm = T)
    Jacvar = var(sharedtypes$Jac, na.rm = T)
    
    Beta.q0.D = Sor + 1

    BetaValues <- cbind(Sor, Sorvar, Jac, Jacvar, Beta.q0.D)
    #add no. of SNPs, no. of NAs
    
    
    
    return(BetaValues)
  }
  
  #Takes Site1 and Site2 from pairwise list and calculates betas using function above
  BetaOutput <- map2((AlProp2[[1]]), (AlProp2[[2]]), calcbeta, data = AlProp, qvalue = qvalue)
  
      }
      if(qvalue == 1){  
        calcbeta <- function(Site1, Site2, data, qvalue) {
    
    data1 <-  filter(data, Site == Site1) #all loci and proportion data for Site 1, as data frame
    data2 <-  filter(data, Site == Site2) #all loci and proportion data for Site 2, as data frame
    
    Alpool <-  0.5*(data1$Prop + data2$Prop) #make vector of average proportions per loci
    

    Shpool <- sapply(Alpool, function (x) ifelse(x==1|x==0, 0, -((x)*log(x))-((1-x)*log(1-x)))) #Vector of pooled Shannon
    Iloc  <- Shpool - 0.5*(data1$Sh+data2$Sh) #Pooled Shannon Vector - 0.5(Site1 Shannon + Site2 Shannon)
    IlocD <- (exp(Shpool)/ #D of pooled I
                (exp(0.5*(data1$Sh +data2$Sh)))) #divided by D of Site 1 and 2 alpha
    Iav = mean(Iloc, na.rm = T)
    Ivar = var(Iloc, na.rm = T)
    Beta.q1.D = mean((IlocD), na.rm = T)
    Beta.q1.D = 1+(log(Beta.q1.D))
    
    BetaValues <- cbind(Iav, Ivar, Beta.q1.D)
  return(BetaValues)
  }
  
  #Takes Site1 and Site2 from pairwise list and calculates betas using function above
  BetaOutput <- map2((AlProp2[[1]]), (AlProp2[[2]]), calcbeta, data = AlProp, qvalue = qvalue)
  }
      if(qvalue == 2){
        calcbeta <- function(Site1, Site2, data, qvalue) {
    
    data1 <-  filter(data, Site == Site1) #all loci and proportion data for Site 1, as data frame
    data2 <-  filter(data, Site == Site2) #all loci and proportion data for Site 2, as data frame
    
    Alpool <-  0.5*(data1$Prop + data2$Prop) #make vector of average proportions per loci

    Hepool <- sapply(Alpool, function (x) 2*x*(1-x)) #Vector of pooled He
    Heloc <- (Hepool - 0.5*(data1$He +data2$He))/Hepool #Pooled He Vector - 0.5(Site1 He + Site2 He)
    
    HelocD <- (1/(1-Hepool))/ #D of pooled He
      (1/(1-(0.5*(data1$He +data2$He)))) #divided by D of Site 1 and 2 alpha
    Heav = mean(Heloc, na.rm = T)
    Hevar = var(Heloc, na.rm = T)
    Beta.q2.D = mean(HelocD, na.rm=T)
    
    BetaValues <- cbind(Heav, Hevar, Beta.q2.D)
    return(BetaValues)
  }
  
  #Takes Site1 and Site2 from pairwise list and calculates betas using function above
  BetaOutput <- map2((AlProp2[[1]]), (AlProp2[[2]]), calcbeta, data = AlProp, qvalue = qvalue)
      }
  
  variables <- length(colnames(BetaOutput[[1]]))
  
  BetaValues  <-  cbind(AlProp2, as.data.frame(t(matrix(unlist(BetaOutput), nrow = variables))))

  colnames(BetaValues)[3:(2+variables)] <- colnames(BetaOutput[[1]])
  
  div_data <- BetaValues

}
   
#Calculate Gamma Diversity (H - Entropy (Diversity), D - Numbers Equivalent Diversity)
if(div.measure == "Gamma"){
 
  if(qvalue == 0){
    
  GammaValues <- AlProp %>%
                 ungroup() %>%
                 summarise(
                 Gamma.q0.H = ((sum((Prop>0&Prop<1), na.rm = T))*2 + #sum in this case is logical (e.g. T =1)
                              sum((Prop==0|Prop==1), na.rm = T))/(sum((Prop>-1), na.rm = T)), #removes NAs from counts
                 Gamma.q0.D = Gamma.q0.H)
  }
  
  if(qvalue == 1){  
      AlProp <- AlProp %>%
                mutate(Sh = ifelse(Prop==1|Prop==0, 0, -((Prop)*log(Prop))-((1-Prop)*log(1-Prop)))) %>%
                ungroup()
    
    
      GammaValues <- AlProp %>%
                     summarise(
              Gamma.q1.H = mean(Sh, na.rm =T), 
              Gamma.q1.Var = var(Sh, na.rm=T),
              Gamma.q1.D = exp(Gamma.q1.H))   #exponetial of averaged Shannon values give 1Da  
  }
  
  if(qvalue == 2){
      AlProp <- AlProp %>%
            mutate(Sh = ifelse(Prop==1|Prop==0, 0, -((Prop)*log(Prop))-((1-Prop)*log(1-Prop))),
            He = 2*Prop*(1-Prop)) %>%
            ungroup()
      
      GammaValues <- AlProp %>%
                     summarise(
              Gamma.q2.H = mean(He, na.rm =T),
              Gamma.q2.Var = var(He, na.rm=T),
              Gamma.q2.D = (1/(1-Gamma.q2.H)))#Calculates 2Da from averaged He values
  
  }
  
  div_data <- GammaValues  
}

  
  return(div_data)}
  else{stop('qvalue must be 0 (Allelic Richness), 1 (Shannon), or 2 (Heterozygosity)')}}
else{stop('Type of diversity (div.measure) must be chosen. Enter "Alpha" for alpha diversity (within populations),  
          "Beta" for beta diversity (between populations), or "Gamma" for total diversity (accross populations)')}
}
