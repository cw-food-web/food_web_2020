#Unified script for MT lake trout food web
#load libraries

library(ggplot2)
library(dplyr)
library(ellipse)
library(gridExtra)
library(rjags)
library(hdrcde)
library(reshape2)
library(tidyr)
library(cowplot)
library(lme4)
library(effects)
library(lmerTest)
library(SIBER)
library(DataCombine)
library(tidyverse)
library(Hmisc)
library(nicheROVER)
library(MixSIAR)
library(cowplot)

#########################################################################################################
#Fish summary table (SI Table 1)
#########################################################################################################
rm(list=ls())
graphics.off()
set.seed(1)

#load data
isotope_dataset <- read.csv("C:/Users/wainr/Documents/nature_08132020/final_datasets/isotope_dataset_MT_lake_trout_food_web_11072020.csv", header=TRUE)
data <- subset(isotope_dataset,Sample_type=="Fish")

#write summary function
#summary function based on: https://github.com/clquezada/tRophicPosition/blob/master/R/summariseIsotopeData.R
summariseIsotopeData <- function (df = NULL, grouping = c("Invasion_category", "Functional_group"),
                                  printSummary = FALSE, ...){
  
  d13C <- d15N <- NULL
  
  if (is.null(checkmate::checkNames(df, c("d13C_corr", "d15N_corr", grouping))))
    stop("Check the grouping variable or the names in your dataframe")
  
  summary <- plyr::ddply(df, grouping, plyr::summarise,
                         n = length(d13C_corr),
                         mean_d13C = mean(d13C_corr),
                         sd_d13C = sd(d13C_corr),
                         se_d13C = sqrt(var(d13C_corr)/length(d13C_corr)),
                         mean_d15N = mean(d15N_corr),
                         sd_d15N = sd(d15N_corr),
                         se_d15N = sqrt(var(d15N_corr)/length(d15N_corr)),
                         meanL = mean(Total_length_mm),
                         seL = sqrt(var(Total_length_mm)/length(Total_length_mm))
                         
  )
  
  if (printSummary)  print(summary)
  
  summary
  return(summary)
}

#create summary dataframe using function and grouping specified above
isotope_summary_df<-summariseIsotopeData(df = subset(data),
                                         grouping = c("Invasion_category", "Functional_group"))

#########################################################################################################
#Lake summary table (SI Table 2)
#########################################################################################################
rm(list=ls())
graphics.off()
set.seed(1)

#load data
isotope_dataset <- read.csv("C:/Users/wainr/Documents/nature_08132020/final_datasets/isotope_dataset_MT_lake_trout_food_web_11072020.csv", header=TRUE)

#write summary function
#summary function based on: https://github.com/clquezada/tRophicPosition/blob/master/R/summariseIsotopeData.R
summariseIsotopeData <- function (df = NULL, grouping = c("Lake_name"),
                                  printSummary = FALSE, ...){
  
  d13C <- d15N <- NULL
  
  if (is.null(checkmate::checkNames(df, c("d13C", "d15N", grouping))))
    stop("Check the grouping variable or the names in your dataframe")
  
  summary <- plyr::ddply(df, grouping, plyr::summarise,
                         Lake_surface_elev_m = mean(Lake_surface_elev_m),
                         Surface_area_ha = mean(Surface_area_ha),
                         Lake_max_depth_m = mean(Lake_max_depth_m),
                         Secchi_m = mean(Secchi_m),
                         Current_conversion = mean(Current_conversion)
                         
  )
  
  if (printSummary)  print(summary)
  
  summary
  return(summary)
}

#create summary dataframe using function and grouping specified above
isotope_summary_df<-summariseIsotopeData(df = subset(isotope_dataset),
                                         grouping = c("Lake_name"))

#use binomial linear regression coefficients to calculate modelled "timestep" for each lake
isotope_summary_df$prop <- isotope_summary_df$Current_conversion/(1-isotope_summary_df$Current_conversion) #logit link function
isotope_summary_df$log <- log10(isotope_summary_df$prop) #logit link function
isotope_summary_df$Time <- (isotope_summary_df$log+1.97331)/.09127 #binomial linear regression coefficients
isotope_summary_df$Model_timestep <- isotope_summary_df$Time+15 # years is "detection period" (lim yhat -> 0) calculated based on binomial linear regression

#########################################################################################################
#SI Table 2, add trophic baselines to lake table
#########################################################################################################
lookup <- data.frame(Lake_name=c("Big Salmon","Bowman","Kintla","Trout","Logging","Quartz","Hungry Horse","McDonald","Swan","Lindbergh"),
                     n_base=rep(NA,10),
                     d13C_baseline=rep(NA,10),
                     d13C_base_se=rep(NA,10),
                     d15N_baseline=rep(NA,10),
                     d15N_base_se=rep(NA,10))

isotope_dataset <- read.csv("C:/Users/wainr/Documents/nature_08132020/final_datasets/isotope_dataset_MT_lake_trout_food_web_11072020.csv", header=TRUE)
#1/10: Big salmon
data <- isotope_dataset
data <- subset(data,data$Lake_name=="Big Salmon")
data <- subset(data,data$real.name%in%c("Baetidae","Leptophlebiidae","Heptageniidae","Ephemerellidae"))
lookup$n_base[1] <- length(data$real.name)
lookup$d13C_baseline[1] <- mean(data$d13C)
lookup$d13C_base_se[1] <- sqrt(var(data$d13C)/length(data$d13C))
lookup$d15N_baseline[1] <- mean(data$d15N)
lookup$d15N_base_se[1] <- sqrt(var(data$d15N)/length(data$d15N))
#2/10: Bowman
data <- isotope_dataset
data <- subset(data,data$Lake_name=="Bowman")
data <- subset(data,data$real.name%in%c("Baetidae","Leptophlebiidae","Heptageniidae","Ephemerellidae"))
lookup$n_base[2] <- length(data$real.name)
lookup$d13C_baseline[2] <- mean(data$d13C)
lookup$d13C_base_se[2] <- sqrt(var(data$d13C)/length(data$d13C))
lookup$d15N_baseline[2] <- mean(data$d15N)
lookup$d15N_base_se[2] <- sqrt(var(data$d15N)/length(data$d15N))
#3/10: Kintla
data <- isotope_dataset
data <- subset(data,data$Lake_name=="Kintla")
data <- subset(data,data$real.name%in%c("Baetidae","Leptophlebiidae","Heptageniidae","Ephemerellidae"))
lookup$n_base[3] <- length(data$real.name)
lookup$d13C_baseline[3] <- mean(data$d13C)
lookup$d13C_base_se[3] <- sqrt(var(data$d13C)/length(data$d13C))
lookup$d15N_baseline[3] <- mean(data$d15N)
lookup$d15N_base_se[3] <- sqrt(var(data$d15N)/length(data$d15N))
#4/10: Trout
data <- isotope_dataset
data <- subset(data,data$Lake_name=="Trout")
data <- subset(data,data$real.name%in%c("Baetidae","Leptophlebiidae","Heptageniidae","Ephemerellidae"))
lookup$n_base[4] <- length(data$real.name)
lookup$d13C_baseline[4] <- mean(data$d13C)
lookup$d13C_base_se[4] <- sqrt(var(data$d13C)/length(data$d13C))
lookup$d15N_baseline[4] <- mean(data$d15N)
lookup$d15N_base_se[4] <- sqrt(var(data$d15N)/length(data$d15N))
#5/10: Logging
data <- isotope_dataset
data <- subset(data,data$Lake_name=="Logging")
data <- subset(data,data$real.name%in%c("Baetidae","Leptophlebiidae","Heptageniidae","Ephemerellidae"))
lookup$n_base[5] <- length(data$real.name)
lookup$d13C_baseline[5] <- mean(data$d13C)
lookup$d13C_base_se[5] <- sqrt(var(data$d13C)/length(data$d13C))
lookup$d15N_baseline[5] <- mean(data$d15N)
lookup$d15N_base_se[5] <- sqrt(var(data$d15N)/length(data$d15N))
#6/10: Quartz
data <- isotope_dataset
data <- subset(data,data$Lake_name=="Quartz")
data <- subset(data,data$real.name%in%c("Baetidae","Leptophlebiidae","Heptageniidae","Ephemerellidae"))
lookup$n_base[6] <- length(data$real.name)
lookup$d13C_baseline[6] <- mean(data$d13C)
lookup$d13C_base_se[6] <- sqrt(var(data$d13C)/length(data$d13C))
lookup$d15N_baseline[6] <- mean(data$d15N)
lookup$d15N_base_se[6] <- sqrt(var(data$d15N)/length(data$d15N))
#7/10: Hungry Horse
data <- isotope_dataset
data <- subset(data,data$Lake_name=="Hungry Horse")
data <- subset(data,data$real.name%in%c("Baetidae","Leptophlebiidae","Heptageniidae","Ephemerellidae"))
lookup$n_base[7] <- length(data$real.name)
lookup$d13C_baseline[7] <- mean(data$d13C)
lookup$d13C_base_se[7] <- sqrt(var(data$d13C)/length(data$d13C))
lookup$d15N_baseline[7] <- mean(data$d15N)
lookup$d15N_base_se[7] <- sqrt(var(data$d15N)/length(data$d15N))
#8/10: McDonald
data <- isotope_dataset
data <- subset(data,data$Lake_name=="McDonald")
data <- subset(data,data$real.name%in%c("Baetidae","Leptophlebiidae","Heptageniidae","Ephemerellidae"))
lookup$n_base[8] <- length(data$real.name)
lookup$d13C_baseline[8] <- mean(data$d13C)
lookup$d13C_base_se[8] <- sqrt(var(data$d13C)/length(data$d13C))
lookup$d15N_baseline[8] <- mean(data$d15N)
lookup$d15N_base_se[8] <- sqrt(var(data$d15N)/length(data$d15N))
#9/10: Swan
data <- isotope_dataset
data <- subset(data,data$Lake_name=="Swan")
data <- subset(data,data$real.name%in%c("Baetidae","Leptophlebiidae","Heptageniidae","Ephemerellidae"))
lookup$n_base[9] <- length(data$real.name)
lookup$d13C_baseline[9] <- mean(data$d13C)
lookup$d13C_base_se[9] <- sqrt(var(data$d13C)/length(data$d13C))
lookup$d15N_baseline[9] <- mean(data$d15N)
lookup$d15N_base_se[9] <- sqrt(var(data$d15N)/length(data$d15N))
#10/10: Lindbergh
data <- isotope_dataset
data <- subset(data,data$Lake_name=="Lindbergh")
data <- subset(data,data$real.name%in%c("Baetidae","Leptophlebiidae","Heptageniidae","Ephemerellidae"))
lookup$n_base[10] <- length(data$real.name)
lookup$d13C_baseline[10] <- mean(data$d13C)
lookup$d13C_base_se[10] <- sqrt(var(data$d13C)/length(data$d13C))
lookup$d15N_baseline[10] <- mean(data$d15N)
lookup$d15N_base_se[10] <- sqrt(var(data$d15N)/length(data$d15N))

#Combine summary tables.  This is SI Table 2.
SI_table_2 <- merge(isotope_summary_df,lookup,by.x="Lake_name",by.y="Lake_name")

#add "corrected" del values in dataset
#isotope_dataset2 <- isotope_dataset
#isotope_dataset2$d15N_baseline <- rep(NA,length(isotope_dataset$d15N))
#isotope_dataset2 <- merge(isotope_dataset2,lookup,by.x="Lake_name",by.y="Lake_name")
#isotope_dataset2$d13C_corr <- isotope_dataset2$d13C-isotope_dataset2$d13C_baseline
#isotope_dataset2$d15N_corr <- isotope_dataset2$d15N-isotope_dataset2$d15N_baseline
#write new csv for dataset
#write.csv(isotope_dataset2,
          #paste0("C:/Users/wainr/Documents/nature_08132020/final_datasets/",
                 #"isotope_dataset_MT_lake_trout_food_web_09102020",
                 #".csv"))

#########################################################################################################
#Figure 1.1: facet wrap plot of d13C and d15N with 95% confidence interval ellipses 
#########################################################################################################
rm(list=ls())
graphics.off()
set.seed(1)

#load data
isotope_dataset <- read.csv("C:/Users/wainr/Documents/nature_08132020/final_datasets/isotope_dataset_MT_lake_trout_food_web_11072020.csv", header=TRUE)
data <- subset(isotope_dataset,Sample_type=="Fish")

#Put variables in correct order
data$Invasion_category <- factor(data$Invasion_category, levels= c("Reference","Mid","Late"))
data$Functional_group <- factor(data$Functional_group, levels= c("Bull trout","Lake trout","Littoral forage fish","Generalist fish","Pelagic forage fish"))

#Plot Figure 1.1
ggplot(data,aes(x=d13C_corr, y=d15N_corr, shape=Functional_group))+
  facet_wrap(~Invasion_category,ncol=5)+
  geom_point(size = 3)+
  stat_ellipse(aes(x=d13C_corr, y=d15N_corr, fill=Functional_group),
               geom="polygon",
               alpha = 0.4)+
  theme_bw()+
  labs(x=expression(~delta^13*'C'~'(\211)'),
       y=expression(~delta^15*'N'~'(\211)'),
       fill = "Species",
       shape = "Species",
       group = "Species")+
  theme(axis.text = element_text(size = 20, color = "black"),
        axis.title.x = element_text(size = 20, vjust=-2),
        axis.title.y = element_text(size = 20, vjust=2),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        strip.text = element_text(size = 20))

#########################################################################################################
#Figure 1.2: SEA.b, ellipse area in SIBER 
#########################################################################################################
#https://cran.r-project.org/web/packages/SIBER/vignettes/Introduction-to-SIBER.html
#https://cran.r-project.org/web/packages/SIBER/vignettes/Ellipse-Overlap.html
rm(list=ls())
graphics.off()
set.seed(1)

#load data
isotope_dataset <- read.csv("C:/Users/wainr/Documents/nature_08132020/final_datasets/isotope_dataset_MT_lake_trout_food_web_11072020.csv", header=TRUE)
data <- subset(isotope_dataset,Sample_type=="Fish")

#prep data to meet siber specs
siber_df <- data[c("d13C_corr", "d15N_corr", "Functional_group", "Invasion_category")] #select siber's columns from full dataset
names(siber_df) <- c("iso1","iso2","group","community") #re-name columns to match siber

siber_object <- createSiberObject(siber_df) #create siber object (the input for the mcmc)
group.ML <- groupMetricsML(siber_object) #calculate TA, SEA, SEA.c

# define siber mcmc model parameters
parms <- list()
parms$n.iter <- 2 * 10^4   # number of iterations to run the model for
parms$n.burnin <- 1 * 10^3 # discard the first set of values
parms$n.thin <- 10     # thin the posterior by this many
parms$n.chains <- 2        # run this many chains

# define siber mcmc model priors
priors <- list()
priors$R <- 1 * diag(2) #uninformative
priors$k <- 2 #uninformative
priors$tau.mu <- 1.0E-3 #uninformative

#run siber mcmc
ellipses.posterior <- siberMVN(siber_object, parms, priors)

#calculate TA, SEA, SEAc for each group
group.ML <- groupMetricsML(siber_object) #ellipse area, NOT ellipse coordinates
#draw posteriors for siber ellipses for each group.community combination
SEA.B <- siberEllipses(ellipses.posterior) #motherlode
#SEA.B credible intervals:
cr.p <- c(0.95, 0.99) # vector of quantiles
#calculate credible intervals of SEA.B
SEA.B.credibles <- lapply(# call to hdrcde:hdr using lapply()
  as.data.frame(SEA.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$hdr},
  prob = cr.p) #vectors of the credible intervals for each group's sea.b

#siber has a built-in SEA.B plot: https://cran.r-project.org/web/packages/SIBER/vignettes/Introduction-to-SIBER.html
#but we'll re-do siberDensityPlot in ggplot

#draw posteriors for siber ellipses area for each group.community combination
SEA.B_2 <- SEA.B # copy siber's list
colnames(SEA.B_2) <- colnames(group.ML) #re-name columns in new list
SEA.B_2 <- as.data.frame(SEA.B_2) #make dataframe from list
SEA.B_2 <- melt(SEA.B_2) #re-arrange variables in list
SEA.B_2$combo <- SEA.B_2$variable #re-name $combo (from melt()) into $variable, as required for separate()
SEA.B_2 <- separate(data = SEA.B_2,
                    col = variable,
                    into = c("invasion.category","Species"),
                    sep = "\\.",
                    extra = "merge") #use separate() to produce two useful columns
SEA.B_2$invasion.category <- factor(SEA.B_2$invasion.category,
                                    levels=c('Reference',
                                             'Mid',
                                             'Late')) #re-order factor levels in dataframe for plotting

SEA.B_2$Species <- factor(SEA.B_2$Species,levels=c("Bull trout","Lake trout","Littoral forage fish","Generalist fish","Pelagic forage fish"))
#save.image(file='siber_outputs.RData') #must save your environment to reproduce plots

#Plot Figure 1.2
ggplot(data = SEA.B_2,
       aes(x=invasion.category, y=value, fill=invasion.category)) +
  geom_boxplot()+
  scale_fill_manual(values=c("white","light gray","grey32"))+
  theme_bw()+
  scale_y_continuous(limits=c(0,25))+
  #scale_fill_discrete(name = "Invasion")+
  labs(y = expression("Isotopic dietary breadth ( \211"^2~")"),
       x = "\nInvasion")+
  facet_wrap(~Species,
             ncol = 5)+
  theme(axis.text.x = element_text(size=20, angle = 90, color = "black", vjust=0.5),
        axis.text.y = element_text(size=20, color = "black"),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20, vjust=3),
        legend.position = "none",
        strip.text = element_text(size=20, color = "black"))

#########################################################################################################
#Figure 1.3 mixed effects linear modelling for d15N
#########################################################################################################
#Figure 1.3, d15N
rm(list=ls())
graphics.off()
set.seed(1)
isotope_dataset <- read.csv("C:/Users/wainr/Documents/nature_08132020/final_datasets/isotope_dataset_MT_lake_trout_food_web_11072020.csv", header=TRUE)
data <- subset(isotope_dataset,Sample_type=="Fish")

#d15N model first part
data$Invasion_category <- factor(data$Invasion_category, levels= c("Reference","Mid","Late"))
data <- subset(data,real.name!="Lake trout")
d15N_lmer<- lmer(d15N_corr ~ Invasion_category*Functional_group+
                   (1|Lake_name),
                 data=data)
summary(d15N_lmer)
invasion<-predictorEffect("Invasion_category", d15N_lmer)
invasion<-as.data.frame(invasion)

#d15N model second part
data <- subset(isotope_dataset,Functional_group%in%c("Lake trout","Bull trout")) #must exclude lake trout because lake trout are not present in all invasion categories
data <- subset(data,Invasion_category != "Reference")
data$Invasion_category <- factor(data$Invasion_category, levels = c("Mid","Late"))
d15N_lmer<- lmer(d15N_corr ~ Invasion_category*Functional_group+
                   (1|Lake_name),
                 data=data)
summary(d15N_lmer)#p-values
invasion2<-predictorEffect("Invasion_category", d15N_lmer)
invasion2<-as.data.frame(invasion2)
invasion2 <- subset(invasion2,Functional_group!="Bull trout")

#combine parts 1 and 2
test <- rbind(invasion,invasion2)
test$Functional_group <- factor(test$Functional_group, levels= c("Bull trout","Lake trout","Littoral forage fish","Generalist fish","Pelagic forage fish"))
test$Invasion_category <- factor(test$Invasion_category, levels= c("Reference","Mid","Late"))

#Plot figure 1.3
ggplot(test, aes (x= Invasion_category, y=fit, group=Functional_group))+
  geom_line()+
  geom_point(size = 4) +
  geom_errorbar(width=.1, aes(ymin=fit-se, ymax=fit+se)) + #you could also change this to fit+/- the se
  labs(y= expression(paste(delta^15, "N"~'(\211)')), x="\nInvasion")+
  facet_wrap(~Functional_group, ncol=5)+
  theme_bw()+
  theme(axis.text.x = element_text(size = 20, angle=90, color = "black", vjust=0.5),
        axis.text.y = element_text(size = 20, color = "black"),
        axis.title.y = element_text(size = 20, vjust = 4),
        axis.title.x = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        strip.text = element_text(size = 20))

#p-values for bull trout d15N
rm(list=ls())
graphics.off()
set.seed(1)
isotope_dataset <- read.csv("C:/Users/wainr/Documents/nature_08132020/final_datasets/isotope_dataset_MT_lake_trout_food_web_11072020.csv", header=TRUE)
data <- subset(isotope_dataset,Sample_type=="Fish")
data$Invasion_category <- factor(data$Invasion_category, levels= c("Reference","Mid","Late"))
data$Functional_group <- factor(data$Functional_group, levels= c("Bull trout","Generalist fish","Pelagic forage fish","Lake trout","Littoral forage fish"))
data <- subset(data,real.name!="Lake trout")
d15N_lmer<- lmer(d15N_corr ~ Invasion_category*Functional_group+
                   (1|Lake_name),
                 data=data)
summary(d15N_lmer) #bull trout d15N p-values (vs ref)

rm(list=ls())
graphics.off()
set.seed(1)
isotope_dataset <- read.csv("C:/Users/wainr/Documents/nature_08132020/final_datasets/isotope_dataset_MT_lake_trout_food_web_11072020.csv", header=TRUE)
data <- subset(isotope_dataset,Sample_type=="Fish")
data$Invasion_category <- factor(data$Invasion_category, levels= c("Late","Reference","Mid"))
data$Functional_group <- factor(data$Functional_group, levels= c("Bull trout","Generalist fish","Pelagic forage fish","Lake trout","Littoral forage fish"))
data <- subset(data,real.name!="Lake trout")
d15N_lmer<- lmer(d15N_corr ~ Invasion_category*Functional_group+
                   (1|Lake_name),
                 data=data)
summary(d15N_lmer) #bull trout d15N p-values (vs late)

#p-values for littoral forage fish d15N
rm(list=ls())
graphics.off()
set.seed(1)
isotope_dataset <- read.csv("C:/Users/wainr/Documents/nature_08132020/final_datasets/isotope_dataset_MT_lake_trout_food_web_11072020.csv", header=TRUE)
data <- subset(isotope_dataset,Sample_type=="Fish")
data$Invasion_category <- factor(data$Invasion_category, levels= c("Reference","Mid","Late"))
data$Functional_group <- factor(data$Functional_group, levels= c("Littoral forage fish","Generalist fish","Pelagic forage fish","Bull trout","Lake trout"))
data <- subset(data,real.name!="Lake trout")
d15N_lmer<- lmer(d15N_corr ~ Invasion_category*Functional_group+
                   (1|Lake_name),
                 data=data)
summary(d15N_lmer) #littoral forage fish d15N p-values (vs ref)

rm(list=ls())
graphics.off()
set.seed(1)
isotope_dataset <- read.csv("C:/Users/wainr/Documents/nature_08132020/final_datasets/isotope_dataset_MT_lake_trout_food_web_11072020.csv", header=TRUE)
data <- subset(isotope_dataset,Sample_type=="Fish")
data$Invasion_category <- factor(data$Invasion_category, levels= c("Late","Reference","Mid"))
data$Functional_group <- factor(data$Functional_group, levels= c("Littoral forage fish","Generalist fish","Pelagic forage fish","Bull trout","Lake trout"))
data <- subset(data,real.name!="Lake trout")
d15N_lmer<- lmer(d15N_corr ~ Invasion_category*Functional_group+
                   (1|Lake_name),
                 data=data)
summary(d15N_lmer) #Littoral forage fish d15N p-values (vs late)

#p-values for generalist fish d15N
rm(list=ls())
graphics.off()
set.seed(1)
isotope_dataset <- read.csv("C:/Users/wainr/Documents/nature_08132020/final_datasets/isotope_dataset_MT_lake_trout_food_web_11072020.csv", header=TRUE)
data <- subset(isotope_dataset,Sample_type=="Fish")
data$Invasion_category <- factor(data$Invasion_category, levels= c("Reference","Mid","Late"))
data$Functional_group <- factor(data$Functional_group, levels= c("Generalist fish","Pelagic forage fish","Bull trout","Lake trout","Littoral forage fish"))
data <- subset(data,real.name!="Lake trout")
d15N_lmer<- lmer(d15N_corr ~ Invasion_category*Functional_group+
                   (1|Lake_name),
                 data=data)
summary(d15N_lmer) #generalist fish d15N p-values (vs ref)

rm(list=ls())
graphics.off()
set.seed(1)
isotope_dataset <- read.csv("C:/Users/wainr/Documents/nature_08132020/final_datasets/isotope_dataset_MT_lake_trout_food_web_11072020.csv", header=TRUE)
data <- subset(isotope_dataset,Sample_type=="Fish")
data$Invasion_category <- factor(data$Invasion_category, levels= c("Late","Reference","Mid"))
data$Functional_group <- factor(data$Functional_group, levels= c("Generalist fish","Pelagic forage fish","Bull trout","Lake trout","Littoral forage fish"))
data <- subset(data,real.name!="Lake trout")
d15N_lmer<- lmer(d15N_corr ~ Invasion_category*Functional_group+
                   (1|Lake_name),
                 data=data)
summary(d15N_lmer) #generalist fish d15N p-values (vs late)

#p-values for Pelagic forage fish d15N
rm(list=ls())
graphics.off()
set.seed(1)
isotope_dataset <- read.csv("C:/Users/wainr/Documents/nature_08132020/final_datasets/isotope_dataset_MT_lake_trout_food_web_11072020.csv", header=TRUE)
data <- subset(isotope_dataset,Sample_type=="Fish")
data$Invasion_category <- factor(data$Invasion_category, levels= c("Reference","Mid","Late"))
data$Functional_group <- factor(data$Functional_group, levels= c("Pelagic forage fish","Bull trout","Lake trout","Littoral forage fish","Generalist fish"))
data <- subset(data,real.name!="Lake trout")
d15N_lmer<- lmer(d15N_corr ~ Invasion_category*Functional_group+
                   (1|Lake_name),
                 data=data)
summary(d15N_lmer) #Pelagic forage fish d15N p-values (vs ref)

rm(list=ls())
graphics.off()
set.seed(1)
isotope_dataset <- read.csv("C:/Users/wainr/Documents/nature_08132020/final_datasets/isotope_dataset_MT_lake_trout_food_web_11072020.csv", header=TRUE)
data <- subset(isotope_dataset,Sample_type=="Fish")
data$Invasion_category <- factor(data$Invasion_category, levels= c("Late","Reference","Mid"))
data$Functional_group <- factor(data$Functional_group, levels= c("Pelagic forage fish","Littoral forage fish","Generalist fish","Bull trout","Lake trout"))
data <- subset(data,real.name!="Lake trout")
d15N_lmer<- lmer(d15N_corr ~ Invasion_category*Functional_group+
                   (1|Lake_name),
                 data=data)
summary(d15N_lmer) #Pelagic forage fish d15N p-values (vs late)

#p-values for Lake trout d15N
rm(list=ls())
graphics.off()
set.seed(1)
isotope_dataset <- read.csv("C:/Users/wainr/Documents/nature_08132020/final_datasets/isotope_dataset_MT_lake_trout_food_web_11072020.csv", header=TRUE)
data <- subset(isotope_dataset,Sample_type=="Fish")
data <- subset(isotope_dataset,Functional_group%in%c("Lake trout","Bull trout")) #must exclude lake trout because lake trout are not present in all invasion categories
data <- subset(data,Invasion_category != "Reference")
data$Invasion_category <- factor(data$Invasion_category, levels = c("Mid","Late"))
data$Functional_group <- factor(data$Functional_group, levels= c("Lake trout","Bull trout"))
d15N_lmer<- lmer(d15N_corr ~ Invasion_category*Functional_group+
                   (1|Lake_name),
                 data=data)
summary(d15N_lmer)#lake trout d15N p-values (vs ref)

rm(list=ls())
graphics.off()
set.seed(1)
isotope_dataset <- read.csv("C:/Users/wainr/Documents/nature_08132020/final_datasets/isotope_dataset_MT_lake_trout_food_web_11072020.csv", header=TRUE)
data <- subset(isotope_dataset,Sample_type=="Fish")
data <- subset(isotope_dataset,Functional_group%in%c("Lake trout","Bull trout")) #must exclude lake trout because lake trout are not present in all invasion categories
data <- subset(data,Invasion_category != "Reference")
data$Invasion_category <- factor(data$Invasion_category, levels = c("Late","Mid"))
data$Functional_group <- factor(data$Functional_group, levels= c("Lake trout","Bull trout"))
d15N_lmer<- lmer(d15N_corr ~ Invasion_category*Functional_group+
                   (1|Lake_name),
                 data=data)
summary(d15N_lmer)#lake trout d15N p-values (vs late)

#########################################################################################################
#Figure 1.4 mixed effects linear modelling for d13C
#########################################################################################################
#Figure 1.4, d13C
rm(list=ls())
graphics.off()
set.seed(1)
isotope_dataset <- read.csv("C:/Users/wainr/Documents/nature_08132020/final_datasets/isotope_dataset_MT_lake_trout_food_web_11072020.csv", header=TRUE)
data <- subset(isotope_dataset,Sample_type=="Fish")

#d13C model first part
data$Invasion_category <- factor(data$Invasion_category, levels= c("Reference","Mid","Late"))
data <- subset(data,real.name!="Lake trout")
d13C_lmer<- lmer(d13C_corr ~ Invasion_category*Functional_group+
                   (1|Lake_name),
                 data=data)
summary(d13C_lmer)
invasion<-predictorEffect("Invasion_category", d13C_lmer)
invasion<-as.data.frame(invasion)

#d13C model second part
data <- subset(isotope_dataset,Functional_group%in%c("Lake trout","Bull trout")) #must exclude lake trout because lake trout are not present in all invasion categories
data <- subset(data,Invasion_category != "Reference")
data$Invasion_category <- factor(data$Invasion_category, levels = c("Mid","Late"))
d13C_lmer<- lmer(d13C_corr ~ Invasion_category*Functional_group+
                   (1|Lake_name),
                 data=data)
summary(d13C_lmer)#p-values
invasion2<-predictorEffect("Invasion_category", d13C_lmer)
invasion2<-as.data.frame(invasion2)
invasion2 <- subset(invasion2,Functional_group!="Bull trout")

#combine parts 1 and 2
test <- rbind(invasion,invasion2)
test$Functional_group <- factor(test$Functional_group, levels= c("Bull trout","Lake trout","Littoral forage fish","Generalist fish","Pelagic forage fish"))
test$Invasion_category <- factor(test$Invasion_category, levels= c("Reference","Mid","Late"))

#Plot figure 1.4
ggplot(test, aes (x= Invasion_category, y=fit, group=Functional_group))+
  geom_line()+
  geom_point(size = 4) +
  geom_errorbar(width=.1, aes(ymin=fit-se, ymax=fit+se)) + #you could also change this to fit+/- the se
  labs(y= expression(paste(delta^13, "C"~'(\211)')), x="\nInvasion")+
  facet_wrap(~Functional_group, ncol=5)+
  theme_bw()+
  theme(axis.text.x = element_text(size = 20, angle=90, color = "black", vjust=0.5),
        axis.text.y = element_text(size = 20, color = "black"),
        axis.title.y = element_text(size = 20, vjust = 4),
        axis.title.x = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        strip.text = element_text(size = 20))

#p-values for bull trout
rm(list=ls())
graphics.off()
set.seed(1)
isotope_dataset <- read.csv("C:/Users/wainr/Documents/nature_08132020/final_datasets/isotope_dataset_MT_lake_trout_food_web_11072020.csv", header=TRUE)
data <- subset(isotope_dataset,Sample_type=="Fish")
data$Invasion_category <- factor(data$Invasion_category, levels= c("Reference","Mid","Late"))
data$Functional_group <- factor(data$Functional_group, levels= c("Bull trout","Generalist fish","Pelagic forage fish","Lake trout","Littoral forage fish"))
data <- subset(data,real.name!="Lake trout")
d13C_lmer<- lmer(d13C_corr ~ Invasion_category*Functional_group+
                   (1|Lake_name),
                 data=data)
summary(d13C_lmer) #bull trout d13C p-values (vs ref)

rm(list=ls())
graphics.off()
set.seed(1)
isotope_dataset <- read.csv("C:/Users/wainr/Documents/nature_08132020/final_datasets/isotope_dataset_MT_lake_trout_food_web_11072020.csv", header=TRUE)
data <- subset(isotope_dataset,Sample_type=="Fish")
data$Invasion_category <- factor(data$Invasion_category, levels= c("Late","Reference","Mid"))
data$Functional_group <- factor(data$Functional_group, levels= c("Bull trout","Generalist fish","Pelagic forage fish","Lake trout","Littoral forage fish"))
data <- subset(data,real.name!="Lake trout")
d13C_lmer<- lmer(d13C_corr ~ Invasion_category*Functional_group+
                   (1|Lake_name),
                 data=data)
summary(d13C_lmer) #bull trout d13C p-values (vs late)

#p-values for littoral forage fish
rm(list=ls())
graphics.off()
set.seed(1)
isotope_dataset <- read.csv("C:/Users/wainr/Documents/nature_08132020/final_datasets/isotope_dataset_MT_lake_trout_food_web_11072020.csv", header=TRUE)
data <- subset(isotope_dataset,Sample_type=="Fish")
data$Invasion_category <- factor(data$Invasion_category, levels= c("Reference","Mid","Late"))
data$Functional_group <- factor(data$Functional_group, levels= c("Littoral forage fish","Generalist fish","Pelagic forage fish","Bull trout","Lake trout"))
data <- subset(data,real.name!="Lake trout")
d13C_lmer<- lmer(d13C_corr ~ Invasion_category*Functional_group+
                   (1|Lake_name),
                 data=data)
summary(d13C_lmer) #littoral forage fish d13C p-values (vs ref)

rm(list=ls())
graphics.off()
set.seed(1)
isotope_dataset <- read.csv("C:/Users/wainr/Documents/nature_08132020/final_datasets/isotope_dataset_MT_lake_trout_food_web_11072020.csv", header=TRUE)
data <- subset(isotope_dataset,Sample_type=="Fish")
data$Invasion_category <- factor(data$Invasion_category, levels= c("Late","Reference","Mid"))
data$Functional_group <- factor(data$Functional_group, levels= c("Littoral forage fish","Generalist fish","Pelagic forage fish","Bull trout","Lake trout"))
data <- subset(data,real.name!="Lake trout")
d13C_lmer<- lmer(d13C_corr ~ Invasion_category*Functional_group+
                   (1|Lake_name),
                 data=data)
summary(d13C_lmer) #littoral forage fish d13C p-values (vs late)

#p-values for generalist fish
rm(list=ls())
graphics.off()
set.seed(1)
isotope_dataset <- read.csv("C:/Users/wainr/Documents/nature_08132020/final_datasets/isotope_dataset_MT_lake_trout_food_web_11072020.csv", header=TRUE)
data <- subset(isotope_dataset,Sample_type=="Fish")
data$Invasion_category <- factor(data$Invasion_category, levels= c("Reference","Mid","Late"))
data$Functional_group <- factor(data$Functional_group, levels= c("Generalist fish","Pelagic forage fish","Bull trout","Lake trout","Littoral forage fish"))
data <- subset(data,real.name!="Lake trout")
d13C_lmer<- lmer(d13C_corr ~ Invasion_category*Functional_group+
                   (1|Lake_name),
                 data=data)
summary(d13C_lmer) #generalist fish d13C p-values (vs ref)

rm(list=ls())
graphics.off()
set.seed(1)
isotope_dataset <- read.csv("C:/Users/wainr/Documents/nature_08132020/final_datasets/isotope_dataset_MT_lake_trout_food_web_11072020.csv", header=TRUE)
data <- subset(isotope_dataset,Sample_type=="Fish")
data$Invasion_category <- factor(data$Invasion_category, levels= c("Late","Reference","Mid"))
data$Functional_group <- factor(data$Functional_group, levels= c("Generalist fish","Pelagic forage fish","Bull trout","Lake trout","Littoral forage fish"))
data <- subset(data,real.name!="Lake trout")
d13C_lmer<- lmer(d13C_corr ~ Invasion_category*Functional_group+
                   (1|Lake_name),
                 data=data)
summary(d13C_lmer) #generalist fish d13C p-values (vs late)

#p-values for Pelagic forage fish
rm(list=ls())
graphics.off()
set.seed(1)
isotope_dataset <- read.csv("C:/Users/wainr/Documents/nature_08132020/final_datasets/isotope_dataset_MT_lake_trout_food_web_11072020.csv", header=TRUE)
data <- subset(isotope_dataset,Sample_type=="Fish")
data$Invasion_category <- factor(data$Invasion_category, levels= c("Reference","Mid","Late"))
data$Functional_group <- factor(data$Functional_group, levels= c("Pelagic forage fish","Bull trout","Lake trout","Littoral forage fish","Generalist fish"))
data <- subset(data,real.name!="Lake trout")
d13C_lmer<- lmer(d13C_corr ~ Invasion_category*Functional_group+
                   (1|Lake_name),
                 data=data)
summary(d13C_lmer) #Pelagic forage fish d13C p-values (vs ref)

rm(list=ls())
graphics.off()
set.seed(1)
isotope_dataset <- read.csv("C:/Users/wainr/Documents/nature_08132020/final_datasets/isotope_dataset_MT_lake_trout_food_web_11072020.csv", header=TRUE)
data <- subset(isotope_dataset,Sample_type=="Fish")
data$Invasion_category <- factor(data$Invasion_category, levels= c("Late","Reference","Mid"))
data$Functional_group <- factor(data$Functional_group, levels= c("Pelagic forage fish","Bull trout","Lake trout","Littoral forage fish","Generalist fish"))
data <- subset(data,real.name!="Lake trout")
d13C_lmer<- lmer(d13C_corr ~ Invasion_category*Functional_group+
                   (1|Lake_name),
                 data=data)
summary(d13C_lmer) #Pelagic forage fish d13C p-values (vs late)

#p-values for Lake trout d13C
rm(list=ls())
graphics.off()
set.seed(1)
isotope_dataset <- read.csv("C:/Users/wainr/Documents/nature_08132020/final_datasets/isotope_dataset_MT_lake_trout_food_web_11072020.csv", header=TRUE)
data <- subset(isotope_dataset,Sample_type=="Fish")
data <- subset(isotope_dataset,Functional_group%in%c("Lake trout","Bull trout")) #must exclude lake trout because lake trout are not present in all invasion categories
data <- subset(data,Invasion_category != "Reference")
data$Invasion_category <- factor(data$Invasion_category, levels = c("Mid","Late"))
data$Functional_group <- factor(data$Functional_group, levels= c("Lake trout","Bull trout"))
d13C_lmer<- lmer(d13C_corr ~ Invasion_category*Functional_group+
                   (1|Lake_name),
                 data=data)
summary(d13C_lmer)#lake trout d13C p-values (vs mid)


#########################################################################################################
#Figure 1.5: niche overlap 
#########################################################################################################
rm(list=ls())
graphics.off()
set.seed(1)

#load data
isotope_dataset <- read.csv("C:/Users/wainr/Documents/nature_08132020/final_datasets/isotope_dataset_MT_lake_trout_food_web_11072020.csv", header=TRUE)
data <- subset(isotope_dataset,Sample_type=="Fish")

#STEP 1: RUN THIS SECTION FOR Invasion_category == Late

#prep data
nicheROVER_df <- subset(data, Invasion_category %in% "Late" & real.name %in% c("Bull trout","Lake trout")) #select Late blt and lkt
#eliminate all extra columns (nicheROVER requires only 3 columns as inputs)
nicheROVER_df <- nicheROVER_df[c("Functional_group","d13C_corr", "d15N_corr")] #columns must be in this order (species, d13C, d15N)
names(nicheROVER_df) <- c("species","D13C","D15N") #re-name columns to match nicherover

#aggregate data as required for nicheROVER
aggregate(nicheROVER_df[2:3], nicheROVER_df[1], mean) #point estimates of mean isotope values

#create fish.par object
nsamples <- 1000
system.time({
  fish.par <- tapply(1:nrow(nicheROVER_df), nicheROVER_df$species, function(ii) niw.post(nsamples = nsamples, 
                                                                                         X = nicheROVER_df[ii, 2:3]))
})
#remove NULLs from fish.par list
#fish.par should be a list of length = number of species (or whatever factor you're aggregating by)
#https://stackoverflow.com/questions/33004238/r-removing-null-elements-from-a-list
fish.par <- fish.par[-which(sapply(fish.par, is.null))] #get rid of NULLs from our fish.par list


#estimate niche overlaps
nsamples <- 1000

# Overlap calculation.  use nsamples = nprob = 10000 (1e4) for higher
# accuracy.  the variable over.stat can be supplied directly to the
# overlap.plot function

over.stat <- overlap(fish.par, nreps = nsamples, nprob = 1000, alpha = c(0.95,.99))

# The mean overlap metrics calculated across iteratations for both niche
# region sizes (alpha = .95 and alpha = .99) can be calculated and displayed
# in an array.
over.mean <- apply(over.stat,
                   c(1:2, 4),
                   mean) * 100
round(over.mean, 2)

over.cred <- apply(over.stat * 100,
                   c(1:2, 4),
                   quantile,
                   prob = c(0.025, 0.975), 
                   na.rm = TRUE)
round(over.cred[, , , 1])  # display alpha = .95 niche region

#generate array of overlap probabilities for species A vs species B using over.mean and and over.cred
over.stat <- overlap(fish.par, nreps = nsamples, nprob = 1000, alpha = 0.95)

#re-do the package's plotting function to prep data for ggplot
over.hist <- apply(over.stat, 1:2, function(x) {
  if (any(is.na(x))) 
    return(NULL)
  else {
    tmp <- x * 100
  }
  tmp
})


#How much do bull trout overlap with lake trout?
data1<-as.data.frame(over.hist[[1,2]])
colnames(data1)<-"Overlap"

#Add column to data1, label rows Late
BLT_Late <- data1
BLT_Late$Species <- "Bull trout"
BLT_Late$Invasion <- "Late"

#lake trout ggplot (Figure 1.5B)
#How much do lake trout overlap bull trout?
data1<-as.data.frame(over.hist[[2,1]])
colnames(data1)<-"Overlap"

#Add column to data1, label rows Late
LKT_Late <- data1
LKT_Late$Species <- "Lake trout"
LKT_Late$Invasion <- "Late"

#STEP 2: RUN THIS SECTION FOR Invasion_category == Mid

#prep data
nicheROVER_df <- subset(data, Invasion_category %in% "Mid" & real.name %in% c("Bull trout","Lake trout")) #select Mid blt and lkt #select blt and lkt
#eliminate all extra columns (nicheROVER requires only 3 columns as inputs)
nicheROVER_df <- nicheROVER_df[c("Functional_group","d13C_corr", "d15N_corr")] #columns must be in this order (species, d13C, d15N)
names(nicheROVER_df) <- c("species","D13C","D15N") #re-name columns to match nicherover

#aggregate data as required for nicheROVER
aggregate(nicheROVER_df[2:3], nicheROVER_df[1], mean) #point estimates of mean isotope values

#create fish.par object
nsamples <- 1000
system.time({
  fish.par <- tapply(1:nrow(nicheROVER_df), nicheROVER_df$species, function(ii) niw.post(nsamples = nsamples, 
                                                                                         X = nicheROVER_df[ii, 2:3]))
})
#remove NULLs from fish.par list
#fish.par should be a list of length = number of species (or whatever factor you're aggregating by)
#https://stackoverflow.com/questions/33004238/r-removing-null-elements-from-a-list
fish.par <- fish.par[-which(sapply(fish.par, is.null))] #get rid of NULLs from our fish.par list


#estimate niche overlaps
nsamples <- 1000

# Overlap calculation.  use nsamples = nprob = 10000 (1e4) for higher
# accuracy.  the variable over.stat can be supplied directly to the
# overlap.plot function

over.stat <- overlap(fish.par, nreps = nsamples, nprob = 1000, alpha = c(0.95,.99))

# The mean overlap metrics calculated across iteratations for both niche
# region sizes (alpha = .95 and alpha = .99) can be calculated and displayed
# in an array.
over.mean <- apply(over.stat,
                   c(1:2, 4),
                   mean) * 100
round(over.mean, 2)

over.cred <- apply(over.stat * 100,
                   c(1:2, 4),
                   quantile,
                   prob = c(0.025, 0.975), 
                   na.rm = TRUE)
round(over.cred[, , , 1])  # display alpha = .95 niche region

#generate array of overlap probabilities for species A vs species B using over.mean and and over.cred
over.stat <- overlap(fish.par, nreps = nsamples, nprob = 1000, alpha = 0.95)

#re-do the package's plotting function to prep data for ggplot
over.hist <- apply(over.stat, 1:2, function(x) {
  if (any(is.na(x))) 
    return(NULL)
  else {
    tmp <- x * 100
  }
  tmp
})


#How much do bull trout overlap with lake trout?
data1<-as.data.frame(over.hist[[1,2]])
colnames(data1)<-"Overlap"

#Add column to data1, label rows Late
BLT_Mid <- data1
BLT_Mid$Species <- "Bull trout"
BLT_Mid$Invasion <- "Mid"

#lake trout ggplot (Figure 1.5B)
#How much do lake trout overlap bull trout?
data1<-as.data.frame(over.hist[[2,1]])
colnames(data1)<-"Overlap"

#Add column to data1, label rows Late
LKT_Mid <- data1
LKT_Mid$Species <- "Lake trout"
LKT_Mid$Invasion <- "Mid"

BLT_both <- rbind(BLT_Mid,BLT_Late) ##combine BLT rows
LKT_both <- rbind(LKT_Mid,LKT_Late) #combine LKT rows

Mid_both <- rbind(BLT_Mid,LKT_Mid) #combine mid rows
#Mid_both$Overlap2 <- Mid_both$Overlap/100 #probability as proportion
Late_both <- rbind(BLT_Late,LKT_Late) #combine late rows
#Late_both$Overlap2 <- Late_both$Overlap/100 #probability as proportion

#BLT_plot
#Mid lines should be red, Late lines should be blue
#BLT_both$Invasion <- factor(BLT_both$Invasion, levels= c("Mid","Late"))
BLT_both$Invasion <- factor(BLT_both$Invasion, levels=c("Mid","Late"))
ggplot(data=BLT_both, aes(x=Overlap,y=..scaled.., fill=Invasion, color = Invasion))+ 
  geom_density(alpha=0.5)+
  scale_fill_manual(values=c("gray","black"))+
  scale_color_manual(values=c("gray","black"))+
  #scale_x_continuous(limits=c(0,1),breaks=c(0,0.5,1))+ #for proportional overlap instead of percent
  scale_x_continuous(limits=c(0,100),breaks=c(0,50,100))+ #for percent instead of proportional overlap
  scale_y_continuous(limits=c(0,1),breaks=c(0,0.5,1))+
  labs(title="A. Bull trout overlap with lake trout",
       x="\nProbability of Overlap (%)",
       y="Scaled Posterior Density\n",
       fill="Invasion")+
  theme_bw()+
  #geom_vline(xintercept=mean(subset(Mid_both,Mid_both$Species=="Lake trout")$Overlap), size=1, color="blue")+ #LKT
  #geom_vline(xintercept = quantile(subset(Mid_both,Mid_both$Species=="Lake trout")$Overlap, prob = 0.025), linetype="dotted", size=2, color="blue")+ #LKT
  #geom_vline(xintercept = quantile(subset(Mid_both,Mid_both$Species=="Lake trout")$Overlap, prob = 0.975), linetype="dotted", size=2, color="blue")+ #LKT
  #geom_vline(xintercept = quantile(subset(Mid_both,Mid_both$Species=="Lake trout")$Overlap2, prob = 0.5), linetype="dotted", size=2, color="blue")+ #LKT
  geom_vline(xintercept = quantile(subset(BLT_both,BLT_both$Invasion=="Late")$Overlap, prob = 0.5), linetype="solid", size=2)+ #Late
  #geom_vline(xintercept=mean(subset(Mid_both,Mid_both$Species=="Bull trout")$Overlap), size=1, color="red")+ #BLT
  #geom_vline(xintercept = quantile(subset(Mid_both,Mid_both$Species=="Bull trout")$Overlap, prob = 0.025), linetype="dotted", size=2, color="red")+ #BLT
  #geom_vline(xintercept = quantile(subset(Mid_both,Mid_both$Species=="Bull trout")$Overlap, prob = 0.975), linetype="dotted", size=2, color="red")+ #BLT
  #geom_vline(xintercept = quantile(subset(Mid_both,Mid_both$Species=="Bull trout")$Overlap2, prob = 0.5), linetype="dotted", size=2, color="red")+ #BLT
  geom_vline(xintercept = quantile(subset(BLT_both,BLT_both$Invasion=="Mid")$Overlap, prob = 0.5), linetype="dotted", size=2)+ #Mid
  theme(#axis.text = element_blank(),
    #axis.title = element_blank()
    #legend.position = "none",
    axis.text = element_text(size = 30, color = "black"),
    axis.title = element_text(size = 30),
    #axis.text.y = element_blank(),
    #axis.title.y = element_blank(),
    title = element_text(size=30),
    legend.text = element_text(size=30))

#LKT_plot
#Mid lines should be blue, Late lines should be red
LKT_both$Invasion <- factor(LKT_both$Invasion, levels=c("Mid","Late"))
ggplot(data=LKT_both, aes(x=Overlap,y=..scaled.., fill=Invasion, color = Invasion))+ 
  geom_density(alpha=0.5)+
  scale_fill_manual(values=c("gray","black"))+
  scale_color_manual(values=c("gray","black"))+
  #scale_x_continuous(limits=c(0,1),breaks=c(0,0.5,1))+ #for proportional overlap instead of percent
  scale_x_continuous(limits=c(0,100),breaks=c(0,50,100))+ #for percent instead of proportional overlap
  scale_y_continuous(limits=c(0,1),breaks=c(0,0.5,1))+
  labs(title="B. Lake trout overlap with bull trout",
       x="\nProbability of Overlap (%)",
       y="Scaled Posterior Density\n",
       fill="Invasion")+
  theme_bw()+
  #geom_vline(xintercept=mean(subset(Mid_both,Mid_both$Species=="Lake trout")$Overlap), size=1, color="blue")+ #LKT
  #geom_vline(xintercept = quantile(subset(Mid_both,Mid_both$Species=="Lake trout")$Overlap, prob = 0.025), linetype="dotted", size=2, color="blue")+ #LKT
  #geom_vline(xintercept = quantile(subset(Mid_both,Mid_both$Species=="Lake trout")$Overlap, prob = 0.975), linetype="dotted", size=2, color="blue")+ #LKT
  #geom_vline(xintercept = quantile(subset(Mid_both,Mid_both$Species=="Lake trout")$Overlap2, prob = 0.5), linetype="dotted", size=2, color="blue")+ #LKT
  geom_vline(xintercept = quantile(subset(LKT_both,LKT_both$Invasion=="Late")$Overlap, prob = 0.5), linetype="solid", size=2)+ #Late
  #geom_vline(xintercept=mean(subset(Mid_both,Mid_both$Species=="Bull trout")$Overlap), size=1, color="red")+ #BLT
  #geom_vline(xintercept = quantile(subset(Mid_both,Mid_both$Species=="Bull trout")$Overlap, prob = 0.025), linetype="dotted", size=2, color="red")+ #BLT
  #geom_vline(xintercept = quantile(subset(Mid_both,Mid_both$Species=="Bull trout")$Overlap, prob = 0.975), linetype="dotted", size=2, color="red")+ #BLT
  #geom_vline(xintercept = quantile(subset(Mid_both,Mid_both$Species=="Bull trout")$Overlap2, prob = 0.5), linetype="dotted", size=2, color="red")+ #BLT
  geom_vline(xintercept = quantile(subset(LKT_both,LKT_both$Invasion=="Mid")$Overlap, prob = 0.5), linetype="dotted", size=2)+ #Mid
  theme(#axis.text = element_blank(),
    #axis.title = element_blank()
    #legend.position = "none",
    axis.text.x = element_text(size = 30, color="black"),
    axis.title = element_text(size = 30),
    axis.text.y = element_text(size = 30, color = "black"),
    axis.title.y = element_text(size = 30, color = "black"),
    title = element_text(size=30),
    legend.text = element_text(size=30))

#Produce SI Table 4
ls_SI_4 <- vector(mode = "list")
#BLT Mid
ls_SI_4$BLT_on_LKT$Mid$CrI_2.5 <- quantile(subset(Mid_both,Mid_both$Species=="Bull trout")$Overlap, prob = 0.025)
ls_SI_4$BLT_on_LKT$Mid$CrI_50 <- quantile(subset(Mid_both,Mid_both$Species=="Bull trout")$Overlap, prob = 0.5)
ls_SI_4$BLT_on_LKT$Mid$Mean <- mean(subset(Mid_both,Mid_both$Species=="Bull trout")$Overlap)
ls_SI_4$BLT_on_LKT$Mid$CrI_97.5 <- quantile(subset(Mid_both,Mid_both$Species=="Bull trout")$Overlap, prob = 0.975)
#BLT Late
ls_SI_4$BLT_on_LKT$Late$CrI_2.5 <- quantile(subset(Late_both,Mid_both$Species=="Bull trout")$Overlap, prob = 0.025)
ls_SI_4$BLT_on_LKT$Late$CrI_50 <- quantile(subset(Late_both,Mid_both$Species=="Bull trout")$Overlap, prob = 0.5)
ls_SI_4$BLT_on_LKT$Late$Mean <- mean(subset(Late_both,Late_both$Species=="Bull trout")$Overlap)
ls_SI_4$BLT_on_LKT$Late$CrI_97.5 <- quantile(subset(Late_both,Late_both$Species=="Bull trout")$Overlap, prob = 0.975)
#LKT Mid
ls_SI_4$LKT_on_BLT$Mid$CrI_2.5 <- quantile(subset(Mid_both,Mid_both$Species=="Lake trout")$Overlap, prob = 0.025)
ls_SI_4$LKT_on_BLT$Mid$CrI_50 <- quantile(subset(Mid_both,Mid_both$Species=="Lake trout")$Overlap, prob = 0.5)
ls_SI_4$LKT_on_BLT$Mid$Mean <- mean(subset(Mid_both,Mid_both$Species=="Lake trout")$Overlap)
ls_SI_4$LKT_on_BLT$Mid$CrI_97.5 <- quantile(subset(Mid_both,Mid_both$Species=="Lake trout")$Overlap, prob = 0.975)
#LKT Late
ls_SI_4$LKT_on_BLT$Late$CrI_2.5 <- quantile(subset(Late_both,Mid_both$Species=="Lake trout")$Overlap, prob = 0.025)
ls_SI_4$LKT_on_BLT$Late$CrI_50 <- quantile(subset(Late_both,Mid_both$Species=="Lake trout")$Overlap, prob = 0.5)
ls_SI_4$LKT_on_BLT$Late$Mean <- mean(subset(Late_both,Late_both$Species=="Lake trout")$Overlap)
ls_SI_4$LKT_on_BLT$Late$CrI_97.5 <- quantile(subset(Late_both,Mid_both$Species=="Lake trout")$Overlap, prob = 0.975)

stats <- c("Cr_2.5_mid","Cr_50_mid","Mean_mid","Cr_97.5_mid","Cr_2.5_late","Cr_50_late","Mean_late","Cr_97.5_late")
overlaps <- c("BLT_in_LKT","LKT_in_BLT")

SI_4<-data.frame(matrix(unlist(ls_SI_4), #make a dataframe out of layman_list
                        nrow=length(overlaps), #one row for each overlap
                        byrow=T),
                 stringsAsFactors=FALSE)
row.names(SI_4) <- overlaps
colnames(SI_4) <- stats
SI_4 <- round(SI_4,1)

#########################################################################################################
#Figure 1.7: binomial linear regression of conversion
#########################################################################################################
rm(list=ls())
graphics.off()
set.seed(1)

data <- read.csv("C:/Users/wainr/Documents/nature_08132020/final_datasets/conversion_09302020_2.csv", header=TRUE)
data$conv <- data$LKT/(data$LKT+data$BLT) #calculate conversion
data$timestep <- data$year-1955 #calculate timestep
conv_lm <- glm(conv ~ timestep, data = data, family = binomial) #lm using timestep
#conv_lm2 <- glm(conv ~ year, data = data, family = binomial) #lm using year
summary(conv_lm) #intercept and estimate
library(rsq)
rsq(conv_lm) #r-squared

#Figure 1.7 (both timestep and year)
ggplot(data = data, aes(x = year,y = conv))+
  geom_jitter(width = 0.2, height = 0, size = 5)+
  stat_smooth(method = "glm",
              method.args = list(family = "binomial"),
              se = TRUE,
              size = 3,
              color = "black",
              fullrange = TRUE)+
  theme_bw()+
  labs(y = "Conversion\n")+
  scale_y_continuous(breaks = c(0,0.5,1))+
  scale_x_continuous(name = "\nSurvey year",
                     breaks = c(1955,1970,1985,2000,2015,2030),
                     limits = c(1954,2031),
                     sec.axis = sec_axis(trans = ~.-1954,
                                         name = "Timestep (years)\n",
                                         breaks = c(-15,0,15,30,45,60,75)))+
  theme(axis.text.x = element_text(size = 20, color = "black"),
        axis.text.y = element_text(size = 20, color = "black"),
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        plot.title = element_text(size = 29))

#Figure 1.7 (timestep only)
ggplot(data = data, aes(x=timestep,y=conv))+
  geom_jitter(width = 0.2, height = 0, size=5)+
  stat_smooth(method = "glm",
              method.args = list(family = "binomial"),
              se = TRUE,
              size=3,
              color="black",
              fullrange=TRUE)+
  theme_bw()+
  labs(x="\nTimestep (years)",
    y="Conversion\n")+
  scale_y_continuous(breaks=c(0,0.5,1))+
  scale_x_continuous(breaks=c(-15,0,15,30,45,60),limits=c(-15,60))+
  theme(axis.text.x = element_text(size=20, color="black"),
        axis.text.y = element_text(size=20, color = "black"),
        axis.title.x = element_text(size=25),
        axis.title.y = element_text(size=25),
        plot.title = element_text(size=29))

#Figure 1.7 (year only)
ggplot(data = data, aes(x=year,y=conv))+
  geom_jitter(width = 0.2, height = 0, size=5)+
  stat_smooth(method = "glm",
              method.args = list(family = "binomial"),
              se = TRUE,
              size=3,
              color="black",
              fullrange=TRUE)+
  theme_bw()+
  labs(x="\nSurvey year",
       y="Conversion\n")+
  scale_y_continuous(breaks=c(0,0.5,1))+
  scale_x_continuous(breaks=c(1955,1970,1985,2000,2015,2030),limits=c(1955,2030))+
  theme(axis.text.x = element_text(size=20, color="black"),
        axis.text.y = element_text(size=20, color = "black"),
        axis.title.x = element_text(size=25),
        axis.title.y = element_text(size=25),
        plot.title = element_text(size=29))
#SI Figure 9
ggplot(data = data, aes(x=timestep,y=conv))+
  stat_smooth(method = "glm",
              method.args = list(family = "binomial"),
              se = TRUE,
              size=3,
              color="black",
              fullrange=TRUE)+
  theme_bw()+
  labs(x="\nTimestep (years)",
       y="Conversion\n")+
  scale_y_continuous(breaks=c(0,0.5,1))+
  scale_x_continuous(breaks=c(-15,0,15,30,45,60),limits=c(-15,60))+
  theme(axis.text.x = element_text(size=20, color="black"),
        axis.text.y = element_text(size=20, color = "black"),
        axis.title.x = element_text(size=25),
        axis.title.y = element_text(size=25),
        plot.title = element_text(size=29))+
  annotate("text", x = 12, y = .4, size=6,
           label = "Quartz", parse = TRUE)+
  annotate("text", x = 14, y = .71, size=6,
           label = "Swan", parse = TRUE)+
  annotate("text", x = 36, y = .78, size=6,
           label = "Logging", parse = TRUE)


#########################################################################################################
#Figure 1.8: littoral macroinvertebrate ordination
#########################################################################################################
rm(list=ls())
graphics.off()
set.seed(1)
library(vegan)
library(remotes)
library(ggplot2)
library(dplyr)
library(metR)
library(directlabels)

#load data
data <- read.csv("C:/Users/wainr/Documents/nature_08132020/final_datasets/transect_data.csv", header=TRUE)
lake <- read.csv("C:/Users/wainr/Documents/nature_08132020/final_datasets/lake_data.csv", header=TRUE)

#community statistics
lake$richness<-rowSums(data[,2:48] > 0) #richness
lake$div<-diversity(data[,2:48]) #diversity
lake$even<-lake$div/log(lake$richness) #evenness
lake$abundance<-rowSums(data[,2:48]) #abundance
species<-data[,2:48]

#nmds scores
nms<-metaMDS(species, k=3)

#permanova
adonis(species~Invasion_category, data = lake)

#stressplot
stressplot(nms)

#prep nmds scores for plotting
#ordisurf:
ordi <- ordisurf(nms ~ lake$Timestep_years) #created the ordisurf object
ordi.grid <- ordi$grid #extracts the ordisurf object
str(ordi.grid) #it's a list though - cannot be plotted as is
ordi.lake <- expand.grid(x = ordi.grid$x, y = ordi.grid$y) #get x and ys
ordi.lake$z <- as.vector(ordi.grid$z) #unravel the matrix for the z scores
ordi.lake.na <- data.frame(na.omit(ordi.lake)) #gets rid of the nas
ordi.lake.na #looks ready for plotting!
lake.NMDS.data <- lake #there are other ways of doing this. But this is the way I do it for ease of plotting
lake.NMDS.data$NMDS1 <- nms$points[ ,1] #this puts the NMDS scores for the plots into a new dataframe. you could put them into an existing one if you preferred.
lake.NMDS.data$NMDS2 <- nms$points[ ,2]

lake.NMDS.data$Invasion_category <- factor(lake.NMDS.data$Invasion_category,
                                   levels=c('Reference',
                                            'Mid',
                                            'Late'))
#color version
ggplot(lake.NMDS.data, aes(x = NMDS1, y = NMDS2))+
  geom_point(aes(shape = Invasion_category, size=Invasion_category), alpha = 0.6)+ #plots the NMDS points, with shape by invasion
  theme_minimal() + #for aesthetics
  stat_ellipse(aes(x=NMDS1, y=NMDS2, group=Invasion_category, fill = Invasion_category),
               type = "norm",
               geom = "polygon",
               alpha = 0.3,
               level = 0.95)+
  labs(shape = "Invasion", size="Invasion", group="Invasion", fill="Invasion")+#another way to set the labels, in this case, for the colour legend
  geom_hline(yintercept=0, linetype="dashed")+
  geom_vline(xintercept=0, linetype="dashed")+
  theme(legend.key = element_blank(),  #removes the box around each legend item
        legend.position = "bottom", #legend at the bottom
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.box.just = "centre",
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        plot.title = element_text(size=30),
        legend.text = element_text(size=20),
        legend.title = element_text(size=25),
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank())

#black and white version
ggplot(lake.NMDS.data, aes(x = NMDS1, y = NMDS2))+
  geom_point(aes(shape = Invasion_category, size=Invasion_category), alpha = 0.6)+ #plots the NMDS points, with shape by invasion
  theme_minimal() + #for aesthetics
  stat_ellipse(aes(x=NMDS1, y=NMDS2, group=Invasion_category, fill = Invasion_category),
               type = "norm",
               geom = "polygon",
               alpha = 0.4,
               level = 0.95)+
  labs(shape = "Invasion", size="Invasion", group="Invasion", fill="Invasion")+#another way to set the labels, in this case, for the colour legend
  geom_hline(yintercept=0, linetype="dashed")+
  geom_vline(xintercept=0, linetype="dashed")+
  scale_fill_manual(values=c("gray","black","dark gray"))+ # get rid of this line to make a color figure
  theme(legend.key = element_blank(),  #removes the box around each legend item
        legend.position = "bottom", #legend at the bottom
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.box.just = "centre",
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        plot.title = element_text(size=30),
        legend.text = element_text(size=20),
        legend.title = element_text(size=25),
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank())


