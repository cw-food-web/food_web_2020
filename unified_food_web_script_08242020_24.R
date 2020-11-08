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
isotope_dataset <- read.csv("C:/Users/wainr/Documents/nature_08132020/final_datasets/isotope_dataset_MT_lake_trout_food_web_10012020.csv", header=TRUE)
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
isotope_dataset <- read.csv("C:/Users/wainr/Documents/nature_08132020/final_datasets/isotope_dataset_MT_lake_trout_food_web_10012020.csv", header=TRUE)

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

isotope_dataset <- read.csv("C:/Users/wainr/Documents/nature_08132020/final_datasets/isotope_dataset_MT_lake_trout_food_web_10012020.csv", header=TRUE)
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
isotope_dataset <- read.csv("C:/Users/wainr/Documents/nature_08132020/final_datasets/isotope_dataset_MT_lake_trout_food_web_10012020.csv", header=TRUE)
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
isotope_dataset <- read.csv("C:/Users/wainr/Documents/nature_08132020/final_datasets/isotope_dataset_MT_lake_trout_food_web_10012020.csv", header=TRUE)
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
#but the default SEA.B plot makes me sad so we'll re-do siberDensityPlot in ggplot

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
isotope_dataset <- read.csv("C:/Users/wainr/Documents/nature_08132020/final_datasets/isotope_dataset_MT_lake_trout_food_web_10012020.csv", header=TRUE)
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
isotope_dataset <- read.csv("C:/Users/wainr/Documents/nature_08132020/final_datasets/isotope_dataset_MT_lake_trout_food_web_10012020.csv", header=TRUE)
data <- subset(isotope_dataset,Sample_type=="Fish")
data$Invasion_category <- factor(data$Invasion_category, levels= c("Reference","Mid","Late"))
data$Functional_group <- factor(data$Functional_group, levels= c("Bull trout","Generalist fish","Pelagic forage fish","Lake trout","Littoral forage fish"))
data <- subset(data,real.name!="Lake trout")
d15N_lmer<- lmer(d15N_corr ~ Invasion_category*Functional_group+
                   (1|Lake_name),
                 data=data)
summary(d15N_lmer) #bull trout d15N p-values

#p-values for littoral forage fish d15N
rm(list=ls())
graphics.off()
set.seed(1)
isotope_dataset <- read.csv("C:/Users/wainr/Documents/nature_08132020/final_datasets/isotope_dataset_MT_lake_trout_food_web_10012020.csv", header=TRUE)
data <- subset(isotope_dataset,Sample_type=="Fish")
data$Invasion_category <- factor(data$Invasion_category, levels= c("Reference","Mid","Late"))
data$Functional_group <- factor(data$Functional_group, levels= c("Littoral forage fish","Generalist fish","Pelagic forage fish","Bull trout","Lake trout"))
data <- subset(data,real.name!="Lake trout")
d15N_lmer<- lmer(d15N_corr ~ Invasion_category*Functional_group+
                   (1|Lake_name),
                 data=data)
summary(d15N_lmer) #littoral forage fish d15N p-values

#p-values for generalist fish d15N
rm(list=ls())
graphics.off()
set.seed(1)
isotope_dataset <- read.csv("C:/Users/wainr/Documents/nature_08132020/final_datasets/isotope_dataset_MT_lake_trout_food_web_10012020.csv", header=TRUE)
data <- subset(isotope_dataset,Sample_type=="Fish")
data$Invasion_category <- factor(data$Invasion_category, levels= c("Reference","Mid","Late"))
data$Functional_group <- factor(data$Functional_group, levels= c("Generalist fish","Pelagic forage fish","Bull trout","Lake trout","Littoral forage fish"))
data <- subset(data,real.name!="Lake trout")
d15N_lmer<- lmer(d15N_corr ~ Invasion_category*Functional_group+
                   (1|Lake_name),
                 data=data)
summary(d15N_lmer) #generalist fish d15N p-values

#p-values for Pelagic forage fish d15N
rm(list=ls())
graphics.off()
set.seed(1)
isotope_dataset <- read.csv("C:/Users/wainr/Documents/nature_08132020/final_datasets/isotope_dataset_MT_lake_trout_food_web_10012020.csv", header=TRUE)
data <- subset(isotope_dataset,Sample_type=="Fish")
data$Invasion_category <- factor(data$Invasion_category, levels= c("Reference","Mid","Late"))
data$Functional_group <- factor(data$Functional_group, levels= c("Pelagic forage fish","Bull trout","Lake trout","Littoral forage fish","Generalist fish"))
data <- subset(data,real.name!="Lake trout")
d15N_lmer<- lmer(d15N_corr ~ Invasion_category*Functional_group+
                   (1|Lake_name),
                 data=data)
summary(d15N_lmer) #Pelagic forage fish d15N p-values

#p-values for Lake trout d15N
rm(list=ls())
graphics.off()
set.seed(1)
isotope_dataset <- read.csv("C:/Users/wainr/Documents/nature_08132020/final_datasets/isotope_dataset_MT_lake_trout_food_web_10012020.csv", header=TRUE)
data <- subset(isotope_dataset,Sample_type=="Fish")
data <- subset(isotope_dataset,Functional_group%in%c("Lake trout","Bull trout")) #must exclude lake trout because lake trout are not present in all invasion categories
data <- subset(data,Invasion_category != "Reference")
data$Invasion_category <- factor(data$Invasion_category, levels = c("Mid","Late"))
data$Functional_group <- factor(data$Functional_group, levels= c("Lake trout","Bull trout"))
d15N_lmer<- lmer(d15N_corr ~ Invasion_category*Functional_group+
                   (1|Lake_name),
                 data=data)
summary(d15N_lmer)#lake trout d15N p-values

#########################################################################################################
#Figure 1.4 mixed effects linear modelling for d13C
#########################################################################################################
#Figure 1.4, d13C
rm(list=ls())
graphics.off()
set.seed(1)
isotope_dataset <- read.csv("C:/Users/wainr/Documents/nature_08132020/final_datasets/isotope_dataset_MT_lake_trout_food_web_10012020.csv", header=TRUE)
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
isotope_dataset <- read.csv("C:/Users/wainr/Documents/nature_08132020/final_datasets/isotope_dataset_MT_lake_trout_food_web_10012020.csv", header=TRUE)
data <- subset(isotope_dataset,Sample_type=="Fish")
data$Invasion_category <- factor(data$Invasion_category, levels= c("Reference","Mid","Late"))
data$Functional_group <- factor(data$Functional_group, levels= c("Bull trout","Generalist fish","Pelagic forage fish","Lake trout","Littoral forage fish"))
data <- subset(data,real.name!="Lake trout")
d13C_lmer<- lmer(d13C_corr ~ Invasion_category*Functional_group+
                   (1|Lake_name),
                 data=data)
summary(d13C_lmer) #bull trout d13C p-values

#p-values for littoral forage fish
rm(list=ls())
graphics.off()
set.seed(1)
isotope_dataset <- read.csv("C:/Users/wainr/Documents/nature_08132020/final_datasets/isotope_dataset_MT_lake_trout_food_web_10012020.csv", header=TRUE)
data <- subset(isotope_dataset,Sample_type=="Fish")
data$Invasion_category <- factor(data$Invasion_category, levels= c("Reference","Mid","Late"))
data$Functional_group <- factor(data$Functional_group, levels= c("Littoral forage fish","Generalist fish","Pelagic forage fish","Bull trout","Lake trout"))
data <- subset(data,real.name!="Lake trout")
d13C_lmer<- lmer(d13C_corr ~ Invasion_category*Functional_group+
                   (1|Lake_name),
                 data=data)
summary(d13C_lmer) #littoral forage fish d13C p-values

#p-values for generalist fish
rm(list=ls())
graphics.off()
set.seed(1)
isotope_dataset <- read.csv("C:/Users/wainr/Documents/nature_08132020/final_datasets/isotope_dataset_MT_lake_trout_food_web_10012020.csv", header=TRUE)
data <- subset(isotope_dataset,Sample_type=="Fish")
data$Invasion_category <- factor(data$Invasion_category, levels= c("Reference","Mid","Late"))
data$Functional_group <- factor(data$Functional_group, levels= c("Generalist fish","Pelagic forage fish","Bull trout","Lake trout","Littoral forage fish"))
data <- subset(data,real.name!="Lake trout")
d13C_lmer<- lmer(d13C_corr ~ Invasion_category*Functional_group+
                   (1|Lake_name),
                 data=data)
summary(d13C_lmer) #generalist fish d13C p-values

#p-values for Pelagic forage fish
rm(list=ls())
graphics.off()
set.seed(1)
isotope_dataset <- read.csv("C:/Users/wainr/Documents/nature_08132020/final_datasets/isotope_dataset_MT_lake_trout_food_web_10012020.csv", header=TRUE)
data <- subset(isotope_dataset,Sample_type=="Fish")
data$Invasion_category <- factor(data$Invasion_category, levels= c("Reference","Mid","Late"))
data$Functional_group <- factor(data$Functional_group, levels= c("Pelagic forage fish","Bull trout","Lake trout","Littoral forage fish","Generalist fish"))
data <- subset(data,real.name!="Lake trout")
d13C_lmer<- lmer(d13C_corr ~ Invasion_category*Functional_group+
                   (1|Lake_name),
                 data=data)
summary(d13C_lmer) #Pelagic forage fish d13C p-values

#p-values for Lake trout d13C
rm(list=ls())
graphics.off()
set.seed(1)
isotope_dataset <- read.csv("C:/Users/wainr/Documents/nature_08132020/final_datasets/isotope_dataset_MT_lake_trout_food_web_10012020.csv", header=TRUE)
data <- subset(isotope_dataset,Sample_type=="Fish")
data <- subset(isotope_dataset,Functional_group%in%c("Lake trout","Bull trout")) #must exclude lake trout because lake trout are not present in all invasion categories
data <- subset(data,Invasion_category != "Reference")
data$Invasion_category <- factor(data$Invasion_category, levels = c("Mid","Late"))
data$Functional_group <- factor(data$Functional_group, levels= c("Lake trout","Bull trout"))
d13C_lmer<- lmer(d13C_corr ~ Invasion_category*Functional_group+
                   (1|Lake_name),
                 data=data)
summary(d13C_lmer)#lake trout d13C p-values


#########################################################################################################
#Figure 1.5: niche overlap 
#########################################################################################################
rm(list=ls())
graphics.off()
set.seed(1)

#load data
isotope_dataset <- read.csv("C:/Users/wainr/Documents/nature_08132020/final_datasets/isotope_dataset_MT_lake_trout_food_web_10012020.csv", header=TRUE)
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
    axis.text.y = element_text(size = 30, color = "white"),
    axis.title.y = element_text(size = 30, color = "white"),
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
#Figure 1.6: diet proportions, MixSIAR
#########################################################################################################
rm(list=ls())
graphics.off()
set.seed(1)

setwd("/Users/wainr/Documents/nature_08132020/final_datasets/")
mix <- load_mix_data(filename="bull_trout_mix2.csv", #d13C and d15N of individual blt
                     iso_names=c("d13C","d15N"), 
                     factors=c("Invasion","Lake"), 
                     fac_random=c(FALSE,TRUE), #Invasion is fixed effect, lake is random effect
                     fac_nested=c(FALSE,TRUE), #Lake is nested under invasion
                     cont_effects=NULL)
source <- load_source_data(file="bull_trout_source2.csv", #mean and sd d13C and d15N of prey
                           source_factors="Invasion", 
                           conc_dep=FALSE, 
                           data_type="means", 
                           mix)
discr <- load_discr_data(filename="bull_trout_trophic_discrimination.csv", #Post 2002 trophic discrimination
                         mix)


model_filename <- "MixSIAR_model.txt"   # Name of the JAGS model file
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(model_filename, resid_err, process_err, mix, source)

jags.1 <- run_model(run="test", mix, source, discr, model_filename, #update from test as necessary
                    alpha.prior = 1, resid_err, process_err)

#save.image(file='mixsiar_outputs.RData') #save your environment to reproduce plots

output_options <- list(summary_save = FALSE,
                       summary_name = "summary_statistics",
                       sup_post = FALSE,
                       plot_post_save_pdf = FALSE,
                       plot_post_name = "posterior_density",
                       sup_pairs = FALSE,
                       plot_pairs_save_pdf = FALSE,
                       plot_pairs_name = "pairs_plot",
                       sup_xy = FALSE,
                       plot_xy_save_pdf = FALSE,
                       plot_xy_name = "xy_plot",
                       gelman = TRUE,
                       heidel = FALSE,
                       geweke = TRUE,
                       diag_save = TRUE,
                       diag_name = "diagnostics",
                       indiv_effect = FALSE,
                       plot_post_save_png = FALSE,
                       plot_pairs_save_png = FALSE,
                       plot_xy_save_png = FALSE)

#output_JAGS (directly below) is the command to generate MixSIAR's default plots
#output_JAGS(jags.1, mix, source, output_options) 
#I needed to customize figures so I got the source code (below)

################################################################################
#Plot Figure 1.6
#Adjusted output_JAGS function
#select all 700+ lines of this and run at once to get the right plots
################################################################################
library(ggplot2)
{
  mcmc.chains <- jags.1$BUGSoutput$n.chains
  N <- mix$N
  n.re <- mix$n.re
  n.effects <- mix$n.effects
  if(n.re==1){
    random_effects <- ifelse(mix$FAC[[1]]$re,mix$FAC[[1]]$name,mix$FAC[[2]]$name)
  }
  if(n.re==2){
    random_effects <- mix$factors
  }
  n.sources <- source$n.sources
  source_names <- source$source_names
  # p.global <- ilr.global <- ilr.fac1 <- ilr.fac2 <- fac1.sig <- fac2.sig <- NULL
  # ind.sig <- ..scaled.. <- p.fac1 <- p.fac2 <- p.ind <- sources <- NULL
  # R2jags::attach.jags(jags.1)
  jags1.mcmc <- coda::as.mcmc(jags.1)
  n.draws <- length(jags.1$BUGSoutput$sims.list$p.global[,1])
  
  # Post-processing for 2 FE or 1FE + 1RE
  #   calculate p.both = ilr.global + ilr.fac1 + ilr.fac2
  if(mix$fere){
    fac2_lookup <- list()
    for(f1 in 1:mix$FAC[[1]]$levels){
      fac2_lookup[[f1]] <- unique(mix$FAC[[2]]$values[which(mix$FAC[[1]]$values==f1)])
    }
    ilr.both <- array(NA,dim=c(n.draws,mix$FAC[[1]]$levels, mix$FAC[[2]]$levels, n.sources-1))
    p.both <- array(NA,dim=c(n.draws,mix$FAC[[1]]$levels, mix$FAC[[2]]$levels, n.sources))
    cross.both <- array(data=NA,dim=c(n.draws,mix$FAC[[1]]$levels, mix$FAC[[2]]$levels,n.sources,n.sources-1))
    e <- matrix(rep(0,n.sources*(n.sources-1)),nrow=n.sources,ncol=(n.sources-1))
    for(i in 1:(n.sources-1)){
      e[,i] <- exp(c(rep(sqrt(1/(i*(i+1))),i),-sqrt(i/(i+1)),rep(0,n.sources-i-1)))
      e[,i] <- e[,i]/sum(e[,i])
    }
    for(i in 1:n.draws){
      for(f1 in 1:mix$FAC[[1]]$levels) {
        for(f2 in fac2_lookup[[f1]]){
          for(src in 1:(n.sources-1)) {
            ilr.both[i,f1,f2,src] <- jags.1$BUGSoutput$sims.list$ilr.global[i,src] + jags.1$BUGSoutput$sims.list$ilr.fac1[i,f1,src] + jags.1$BUGSoutput$sims.list$ilr.fac2[i,f2,src];
            cross.both[i,f1,f2,,src] <- (e[,src]^ilr.both[i,f1,f2,src])/sum(e[,src]^ilr.both[i,f1,f2,src]);
            # ilr.both[,f1,f2,src] <- ilr.global[,src] + ilr.fac1[,f1,src] + ilr.fac2[,f2,src];
          }
          for(src in 1:n.sources) {
            p.both[i,f1,f2,src] <- prod(cross.both[i,f1,f2,src,]);
          }
          p.both[i,f1,f2,] <- p.both[i,f1,f2,]/sum(p.both[i,f1,f2,]);
        } # f2
      } # f1
    }
  } # end fere
  
  ###########################################################################################
  # XY/Trace Plots
  ###########################################################################################
  
  # XY plots for p.global and factor SD's
  if(!output_options[[9]]){  # if 'suppress XY plot' is NOT checked
    # XY plot for p.global
    dev.new()
    print(lattice::xyplot(coda::as.mcmc(jags.1$BUGSoutput$sims.list$p.global),strip=lattice::strip.custom(factor.levels=source_names)))
    
    # Save the xy p.global plot to file
    if(output_options[[10]]){ # svalue(plot_xy_save_pdf)
      mypath <- file.path(paste(getwd(),"/",output_options[[11]],"_diet_p.pdf",sep=""))  # svalue(plot_xy_name)
      dev.copy2pdf(file=mypath)
    }
    if(output_options[[20]]){ # svalue(plot_xy_save_png)
      mypath <- file.path(paste(getwd(),"/",output_options[[11]],"_diet_p.png",sep=""))  # svalue(plot_xy_name)
      dev.copy(png,mypath)
    }
    
    # XY plot for the factor SDs
    if(output_options[[17]]){ # include_indiv ('ind.sig' is in the model)
      dev.new()
      traceplot_labels <- rep("",length(random_effects)+1)  # +1 because we need to add "Individual SD"
      if(n.re > 0){
        for(i in 1:length(random_effects)){
          traceplot_labels[i] <- paste(random_effects[i]," SD",sep="")
        }
      }
      traceplot_labels[length(random_effects)+1] <- "Individual SD"
      if(n.re==2) print(lattice::xyplot(coda::as.mcmc(cbind(jags.1$BUGSoutput$sims.list$fac1.sig,jags.1$BUGSoutput$sims.list$fac2.sig,jags.1$BUGSoutput$sims.list$ind.sig)),strip=lattice::strip.custom(factor.levels=traceplot_labels)))
      if(n.re==1){
        if(mix$FAC[[1]]$re){
          print(lattice::xyplot(coda::as.mcmc(cbind(jags.1$BUGSoutput$sims.list$fac1.sig,jags.1$BUGSoutput$sims.list$ind.sig)),strip=lattice::strip.custom(factor.levels=traceplot_labels)))
        } else { # FAC 2 is the 1 random effect
          print(lattice::xyplot(coda::as.mcmc(cbind(jags.1$BUGSoutput$sims.list$fac2.sig,jags.1$BUGSoutput$sims.list$ind.sig)),strip=lattice::strip.custom(factor.levels=traceplot_labels)))
        }
      }
      if(n.re==0) print(lattice::xyplot(coda::as.mcmc(jags.1$BUGSoutput$sims.list$ind.sig),strip=lattice::strip.custom(factor.levels=traceplot_labels)))
    } else { # Individual SD is not in the model (no 'ind.sig')
      if(n.re > 0){
        dev.new()
        traceplot_labels <- rep("",length(random_effects))
        for(i in 1:length(random_effects)) { traceplot_labels[i] <- paste(random_effects[i]," SD",sep="") }
        if(n.re==2) print(lattice::xyplot(coda::as.mcmc(cbind(jags.1$BUGSoutput$sims.list$fac1.sig,jags.1$BUGSoutput$sims.list$fac2.sig)),strip=lattice::strip.custom(factor.levels=traceplot_labels)))
        if(n.re==1){
          if(mix$FAC[[1]]$re){
            print(lattice::xyplot(coda::as.mcmc(cbind(jags.1$BUGSoutput$sims.list$fac1.sig)),strip=lattice::strip.custom(factor.levels=traceplot_labels)))
          } else { # FAC 2 is the 1 random effect
            print(lattice::xyplot(coda::as.mcmc(cbind(jags.1$BUGSoutput$sims.list$fac2.sig)),strip=lattice::strip.custom(factor.levels=traceplot_labels)))
          }
        }
      }
    }
    # Save the xy factor SD plot to file
    if(output_options[[10]]){ # svalue(plot_xy_save_pdf)
      mypath <- file.path(paste(getwd(),"/",output_options[[11]],"_SD.pdf",sep=""))  # svalue(plot_xy_name)
      dev.copy2pdf(file=mypath)
    }
    if(output_options[[20]]){ # svalue(plot_xy_save_png)
      mypath <- file.path(paste(getwd(),"/",output_options[[11]],"_SD.png",sep=""))  # svalue(plot_xy_name)
      dev.copy(png,mypath)
    }
  }
  
  # Fancy pairs plot of p.global
  # Contour plots in the upper right, histograms on the diagonal, correlation coefficients in the lower left
  if(!output_options[[6]]){   # if 'suppress pairs plot' is NOT checked
    dev.new()
    # Function: panel.hist (from ?pairs)
    # Purpose: creates histograms on the diagonal of the pairs plot matrix
    panel.hist <- function(x, ...){
      usr <- par("usr"); on.exit(par(usr), add=TRUE)
      par(usr = c(usr[1:2], 0, 1.5) )
      h <- hist(x, plot = FALSE)
      breaks <- h$breaks; nB <- length(breaks)
      y <- h$counts; y <- y/max(y)
      rect(breaks[-nB], 0, breaks[-1], y, col='blue', xlim=c(0,1),...)
    }
    # Function: panel.cor (from http://personality-project.org/r/r.graphics.html)
    # Purpose: prints correlation coefficients in the lower panel,
    #          scales text sizes to the correlation coefficient magnitudes
    panel.cor <- function(x, y, digits=2, prefix="", cex.cor){
      usr <- par("usr"); on.exit(par(usr), add=TRUE)
      par(usr = c(0, 1, 0, 1))
      r = (cor(x, y,use="pairwise"))
      txt <- format(c(r, 0.123456789), digits=digits)[1]
      txt <- paste(prefix, txt, sep="")
      if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
      text(0.5, 0.5, txt, cex = cex * abs(r))
    }
    # Function: panel.contour (inspired by http://stats.stackexchange.com/questions/31726/scatterplot-with-contour-heat-overlay)
    # Purpose: replaces scatterplots with colored contour plots
    panel.contour <- function(x,y){
      n.lines <- 4  # number of contour lines
      my.cols <- rev(RColorBrewer::brewer.pal(n.lines, "RdYlBu"))   # gets some pretty colors
      z <- MASS::kde2d(x,y)   # calculates the 2D kernel density that the contour function needs
      contour(z, drawlabels=FALSE, nlevels=n.lines, col=my.cols, add=TRUE)
    }
    pairs(jags.1$BUGSoutput$sims.list$p.global, labels=source_names, diag.panel=panel.hist, lower.panel=panel.cor, upper.panel=panel.contour)
    
    # Save the plot to file
    if(output_options[[7]]){ # svalue(plot_pairs_save_pdf)
      mypath <- file.path(paste(getwd(),"/",output_options[[8]],".pdf",sep=""))  # svalue(plot_pairs_name)
      dev.copy2pdf(file=mypath)
    }
    if(output_options[[19]]){ # svalue(plot_pairs_save_png)
      mypath <- file.path(paste(getwd(),"/",output_options[[8]],".png",sep=""))  # svalue(plot_pairs_name)
      dev.copy(png,mypath)
    }
  }
  
  ######################################################################
  # Posterior density plots
  ######################################################################
  if(!output_options[[3]]){   # if 'suppress posterior plots' is NOT checked
    n.draws <- length(jags.1$BUGSoutput$sims.list$p.global[,1])   # number of posterior draws
    if(mix$n.fe == 0){ # only if there are no fixed effects, otherwise p.global is meaningless
      # Posterior density plot for p.global
      dev.new()
      df <- data.frame(sources=rep(NA,n.draws*n.sources), x=rep(NA,n.draws*n.sources))  # create empty data frame
      for(i in 1:n.sources){
        df$x[seq(1+n.draws*(i-1),i*n.draws)] <- as.matrix(jags.1$BUGSoutput$sims.list$p.global[,i]) # fill in the p.global[i] values
        df$sources[seq(1+n.draws*(i-1),i*n.draws)] <- rep(source_names[i],n.draws)  # fill in the source names
      }
      my.title <- "Overall Population"
      print(ggplot2::ggplot(df, ggplot2::aes(x=x, fill=sources, colour=sources)) +
              ggplot2::geom_density(alpha=.3, ggplot2::aes(y=..scaled..)) +
              ggplot2::theme_bw() +
              ggplot2::xlab("Proportion of Diet") +
              ggplot2::ylab("Scaled Posterior Density") +
              #ggplot2::xlim(0,1) +
              ggplot2::scale_y_continuous(breaks=c(0,0.5,1)) +
              ggplot2::labs(title = my.title) +
              ggplot2::theme(legend.position="right"))
      
      # Save the plot to file
      if(output_options[[4]]){ # svalue(plot_post_save_pdf)
        mypath <- file.path(paste(getwd(),"/",output_options[[5]],"_diet_p_global.pdf",sep=""))  # svalue(plot_post_name)
        dev.copy2pdf(file=mypath)
      }
      if(output_options[[18]]){ # svalue(plot_post_save_png)
        mypath <- file.path(paste(getwd(),"/",output_options[[5]],"_diet_p_global.png",sep=""))  # svalue(plot_post_name)
        dev.copy(png,mypath)
      }
    }
    
    if(n.effects >= 1 & mix$n.fe != 2){
      # Posterior density plots for p.fac1's
      for(f1 in 1:mix$FAC[[1]]$levels){    # formerly factor1_levels
        dev.new()
        df <- data.frame(sources=rep(NA,n.draws*n.sources), x=rep(NA,n.draws*n.sources))  # create empty data frame
        for(src in 1:n.sources){
          df$x[seq(1+n.draws*(src-1),src*n.draws)] <- as.matrix(jags.1$BUGSoutput$sims.list$p.fac1[,f1,src]) # fill in the p.fac1[f1] values
          df$sources[seq(1+n.draws*(src-1),src*n.draws)] <- rep(source_names[src],n.draws)  # fill in the source names
        }
        my.title <- mix$FAC[[1]]$labels[f1]  # formerly factor1_names
        print(ggplot2::ggplot(df, ggplot2::aes(x=x, fill=sources, colour=sources)) +
                #        geom_density(alpha=.3) +
                ggplot2::geom_density(alpha=.3, ggplot2::aes(y=..scaled..)) +
                #ggplot2::xlim(0,1) +
                ggplot2::scale_y_continuous(limits=c(0,1),breaks=c(0,0.5,1)) +
                ggplot2::scale_x_continuous(limits=c(0,1),breaks=c(0,0.5,1)) +
                ggplot2::theme_bw() +
                ggplot2::xlab("Proportion of Diet") +
                ggplot2::ylab("Scaled Posterior Density") +
                ggplot2::labs(title = my.title,
                              color = "Diet item",
                              fill = "Diet item") +
                ggplot2::theme(legend.position="right",axis.text = element_text(size = 20),
                               axis.title = element_text(size = 20),
                               legend.text = element_text(size = 15),
                               legend.title = element_text(size = 15),
                               plot.title = element_text(size = 20)))
        
        # Save the plot to file
        if(output_options[[4]]){ # svalue(plot_post_save_pdf)
          mypath <- file.path(paste(getwd(),"/",output_options[[5]],"_diet_p_",mix$FAC[[1]]$labels[f1],".pdf",sep=""))  # svalue(plot_post_name), factor1_names
          dev.copy2pdf(file=mypath)
        }
        if(output_options[[18]]){ # svalue(plot_post_save_png)
          mypath <- file.path(paste(getwd(),"/",output_options[[5]],"_diet_p_",mix$FAC[[1]]$labels[f1],".png",sep=""))  # svalue(plot_post_name), factor1_names
          dev.copy(png,mypath)
        }
      } # end p.fac1 posterior plots
      
      if(n.re==2){
        # Posterior density plots for p.fac2's
        for(f2 in 1:mix$FAC[[2]]$levels){  # formerly factor2_levels
          dev.new()
          df <- data.frame(sources=rep(NA,n.draws*n.sources), x=rep(NA,n.draws*n.sources))  # create empty data frame
          for(src in 1:n.sources){
            df$x[seq(1+n.draws*(src-1),src*n.draws)] <- as.matrix(jags.1$BUGSoutput$sims.list$p.fac2[,f2,src]) # fill in the p.fac2 values
            df$sources[seq(1+n.draws*(src-1),src*n.draws)] <- rep(source_names[src],n.draws)  # fill in the source names
          }
          #my.title <- paste(factor1_names[f1],", ", random_effects[2]," ",f2,sep="") # plot title (ex. "Region 1, Pack 3")
          my.title <- mix$FAC[[2]]$labels[f2] # formerly factor2_names
          print(ggplot2::ggplot(df, ggplot2::aes(x=x, fill=sources, colour=sources)) +
                  #             geom_density(alpha=.3) +
                  ggplot2::geom_density(alpha=.3, ggplot2::aes(y=..scaled..)) +
                  ggplot2::theme_bw() +
                  ggplot2::xlim(0,1) +
                  ggplot2::xlab("Proportion of Diet") +
                  ggplot2::ylab("Scaled Posterior Density") +
                  ggplot2::labs(title = my.title) +
                  ggplot2::theme(legend.position="none",
                                 axis.text = element_text(size = 20),
                                 axis.title = element_text(size = 20),
                                 legend.text = element_text(size = 15),
                                 legend.title = element_text(size = 15)))
          
          # Save the plot as a pdf file
          if(output_options[[4]]){ # svalue(plot_post_save_pdf)
            mypath <- file.path(paste(getwd(),"/",output_options[[5]],"_diet_p_",mix$FAC[[2]]$labels[f2],".pdf",sep="")) #  svalue(plot_post_name), factor2_names
            dev.copy2pdf(file=mypath)
          }
          if(output_options[[18]]){  # svalue(plot_post_save_png)
            mypath <- file.path(paste(getwd(),"/",output_options[[5]],"_diet_p_",mix$FAC[[2]]$labels[f2],".png",sep="")) #  svalue(plot_post_name), factor2_names
            dev.copy(png,mypath)
          }
        }# end p.fac2 posterior plots
      } # end if(n.re==2)
    } # end if(n.effects >=1 & n.fe != 2)
    
    # Posterior density plots for p.both (when 2 FE or 1FE + 1RE)
    if(mix$fere){
      for(f1 in 1:mix$FAC[[1]]$levels) {
        for(f2 in fac2_lookup[[f1]]){
          dev.new()
          df <- data.frame(sources=rep(NA,n.draws*n.sources), x=rep(NA,n.draws*n.sources))  # create empty data frame
          for(src in 1:n.sources){
            df$x[seq(1+n.draws*(src-1),src*n.draws)] <- as.matrix(p.both[,f1,f2,src]) # fill in the p.both values
            df$sources[seq(1+n.draws*(src-1),src*n.draws)] <- rep(source_names[src],n.draws)  # fill in the source names
          }
          my.title <- paste(mix$FAC[[1]]$labels[f1],mix$FAC[[2]]$labels[f2],sep=" ") # formerly factor2_names
          print(ggplot2::ggplot(df, ggplot2::aes(x=x, fill=sources, colour=sources)) +
                  ggplot2::geom_density(alpha=.3, ggplot2::aes(y=..scaled..)) +
                  ggplot2::theme_bw() +
                  ggplot2::xlim(0,1) +
                  ggplot2::xlab("Proportion of Diet") +
                  ggplot2::ylab("Scaled Posterior Density") +
                  ggplot2::labs(title = my.title) +
                  ggplot2::theme(legend.position="none"))
          
          # Save the plot as a pdf file
          if(output_options[[4]]){ # svalue(plot_post_save_pdf)
            mypath <- file.path(paste(getwd(),"/",output_options[[5]],"_diet_p_",mix$FAC[[1]]$labels[f1],"_",mix$FAC[[2]]$labels[f2],".pdf",sep="")) #  svalue(plot_post_name), factor2_names
            dev.copy2pdf(file=mypath)
          }
          if(output_options[[18]]){  # svalue(plot_post_save_png)
            mypath <- file.path(paste(getwd(),"/",output_options[[5]],"_diet_p_",mix$FAC[[1]]$labels[f1],"_",mix$FAC[[2]]$labels[f2],".png",sep="")) #  svalue(plot_post_name), factor2_names
            dev.copy(png,mypath)
          }
        } # f2
      } # f1
    }
    
    # Posterior density plot for fac1.sig, fac2.sig, and ind.sig
    if(n.re > 0 || output_options[[17]]){ # only have an SD posterior plot if we have Individual, Factor1, or Factor2 random effects)
      dev.new()
      n.re_ind <- n.re + as.numeric(output_options[[17]]) # this*n.draws will be the length of the plot data frame
      level <- c()
      x <- c()
      if(output_options[[17]]){ # if Individual is in the model, add ind.sig to the SD plot
        level <- c(level,rep("Individual SD",n.draws))
        x <- c(x,jags.1$BUGSoutput$sims.list$ind.sig)
      }
      if(n.re==1){ # if Factor.1 is in the model, add fac1.sig to the SD plot
        if(mix$FAC[[1]]$re){
          level <- c(level,rep(paste(mix$FAC[[1]]$name," SD",sep=""),n.draws))
          x <- c(x,jags.1$BUGSoutput$sims.list$fac1.sig)
        } else { # FAC 2 is the random effect
          level <- c(level,rep(paste(mix$FAC[[2]]$name," SD",sep=""),n.draws))
          x <- c(x,jags.1$BUGSoutput$sims.list$fac2.sig)
        }
      }
      if(n.re==2){ # if Factor.2 is in the model, add fac1.sig and fac2.sig to the SD plot
        level <- c(level,rep(paste(random_effects[1]," SD",sep=""),n.draws), rep(paste(random_effects[2]," SD",sep=""),n.draws))
        x <- c(x,jags.1$BUGSoutput$sims.list$fac1.sig,jags.1$BUGSoutput$sims.list$fac2.sig)
      }
      df2 <- data.frame(level=level, x=x) # create the SD plot data frame
      
      print(ggplot2::ggplot(df2, ggplot2::aes(x=x, fill=level, colour=level)) +
              #        geom_density(alpha=.3) +
              ggplot2::geom_density(alpha=.3) +
              ggplot2::theme_bw() +
              ggplot2::xlab(expression(sigma)) +
              ggplot2::ylab("Test line 334") +
              ggplot2::theme(legend.position="none"))   # + xlim(0,2)
      
      # Save the plot to file
      if(output_options[[4]]){ # svalue(plot_post_save_pdf)
        mypath <- file.path(paste(getwd(),"/",output_options[[5]],"_SD.pdf",sep=""))  # svalue(plot_post_name)
        dev.copy2pdf(file=mypath)
      }
      if(output_options[[18]]){ # svalue(plot_post_save_png)
        mypath <- file.path(paste(getwd(),"/",output_options[[5]],"_SD.png",sep=""))  # svalue(plot_post_name)
        dev.copy(png,mypath)
      }
    }
  }
  
  # Calculate the summary statistics for the variables we're interested in (p.global's and factor SD's, maybe p.ind's)
  # We print them out later, at the very bottom
  sig_labels <- NULL; ind_labels <- NULL; fac1_labels <- NULL; fac2_labels <- NULL; sig_stats <- NULL;
  getQuant <- function(x) quantile(x,probs=c(.025,.05,.25,.5,.75,.95,.975))
  getMeanSD <- function(x) cbind(round(apply(x,2,mean),3),round(apply(x,2,sd),3))
  
  stats <- NULL
  sig_stats <- NULL
  sig_labels <- NULL
  eps_stats <- NULL
  eps_labels <- NULL
  # print(mix)
  # print(mix$n.fe)
  if(mix$n.fe == 0){
    global_quants <- t(round(apply(jags.1$BUGSoutput$sims.list$p.global,2,getQuant),3))
    global_means <- getMeanSD(jags.1$BUGSoutput$sims.list$p.global)
    stats <- cbind(global_means, global_quants)
    global_labels <- rep(NA,n.sources)
    for(src in 1:n.sources){
      global_labels[src] <- paste("p.global.",source_names[src],sep="")
    }
    rownames(stats) <- global_labels
  }
  if(n.effects > 0 & mix$n.fe != 2){
    fac1_quants <- as.matrix(reshape::cast(reshape2::melt(round(apply(jags.1$BUGSoutput$sims.list$p.fac1,c(2,3),getQuant),3)),Var3+Var2~Var1)[,-c(1,2)])
    fac1_quants <- t(apply(fac1_quants,1,sort)) # BUG FIX 10/28/14, quantiles were out of order from cast/melt (thanks to Jason Waite)
    fac1_means <- cbind(reshape2::melt(round(apply(jags.1$BUGSoutput$sims.list$p.fac1,c(2,3),mean),3))$value, reshape2::melt(round(apply(jags.1$BUGSoutput$sims.list$p.fac1,c(2,3),sd),3))$value)
    fac1_stats <- cbind(fac1_means,fac1_quants)
    fac1_labels <- rep(NA,mix$FAC[[1]]$levels*n.sources)
    for(src in 1:n.sources){
      for(f1 in 1:mix$FAC[[1]]$levels){
        fac1_labels[mix$FAC[[1]]$levels*(src-1)+f1] <- paste("p.",mix$FAC[[1]]$labels[f1],".",source_names[src],sep="")
      }
    }
    rownames(fac1_stats) <- fac1_labels
    stats <- rbind(stats,fac1_stats)
    if(mix$FAC[[1]]$re){
      sig_stats <- cbind(getMeanSD(jags.1$BUGSoutput$sims.list$fac1.sig),t(round(apply(jags.1$BUGSoutput$sims.list$fac1.sig,2,getQuant),3)))
      sig_labels <- paste(mix$FAC[[1]]$name,".SD",sep="")
    }
  }
  if(n.re==2){
    fac2_quants <- as.matrix(reshape::cast(reshape2::melt(round(apply(jags.1$BUGSoutput$sims.list$p.fac2,c(2,3),getQuant),3)),Var3+Var2~Var1)[,-c(1,2)])
    fac2_quants <- t(apply(fac2_quants,1,sort)) # BUG FIX 10/28/14, quantiles were out of order from cast/melt (thanks to Jason Waite)
    fac2_means <- cbind(reshape2::melt(round(apply(jags.1$BUGSoutput$sims.list$p.fac2,c(2,3),mean),3))$value, reshape2::melt(round(apply(jags.1$BUGSoutput$sims.list$p.fac2,c(2,3),sd),3))$value)
    fac2_stats <- cbind(fac2_means,fac2_quants)
    fac2_labels <- rep(NA,mix$FAC[[2]]$levels*n.sources)
    for(src in 1:n.sources){
      for(f2 in 1:mix$FAC[[2]]$levels){
        fac2_labels[mix$FAC[[2]]$levels*(src-1)+f2] <- paste("p.",mix$FAC[[2]]$labels[f2],".",source_names[src],sep="")
      }
    }
    rownames(fac2_stats) <- fac2_labels
    stats <- rbind(stats,fac2_stats)
    if(mix$FAC[[2]]$re){
      sig_stats <- rbind(sig_stats,cbind(getMeanSD(jags.1$BUGSoutput$sims.list$fac2.sig),t(round(apply(jags.1$BUGSoutput$sims.list$fac2.sig,2,getQuant),3))))
      sig_labels <- c(sig_labels,paste(mix$FAC[[2]]$name,".SD",sep=""))
    }
  }
  if(mix$fere){
    fac2_quants <- matrix(NA,nrow=n.sources*length(unlist(fac2_lookup)),ncol=7)
    fac2_means <- matrix(NA,nrow=n.sources*length(unlist(fac2_lookup)),ncol=2)
    fac2_labels <- rep(NA,n.sources*length(unlist(fac2_lookup)))
    i <- 1
    for(f1 in 1:mix$FAC[[1]]$levels) {
      for(f2 in fac2_lookup[[f1]]){
        for(src in 1:n.sources){
          fac2_quants[i,] <- getQuant(p.both[,f1,f2,src])
          fac2_means[i,] <- c(mean(p.both[,f1,f2,src]),sd(p.both[,f1,f2,src]))
          fac2_labels[i] <- paste("p",mix$FAC[[1]]$labels[f1],mix$FAC[[2]]$labels[f2],source_names[src],sep=".")
          i <- i+1
        }
      }
    }
    # fac2_quants <- as.matrix(cast(melt(round(apply(p.both,c(2,3,4),getQuant,na.rm=TRUE),3)),X4+X3+X2~X1)[,-c(1,2)])
    # fac2_quants <- t(apply(fac2_quants,1,sort)) # BUG FIX 10/28/14, quantiles were out of order from cast/melt (thanks to Jason Waite)
    # fac2_means <- cbind(melt(round(apply(p.fac2,c(2,3),mean),3))$value, melt(round(apply(p.fac2,c(2,3),sd),3))$value)
    fac2_stats <- round(cbind(fac2_means,fac2_quants),3)
    rownames(fac2_stats) <- fac2_labels
    stats <- rbind(stats,fac2_stats)
    if(mix$FAC[[2]]$re){
      sig_stats <- rbind(sig_stats,cbind(getMeanSD(jags.1$BUGSoutput$sims.list$fac2.sig),t(round(apply(jags.1$BUGSoutput$sims.list$fac2.sig,2,getQuant),3))))
      sig_labels <- c(sig_labels,paste(mix$FAC[[2]]$name,".SD",sep=""))
    }
  }
  
  if(output_options[[17]]){ # include_indiv (if Individual is in the model)
    ind_quants <- as.matrix(reshape::cast(reshape2::melt(round(apply(p.ind,c(2,3),getQuant),3)),X3+X2~X1)[,-c(1,2)])
    ind_quants <- t(apply(ind_quants,1,sort)) # BUG FIX 10/28/14, quantiles were out of order from cast/melt (thanks to Jason Waite)
    ind_means <- cbind(reshape2::melt(round(apply(p.ind,c(2,3),mean),3))$value, reshape2::melt(round(apply(p.ind,c(2,3),sd),3))$value)
    ind_stats <- cbind(ind_means,ind_quants)
    ind_labels <- rep(NA,N*n.sources)
    for(src in 1:n.sources){
      for(j in 1:N){
        ind_labels[N*(src-1)+j] <- paste("p.Ind ",j,".",source_names[src],sep="")
      }
    }
    sig_stats <- rbind(sig_stats,cbind(getMeanSD(jags.1$BUGSoutput$sims.list$ind.sig),t(round(apply(jags.1$BUGSoutput$sims.list$ind.sig,2,getQuant),3))))
    sig_labels <- c(sig_labels,"Individual.SD")
    rownames(ind_stats) <- ind_labels
    stats <- rbind(stats, ind_stats)
  }
  
  # Add SD stats to the top of the summary
  rownames(sig_stats) <- sig_labels
  stats <- rbind(sig_stats,stats)
  
  # Add epsilon (multiplicative error term) to stat summary
  # Also plot posterior density
  epsTF <- "resid.prop" %in% names(jags.1$BUGSoutput$sims.list)
  if(epsTF){
    eps_stats <- cbind(getMeanSD(jags.1$BUGSoutput$sims.list$resid.prop),t(round(apply(jags.1$BUGSoutput$sims.list$resid.prop,2,getQuant),3)))
    eps_labels <- paste0("Epsilon.", 1:mix$n.iso)
    rownames(eps_stats) <- eps_labels
    stats <- rbind(eps_stats,stats)
    
    # posterior plot
    level <- c()
    x <- c()
    for(j in 1:mix$n.iso){
      level <- c(level,rep(eps_labels[j], n.draws))
      x <- c(x, jags.1$BUGSoutput$sims.list$resid.prop[,j])    
    }
    df2 <- data.frame(level=level, x=x) 
    
    dev.new()
    print(ggplot2::ggplot(df2, ggplot2::aes(x=x, fill=level, colour=level)) +
            ggplot2::geom_density(alpha=.3) +
            ggplot2::theme_bw() +
            ggplot2::xlab(expression(epsilon)) +
            ggplot2::ylab("Posterior Density") +
            ggplot2::theme(legend.position=c(.95,.95), legend.justification=c(1,1), legend.title=ggplot2::element_blank()))   # + xlim(0,2)
    
  }
  colnames(stats) <- c("Mean","SD","2.5%","5%","25%","50%","75%","95%","97.5%")
  
  # Pack 1 stats only
  #stats[grep("Pack 1",rownames(stats)),]
  
  # Region stats only
  #stats[grep("Region",rownames(stats)),]
  
  # Region stats, by Region
  # byVec <- function(x){ind <- NULL; for(i in 1:length(x)){ ind <- c(ind,grep(x[i],rownames(stats)))}; return(ind)}
  # stats[byVec(mix$RE[[1]]$labels),]
  
  # All means
  # stats[,"Mean"]
  
  # Region means only
  # stats[byVec(mix$RE[[1]]$labels),"Mean"]
  
  ################################################################################
  # Calulate diagnostics
  ################################################################################
  # Get number of variables in the model
  n.var <- coda::nvar(jags1.mcmc)
  # Gelman-Rubin diagnostic
  if(output_options[[12]]){  # if Gelman is checked
    if(mcmc.chains == 1){
      gelman <- "*** Error: Gelman diagnostic requires more than one chain ***"
    }
    if(mcmc.chains > 1){    # Gelman diagnostic requires more than one chain
      # Gelman diagnostic, for when the multivariate Gelman fails (matrix not positive definite)
      # Remove the test results for dummy/empty variables
      gelman <- matrix(NA, nrow=n.var, ncol=2)
      for (v in 1:coda::nvar(jags1.mcmc)) {
        gelman[v,] <- coda::gelman.diag(jags1.mcmc[,v])$psrf
      }
      #gelman <- gelman[ind,]
      colnames(gelman) <- c("Point est.","Upper C.I.")
      rownames(gelman) <- coda::varnames(jags1.mcmc)
      #rownames(gelman) <- c(sig_labels,global_labels,fac1_labels,fac2_labels,ind_labels)
      gelman.all <- gelman[which(!is.nan(gelman[,1])),] # Remove dummy variables (show up as NA)
      gelman_short <- gelman[order(gelman[,1],decreasing=T),]
      if(n.var>10) gelman_short <- gelman_short[1:10,]
      gelman_fail <- c(length(which(gelman[,1]>1.01)), length(which(gelman[,1]>1.05)), length(which(gelman[,1]>1.1)))
    }
  }
  
  # Heidelberger and Welch's diagnostic
  # Remove the test results for dummy/empty variables
  if(output_options[[13]]){   # if Heidel is checked
    heidel <- coda::heidel.diag(jags1.mcmc)
    w <- which(!is.na(heidel[[1]][,"pvalue"]))  # find all the non-dummy variables
    heidel.all <- data.frame(matrix(NA,nrow=length(w),ncol=3*mcmc.chains))  # create empty data frame
    colstring <- rep(NA,mcmc.chains*3)  # vector of column names
    for(i in 1:mcmc.chains){
      heidel.tmp <- as.data.frame(heidel[[i]][w,c("stest","pvalue","htest")]) # stest, pvalue, and htest are the relevant statistics - get them
      heidel.all[,(3*i-2):(3*i)] <- heidel.tmp
      colstring[(3*i-2):(3*i)] <- c(paste("stest.",i,sep=""), paste("pval.",i,sep=""), paste("hwtest.",i,sep="")) # create the appropriate column names
    }
    #heidel.all <- heidel.all[ind,]
    #rownames(heidel.all) <- c(sig_labels,global_labels,fac1_labels,fac2_labels,ind_labels)
    rownames(heidel.all) <- coda::varnames(jags1.mcmc)[w]
    colnames(heidel.all) <- colstring
    heidel.all <- round(heidel.all,3)
    heidel.all <- replace(heidel.all,heidel.all==0,"fail")  # A normal call to 'heidel.diag' prints "fail" and "pass", for some reason they turn to 0's and 1's
    heidel.all <- replace(heidel.all,heidel.all==1,"pass")  # when you access the statistics directly.  Here we turn the 0's and 1's back into "fail" and "pass"
    # When the stationarity test fails, hwtest returns <NA>...change these NAs to 'fail'
    heidel.all <- replace(heidel.all,is.na(heidel.all),"fail")
    # Count the number of failures (2 tests per chain - 'stationarity' and 'half-width')
    stest_fail <- rep(NA,mcmc.chains); hwtest_fail <- rep(NA,mcmc.chains)
    for(i in 1:mcmc.chains){
      stest_fail[i] <- sum(heidel.all[,3*i-2]=="fail")
      hwtest_fail[i] <- sum(heidel.all[,3*i]=="fail")
    }
    heidel_fail <- rbind(stest_fail,hwtest_fail)
    rownames(heidel_fail) <- c("Stationarity","Half-width")
    colnames(heidel_fail) <- paste("Chain",1:mcmc.chains)
  }
  
  # Geweke diagnostic
  # Remove the test results for dummy/empty variables
  if(output_options[[14]]){ # if Geweke is checked
    geweke <- coda::geweke.diag(jags1.mcmc)
    geweke.all <- data.frame(matrix(NA,nrow=n.var,ncol=mcmc.chains))    # create empty data frame
    colstring <- rep(NA,mcmc.chains)    # vector of column names
    for(i in 1:mcmc.chains){
      geweke.tmp <- as.data.frame(geweke[[i]]$z) # get the relevant geweke statistics
      geweke.all[,i] <- geweke.tmp
      colstring[i] <- c(paste("chain",i,sep=""))  # create the column names "chain1", "chain2", etc.
    }
    #geweke.all <- geweke.all[ind,]
    #rownames(geweke.all) <- c(sig_labels,global_labels,fac1_labels,fac2_labels,ind_labels)
    rownames(geweke.all) <- coda::varnames(jags1.mcmc)
    colnames(geweke.all) <- colstring
    geweke.all <- round(geweke.all,3)
    w <- which(!is.nan(geweke[[1]]$z))  # find all the non-dummy variables
    geweke.all <- geweke.all[w,]
    geweke_fail <- matrix(NA,nrow=1,ncol=mcmc.chains)
    for(i in 1:mcmc.chains){
      geweke_fail[1,i] <- sum(abs(geweke.all[,i])>1.96)
    }
    colnames(geweke_fail) <- paste("Chain",1:mcmc.chains)
    rownames(geweke_fail) <- "Geweke"
  }
  
  ################################################################################
  # Print diagnostics
  ################################################################################
  
  if(output_options[[12]]){  # svalue(gelman)
    cat("
        ################################################################################
        # Gelman-Rubin Diagnostic
        ################################################################################
        Generally the Gelman diagnostic should be < 1.05
        ",paste("Out of ",n.var," variables: ",gelman_fail[1]," > 1.01",sep=""),"
        ",paste(gelman_fail[2]," > 1.05",sep=""),"
        ",paste(gelman_fail[3]," > 1.1",sep=""),"
        The worst variables are:
        ",sep="")
    print(gelman_short)
    
    #print(gelman)
    
    if(output_options[[15]]){  # svalue(diag_save)
      mypath <- file.path(paste(getwd(),"/",output_options[[16]],".txt",sep=""))  # svalue(diag_name)
      out <- capture.output(gelman)
      out2 <- capture.output(gelman_short)
      cat("
          ################################################################################
          # Gelman-Rubin Diagnostic
          ################################################################################
          Generally the Gelman diagnostic should be < 1.05
          ",paste("Out of ",n.var," variables: ",gelman_fail[1]," > 1.01",sep=""),"
          ",paste(gelman_fail[2]," > 1.05",sep=""),"
          ",paste(gelman_fail[3]," > 1.1",sep=""),"
          The worst variables are:
          ",out2,"
          And here are the Gelman diagnostics for all variables:
          ",out,sep="\n", file=mypath, append=FALSE)
    } # end save Gelman
  } # end Gelman printout
  
  if(output_options[[13]]){  # svalue(heidel)
    cat("
        ################################################################################
        # Heidelberger and Welch Diagnostic
        ################################################################################
        A few failures is normal and acceptable...
        Number of failures in each chain (out of ",n.var," variables):
        ",sep="")
    print(heidel_fail)
    #print(heidel.all)
    
    if(output_options[[15]]){  # svalue(diag_save)
      mypath <- file.path(paste(getwd(),"/",output_options[[16]],".txt",sep=""))  # svalue(diag_name)
      out <- capture.output(heidel.all)
      out2 <- capture.output(heidel_fail)
      cat("
          ################################################################################
          # Heidelberger and Welch Diagnostic
          ################################################################################
          A few failures is normal and acceptable...
          Number of failures in each chain (out of ",n.var," variables):
          ",out2,"
          And here are the Heidelberger-Welch diagnostics for all variables:
          ",out,sep="\n", file=mypath, append=output_options[[12]]) # svalue(gelman)
    } # end save Heidel
  } # end Heidel printout
  
  if(output_options[[14]]){ # svalue(geweke)
    cat("
        ################################################################################
        # Geweke Diagnostic
        ################################################################################
        The Geweke diagnostic is a standard z-score, so we'd expect 5% to be outside +/-1.96
        Number of variables outside +/-1.96 in each chain (out of ",n.var,"):
        ",sep="")
    print(geweke_fail)
    #print(geweke.all)
    
    if(output_options[[15]]){  # svalue(diag_save)
      mypath <- file.path(paste(getwd(),"/",output_options[[16]],".txt",sep=""))  # svalue(diag_name)
      out <- capture.output(geweke.all)
      out2 <- capture.output(geweke_fail)
      cat("
          ################################################################################
          # Geweke Diagnostic
          ################################################################################
          The Geweke diagnostic is a standard z-score, so we'd expect 5% to be outside +/-1.96
          Number of variables outside +/-1.96 in each chain (out of ",n.var,"):
          ",out2,"
          And here are the Geweke diagnostics for all variables:
          ",out,sep="\n", file=mypath, append=output_options[[12]]||output_options[[13]]) # svalue(gelman) || svalue(heidel)
    } # end Geweke save
  } # end Geweke printout
  
  DIC <- jags.1$BUGSoutput$DIC
  cat("
      ################################################################################
      # Summary Statistics
      ################################################################################
      DIC = ",DIC,sep="")
  out1 <- capture.output(stats)
  cat("
      ",out1,sep="\n")
  
  if(output_options[[1]]){  # svalue(summary_save)
    mypath <- file.path(paste(getwd(),"/",output_options[[2]],".txt",sep=""))  # svalue(summary_name)
    cat("
        #################################################################
        # Summary Statistics
        #################################################################
        DIC = ",DIC,sep="", file=mypath, append=FALSE)
    cat("
        ",out1,sep="\n", file=mypath, append=TRUE)
  }
  
  # Plot any continuous effects
  if(mix$n.ce > 0){
    plot_continuous_var(jags.1,mix,source,output_options)
  }
  
  # Use ggmcmc package to create diagnostic plots
  diag_filename <- paste(getwd(),"/",output_options[[16]],".pdf",sep="")
  ggmcmc::ggmcmc(ggmcmc::ggs(jags1.mcmc),file=diag_filename,plot=c("Rhat","geweke","density","traceplot","running","autocorrelation","crosscorrelation"))
  
  # Return p.both if 2 FE or 1FE + 1RE
  if(mix$fere){
    return(p.both)
  } else return(NULL) # otherwise return nothing
  
} # end function output_JAGS



#########################################################################################################
#Figure 1.7: binomial linear regression of conversion
#########################################################################################################
rm(list=ls())
graphics.off()
set.seed(1)

data <- read.csv("C:/Users/wainr/Documents/nature_08132020/final_datasets/conversion_09302020_2.csv", header=TRUE)
data$conv <- data$LKT/(data$LKT+data$BLT) #calculate conversion
data$timestep <- data$year-1969 #calculate timestep
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
                     sec.axis = sec_axis(trans = ~.-1969,
                                       name = "Timestep (years)\n",
                                       breaks = c(-15,0,15,30,45,60)))+
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






