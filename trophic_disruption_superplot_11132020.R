rm(list=ls())
graphics.off()
set.seed(1)

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

#script to make trophic disruption super-plot
#Make objects, then plots

######################################
#Make the SIBER output
######################################

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

######################################
#Make the d15N output
######################################
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
d15N_output <- rbind(invasion,invasion2)
d15N_output$Functional_group <- factor(d15N_output$Functional_group, levels= c("Bull trout","Lake trout","Littoral forage fish","Generalist fish","Pelagic forage fish"))
d15N_output$Invasion_category <- factor(d15N_output$Invasion_category, levels= c("Reference","Mid","Late"))

######################################
#Make the d13C output
######################################
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
d13C_output <- rbind(invasion,invasion2)
d13C_output$Functional_group <- factor(d13C_output$Functional_group, levels= c("Bull trout","Lake trout","Littoral forage fish","Generalist fish","Pelagic forage fish"))
d13C_output$Invasion_category <- factor(d13C_output$Invasion_category, levels= c("Reference","Mid","Late"))

######################################
#Make the individual plots
######################################

A <-
  ggplot(data = SEA.B_2,
       aes(x=invasion.category, y=value, fill=invasion.category)) +
  geom_boxplot()+
  scale_fill_manual(values=c("white","light gray","grey32"))+
  theme_bw()+
  scale_y_continuous(limits=c(0,25))+
  labs(y = expression("Diet breadth ( \211"^2~")"),
       x = "\nInvasion")+
  facet_wrap(~Species,
             ncol = 5)+
  xlab(NULL)+
  theme(legend.position="none",
        strip.text = element_text(size = 15),
        axis.text.x = element_text(size = 15, color = "black"),
        axis.text.y = element_text(size = 15, color = "black"),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        plot.title = element_text(size = 15))


B <-
  ggplot(d15N_output, aes (x= Invasion_category, y=fit, group=Functional_group))+
  geom_line()+
  geom_point(size = 2) +
  geom_errorbar(width=.1, aes(ymin=fit-se, ymax=fit+se)) + #you could also change this to fit+/- the se
  labs(y= expression(paste(delta^15, "N"~'( \211 )')), x="\nInvasion")+
  facet_wrap(~Functional_group, ncol=5)+
  xlab(NULL)+
  theme_bw()+
  theme(strip.text = element_text(size = 15),
        axis.text.x = element_text(size = 15, color = "black"),
        axis.text.y = element_text(size = 15, color = "black"),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        plot.title = element_text(size = 15))

C <- 
  ggplot(d13C_output, aes (x= Invasion_category, y=fit, group=Functional_group))+
  geom_line()+
  geom_point(size = 2) +
  geom_errorbar(width=.1, aes(ymin=fit-se, ymax=fit+se)) + #you could also change this to fit+/- the se
  labs(y= expression(paste(delta^13, "C"~'( \211 )')), x="\nInvasion")+
  facet_wrap(~Functional_group, ncol=5)+
  theme_bw()+
  theme(strip.text = element_text(size = 15),
        axis.text.x = element_text(size = 15, color = "black"),
        axis.text.y = element_text(size = 15, color = "black"),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        plot.title = element_text(size = 15))

######################################
#Make the Figure 2 superplot
######################################
plot_grid(NULL, A, NULL, B, NULL,C,
          rel_widths = c(0.1,3,0.1,3,0.1,3),
          align = 'vh',
          labels = c('2A','','2B','','2C',''),
          nrow = 3,
          ncol = 2)

