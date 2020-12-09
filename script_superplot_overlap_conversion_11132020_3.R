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

#script to make conversion and overlap superplot
#Make objects, then plots

######################################
#Overlap objects
######################################
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

LKT_both$Invasion <- factor(LKT_both$Invasion, levels=c("Mid","Late"))
BLT_both$Invasion <- factor(BLT_both$Invasion, levels=c("Mid","Late"))
######################################
#Make the conversion objects
######################################
data <- read.csv("C:/Users/wainr/Documents/nature_08132020/final_datasets/conversion_09302020_2.csv", header=TRUE)
data$conv <- data$LKT/(data$LKT+data$BLT) #calculate conversion
data$timestep <- data$year-1955 #calculate timestep


######################################
#Then make the plots
######################################
p1 <- 
  ggplot(data=BLT_both, aes(x=Overlap,y=..scaled.., fill=Invasion, color = Invasion))+ 
  geom_density(alpha=0.5)+
  scale_fill_manual(values=c("gray","black"))+
  scale_color_manual(values=c("gray","black"))+
  scale_x_continuous(limits=c(0,100),breaks=c(0,50,100))+ #for percent instead of proportional overlap
  scale_y_continuous(limits=c(0,1),breaks=c(0,0.5,1))+
  labs(title="Bull trout overlap with lake trout",
       x="\nDiet Overlap (%)\n",
       y="Scaled Density\n",
       fill="Invasion")+
  theme_bw()+
  geom_vline(xintercept = quantile(subset(BLT_both,BLT_both$Invasion=="Late")$Overlap, prob = 0.5), linetype="solid", size=2)+ #Late
  geom_vline(xintercept = quantile(subset(BLT_both,BLT_both$Invasion=="Mid")$Overlap, prob = 0.5), linetype="dotted", size=2)+ #Mid
  theme(axis.text.x = element_text(size = 15, color = "black"),
        axis.text.y = element_text(size = 15, color = "black"),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        plot.title = element_text(size = 15))

p2 <-
  ggplot(data=LKT_both, aes(x=Overlap,y=..scaled.., fill=Invasion, color = Invasion))+ 
  geom_density(alpha=0.5)+
  scale_fill_manual(values=c("gray","black"))+
  scale_color_manual(values=c("gray","black"))+
  scale_x_continuous(limits=c(0,100),breaks=c(0,50,100))+ #for percent instead of proportional overlap
  scale_y_continuous(limits=c(0,1),breaks=c(0,0.5,1))+
  labs(title="Lake trout overlap with bull trout",
       x="\nDiet Overlap (%)\n",
       y="Scaled Density\n",
       fill="Invasion")+
  theme_bw()+
  geom_vline(xintercept = quantile(subset(LKT_both,LKT_both$Invasion=="Late")$Overlap, prob = 0.5), linetype="solid", size=2)+ #Late
  geom_vline(xintercept = quantile(subset(LKT_both,LKT_both$Invasion=="Mid")$Overlap, prob = 0.5), linetype="dotted", size=2)+ #Mid
  ylab(NULL)+
  theme(axis.text.x = element_text(size = 15, color = "black"),
        axis.text.y = element_text(size = 15, color = "black"),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        plot.title = element_text(size = 15))

prow <- plot_grid(
  p1 + theme(legend.position="none"),
  p2 + theme(legend.position="none"),
  align = 'vh',
  labels = c("3A", "3B"),
  hjust = -1,
  nrow = 1
)

legend <- get_legend(
  # create some space to the left of the legend
  p1 + theme(legend.box.margin = margin(0, 0, 0, 12),
             legend.text = element_text(size=15),
             legend.title = element_text(size=15))
)

top <- plot_grid(prow, legend, rel_widths = c(3, .4))

C <-
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
                     breaks = c(1955,1970,1985,2000,2015,2030,2045),
                     limits = c(1955,2045),
                     sec.axis = sec_axis(trans = ~.-1955,
                                         name = "Invasion timeline (years)\n",
                                         breaks = c(0,15,30,45,60,75,90)))+
  theme(axis.text.x = element_text(size = 15, color = "black"),
        axis.text.y = element_text(size = 15, color = "black"),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        plot.title = element_text(size = 15))

bottom <- plot_grid(
  C,NULL,
  rel_widths = c(3, .4),
  align = 'vh',
  labels = c("3C",''),
  hjust = -1,
  nrow = 1
)
######################################
#Figure 3 superplot (conversion and BLT/LKT overlap)
######################################
plot_grid(top, bottom,
          nrow=2)

