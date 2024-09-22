#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#Seasonal SWC thresholds ####
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#Select and check data ####
file <- file.choose(new = FALSE)#NLFA.csv
dat <- read.csv(file,header=T,dec=".",sep=",", row.names=1)

## install packages ####
library(ggplot2)                  
library(devtools)                
library(ggthemes)                 
library(scales)
library(RColorBrewer)
library(ggpmisc)
library(wesanderson)
library(plyr)
library(dplyr)
library(lubridate)
library(nlme)
library(MetBrewer)
#install.packages("ggpmisc")
## setting the seed ####
#  makes output reproducible 
set.seed(42)


## Shortcut functions ####

g <- glue::glue
r <- rex::rex
write_tsv <- function(...) data.table::fwrite(..., sep = "\t")

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=TRUE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)}
  
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm))},measurevar)
  datac <- rename(datac, c("mean" = measurevar))
  datac$se <- datac$sd / sqrt(datac$N) 
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}


my.graph.theme<-
  theme_bw()+
  theme(plot.title= element_text(face="bold", colour = "black", size=12),
        axis.title.x= element_text(face="plain", color= "black", size=11),
        axis.title.y= element_text(face="plain", color= "black", size=11),
        axis.text.x= element_text(face="plain", colour="black", size=11),
        axis.text.y= element_text(face="plain", colour="black", size=11),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())

names(dat)

dat$date <- as.Date(dat$date, "%Y/%m/%d")
ag <- aggregate(Theta ~ date+Trt+Measurement.Depth, dat, mean) #

ag1<-ag%>% mutate(Treatment = case_when(Trt == "C0T0" ~"Ambient",
                                Trt == "C0T0D"~"Drought",
                                Trt == "C2T2" ~"Future Climate",
                                Trt == "C2T2D"~"Future Climate + Drought"))
names(ag1)
ag1$Theta100<-ag1$Theta*100

SWC<-ggplot(ag1, aes(x=date, y=Theta100, group= Treatment, color=Treatment)) +
  scale_color_manual(values=c("#5189b8","#f4e093" , "#943124", "#e77e39"))+
  geom_line(size=0.5)+
  #scale_x_date(date_breaks = "3 month",date_labels = "%Y %b")+
  scale_x_date(limits=as.Date(c("2020-06-01","2020-08-15")))+
  ylim(0, 50) +
  facet_grid(Measurement.Depth ~., scales = "free",space= "free")+
  #ggtitle("Soil water content - daily averages")+
  ylab(expression(paste("Vol. soil moisture (%)"))) +
  my.graph.theme
SWC

ggsave(plot=SWC,"Fig. S1.png", width = 6, height = 4)

