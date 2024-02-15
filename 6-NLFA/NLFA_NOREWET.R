###Load packages
library(ggplot2)
library(ggpubr)
library(reshape2)
library(forcats)
library(ggResidpanel)
library(car)

## Load data

file <- file.choose(new = FALSE)#NLFA.csv
data <- read.csv(file,header=T,dec=".",sep=",", row.names=1)
data<-subset(data, Time!="Rewet")


#run statistic for each group
library(nlme)
#column<-sub(" rewet.*", "", data_ind[,7])
data_stat<-data
#data_stat$Treatment2<-column
drought <- c("Drought","Future Climate + Drought")
notdrought<- c("Ambient","Future Climate")
future <- c("Future Climate","Future Climate + Drought")
ambient <- c("Ambient","Drought")
#rewet <- c("Rewet")
#data_stat$varia<-row.names(data_stat)
#rewet<-c("aD1.rL","aD2.rL","aD3.rL","aD4.rL","aFD1.rL","aFD2.rL","aFD3.rL","aFD4.rL")
#notrewet<-c("a1L","a2L","a3L","a4L","aF1L","aF2L","aF3L","aF4L")
#fut<-c("aF1L","aF2L","aF3L","aF4L","aFD1.rL","aFD2.rL","aFD3.rL","aFD4.rL")
#notfut<-c("aD1.rL","aD2.rL","aD3.rL","aD4.rL","a1L","a2L","a3L","a4L")

data_stat$plot<-c("a1","a2","a3","a4","aD1","aD2","aD3","aD4","aF1","aF2","aF3","aF4","aFD1","aFD2","aFD3","aFD4","a1","a2","a3","a4","aD1","aD2","aD3","aD4","aF1","aF2","aF3","aF4","aFD1","aFD2","aFD3","aFD4")

#data_stat$Drought <- ifelse(data_stat$Treatment %in% drought, "drought", "no")
data_stat$Drought[data_stat$Treatment %in% drought] <- 'drought'
data_stat$Drought[data_stat$Treatment %in% notdrought] <- 'not'
#data_stat$Future <- ifelse(data_stat$Treatment %in% future, "future", "ambient")

data_stat$Future[data_stat$Treatment %in% future] <- 'future'
data_stat$Future[data_stat$Treatment %in% ambient] <- 'ambient'

#data_stat$Rewet[data_stat$varia %in% rewet] <- 'rewet'
#data_stat$Rewet[data_stat$varia %in% notrewet] <- 'not'

#data_stat$Fut2[data_stat$varia %in% fut] <- 'fut'
#data_stat$Fut2[data_stat$varia %in% notfut] <- 'amb'

# names of groups
group.names <- colnames(data_stat)[3:10]
no.groups <- length(group.names)

# create a named list to hold the fitted models
fitlist <- as.list(1:no.groups)
names(fitlist) <- group.names
anovalist<- as.list(1:no.groups)
names(anovalist) <- group.names
shapirolist<- as.list(1:no.groups)
names(shapirolist) <- group.names
levenelist<- as.list(1:no.groups)
names(levenelist) <- group.names
# loop over lipid names
for(i in group.names){ 
  
  # print status
  print(paste("Running entity:", i, "which is", which(group.names==i), "out of", no.groups))
  
  # create temporary data matrix and model formula
  tmp <- subset(data_stat[, c(i,"Drought","Future","Time")], Time=="Drought")
  
  #start with untransformed variables and use log if data is not normally distributed
  #fml <- as.formula( paste("log(", i,")", "~", paste(c("Drought","Future"), collapse="*") ) )
  fml <- as.formula( paste(i, "~", paste(c("Drought","Future"), collapse="*") ) )#use if not normally distribute

  #develop formula for levenetest
  lvn <-as.formula( paste( "residua", "~", paste(c("Drought","Future"), collapse="*") ) )
  
  # assign fit to list by name
  lme_fit <- lme(fml, random=~1|plot, data=subset(data_stat, Time=="Drought"))
  residua<-resid(lme_fit)
  anovalist[[i]] <- anova(lme_fit)
  shapirolist[[i]] <- shapiro.test(resid(lme_fit))
  levenelist[[i]] <- leveneTest(lvn, data=subset(data_stat, Time=="Drought"))
  print(resid_panel(lme_fit))
  #add weights if leven test fails weights = varIdent(form =~ 1 | Drought *Future),
}
fitlist
anovalist
shapirolist
levenelist



# create a named list1 to hold the fitted models
fitlist1 <- as.list(1:no.groups)
names(fitlist1) <- group.names
anovalist1<- as.list(1:no.groups)
names(anovalist1) <- group.names
shapirolist1<- as.list(1:no.groups)
names(shapirolist1) <- group.names
levenelist1<- as.list(1:no.groups)
names(levenelist1) <- group.names
# loop over lipid names

for(i in group.names){ 
  
  # print status
  print(paste("Running entity:", i, "which is", which(group.names==i), "out of", no.groups))
  
  # create temporary data matrix and model formula
  tmp <- subset(data_stat[, c(i,"Fut2","Rewet")], Rewet=="rewet" |Rewet=="not")
  fml <- as.formula( paste( i, "~", paste(c("Rewet","Fut2"), collapse="*") ) )
  
  #develop formula for levenetest
  lvn <-as.formula( paste( "residua", "~", paste(c("Rewet","Fut2"), collapse="*") ) )
  
  # assign fit to list1 by name
  lme_fit <- lme(fml, random=~1|plot, data=subset(data_stat, Rewet=="rewet" |Rewet=="not"))
  residua<-resid(lme_fit)
  anovalist[[i]] <- anova(lme_fit)
  shapirolist[[i]] <- shapiro.test(resid(lme_fit))
  levenelist[[i]] <- leveneTest(lvn, data=subset(data_stat, Time=="Drought"))
  print(resid_panel(lme(fml, random=~1|plot, data=subset(data_stat, Rewet=="rewet" |Rewet=="not"))))
}
fitlist1
anovalist1
shapirolist1
levenelist1
#Not normal or heterogeneous
#total (shapiro fail)
#fungi (shapiro fail)
#general (shapiro fail)
#percentofPLFAamf (leven fail)

# create a named list2 to hold the fitted models
fitlist2 <- as.list(1:no.groups)
names(fitlist2) <- group.names
anovalist2<- as.list(1:no.groups)
names(anovalist2) <- group.names
shapirolist2<- as.list(1:no.groups)
names(shapirolist2) <- group.names
levenelist2<- as.list(1:no.groups)
names(levenelist2) <- group.names
# loop over lipid names
for(i in group.names){ 
  
  # print status
  print(paste("Running entity:", i, "which is", which(group.names==i), "out of", no.groups))
  
  # create temporary data matrix and model formula
  tmp <- subset(data_stat[, c(i,"Drought","Future","Time")], Time=="Recovery")
  fml <- as.formula( paste( i, "~", paste(c("Drought","Future"), collapse="*") ) )
  
  #develop formula for levenetest
  lvn <-as.formula( paste( "residua", "~", paste(c("Drought","Future"), collapse="*") ) )
  
  # assign fit to list2 by name
  lme_fit <- lme(fml, random=~1|plot, data=subset(data_stat, Time=="Recovery"))
  residua<-resid(lme_fit)
  anovalist2[[i]] <- anova(lme_fit)
  shapirolist2[[i]] <- shapiro.test(resid(lme_fit))
  levenelist2[[i]] <- leveneTest(lvn, data=subset(data_stat, Time=="Recovery"))
  print(resid_panel(lme(fml, random=~1|plot, data=subset(data_stat, Time=="Recovery"))))
}
fitlist2
anovalist2
shapirolist2
levenelist2
#Not normal or heterogeneous

#fungi (shapiro fail)

#percentofPLFAamf (shapiro fail)
#percentofPLFAgramneg (leven fail)


#run model for non normal or heterogeneus variables
#FB ratio during recovery
m3<-lme(log(percentofPLFA)~Drought*Future,random=~1|plot, data=subset(data_stat, Time=="Drought"))
m3<-lme((percentofPLFA)~Rewet*Fut2,random=~1|plot, data=subset(data_stat, Rewet=="rewet" |Rewet=="not"))
m3<-lme((percentofPLFA)~Drought*Future,random=~1|plot, data=subset(data_stat, Time=="Recovery"))

summary(m3)
anova(m3)
plot(m3,resid(.,scaled=TRUE)~fitted(.),abline=0)
plot(ranef(m3))
qqnorm(m3)
shapiro.test(resid(m3))
leveneTest(resid(m3)~Drought*Future,data=subset(data_stat, Time=="Drought"))
leveneTest(resid(m3)~Rewet*Fut2,data=subset(data_stat, Rewet=="rewet" |Rewet=="not"))
leveneTest(resid(m3)~Drought*Future,data=subset(data_stat, Time=="Recovery"))






#manually summarize the statistic for the main figure
statistic.nlfa <- data.frame(x = 2.5, y =c(3.8,3.5), 
                                Treatment=c("Ambient","Drought", "Future Climate", "Future Climate + Drought"),
                                Time ='Drought', 
                                labl = c("Drought: <.001","","",""))
statistic.percentage <- data.frame(x = 2.5, y =450, 
                                Treatment=c("Ambient","Drought", "Future Climate", "Future Climate + Drought"),
                                Time ='Drought', 
                                labl = c("Drought: <.001","","",""))



color<-c( "#5189b8","#f4e093" , "#943124", "#e77e39")
## create figure 4: NLFA fungi

p1<-ggplot(data, aes(x = Treatment, y=fungi.ms, fill=Treatment)) +
  scale_fill_manual(values=color)+
  geom_boxplot(alpha=0.8) + 
  geom_point(shape = 21,
             color = "black", size = 2)+
  facet_grid(. ~ fct_relevel(Time,'Drought','Rewet','Recovery'), scales = "free", space = "free")+
   theme_bw()+
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x =element_blank(),
    panel.grid = element_blank(),
    plot.margin=unit(c(2,2,15,2), "mm"))+ylim(0,4)+
  ylab(expression(atop("Fungi NLFA production rates", paste(" ("*mg*" " * C*" "*g* " "* C^-1 * h^-1 * ")")))) +
geom_text(mapping= aes(x=x, y=y, label=labl), data=statistic.nlfa)
p1 

#ggsave(file="CUE.png", width=5, height=4)
## PLFA CUE


p2<-ggplot(data, aes(x = Treatment, y=percentofPLFA, fill=Treatment)) +
  scale_fill_manual(values=color)+
  geom_boxplot(alpha=0.8) + 
  geom_point(shape = 21,
             color = "black", size = 2)+
  facet_grid(. ~ fct_relevel(Time,'Drought','Rewet','Recovery'), scales = "free", space = "free")+
  theme_bw()+
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x =element_blank(),
    panel.grid = element_blank(),
    plot.margin=unit(c(2,2,15,2), "mm"))+
  ylab(expression(paste("Fungal NLFA to PLFA ratio (%)")))+
  geom_text(mapping= aes(x=x, y=y, label=labl), data=statistic.percentage)
p2

#ggsave(file="CUE_PLFA.png", width=5, height=4)

#Plot all together
cue<-ggarrange(
  p1, 
  p2 , 
  common.legend = TRUE, legend="bottom",
  ncol=2,
  nrow=1,
  widths=c(1,1),
  labels = c("a)", "b)")  
)
cue
ggsave(plot=cue,"Fig. 4.png", width = 8, height = 4)





#manually summarize the statistic for supplementary figure
statistic.gramneg <- data.frame(x = 2.5, y =c(0.5,0.45), 
                             Treatment=c("Ambient","Drought", "Future Climate", "Future Climate + Drought"),
                             Time ='Drought', 
                             labl = c("Drought: <.01","","",""))
statistic.percgramneg <- data.frame(x = 2.5, y =12, 
                                   Treatment=c("Ambient","Drought", "Future Climate", "Future Climate + Drought"),
                                   Time ='Drought', 
                                   labl = c("Drought: <.001","","",""))



color<-c( "#5189b8","#f4e093" , "#943124", "#e77e39")
## create figure S7: NLFA gram neg

p3<-ggplot(data, aes(x = Treatment, y=gram.neg.ms, fill=Treatment)) +
  scale_fill_manual(values=color)+
  geom_boxplot(alpha=0.8) + 
  geom_point(shape = 21,
             color = "black", size = 2)+
  facet_grid(. ~ fct_relevel(Time,'Drought','Rewet','Recovery'), scales = "free", space = "free")+
  theme_bw()+
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x =element_blank(),
    panel.grid = element_blank(),
    plot.margin=unit(c(2,2,15,2), "mm"))+
  ylab(expression(atop("Gram(-) NLFA production rates", paste(" ("*mg*" " * C*" "*g* " "* C^-1 * h^-1 * ")"))))  +
  geom_text(mapping= aes(x=x, y=y, label=labl), data=statistic.gramneg)
p3 

#ggsave(file="CUE.png", width=5, height=4)
## PLFA CUE


p4<-ggplot(data, aes(x = Treatment, y=percentofPLFAgramneg, fill=Treatment)) +
  scale_fill_manual(values=color)+
  geom_boxplot(alpha=0.8) + 
  geom_point(shape = 21,
             color = "black", size = 2)+
  facet_grid(. ~ fct_relevel(Time,'Drought','Rewet','Recovery'), scales = "free", space = "free")+
  theme_bw()+
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x =element_blank(),
    panel.grid = element_blank(),
    plot.margin=unit(c(2,2,15,2), "mm"))+
  ylab(expression(paste("Gram(-) NLFA as % of PLFA")))+
  geom_text(mapping= aes(x=x, y=y, label=labl), data=statistic.percgramneg)
p4

#ggsave(file="CUE_PLFA.png", width=5, height=4)
#manually summarize the statistic for supplementary figure
statistic.amf <- data.frame(x = 2.5, y =c(0.5,0.45), 
                                Treatment=c("Ambient","Drought", "Future Climate", "Future Climate + Drought"),
                                Time ='Drought', 
                                labl = c("Drought: <.001","","",""))
statistic.percamf <- data.frame(x = 2.5, y =12, 
                                    Treatment=c("Ambient","Drought", "Future Climate", "Future Climate + Drought"),
                                    Time ='Drought', 
                                    labl = c("Drought: <.001","","",""))



color<-c( "#5189b8","#f4e093" , "#943124", "#e77e39")
## create figure S7: NLFA gram neg

p5<-ggplot(data, aes(x = Treatment, y=amf.ms, fill=Treatment)) +
  scale_fill_manual(values=color)+
  geom_boxplot(alpha=0.8) + 
  geom_point(shape = 21,
             color = "black", size = 2)+
  facet_grid(. ~ fct_relevel(Time,'Drought','Rewet','Recovery'), scales = "free", space = "free")+
  theme_bw()+
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x =element_blank(),
    panel.grid = element_blank(),
    plot.margin=unit(c(2,2,15,2), "mm"))+
  ylab(expression(atop("Fungi NLFA production rates", paste(" ("*mg*" " * C*" "*g* " "* C^-1 * h^-1 * ")")))) +
  geom_text(mapping= aes(x=x, y=y, label=labl), data=statistic.amf)
p5

#ggsave(file="CUE.png", width=5, height=4)
## PLFA CUE


p6<-ggplot(data, aes(x = Treatment, y=percentofPLFAamf, fill=Treatment)) +
  scale_fill_manual(values=color)+
  geom_boxplot(alpha=0.8) + 
  geom_point(shape = 21,
             color = "black", size = 2)+
  facet_grid(. ~ fct_relevel(Time,'Drought','Rewet','Recovery'), scales = "free", space = "free")+
  theme_bw()+
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x =element_blank(),
    panel.grid = element_blank(),
    plot.margin=unit(c(2,2,15,2), "mm"))+
  ylab(expression(paste("Fungal NLFA as % of fungal PLFA")))+
  geom_text(mapping= aes(x=x, y=y, label=labl), data=statistic.percamf)
p6


#manually summarize the statistic for supplementary figure
statistic.fung <- data.frame(x = 2.5, y =c(2000,1900), 
                            Treatment=c("Ambient","Drought", "Future Climate", "Future Climate + Drought"),
                            Time ='Drought', 
                            labl = c("Drought: =0.002","Interaction =0.036","",""))
statistic.tot <- data.frame(x = 2.5, y =c(12000,0.45), 
                                Treatment=c("Ambient","Drought", "Future Climate", "Future Climate + Drought"),
                                Time ='Drought', 
                                labl = c("Interaction =0.005","","",""))



color<-c( "#5189b8","#f4e093" , "#943124", "#e77e39")
## create figure S7: NLFA gram neg

p7<-ggplot(data, aes(x = Treatment, y=fungi.mass, fill=Treatment)) +
  scale_fill_manual(values=color)+
  geom_boxplot(alpha=0.8) + 
  geom_point(shape = 21,
             color = "black", size = 2)+
  facet_grid(. ~ fct_relevel(Time,'Drought','Rewet','Recovery'), scales = "free", space = "free")+
  theme_bw()+
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x =element_blank(),
    panel.grid = element_blank(),
    plot.margin=unit(c(2,2,15,2), "mm"))+
  ylab(expression(atop("Fungi NLFA", paste(" ("*noml*" " * C*" "*g* " "* soil^-1* ")")))) +
  geom_text(mapping= aes(x=x, y=y, label=labl), data=statistic.fung)
p7

#ggsave(file="CUE.png", width=5, height=4)
## PLFA CUE


p8<-ggplot(data, aes(x = Treatment, y=tot.mass, fill=Treatment)) +
  scale_fill_manual(values=color)+
  geom_boxplot(alpha=0.8) + 
  geom_point(shape = 21,
             color = "black", size = 2)+
  facet_grid(. ~ fct_relevel(Time,'Drought','Rewet','Recovery'), scales = "free", space = "free")+
  theme_bw()+
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x =element_blank(),
    panel.grid = element_blank(),
    plot.margin=unit(c(2,2,15,2), "mm"))+
  ylab(expression(atop("Total NLFA", paste(" ("*noml*" " * C*" "*g* " "* soil^-1* ")")))) +
  geom_text(mapping= aes(x=x, y=y, label=labl), data=statistic.tot)
p8
#ggsave(file="CUE_PLFA.png", width=5, height=4)

#Plot all together
cue<-ggarrange(
  p3, 
  p4 ,
  p7,
  p8,
  common.legend = TRUE, legend="bottom",
  ncol=2,
  nrow=2,
  widths=c(1,1),
  labels = c("a)", "b)","c)", "d)")  
)
cue
ggsave(plot=cue,"Fig. S11.png", width = 8, height = 8)



