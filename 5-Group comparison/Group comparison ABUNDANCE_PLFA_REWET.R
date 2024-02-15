###Load packages
library(ggplot2)
library(ggpubr)
library(reshape2)
library(forcats)
library(tidyverse)
library(ggResidpanel)
library (car)
## Load data

file <- file.choose(new = FALSE)#MultivariateABUNDANCE_PLFA.csv
data <- read.csv(file,header=T,dec=".",sep=",", row.names=1)
data$Treatment <- gsub(' - Rewet', '', data$Treatment)
data$Treatment <- gsub('Future', 'Future Climate', data$Treatment)

#restructure data set
z<-data[,3:20]
colnames(z) <- c("gram.pos","gram.pos",  "gram.pos",	"gram.neg", "amf",	"general",		"actino", "gram.pos",	"gram.pos",	"gram.neg",	"actino",	 "fungi",	"fungi",	"general",	"general",	"general","actino",	"gram.neg")
data_group<- as.data.frame(t(rowsum(t(z), colnames(z))))
data_group$fbratio<-data_group$fungi/(data_group$gram.pos+data_group$gram.neg+data_group$actino)
data_group$Treatment<-data$Treatment
data_group$Time<-data$Time
data_group<-subset(data_group, Time!="Rewet")


#run statistic for each group
library(nlme)
#column<-sub(" rewet.*", "", data_ind[,7])
data_stat<-data_group
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
group.names <- colnames(data_stat)[1:7]
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
# loop over lipid names
for(i in group.names){ 
  
  # print status
  print(paste("Running entity:", i, "which is", which(group.names==i), "out of", no.groups))
  
  # create temporary data matrix and model formula
  tmp <- subset(data_stat[, c(i,"Drought","Future","Time")], Time=="Drought")
  
  #start with untransformed variables and use log if data is not normally distributed
  #fml <- as.formula( paste("log(", i,")", "~", paste(c("Drought","Future"), collapse="*") ) )
  fml <- as.formula( paste("log(", i,")", "~", paste(c("Drought","Future"), collapse="*") ) )#use if not normally distribute
  
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
  anovalist1[[i]] <- anova(lme_fit)
  shapirolist1[[i]] <- shapiro.test(resid(lme_fit))
  levenelist1[[i]] <- leveneTest(lvn, data=subset(data_stat, Time=="Drought"))
  print(resid_panel(lme(fml, random=~1|plot, data=subset(data_stat, Rewet=="rewet" |Rewet=="not"))))
}
fitlist1
anovalist1
shapirolist1
levenelist1



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

#run model for non normal or heterogeneus variables
#FB ratio during recovery
m3<-lme(log(fbratio)~Drought*Future,random=~1|plot, data=subset(data_stat, Time=="Recovery"))

summary(m3)
anova(m3)
plot(m3,resid(.,scaled=TRUE)~fitted(.),abline=0)
plot(ranef(m3))
qqnorm(m3)
shapiro.test(resid(m3))


#manually summarize the statistic for the three plotted groups
statistic.grampos <- data.frame(x = 2.5, y =1250, 
                        Treatment=c("Ambient","Drought", "Future Climate", "Future Climate + Drought"),
                        Time ='Recovery', 
                        labl = c("Future: <.05","","",""))
statistic.fbratio <- data.frame(x = 2.5, y =c(0.52,0.5),
                             Treatment=c("Ambient","Drought", "Future Climate", "Future Climate + Drought"),
                             Time ='Drought', 
                             labl = c("Drought: <.01","","",""))

#Plot for Fig. 3
color<-c( "#5189b8","#f4e093" , "#943124", "#e77e39")

#plot for gram positive
data_group$y<-(data_group$gram.pos+data_group$actino)

p1<-ggplot(data_group, aes(x = Treatment, y=y, fill=Treatment)) +
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
  ylab(expression(atop("Gram(+) abundance", paste(" ("*nmol*" " *C*" " * g*" "* soil^-1 * ")")))) +
  geom_text(mapping= aes(x=x, y=y, label=labl), data=statistic.grampos)
p1 



#plot for gram negative


p2<-ggplot(data_group, aes(x = Treatment, y=gram.neg, fill=Treatment)) +
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
  ylab(expression(atop("Gram(-) abundance", paste(" ("*nmol*" " *C*" " * g*" "* soil^-1 * ")"))))  

p2 


#plot for fungi


p3<-ggplot(data_group, aes(x = Treatment, y=fungi, fill=Treatment)) +
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
  ylab(expression(atop("Fungi abundance", paste(" ("*nmol*" " *C*" " * g*" "* soil^-1 * ")"))))  
p3

p7<-ggplot(data_group, aes(x = Treatment, y=fbratio, fill=Treatment)) +
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
  ylab(expression("F:B ratio")) +
  geom_text(mapping= aes(x=x, y=y, label=labl), data=statistic.fbratio)
p7


#Plot all together
all<-ggarrange(
  p1, 
  p2,
  p3 , 
  p7,
  common.legend = TRUE, legend="bottom",
  ncol=2,
  nrow=2,
  widths = c(1,1),
  labels = c("a)", "b)", "c)","d)")
)
all
ggsave(plot=all,"Fig. S9.png", width = 8, height = 8)


#plot the correlation plots
#load the file with respiration data: Respiration comparison.csv
file <- file.choose(new = FALSE)
resp.data <- read.csv(file,header=T,dec=".",sep=",", row.names=1)
resp.data<-subset(resp.data, Time!="Rewet")
data_group$resp<-resp.data$Resp_label



p.names <- colnames(data_stat)[c(1,3:7)]
plot_list = list()
for (i in p.names) {
  p = ggplot(data_group, aes_string(x=i, y=data_group$resp)) +
   
    geom_point(aes( fill=Treatment),shape = 21,
               color = "black", size = 3)+
    scale_fill_manual(values=color)+
    geom_smooth(method=lm, na.rm = TRUE, fullrange= TRUE,
                aes(group=1),colour="black")+
    stat_cor(aes(group=1),method = "pearson", label.x = 0.1, label.y = 2500)+
    theme_bw()+
    theme(
      panel.grid = element_blank(),
      plot.margin=unit(c(2,2,15,2), "mm"))+
    ylab(expression(atop("Respiration rates", paste(" ("*ng*C*" " * g*" "*soil^-1 * h^-1 * ")"))))
  plot_list[[i]] = p
}

p4<-plot_list$gram.pos+
  xlab(expression(atop("Gram(+)  abundance", paste(" ("*nmol*" " *C*" " * g*" "* soil^-1 * ")"))))
p5<-plot_list$gram.neg+
  xlab(expression(atop("Gram(-)  abundance", paste(" ("*nmol*" " *C*" " * g*" "* soil^-1 * ")"))))
p6<-plot_list$fungi+
  xlab(expression(atop("Fungi  abundance", paste(" ("*nmol*" " *C*" " * g*" "* soil^-1 * ")"))))
p8<-plot_list$fbratio+
  xlab(expression("F:B ratio"))
p9<-plot_list$actino+
  xlab(expression(atop("Actinobacteria  abundance", paste(" ("*nmol*" " *C*" " * g*" "* soil^-1 * ")"))))
p10<-plot_list$general+
  xlab(expression(atop("General  abundance", paste(" ("*nmol*" " *C*" " * g*" "* soil^-1 * ")"))))

p4
p5
p6
p8
p9
p10

#Plot all together
all<-ggarrange(
  p4, 
  p5,
  p6 , 
  p8, 
  p9,
  p10,common.legend = TRUE, legend="bottom",
  ncol=3,
  nrow=2,
  widths = c(1,1),
  labels = c("a)", "b)", "c)","d)", "e)","f)")
)
all
ggsave(plot=all,"Fig. S10.png", width = 8, height = 8)





statistic.actino <- data.frame(x = 2.5, y =c(3,2.8), 
                                Treatment=c("Ambient","Drought", "Future Climate", "Future Climate + Drought"),
                                Time ='Drought', 
                                labl = c("Drought: <.001","","",""))
#Plot for Supplementary Fig
color<-c( "#5189b8","#f4e093" , "#943124", "#e77e39")

#plot for gram positive

f1<-ggplot(data_group, aes(x = Treatment, y=actino, fill=Treatment)) +
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
  ylab(expression(atop("Actinobacteria growth rates", paste(" ("*mgC*" " * g*" "*mic* " "* C^-1 * h^-1 * ")"))))
f1 

f2<-ggplot(data_group, aes(x = Treatment, y=general, fill=Treatment)) +
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
  ylab(expression(atop("General biomarkers growth rates", paste(" ("*mgC*" " * g*" "*mic* " "* C^-1 * h^-1 * ")"))))
f2 


#Plot all together
al<-ggarrange(
  f1, 
  f2,
  common.legend = TRUE, legend="bottom",
  ncol=2,
  nrow=1,
  widths = c(1,1),
  labels = c("a)", "b)")
)
al
ggsave(plot=al,"Fig. S7.png", width = 8, height = 4)

