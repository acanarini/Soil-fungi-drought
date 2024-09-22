###Load packages
library(ggplot2)
library(ggpubr)
library(reshape2)
library(forcats)
library(tidyverse)
library(ggResidpanel)
library (car)
## Load data

file <- file.choose(new = FALSE)
data <- read.csv(file,header=T,dec=".",sep=",", row.names=1)
data$Treatment <- gsub('Future', 'Future Climate', data$Treatment)

#restructure data set
z<-data[,3:20]
colnames(z) <- c("gram.pos","gram.pos",  "gram.pos",	"gram.neg", "amf",	"general",		"actino", "gram.pos",	"gram.pos",	"gram.neg",	"actino",	 "fungi",	"fungi",	"general",	"general",	"general","actino",	"gram.neg")
data_group<- as.data.frame(t(rowsum(t(z), colnames(z))))
data_group$fbratio<-data_group$fungi/(data_group$gram.pos+data_group$gram.neg+data_group$actino)
data_group$Treatment<-data$Treatment
data_group$Time<-data$Time


#run statistic for each group
library(nlme)
data_stat<-data_group
drought <- c("Drought","Future Climate + Drought")
notdrought<- c("Ambient","Future Climate")
future <- c("Future Climate","Future Climate + Drought")
ambient <- c("Ambient","Drought")

data_stat$plot<-c("a1","a2","a3","a4","aD1","aD2","aD3","aD4","aF1","aF2","aF3","aF4","aFD1","aFD2","aFD3","aFD4","a1","a2","a3","a4","aD1","aD2","aD3","aD4","aF1","aF2","aF3","aF4","aFD1","aFD2","aFD3","aFD4")

data_stat$Drought[data_stat$Treatment %in% drought] <- 'drought'
data_stat$Drought[data_stat$Treatment %in% notdrought] <- 'not'

data_stat$Future[data_stat$Treatment %in% future] <- 'future'
data_stat$Future[data_stat$Treatment %in% ambient] <- 'ambient'


#run model for each variable (transform if not normal or heterogeneous and test)
#example
m3<-lme(fbratio~Drought*Future,random=~1|plot, data=subset(data_stat, Time=="Recovery"))

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
                                labl = c("Drought: <.01","Rewet: <.05","",""))

#Plot for Fig. 3
color<-c( "#5189b8","#f4e093" , "#943124", "#e77e39")

#plot for gram positive
data_group$y<-(data_group$gram.pos+data_group$actino)

p1<-ggplot(data_group, aes(x = Treatment, y=y, fill=Treatment)) +
  scale_fill_manual(values=color)+
  geom_boxplot(alpha=0.8) + 
  geom_point(shape = 21,
             color = "black", size = 2)+
  facet_grid(. ~ fct_relevel(Time,'Drought','Recovery'), scales = "free", space = "free")+
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
  facet_grid(. ~ fct_relevel(Time,'Drought','Recovery'), scales = "free", space = "free")+
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
  facet_grid(. ~ fct_relevel(Time,'Drought','Recovery'), scales = "free", space = "free")+
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
  facet_grid(. ~ fct_relevel(Time,'Drought','Recovery'), scales = "free", space = "free")+
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
ggsave(plot=all,"Fig. S7.png", width = 8, height = 8)


#plot the correlation plots
#load the file with respiration data: Respiration comparison.csv
file <- file.choose(new = FALSE)
resp.data <- read.csv(file,header=T,dec=".",sep=",", row.names=1)
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
ggsave(plot=all,"Fig. S9.png", width = 8, height = 8)



