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


#run statistic for each group
library(nlme)
data_stat<-data
drought <- c("Drought","Future Climate + Drought")
notdrought<- c("Ambient","Future Climate")
future <- c("Future Climate","Future Climate + Drought")
ambient <- c("Ambient","Drought")

data_stat$plot<-c("a1","a2","a3","a4","aD1","aD2","aD3","aD4","aF1","aF2","aF3","aF4","aFD1","aFD2","aFD3","aFD4","a1","a2","a3","a4","aD1","aD2","aD3","aD4","aF1","aF2","aF3","aF4","aFD1","aFD2","aFD3","aFD4")

data_stat$Drought[data_stat$Treatment %in% drought] <- 'drought'
data_stat$Drought[data_stat$Treatment %in% notdrought] <- 'not'

data_stat$Future[data_stat$Treatment %in% future] <- 'future'
data_stat$Future[data_stat$Treatment %in% ambient] <- 'ambient'

#run model for each variable
#For example
m3<-lme(log(percentofPLFA)~Drought*Future,random=~1|plot, data=subset(data_stat, Time=="Drought"))

summary(m3)
anova(m3)
plot(m3,resid(.,scaled=TRUE)~fitted(.),abline=0)
plot(ranef(m3))
qqnorm(m3)
shapiro.test(resid(m3))
leveneTest(resid(m3)~Drought*Future,data=subset(data_stat, Time=="Drought"))


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
  facet_grid(. ~ fct_relevel(Time,'Drought','Recovery'), scales = "free", space = "free")+
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



p2<-ggplot(data, aes(x = Treatment, y=percentofPLFA, fill=Treatment)) +
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
  ylab(expression(paste("Fungal NLFA to PLFA ratio (%)")))+
  geom_text(mapping= aes(x=x, y=y, label=labl), data=statistic.percentage)
p2


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
## create figure 

p3<-ggplot(data, aes(x = Treatment, y=gram.neg.ms, fill=Treatment)) +
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
  ylab(expression(atop("Gram(-) NLFA production rates", paste(" ("*mg*" " * C*" "*g* " "* C^-1 * h^-1 * ")"))))  +
  geom_text(mapping= aes(x=x, y=y, label=labl), data=statistic.gramneg)
p3 


p4<-ggplot(data, aes(x = Treatment, y=percentofPLFAgramneg, fill=Treatment)) +
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
  ylab(expression(paste("Gram(-) NLFA as % of PLFA")))+
  geom_text(mapping= aes(x=x, y=y, label=labl), data=statistic.percgramneg)
p4


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


p7<-ggplot(data, aes(x = Treatment, y=fungi.mass, fill=Treatment)) +
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
  ylab(expression(atop("Fungi NLFA", paste(" ("*noml*" " * C*" "*g* " "* soil^-1* ")")))) +
  geom_text(mapping= aes(x=x, y=y, label=labl), data=statistic.fung)
p7



p8<-ggplot(data, aes(x = Treatment, y=tot.mass, fill=Treatment)) +
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
  ylab(expression(atop("Total NLFA", paste(" ("*noml*" " * C*" "*g* " "* soil^-1* ")")))) +
  geom_text(mapping= aes(x=x, y=y, label=labl), data=statistic.tot)
p8

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
ggsave(plot=cue,"Fig. S10.png", width = 8, height = 8)

