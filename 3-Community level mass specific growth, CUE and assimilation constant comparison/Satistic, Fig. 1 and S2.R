###Load packages
library(ggplot2)
library(ggpubr)
library(reshape2)
library(forcats)

######## Load data

file <- file.choose(new = FALSE)#Growth and respiration.csv
cue <- read.csv(file,header=T,dec=".",sep=",", row.names=1)


############run statistic for CUE
library(nlme)
data_stat<-cue

drought <- c("Drought","Future Climate + Drought")
notdrought<- c("Ambient","Future Climate")
future <- c("Future Climate","Future Climate + Drought")
ambient <- c("Ambient","Drought")

data_stat$varia<-row.names(data_stat)

data_stat$plot<-c("a1","a2","a3","a4","aD1","aD2","aD3","aD4","aF1","aF2","aF3","aF4","aFD1","aFD2","aFD3","aFD4","a1","a2","a3","a4","aD1","aD2","aD3","aD4","aF1","aF2","aF3","aF4","aFD1","aFD2","aFD3","aFD4")

data_stat$Drought[data_stat$Treatment %in% drought] <- 'drought'
data_stat$Drought[data_stat$Treatment %in% notdrought] <- 'not'

data_stat$Future[data_stat$Treatment %in% future] <- 'future'
data_stat$Future[data_stat$Treatment %in% ambient] <- 'ambient'

#TEST CUE
data_stat$y<-data_stat$Msgrowth #select the variable to test

#run model for DROUGHT
m1<-lme(y~Drought*Future,random=~1|plot, data=subset(data_stat, Time=="Drought"))

summary(m1)
anova(m1)

plot(m1,resid(.,scaled=TRUE)~fitted(.),abline=0)
plot(ranef(m1))
qqnorm(m1)
shapiro.test(resid(m1))
leveneTest(resid(m1)~Drought*Future,data=subset(data_stat, Time=="Drought"))

#run model for RECOVERY
m3<-lme(y~Drought*Future,random=~1|plot, data=subset(data_stat, Time=="Recovery"))

summary(m3)
anova(m3)
plot(m3,resid(.,scaled=TRUE)~fitted(.),abline=0)
plot(ranef(m3))
qqnorm(m3)
shapiro.test(resid(m3))
leveneTest(resid(m3)~Drought*Future,data=subset(data_stat, Time=="Recovery"))


#####PLOT Fig. 1

statistic <- data.frame(x = 2.5, y =0.02, Treatment=c("Ambient","Drought", "Future Climate", "Future Climate + Drought"),Time ='Drought', labl = c("Drought:=0.002","","",""))

color<-c( "#5189b8","#f4e093" , "#943124", "#e77e39")
## create figure 1: Growth, respiration, CUE
##  PLFA Growth

p1<-ggplot(cue, aes(x = Treatment, y=PLFA_day, fill=Treatment)) +
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
  ylab(expression(atop("Soil microbial growth rates", paste(" (" * day^-1 *")")))) +
  geom_text(mapping= aes(x=x, y=y, label=labl), data=statistic)+theme(legend.position = "none")
p1 

## PLFA respiration
statistic <- data.frame(x = 2.5, y =5.5, Treatment=c("Ambient","Drought", "Future Climate", "Future Climate + Drought"),Time ='Drought', labl = c("Drought: =0.016","","",""))
statistic2 <- data.frame(x = 2.5, y =5.5, Treatment=c("Ambient","Drought", "Future Climate", "Future Climate + Drought"),Time ='Recovery', labl = c("Future: <.001","","",""))

p2<-ggplot(cue, aes(x = Treatment, y=Msresp_PLFA, fill=Treatment)) +
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
    plot.margin=unit(c(2,2,15,2), "mm"))+ylim(0,6)+
  ylab(expression(atop("Mass specific respiration rates", paste(" ("*mg*" "*C*" " * g*" "*mic* " "* C^-1 * h^-1 * ")")))) +
  geom_text(mapping= aes(x=x, y=y, label=labl), data=statistic)+
  geom_text(mapping= aes(x=x, y=y, label=labl), data=statistic2)+
  theme(legend.position = "none")
p2



## PLFA CUE
statistic <- data.frame(x = 2.5, y = 0.35, Treatment=c("Ambient","Drought", "Future Climate", "Future Climate + Drought"),Time ='Recovery', labl = c("","","",'Future: <.001'))
p3<-ggplot(cue, aes(x = Treatment, y=CUE_PLFA, fill=Treatment)) +
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
    plot.margin=unit(c(2,2,15,2), "mm"))+ylim(0,0.4)+
  ylab(expression(atop("      ", paste(CUE)))) +
  geom_text(mapping= aes(x=x, y=y, label=labl), data=statistic)
p3 

#extract legend
leg1 <- get_legend(p3)
# Convert to a ggplot and print
leg1<-as_ggplot(leg1)
leg1

p3<-p3+theme(legend.position = "none")
p3


#plot all together for Figure 1
plot<-ggarrange(
  p1, 
  p2 ,
  p3,
  common.legend = TRUE, 
  legend="bottom",
  ncol=3,
  nrow=1,
  widths=c(1,1),
  labels = c("a)", "b)","c)","")  
)
plot
ggsave(plot=plot,"Fig. 1.png", width = 10, height = 4)




#####Fig. S2
statistic <- data.frame(x = 2.5, y =0.04, Treatment=c("Ambient","Drought", "Future Climate", "Future Climate + Drought"),Time ='Drought', labl = c("Drought:<.001","","",""))

p1<-ggplot(cue, aes(x = Treatment, y=DNA_day, fill=Treatment)) +
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
  ylab(expression(atop("Soil microbial growth rates", paste(" ("* day^-1 *  ")")))) +
  geom_text(mapping= aes(x=x, y=y, label=labl), data=statistic)+theme(legend.position = "none")
p1



statistic <- data.frame(x = 2.5, y =7.5, Treatment=c("Ambient","Drought", "Future Climate", "Future Climate + Drought"),Time ='Drought', labl = c("Drought: =0.004","","",""))
statistic2 <- data.frame(x = 2.5, y =7.5, Treatment=c("Ambient","Drought", "Future Climate", "Future Climate + Drought"),Time ='Recovery', labl = c("Future: <.001","","",""))
p2<-ggplot(cue, aes(x = Treatment, y=Msresp, fill=Treatment)) +
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
    plot.margin=unit(c(2,2,15,2), "mm")) +ylim(0,8)+
  ylab(expression(atop("Mass specific respiration rates", paste(" ("*mg*" "*C*" " * g*" "*mic* " "* C^-1 * h^-1 * ")")))) +
  geom_text(mapping= aes(x=x, y=y, label=labl), data=statistic)+
  geom_text(mapping= aes(x=x, y=y, label=labl), data=statistic2)+theme(legend.position = "none")
p2


statistic2 <- data.frame(x = 2.5, y =0.48, Treatment=c("Ambient","Drought", "Future Climate", "Future Climate + Drought"),Time ='Recovery', labl = c("Future: <.001","","",""))
p3<-ggplot(cue, aes(x = Treatment, y=CUE, fill=Treatment)) +
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
    plot.margin=unit(c(2,2,15,2), "mm"))+ylim(0,0.5)+
  ylab(bquote(CUE))+
  geom_text(mapping= aes(x=x, y=y, label=labl), data=statistic2)
p3



## Compare 18O and PLFA CUE



p4<-ggplot(cue, aes(x = CUE, y=CUE_PLFA, fill=Treatment)) +
  scale_fill_manual(values=color)+
  geom_point(shape = 21,
             color = "black", size = 3)+
  geom_smooth(method=lm, na.rm = TRUE, fullrange= TRUE,
              aes(group=1),colour="black")+
  stat_cor(aes(group=1),method = "pearson", label.x = 0.1, label.y = 0.4)+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    plot.margin=unit(c(2,2,15,2), "mm"))+
  ylab(bquote(CUE[PLFAs]))+xlab(bquote(CUE[DNA]))+theme(legend.position = "none")

p4


## Compare 18O and PLFA MSgrowth



p5<-ggplot(cue, aes(x = Msgrowth, y=Msgrowth_PLFA, fill=Treatment)) +
  scale_fill_manual(values=color)+
  geom_point(shape = 21,
             color = "black", size = 3)+
  geom_smooth(method=lm, na.rm = TRUE, fullrange= TRUE,
              aes(group=1),colour="black")+
  stat_cor(aes(group=1),method = "pearson", label.x = 0.3, label.y = 1)+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    plot.margin=unit(c(2,2,15,2), "mm"))+
  ylab(expression(atop("PLFA growth rates", paste(" ("*mg*" "*C*" " * g*" "*mic* " "* C^-1 * h^-1 * ")")))) +
  xlab(expression(paste("DNA growth rates", " ("*mg*" "*C*" " * g*" "*mic* " "* C^-1 * h^-1 * ")")))+theme(legend.position = "none")
p5



statistic <- data.frame(x = 2.5, y =0.9, Treatment=c("Ambient","Drought", "Future Climate", "Future Climate + Drought"),Time ='Drought', labl = c("Drought:=0.002","","",""))


color<-c( "#5189b8","#f4e093" , "#943124", "#e77e39")

##  PLFA Growth

p6<-ggplot(cue, aes(x = Treatment, y=Msgrowth_PLFA, fill=Treatment)) +
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
    plot.margin=unit(c(2,2,15,2), "mm"))+ylim(0,1)+
  ylab(expression(atop("Mass specific growth rates PLFA", paste(" ("*mg*" "*C*" " * g*" "*mic* " "* C^-1 * h^-1 * ")")))) +
  geom_text(mapping= aes(x=x, y=y, label=labl), data=statistic)+theme(legend.position = "none")
p6 


statistic <- data.frame(x = 2.5, y =2, Treatment=c("Ambient","Drought", "Future Climate", "Future Climate + Drought"),Time ='Drought', labl = c("Drought:<.001","","",""))

p7<-ggplot(cue, aes(x = Treatment, y=Msgrowth, fill=Treatment)) +
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
    plot.margin=unit(c(2,2,15,2), "mm"))+ylim(0,2)+
  ylab(expression(atop("Mass specific growth rates DNA", paste(" ("*mg*" "*C*" " * g*" "*mic* " "* C^-1 * h^-1 * ")")))) +
  geom_text(mapping= aes(x=x, y=y, label=labl), data=statistic)+theme(legend.position = "none")
p7


#extract legend
leg1 <- get_legend(p3)
# Convert to a ggplot and print
leg1<-as_ggplot(leg1)
leg1

p3<-p3+theme(legend.position = "none")
p3


#plot all together for Figure S5
plot1<-ggarrange(
  p1, 
  p2 ,
  p3,
  p4,
  p5,
  p6,
  p7,
  leg1,
  common.legend = FALSE, 
  ncol=3,
  nrow=3,
  widths=c(1,1),
  labels = c("a)", "b)","c)","d)","e)","f)","g)","")  
)
plot1
ggsave(plot=plot1,"Fig. S5.png", width = 12, height = 10)
