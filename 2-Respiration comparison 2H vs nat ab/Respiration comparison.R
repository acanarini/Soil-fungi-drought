library(ggplot2)
library(ggpubr)
library(reshape2)
library(forcats)

#load data: meta_data.csv
file <- file.choose(new = FALSE)
data <- read.csv(file,header=T,dec=".",sep=",")


#create variable for plotting

dat.m <- melt(data,id.vars=c('Treatment','Time'), measure.vars=c('Resp_label','Resp_na'))


p<-ggplot(dat.m, aes(x = Treatment, y=value, fill=variable)) +
  
  geom_boxplot(alpha=0.8) + 
  scale_fill_manual(labels = c("Deuterium", "Natural abundance"),values=c( "#00AFBB", "#E7B800"))+
  geom_point(shape = 21,
             color = "black", size = 3,position = position_jitterdodge())+
  facet_grid(. ~ fct_relevel(Time,'Drought','Rewet'), scales = "free", space = "free")+
  theme_bw()+
  theme(legend.title = element_blank(),
        panel.grid = element_blank(),
        plot.margin=unit(c(2,2,15,2), "mm"))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  theme(legend.position = c(0.20, 0.85), legend.background = element_rect(fill = "white", color = "white"),
        axis.title.x = element_blank()
  )+
ylab(expression(paste("Respiration rates", " ("*ng*C*" " * g*" "*soil^-1 * h^-1 * ")")))
p 

#ggsave(file="Resp rates comparison.png", width=6, height=4)


## Compare respiration of labelled and nat ab with regression

p1<-ggplot(data, aes(x = Resp_na, y=Resp_label)) +
  geom_point(shape = 21,
             color = "black", size = 3, fill="grey")+
  geom_abline(slope = 1,intercept = 0)+
  geom_text(x=2300, y=2700, label = "1:1 Line", color = "black")+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    plot.margin=unit(c(2,2,15,2), "mm")) +
ylab(expression(paste("Respiration rates (Natural abundance)")))+xlab(expression(paste("Respiration rates (Deuterium)")))+
  xlim(0, 3000)+
  ylim(0, 3000)

p1
#ggsave(file="Resp rates comparison line.png", width=6, height=4)

#Plot all together
rest<-ggarrange(
  p, 
  p1 , 
  common.legend = FALSE,
  ncol=2,
  nrow=1,
  widths=c(1,0.8),
  labels = c("a)", "b)")
  
)
rest
ggsave(plot=rest,"Fig. S4.png", width = 10, height = 5)

## Compare respiration of labelled and nat ab pair t test

library("gmodels")

t.test(log(data$Resp_na), log(data$Resp_label), paired = TRUE, alternative = "two.sided")

shapiro.test(log(data$Resp_na)) 
shapiro.test(log(data$Resp_label)) 

########Results
# t.test(log(data$Resp_na), log(data$Resp_label), paired = TRUE, alternative = "two.sided")

#Paired t-test

#data:  log(data$Resp_na) and log(data$Resp_label)
#t = 1.7867, df = 31, p-value = 0.08376
#alternative hypothesis: true mean difference is not equal to 0
#95 percent confidence interval:
#  -0.01613339  0.24421948
#sample estimates:
#  mean difference 
#0.114043 

#> shapiro.test(log(data$Resp_na)) 

#Shapiro-Wilk normality test

#data:  log(data$Resp_na)
#W = 0.97089, p-value = 0.5242

#> shapiro.test(log(data$Resp_label)) 

#Shapiro-Wilk normality test

#data:  log(data$Resp_label)
#W = 0.97365, p-value = 0.6056