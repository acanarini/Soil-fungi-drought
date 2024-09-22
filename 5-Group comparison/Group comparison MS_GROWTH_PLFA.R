###Load packages
library(ggplot2)
library(ggpubr)
library(reshape2)
library(forcats)
library(tidyverse)
library(ggResidpanel)
library (car)
## Load data

file <- file.choose(new = FALSE)#MultivariateGROWTH_PLFA.csv
data <- read.csv(file,header=T,dec=".",sep=",", row.names=1)
data$Treatment <- gsub('Future', 'Future Climate', data$Treatment)

#restructure data set
z<-data[,3:20]
colnames(z) <- c("gram.pos","gram.pos",  "gram.pos",	"gram.neg", "amf",	"general",		"actino", "gram.pos",	"gram.pos",	"gram.neg",	"actino",	 "fungi",	"fungi",	"general",	"general",	"general","actino",	"gram.neg")
data_group<- as.data.frame(t(rowsum(t(z), colnames(z))))
data_group$fbratio<-data_group$fungi/(data_group$gram.pos+data_group$gram.neg+data_group$actino)
data_group$gram.pos<-data_group$gram.pos+data_group$actino
data_group$Treatment<-data$Treatment
data_group$Time<-data$Time
#colnames(data) <- c( "Treatment", "Time", "i15:0","a15:0" ,"i16:0","16:1w7","16:1w5", "16:0" ,"10Me16:0","i17:0", "a17:0","cy17:0", "10Me17:0", "18:2w6:9","cis 18:1w9", "trans 18:1w9", "18:1a","18:0", "10Me18:0","cy19:0")


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

#run model for each variable
#example
m3<-lme(gram.pos~Drought*Future,random=~1|plot, data=subset(data_stat, Time=="Drought"))

summary(m3)
anova(m3)
plot(m3,resid(.,scaled=TRUE)~fitted(.),abline=0)
plot(ranef(m3))
qqnorm(m3)
shapiro.test(resid(m3))
leveneTest(resid(m3)~Drought*Future,data=subset(data_stat, Time=="Drought"))


#manually summarize the statistic for the three plotted groups
statistic.grampos <- data.frame(x = 2.5, y =c(6,5.6), 
                                Treatment=c("Ambient","Drought", "Future Climate", "Future Climate + Drought"),
                                Time ='Drought', 
                                labl = c("Drought: <.001","","",""))
statistic.gramneg <- data.frame(x = 2.5, y =c(2.2,2.05), 
                                Treatment=c("Ambient","Drought", "Future Climate", "Future Climate + Drought"),
                                Time ='Drought', 
                                labl = c("Drought: <.001","","",""))


statistic.fbratio <- data.frame(x = 2.5, y =1, 
                                Treatment=c("Ambient","Drought", "Future Climate", "Future Climate + Drought"),
                                Time ='Drought', 
                                labl = c("Drought: <.001","","",""))

#Plot for Fig. 3
color<-c( "#5189b8","#f4e093" , "#943124", "#e77e39")

#plot for gram positive


p1<-ggplot(data_group, aes(x = Treatment, y=gram.pos, fill=Treatment)) +
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
  ylab(expression(atop("Gram(+) growth rates", paste(" ("*mg*" " *C*" " * g*" "* C^-1 * h^-1 * ")")))) +
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
  ylab(expression(atop("Gram(-) growth rates", paste(" ("*mg*" " *C*" " * g*"  "* C^-1 * h^-1 * ")")))) +
  geom_text(mapping= aes(x=x, y=y, label=labl), data=statistic.gramneg)

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
  ylab(expression(atop("Fungi growth rates", paste(" ("*mg*" " *C*" " * g*"  "* C^-1 * h^-1 * ")"))))
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
  ylab(expression("F:B ratio of growth rates")) +
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
ggsave(plot=all,"Fig. 3.png", width = 8, height = 8)


###############plot the correlation plots
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
  xlab(expression(atop("Gram(+) growth rates", paste(" ("*mgC*" " * g*"  "* C^-1 * h^-1 * ")"))))
p5<-plot_list$gram.neg+
  xlab(expression(atop("Gram(-) growth rates", paste(" ("*mgC*" " * g*"  "* C^-1 * h^-1 * ")"))))
p6<-plot_list$fungi+
  xlab(expression(atop("Fungi growth rates", paste(" ("*mgC*" " * g*"  "* C^-1 * h^-1 * ")"))))
p8<-plot_list$fbratio+
  xlab(expression("F:B ratio"))
p9<-plot_list$actino+
  xlab(expression(atop("Actinobacteria growth rates", paste(" ("*mgC*" " * g*"  "* C^-1 * h^-1 * ")"))))
p10<-plot_list$general+
  xlab(expression(atop("General growth rates", paste(" ("*mgC*" " * g*"  "* C^-1 * h^-1 * ")"))))

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
ggsave(plot=all,"Fig. S6.png", width = 8, height = 8)


color<-c( "#5189b8","#f4e093" , "#943124", "#e77e39")

#plot for actino

f1<-ggplot(data_group, aes(x = Treatment, y=actino, fill=Treatment)) +
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
  ylab(expression(atop("Actinobacteria growth rates", paste(" ("*mgC*" " * g*"  "* C^-1 * h^-1 * ")"))))
f1 

f2<-ggplot(data_group, aes(x = Treatment, y=general, fill=Treatment)) +
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
  ylab(expression(atop("General biomarkers growth rates", paste(" ("*mgC*" " * g*"  "* C^-1 * h^-1 * ")"))))
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
ggsave(plot=al,"Fig. S5.png", width = 8, height = 4)


######PLOT ABUNDANCE VS GROWTH#####
## Load data

file <- file.choose(new = FALSE)#MultivariateABUNDANCE_PLFA.csv
data_ab <- read.csv(file,header=T,dec=".",sep=",", row.names=1)
data_ab$Treatment <- gsub('Future', 'Future Climate', data_ab$Treatment)
colnames(data_ab) <- c( "Treatment", "Time", "i15:0","a15:0" ,"i16:0","16:1w7","16:1w5", "16:0" ,"10Me16:0","i17:0", "a17:0","cy17:0", "10Me17:0", "18:2w6:9","cis 18:1w9", "trans 18:1w9", "18:1a","18:0", "10Me18:0","cy19:0")

library(tidyverse)

# Condense multiple columns into one column called "values"
df_long <- data %>%
  pivot_longer(cols = -c(Treatment, Time), # Exclude the Geneid column from being pivoted
               names_to = "biomarker",
               values_to = "values")

# View the transformed dataframe
print(df_long)

# Condense multiple columns into one column called "values"
df_long_ab <- data_ab %>%
  pivot_longer(cols = -c(Treatment, Time), # Exclude the Geneid column from being pivoted
               names_to = "biomarker",
               values_to = "values")

# View the transformed dataframe
print(df_long_ab)

df_long$abundance<-df_long_ab$values



### summarize without treatment
# Calculate mean and standard error for values and abundance
data_summary <- df_long %>%
  group_by(Time, biomarker) %>%
  summarize(
    mean_values = mean(values),
    se_values = sd(values) / sqrt(n()),
    mean_abundance = mean(abundance),
    se_abundance = sd(abundance) / sqrt(n())
  )


# Order biomarkers by highest abundance
biomarker_order <- data_summary %>%
  group_by(biomarker) %>%
  summarize(max_abundance = max(mean_abundance)) %>%
  arrange(desc(max_abundance)) %>%
  pull(biomarker)

# Convert biomarkers to a factor with levels in the desired order
data_summary$biomarker <- factor(data_summary$biomarker, levels = biomarker_order)

# Create the plot
p1 <- ggplot(data_summary, aes(x = biomarker, group = Time)) +
  geom_point(aes(y = mean_values, color = "PLFA production rates"), position = position_dodge(width = 0.75)) +
  geom_errorbar(aes(ymin = mean_values - se_values, ymax = mean_values + se_values, color = "PLFA production rates"), width = 0.2, position = position_dodge(width = 0.75)) +
  geom_line(aes(y = mean_values, color = "PLFA production rates"), position = position_dodge(width = 0.75), linetype = "dashed") +
  geom_point(aes(y = mean_abundance / 300, color = "PLFA abundance"), position = position_dodge(width = 0.75)) +
  geom_errorbar(aes(ymin = (mean_abundance - se_abundance) / 300, ymax = (mean_abundance + se_abundance) / 300, color = "PLFA abundance"), width = 0.2, position = position_dodge(width = 0.75)) +
  geom_line(aes(y = mean_abundance / 300, color = "PLFA abundance"), position = position_dodge(width = 0.75), linetype = "dashed") +
  scale_y_continuous(
    name = expression(atop("Production rates", paste(" ("*mg*" " *C*" " * g*" "* C^-1 * h^-1 * ")"))),
    sec.axis = sec_axis(~ . * 300, name = expression(atop("Abundance", paste(" ("*nmol*" " *C*" " * g*" "* soil^-1 * ")"))))
  ) +
  scale_color_manual(
    name = "Legend",
    values = c("PLFA production rates" = "red", "PLFA abundance" = "grey"),
    breaks = c("PLFA production rates", "PLFA abundance"),  # Ensure the correct order
    labels = c("PLFA production rates", "PLFA abundance")
  ) +
  facet_wrap(~ Time) +
  labs(x = "Biomarker", color = "Legend") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")  # Rotate x-axis labels for better readability

# Display the plot
print(p1)


ggsave(plot=p1,"Fig. S4.png", width = 8, height = 4)


