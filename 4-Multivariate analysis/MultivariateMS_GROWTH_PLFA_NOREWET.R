###Load packages
library(ggplot2)
library(ggpubr)
library(reshape2)
library(forcats)
library(tidyverse)
library(dplyr)
## Load data

file <- file.choose(new = FALSE)
data <- read.csv(file,header=T,dec=".",sep=",", row.names=1)



#first plot data for Drought and Rewet time points
data2<-data[ c(1:16),c(3:20) ]
drought_label<-data[ c(1:16), ]$Treatment

# Transform this data in %
data_percentage_d <- t(apply(data2, 1, function(x){x*100/sum(x,na.rm=T)}))
data_percentage_d<-as.data.frame(data_percentage_d)
colnames(data_percentage_d) <- c( "i15:0","a15:0" ,"i16:0","16:1w7","16:1w5", "16:0" ,"10Me16:0","i17:0", "a17:0","cy17:0", "10Me17:0", "18:2w6:9","cis 18:1w9", "trans 18:1w9", "18:1a","18:0", "10Me18:0","cy19:0")
data_percentage_d$Treatment<-data[ c(1:16),]$Treatment
data_percentage_d$Time<-data[ c(1:16),]$Time

######plot dfunction
element_textbox <- function(...) {
  el <- element_text(...)
  class(el) <- c("element_textbox", class(el))
  el
}

element_grob.element_textbox <- function(element, ...) {
  text_grob <- NextMethod()
  rect_grob <- element_grob(calc_element("strip.background", theme_bw()))
  
  ggplot2:::absoluteGrob(
    grid::gList(
      element_grob(calc_element("strip.background", theme_bw())),
      text_grob
    ),
    height = grid::grobHeight(text_grob), 
    width = grid::unit(1, "npc")
  )
}


##############dca test##############
library(vegan) 
ord <- decorana(data_percentage_d[,c(1:18) ]) #DCA test: The length of the first DCA axis above 4 you can use CCA, below 3 you can use PCA and RDA, between 3 and 4 you can use both
ord

#PCA and plot
library("FactoMineR")
library("corrplot")
library("factoextra")
res.pca <- PCA(data_percentage_d, quali.sup=19:20, graph = FALSE,scale.unit = TRUE)#run PCA
fviz_eig(res.pca, addlabels = TRUE)#show explained variance by axis
var <- get_pca_var(res.pca)#get variables out of PCA
#### create groups for the PLFA peaks
res.km <- kmeans(var$coord, centers = 3, nstart = 25)
grp <- as.character(res.km$cluster)
colors <- as.character(res.km$cluster)
grp[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18)] <- c("gram(+)","gram(+)",  "gram(+)",	"gram(-)", "fungi",	"general",		"actino", "gram(+)",	"gram(+)",	"gram(-)",	"actino",	 "fungi",	"fungi",	"general",	"general",	"general","actino",	"gram(-)")
group.colors <- c("gram(+)" = "#3B8BB3", "gram(-)" = "#3BABB3", "general" ="grey", "fungi" = "#D78962", "actino" = "white")

# Contributions of variables to PC1 and PC2
var1<-as.data.frame(var$contrib)
var2<-as.data.frame(var$coord)
var1$Group<-grp
var1$ID<-(row.names(var1))
var1$DIR1<-var2$Dim.1
var1$DIR2<-var2$Dim.2
var1$PCA1<-(var1$Dim.1*(var1$DIR1/abs(var1$DIR1)))
var1$PCA2<-(var1$Dim.2*(var1$DIR2/abs(var1$DIR2)))


cont1<-var1 %>% 
  top_n(7, Dim.1) %>% 
  mutate(
    ID = factor(ID, levels = ID[order(Dim.1, decreasing = TRUE)]),
    label_hjust = ifelse(PCA1 < 6, -0.01, 1.1)
  ) %>%
  ggplot(., aes(x=reorder(ID, -PCA1), y=PCA1, fill=Group))+
  geom_bar(stat='identity', color="black", alpha = 0.5)+
  geom_text(
    aes(label = ID, hjust = label_hjust)
  )+
  coord_flip()+ scale_x_discrete(limits=rev,position = 'top')+
  scale_fill_manual(values=group.colors)+
  theme(panel.background = element_blank(),
        legend.position = "none",
        axis.text.y = element_blank(),
        axis.ticks.y  = element_blank(),
        plot.title = element_blank(),
        axis.title.y = element_blank(),
        axis.line.x=element_line())+
  ylab("Contribution to PCA 1 (%)")#+labs(tag = "c)")
cont1



# Contributions of variables to PC2
cont2<-var1 %>% 
  top_n(7, Dim.2) %>% 
  mutate(
    ID = factor(ID, levels = ID[order(Dim.2, decreasing = TRUE)]),
    label_hjust = ifelse(PCA2 < 5, -1.1, 1.1)
  ) %>%
  ggplot(., aes(x=reorder(ID, -PCA2), y=PCA2, fill=Group))+
  geom_bar(stat='identity', color="black", alpha = 0.5)+
  geom_text(
    aes(label = ID, hjust = label_hjust)
  )+
  coord_flip()+ scale_x_discrete(limits=rev,position = 'top')+
  scale_fill_manual(values=group.colors)+
  theme(panel.background = element_blank(),
        legend.position = "none",
        axis.text.y = element_blank(),
        axis.ticks.y  = element_blank(),
        plot.title = element_blank(),
        axis.title.y = element_blank(),
        axis.line.x=element_line())+
  ylab("Contribution to PCA 2 (%)")#+labs(tag = "d)")
cont2


#create graph for samples
sam1 <-fviz_pca_ind(res.pca, 
                    # Fill individuals by groups
                    geom.ind = "point",
                    pointshape = 21,
                    pointsize = 3,
                    fill.ind = drought_label,
                    col.ind = "black",
                    addEllipses = TRUE,
                    ellipse.type = "confidence",
                    repel = TRUE,
                    mean.point=FALSE,
                    palette = c("#5189b8","#f4e093" ,  "#943124", "#e77e39"))+      # Indiviual fill color 
  theme_classic()+
  theme(plot.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position = "none")+
  labs(x = "PC1 (43%)",y="PC2 (30.3%)")+ labs(title = "Drought") +
  theme(
    plot.title = element_textbox(
      hjust = 0.5, margin = margin(t = 5, b = 5)
    )
  )

sam1


#second plot data for Recovery time point

data3<-data[ c(25:40),c(3:20) ]
data$Treatment <- gsub(' - Rewet', '', data$Treatment)
data$Treatment <- gsub('Future', 'Future Climate', data$Treatment)

#data3<-data3[ ,-11 ]
recovery_label<-data[ c(25:40), ]$Treatment
# Transform this data in %
data_percentage_r <- t(apply(data3, 1, function(x){x*100/sum(x,na.rm=T)}))
data_percentage_r<-as.data.frame(data_percentage_r)
colnames(data_percentage_r) <- c( "i15:0","a15:0" ,"i16:0","16:1w7","16:1w5", "16:0" ,"10Me16:0","i17:0", "a17:0","cy17:0", "10Me17:0", "18:2w6:9","cis 18:1w9", "trans 18:1w9", "18:1a","18:0", "10Me18:0","cy19:0")
data_percentage_r$Treatment<-data[ c(25:40),]$Treatment
data_percentage_r$Time<-data[ c(25:40),]$Time

##############dca test##############
ord <- decorana(data_percentage_r[,c(1:18) ]) #DCA test: The length of the first DCA axis above 4 you can use CCA, below 3 you can use PCA and RDA, between 3 and 4 you can use both
ord

#PCA and plots
res.pca1 <- PCA(data_percentage_r,quali.sup=19:20, graph = FALSE,scale.unit = TRUE)#run PCA
fviz_eig(res.pca1, addlabels = TRUE)#show explained variance by axis
var4<- get_pca_var(res.pca1)#get variables out of PCA
#### create groups for the PLFA peaks
res.km1 <- kmeans(var4$coord, centers = 3, nstart = 25)
grp1 <- as.character(res.km1$cluster)
grp1[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18)] <- c("gram(+)","gram(+)",  "gram(+)",	"gram(-)", "fungi",	"general",		"actino", "gram(+)",	"gram(+)",	"gram(-)",	"actino",	 "fungi",	"fungi",	"general",	"general",	"general","actino",	"gram(-)")
group.colors <- c("gram(+)" = "#3B8BB3", "gram(-)" = "#3BABB3", "general" ="grey", "fungi" = "#D78962", "actino" = "white")

# Contributions of variables to PC1 and PC2
var5<-as.data.frame(var4$contrib)
var6<-as.data.frame(var4$coord)
var5$Group<-grp
var5$ID<-(row.names(var5))
var5$DIR1<-var6$Dim.1
var5$DIR2<-var6$Dim.2
var5$PCA1<-(var5$Dim.1*(var5$DIR1/abs(var5$DIR1)))
var5$PCA2<-(var5$Dim.2*(var5$DIR2/abs(var5$DIR2)))


cont3<-var5 %>% 
  top_n(7, Dim.1) %>% 
  mutate(
    ID = factor(ID, levels = ID[order(Dim.1, decreasing = TRUE)]),
    label_hjust = ifelse(PCA1 > 1, 1.1, ifelse(PCA1 < -8.9, -0.1, -1))
  ) %>%
  ggplot(., aes(x=reorder(ID, -PCA1), y=PCA1, fill=Group))+
  geom_bar(stat='identity', color="black", alpha = 0.5)+
  geom_text(
    aes(label = ID, hjust = label_hjust)
  )+
  coord_flip()+ scale_x_discrete(limits=rev,position = 'top')+
  scale_fill_manual(values=group.colors)+
  theme(panel.background = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none",
        axis.ticks.y  = element_blank(),
        plot.title = element_blank(),
        axis.title.y = element_blank(),
        axis.line.x=element_line())+
  ylab("Contribution to PCA 1 (%)")#+labs(tag = "e)")
cont3

# Contributions of variables to PC2
cont4<-var5 %>% 
  top_n(7, Dim.2) %>% 
  mutate(
    ID = factor(ID, levels = ID[order(Dim.2, decreasing = TRUE)]),
    label_hjust = ifelse(PCA2 < 8, -0.1, 1.1)
  ) %>%
  ggplot(., aes(x=reorder(ID, -PCA2), y=PCA2, fill=Group))+
  geom_bar(stat='identity', color="black", alpha = 0.5)+
  geom_text(
    aes(label = ID, hjust = label_hjust)
  )+
  coord_flip()+ scale_x_discrete(limits=rev,position = 'top')+
  scale_fill_manual(values=group.colors)+
  theme(panel.background = element_blank(),
        #legend.position = "none",
        axis.text.y = element_blank(),
        axis.ticks.y  = element_blank(),
        plot.title = element_blank(),
        axis.title.y = element_blank(),
        axis.line.x=element_line())+
  ylab("Contribution to PCA 2 (%)")#+labs(tag = "f)")
cont4

#extract legend
leg1 <- get_legend(cont4)
# Convert to a ggplot and print
leg1<-as_ggplot(leg1)
leg1

cont4<-cont4+theme(legend.position = "none")
cont4

#create graph for samples
sam2 <-fviz_pca_ind(res.pca1, 
                    # Fill individuals by groups
                    geom.ind = "point",
                    pointshape = 21,
                    pointsize = 3,
                    fill.ind = recovery_label,
                    col.ind = "black",
                    addEllipses = TRUE,
                    ellipse.type = "confidence",
                    repel = TRUE,
                    mean.point=FALSE,
                    legend.title = "Treatment",
                    palette = c("#5189b8","#f4e093" ,  "#943124", "#e77e39"))+      # Indiviual fill color 
  theme_classic()+
  theme(plot.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))+
  labs(x = "PC1 (53%)",y="PC2 (21.3%)")+ labs(title = "Recovery") +
  theme(
    plot.title = element_textbox(
      hjust = 0.5, margin = margin(t = 5, b = 5)
    )
  )

#extract legend
leg2 <- get_legend(sam2)
# Convert to a ggplot and print
leg2<-as_ggplot(leg2)
leg2
sam2<-sam2+theme(legend.position = "none")
sam2

#plot all together for Figure S6
library(patchwork)

layout <- c(
  area(1, 1,2,2),
  area(3, 1, 3,1),
  area(3, 2, 3, 2),
  area(1, 3,2,4),
  area(3, 3, 3,3),
  area(3, 4, 3, 4),
  area(1, 5,2,5),
  area(3, 5, 3,5)
)

plot(layout)

all<- sam1 + cont1 + cont2+sam2 + cont3 + cont4+ leg2 + leg1+
  plot_layout(design = layout)

all

ggsave(plot=all,"Fig. 2.png", width = 12, height = 6)






#run Permanovas
#Drought and rewet
data_percentage_d$drought_label<-drought_label

drought<-c("Drought","Future + Drought")
ambient<-c("Ambient","Future")
data_percentage_d$Drought[data_percentage_d$drought_label %in% drought] <- 'drought'
data_percentage_d$Drought[data_percentage_d$drought_label %in% ambient] <- 'ambient'

future<-c("Future","Future + Drought")
control<-c("Ambient","Drought")
data_percentage_d$Future[data_percentage_d$drought_label %in% future] <- 'future'
data_percentage_d$Future[data_percentage_d$drought_label %in% control] <- 'control'


model1=adonis2(scale(data_percentage_d[c(1:16),c(1:18) ]) ~ Drought*Future,
               method = "euclidean",  data=data_percentage_d, permutations=999)
summary(model1)
model1


#Recovery scaled data
data_percentage_r$drought_label<-recovery_label

drought<-c("Drought","Future Climate + Drought")
ambient<-c("Ambient","Future Climate")
data_percentage_r$Drought[data_percentage_r$drought_label %in% drought] <- 'drought'
data_percentage_r$Drought[data_percentage_r$drought_label %in% ambient] <- 'ambient'

future<-c("Future Climate","Future Climate + Drought")
control<-c("Ambient","Drought")
data_percentage_r$Future[data_percentage_r$drought_label %in% future] <- 'future'
data_percentage_r$Future[data_percentage_r$drought_label %in% control] <- 'control'



model3=adonis2(scale(data_percentage_r[,c(1:18) ]) ~ Drought*Future,method = "euclidean",  data=data_percentage_r, permutations=999)
summary(model3)
model3







