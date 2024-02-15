library(ggplot2)
library(ggpubr)
library(ggtext)


#load data: meta_data.csv
file <- file.choose(new = FALSE)
data <- read.csv(file,header=T,dec=".",sep=",")

#Data description
#Treatments: A=Ambient; D=Drought; F=Future climate; FD=Future climate and Drought
#time= hours
#values = atom %

#select data sets
data_deut<-subset(data, type=="Deuterium")
data_18O<-subset(data, type=="18O")


##INDIVIDUAL ANALYSIS FOR EACH TREATMENT
##make 4 different graphs 

#subset data
data_A<-subset(data_deut, Treatment=="Ambient")
#get atom% and time
SE<-data_A$time
SA1<-SE
vA1<-data_A$values
#plot
plot (SA1,vA1, 
      xlab="Time (hours)", 
      ylab="Atom % ", 
      main="Ambient", 
      pch=17, col="red")

#model from paper in GCB
MMcurveA1<-formula(vA1~Vmax+((Vin-Vmax)*(exp(-Km*SA1))))
# We can also build a simple data frame from (S,v) data
kinDataA1 <- data.frame(SA1,vA1)

bestfitA1 <- nls(MMcurveA1, kinDataA1, start=list(Vin=70,Vmax=30,Km=0.05))
summary(bestfitA1)
# Build a theoretical line defining the best fit curve
# First, make a finely detailed set of points between 0 and 48, at 0.1 intervals
# These will be the atm% that are used to calculate 
# the predicted atm%
SconcRange <- seq(0,48,0.1)

# Then, calculate the predicted atm% using the predict function
theorLineA1 <- predict(bestfitA1,list(SA1=SconcRange))

# Best fit values of Km and Vmax obtained by coef function, stored in bestFitVals
bestFitValsA1 <- coef(bestfitA1)

# Now plot the data, the best fit line, and put the best fit coefficients in the plot
plot (kinDataA1, 
      xlab="Time (hours)", 
      ylab="Atom %", 
      title(main="Fitted Ambient data"))

# Draw the line
# lines() function adds to the existing plot as does points()
# here, we want a line of best fit, not points, so we use lines()
lines(SconcRange,theorLineA1,col="red")

#####model for A####
#get Vmax and Km
Ain<-bestFitValsA1[1]
Amax<-bestFitValsA1[2]
Ak<-bestFitValsA1[3]

####generate observed model (eq2) and predicted (eq1)####

#A
eq1 <- function(x){Amax  +(0.0156-Amax  )*(exp(- Ak*x))}
eq2 <- function(x){Amax  +(Ain-Amax  )*(exp(- Ak*x))}

#plot graph for A
A = ggplot(data = data_A, aes(x = time, y = values))+
  geom_point(fill = "#5189b8", shape = 21, colour = "black", size = 2, stroke = 1)+
  
  stat_function(fun=eq1, geom="line",colour = "#5189b8", size = 0.5, linetype = 2)+
  
  stat_function(fun=eq2, geom="line",colour = "#5189b8", size = 0.5, linetype = 1)+
  
  theme(
    panel.background = element_rect(fill = "white", colour = "black",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  )+
  labs(title = "Ambient") + ylab(bquote(' '^2*'H atom%'))+xlab('Time (hours)')+
  theme(
    plot.title = element_textbox(
      hjust = 0.5, margin = margin(t = 5, b = 5)
    )
  ) +ylim(0,75)
A








#subset data
data_D<-subset(data_deut, Treatment=="Drought")
#get atom% and time
SE<-data_D$time
SD1<-SE
vD1<-data_D$values
#plot
plot (SD1,vD1, 
      xlab="Time (hours)", 
      ylab="Atom % ", 
      main="Ambient", 
      pch=17, col="red")

#model from paper in GCB
MMcurveD1<-formula(vD1~Vmax+((Vin-Vmax)*(exp(-Km*SD1))))
# We can also build a simple data frame from (S,v) data
kinDataD1 <- data.frame(SD1,vD1)

bestfitD1 <- nls(MMcurveD1, kinDataD1, start=list(Vin=30,Vmax=20,Km=0.01))
summary(bestfitD1)
# Build a theoretical line defining the best fit curve
# First, make a finely detailed set of points between 0 and 48, at 0.1 intervals
# These will be the atm% that are used to calculate 
# the predicted atm%
SconcRange <- seq(0,48,0.1)

# Then, calculate the predicted atm% using the predict function
theorLineD1 <- predict(bestfitD1,list(SD1=SconcRange))

# Best fit values of Km and Vmax obtained by coef function, stored in bestFitVals
bestFitValSD1 <- coef(bestfitD1)

# Now plot the data, the best fit line, and put the best fit coefficients in the plot
plot (kinDataD1, 
      xlab="Time (hours)", 
      ylab="Atom %", 
      title(main="Fitted Ambient data"))

# Draw the line
# lines() function adds to the existing plot as does points()
# here, we want a line of best fit, not points, so we use lines()
lines(SconcRange,theorLineD1,col="red")


#####model for D####
#get Vmax and Km
Din<-bestFitValSD1[1]
Dmax<-bestFitValSD1[2]
Dk<-bestFitValSD1[3]


#D
eq3 <- function(x){Dmax  +(0.0156-Dmax  )*(exp(- Dk*x))}
eq4 <- function(x){Dmax  +(Din-Dmax  )*(exp(- Dk*x))}
#plot graph for D
D = ggplot(data = data_D, aes(x = time, y = values))+
  geom_point(fill = "#f4e093", shape = 21, colour = "black", size = 2, stroke = 1)+
  
  stat_function(fun=eq3, geom="line",colour = "#f4e093", size = 0.5, linetype = 2)+
  
  stat_function(fun=eq4, geom="line",colour = "#f4e093", size = 0.5, linetype = 1)+
  theme(
    panel.background = element_rect(fill = "white", colour = "black",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  )+
  labs(title = "Drought") +ylab(bquote(' '^2*'H atom%'))+xlab('Time (hours)')+
  theme(
    plot.title = element_textbox(
      hjust = 0.5, margin = margin(t = 5, b = 5)
    )
  ) +ylim(0,35)
D











#subset data
data_F<-subset(data_deut, Treatment=="Future Climate")
#get atom% and time
SE<-data_F$time
SF1<-SE
vF1<-data_F$values
#plot
plot (SF1,vF1, 
      xlab="Time (hours)", 
      ylab="Atom % ", 
      main="Ambient", 
      pch=17, col="red")

#model from paper in GCB
MMcurveF1<-formula(vF1~Vmax+((Vin-Vmax)*(exp(-Km*SF1))))
# We can also build a simple data frame from (S,v) data
kinDataF1 <- data.frame(SF1,vF1)

bestfitF1 <- nls(MMcurveF1, kinDataF1, start=list(Vin=30,Vmax=20,Km=0.1))
summary(bestfitF1)
# Build a theoretical line defining the best fit curve
# First, make a finely detailed set of points between 0 and 48, at 0.1 intervals
# These will be the atm% that are used to calculate 
# the predicted atm%
SconcRange <- seq(0,48,0.1)

# Then, calculate the predicted atm% using the predict function
theorLineF1 <- predict(bestfitF1,list(SF1=SconcRange))

# Best fit values of Km and Vmax obtained by coef function, stored in bestFitVals
bestFitValSF1 <- coef(bestfitF1)

# Now plot the data, the best fit line, and put the best fit coefficients in the plot
plot (kinDataF1, 
      xlab="Time (hours)", 
      ylab="Atom %", 
      title(main="Fitted Ambient data"))

# Draw the line
# lines() function adds to the existing plot as does points()
# here, we want a line of best fit, not points, so we use lines()
lines(SconcRange,theorLineF1,col="red")



#####model for F####
#get Vmax and Km
Fin<-bestFitValSF1[1]
Fmax<-bestFitValSF1[2]
Fk<-bestFitValSF1[3]
#F
eq5 <- function(x){Fmax  +(0.0156-Fmax  )*(exp(- Fk*x))}
eq6 <- function(x){Fmax  +(Fin-Fmax  )*(exp(- Fk*x))}
#plot graph for F
F = ggplot(data = data_F, aes(x = time, y = values))+
  geom_point(fill = "#943124", shape = 21, colour = "black", size = 2, stroke = 1)+
  
  stat_function(fun=eq5, geom="line",colour = "#943124", size = 0.5, linetype = 2)+
  
  stat_function(fun=eq6, geom="line",colour = "#943124", size = 0.5, linetype = 1)+
  
  theme(
    panel.background = element_rect(fill = "white", colour = "black",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  )+
  labs(title = "Future Climate") +ylab(bquote(' '^2*'H atom%'))+xlab('Time (hours)')+
  theme(
    plot.title = element_textbox(
      hjust = 0.5, margin = margin(t = 5, b = 5)
    )
  ) +ylim(0,75)
F








#subset data
data_FD<-subset(data_deut, Treatment=="Future Climate + Drought")
#get atom% and time
SE<-data_FD$time
SFD1<-SE
vFD1<-data_FD$values
#plot
plot (SFD1,vFD1, 
      xlab="Time (hours)", 
      ylab="Atom % ", 
      main="Ambient", 
      pch=17, col="red")

#model from paper in GCB
MMcurveFD1<-formula(vFD1~Vmax+((Vin-Vmax)*(exp(-Km*SFD1))))
# We can also build a simple data frame from (S,v) data
kinDataFD1 <- data.frame(SFD1,vFD1)

bestfitFD1 <- nls(MMcurveFD1, kinDataFD1, start=list(Vin=30,Vmax=20,Km=0.01))
summary(bestfitFD1)
# Build a theoretical line defining the best fit curve
# First, make a finely detailed set of points between 0 and 48, at 0.1 intervals
# These will be the atm% that are used to calculate 
# the predicted atm%
SconcRange <- seq(0,48,0.1)

# Then, calculate the predicted atm% using the predict function
theorLineFD1 <- predict(bestfitFD1,list(SFD1=SconcRange))

# Best fit values of Km and Vmax obtained by coef function, stored in bestFitVals
bestFitValSFD1 <- coef(bestfitFD1)

# Now plot the data, the best fit line, and put the best fit coefficients in the plot
plot (kinDataFD1, 
      xlab="Time (hours)", 
      ylab="Atom %", 
      title(main="Fitted Ambient data"))

# Draw the line
# lines() function adds to the existing plot as does points()
# here, we want a line of best fit, not points, so we use lines()
lines(SconcRange,theorLineFD1,col="red")

#####model for FD####
#get Vmax and Km
FDin<-bestFitValSFD1[1]
FDmax<-bestFitValSFD1[2]
FDk<-bestFitValSFD1[3]
#FD
eq7 <- function(x){FDmax  +(0.0156-FDmax  )*(exp(- FDk*x))}
eq8 <- function(x){FDmax  +(FDin-FDmax  )*(exp(- FDk*x))}

#plot graph for FD
FD = ggplot(data = data_FD, aes(x = time, y = values))+
  geom_point(fill = "#e77e39", shape = 21, colour = "black", size = 2, stroke = 1)+
  
  stat_function(fun=eq7, geom="line",colour = "#e77e39", size = 0.5, linetype = 2)+
  
  stat_function(fun=eq8, geom="line",colour = "#e77e39", size = 0.5, linetype = 1)+
  
  theme(
    panel.background = element_rect(fill = "white", colour = "black",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  )+
  labs(title = "Future Climate + Drought") +ylab(bquote(' '^2*'H atom%'))+xlab('Time (hours)')+
  theme(
    plot.title = element_textbox(
      hjust = 0.5, margin = margin(t = 5, b = 5)
    )
  ) +ylim(0,35)
FD





#Plot all together for Fi S1
rest<-ggarrange(
   A, 
   D + 
     theme(axis.title.y = element_blank()), 
   F, 
   FD+ 
     theme(axis.title.y = element_blank()), 
   common.legend = FALSE,
   ncol=2,
   nrow=2,
   heights=c(1,1.25)
   
)
rest
ggsave(plot=rest,"Fig. S2.png", width = 6, height = 4)


#Calculate integrals
###Integral
#get integral of the curve
require(pracma)
valueA<-integrate(eq1, lower = 0, upper = 48)$value
valueD<-integrate(eq3, lower = 0, upper = 48)$value
valueF<-integrate(eq5, lower = 0, upper = 48)$value
valueFD<-integrate(eq7, lower = 0, upper = 48)$value

#get the atm%
totA<-valueA/48
totD<-valueD/48
totF<-valueF/48
totFD<-valueFD/48
totA
totD
totF
totFD

######Calculate atom% after 48hours
SconcRange <- seq(0,48,0.1)

#Treatment A
#model from paper in GCB
yA<-Amax+((0.0156-Amax)*(exp(-Ak*SconcRange)))
yA
yD<-Dmax+((0.0156-Dmax)*(exp(-Dk*SconcRange)))
yD
yF<-Fmax+((0.0156-Fmax)*(exp(-Fk*SconcRange)))
yF
yFD<-FDmax+((0.0156-FDmax)*(exp(-FDk*SconcRange)))
yFD
