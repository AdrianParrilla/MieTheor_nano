library(ggplot2)
library(readxl)
library(flowCore)
library(RColorBrewer)
library(viridis)

# Calibration plot for the beads and EV model with 3 different refractive indexes

mie <- read_excel("path/to/data/EV_results.xlsx")  
# Select Ev columns and transform intensitites to log
mieEV <- mie[,c(1,2,10,11)]
colnames(mieEV) <- c("size","mietheory","miemodel","Evmodel")
mieEV$mietheory <- log10(mieEV$mietheory)
mieEV <- as.data.frame(mieEV)

# Select beads columns 
beads <- mie[c(8:13),c(5,7)]
colnames(beads) <- c("size","log10median")
beads$size<- as.numeric(beads$size)
beads$log10median<- as.numeric(beads$log10median)
beads <- as.data.frame(beads)

plt <- ggplot(mieEV,aes(x= size)) 

plt + geom_line(aes(y= miemodel, color= "Mie theory"))+
  geom_jitter(data = beads, aes(x= size, y= log10median, color = "Gigamix beads"),size= 2)+
  geom_line(aes(y= Evmodel, color= "EV model"))+
  scale_x_continuous("Particle size (nm)", limits = c(0,2000),breaks = seq(0,2000,200))+
  scale_y_continuous("Scatter intensity (a.u.)" , limits = c(1,7), breaks = c(1:7),labels = scales::math_format(10^.x),
                     sec.axis = sec_axis(~., name="Scattering Cross-Section (nm^2)", breaks = c(1:7),labels = scales::math_format(10^.x)))+
  annotation_logticks(sides = "lr",long = unit(2,"mm"),mid = unit(1.5,"mm"),short = unit(.5,"mm"))+
  ggtitle("EV scattering model")+
  scale_color_manual(values = c("Gigamix beads"="red","Mie theory"="black", "EV model" = "chartreuse4"))+
  theme_classic()+
  theme(
    panel.grid.major = element_line(color = "gray89", size = 0.5),
    panel.grid.minor.y  = element_line(color = "gray89", size = 0.25),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5,size = 11),
    legend.background = element_blank(),
    legend.position = "none",
    #legend.margin = margin(2, 5, 2, 5),
    axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0),size = 11),
    axis.title.x =  element_text(margin = margin(t = 10, r = 0, b = 0, l = 0),size = 11),
    axis.title.y.right = element_text(margin = margin(t = 0, r = 0, b = 0, l = 15)))

# regresion plot

reg <- mie[c(_:_),c(_,_)]
colnames(reg) <- c("logbeads","logmodel")
reg$logbeads <- as.numeric(reg$logbeads)
reg$logmodel <- as.numeric(reg$logmodel)
reg <- as.data.frame(reg)

plot <- ggplot(reg, aes(x = logmodel, y = logbeads))
plot + geom_line()+geom_point()
