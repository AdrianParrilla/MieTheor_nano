library(ggplot2)
library(readxl)
library(flowCore)
library(RColorBrewer)
library(viridis)

#-------------------EV model------------
# Select EV columns and transform intensitites to log
mie <- read_excel("path/to/data/EV_results.xlsx")  
mieEV <- mie[,c(1,2,10,11)]
mieEV[5] <- round(mieEV$`Log EV *F`,2)
colnames(mieEV) <- c("size","mietheory","miemodel","Evmodel","round Ev model")
mieEV <- as.data.frame(mieEV)

# Open fcs file, select desired columns and log transform
fcs <-  read.FCS("B03 GFP + BAF.fcs",truncate_max_range = FALSE)  #change working directory if needed
fcs <- as.data.frame( exprs(fcs) )
fcs_sub <- fcs[,c(5,8)]
fcs_sub[3] <- log10(fcs_sub$`SSC-H`)
fcs_sub[4] <- log10(fcs_sub$`GFP-H`)
colnames(fcs_sub) <- c("SSCH","GFP-H","log SSCh","log GFP-H")
fcs_sub[is.na(fcs_sub)] <- 0

# Plot events
dotplot <- ggplot(fcs_sub, aes(x=`log GFP-H`, y=`log SSCh`))
dotplot + geom_bin2d(bins = 300) +
  scale_fill_viridis_c(option = "turbo",begin = 0.03,end = 1)+ 
  scale_x_continuous("GFP-H",limits = c(1,5),labels = scales::math_format(10^.x))+
  scale_y_continuous("SSCH",limits = c(3,6),labels = scales::math_format(10^.x)) +
  annotation_logticks(sides = "bl",long = unit(2,"mm"),mid = unit(1.5,"mm"),short = unit(.5,"mm"))+
  theme_classic()+
  theme(legend.position = "none",
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) 

sschisto <- ggplot(fcs_sub, aes(x=`log SSCh`))
sschisto + geom_area(stat = "bin",bins = 200,color="black",fill="gray78")+
  scale_x_continuous("Log SSC-H",limits = c(2,6)) +
  theme_classic()

# get the size and side scatter columns
sct <- data.frame(fcs_sub[,c(1,3)])
sct[3] <- round(sct[2],2)
sct[4] <- fcs_sub$`log GFP-H`
colnames(sct) <- c("sct","log sct", "round log sct","log GFP")


#Looks what values in the EV model are equal to the value in the row j of the SSC-H list
evsizedistrib <- c()  
scts <- c()
GFP <- c()

for (j in 1: length(sct$`round log sct`) ) {
  
  res <- mieEV[which(mieEV$`round Ev model` %in% sct$`round log sct`[j]),1]  
  
  if (length(res) == 0){
    next               # if the value is not in the curve, skip it anlend continue.
  }
  else if (length(res) == 1){
    evsizedistrib <- c(evsizedistrib,res)
    scts <- c(scts, sct$`log sct`[j])
    GFP <- c(GFP, sct$`log GFP`[j])
  }
  else if (length(res) > 1) {
    evsizedistrib <- c(evsizedistrib,sample(res,1))
    scts <- c(scts, sct$`log sct`[j])
    GFP <- c(GFP, sct$`log GFP`[j])  
  }
  
}

mix <- data.frame(evsizedistrib,scts,GFP)

#---- get mode---

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

# Plot size distribution

evsizhist <- ggplot(mix, aes(x=evsizedistrib))
evsizhist + geom_area(stat = "bin",bins = 250,color="black",fill="gray78", alpha = 0.3)+
  scale_x_continuous("Particle size (nm)",limits = c(0,2600), breaks = seq(0,2600,400 )) +
  scale_y_continuous("Count",limits = c(0,40000)) +
  annotate("text", x = getmode(evsizedistrib), y= 38000, label = getmode(evsizedistrib), size = 3.5)+
  theme_classic()+
  theme(
    axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
    axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))
  )

# Plot size calibrated dot plot
graph <- ggplot(mix, aes(x=GFP, y=evsizedistrib))
graph + geom_bin2d(bins = 250) +
  scale_fill_viridis_c(option = "turbo",begin = 0.03,end = 1)+ 
  scale_x_continuous("GFP-H",limits = c(0,5),labels = scales::math_format(10^.x))+
  scale_y_continuous("Particle size (nm)",limits = c(400,2600),breaks = seq(400,2600,400)) +
  annotation_logticks(sides = "b",long = unit(2,"mm"),mid = unit(1.5,"mm"),short = unit(.5,"mm"))+
  theme_classic()+
  theme(legend.position = "none",
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=12),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))

