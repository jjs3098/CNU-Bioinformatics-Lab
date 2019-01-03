##figure 5############################################
#source('MetPC.R') #our R package
library(ggplot2)
library(gridExtra)

sdata <- read.csv("Standard_1.csv",header=T)
ldata <- read.csv("Standard_2.csv",header=T)

sdata <- data.frame(Name=sdata$Name,Area=sdata$Area,Spectra=sdata$Spectra)
ldata <- data.frame(Name=ldata$Name,Area=ldata$Area,Spectra=ldata$Spectra)

#peak merging
psdata <- pmerge(sdata)
pldata <- pmerge(ldata)

#calculate observed data
obdata <- cal_ob(psdata$Spectra,pldata$Spectra)

#gernerate metabolite identification table
result <- MetID(500,pldata,0.8,obdata,muT = 2,muF = 15,muF2=50,sigmaT = 3,sigmaF = 15,sigmaF2 = 50)
write.csv(result, file = "supplimentary material2.csv")

##figure 6#############################################
pars <- estpar_tff(500,obdata,muT=2,muF1=15,muF2=50,sigmaT=3,sigmaF1=15,sigmaF2=50)

plot_pars <- function(pars,iter){
  plot1 <- ggplot(data.frame(Iteration=1:iter,rho=pars$rho),aes(x=Iteration, y=rho))+geom_line()
  plot2 <- ggplot(data.frame(Iteration=1:iter,tau=pars$tau),aes(x=Iteration, y=tau))+geom_line()
  plot3 <- ggplot(data.frame(Iteration=1:iter,muT=pars$muT),aes(x=Iteration, y=muT))+geom_line()
  plot4 <- ggplot(data.frame(Iteration=1:iter,sigmaT=pars$sigmaT),aes(x=Iteration, y=sigmaT))+geom_line()
  grid.arrange(plot1,plot2,plot3,plot4,nrow=2,ncol=2)
}

plot_pars(pars,500)


##figure 7#############################################
fT <- pars$pis[500,1]*dnorm(seq(0,90,0.1),pars$muT[500],sqrt(pars$sigmaT[500]))
fF1 <- pars$pis[500,2]*dnorm(seq(0,90,0.1),pars$muF1[500],sqrt(pars$sigmaF1[500]))
fF2 <- (1-pars$pis[500,1]-pars$pis[500,2])*dnorm(seq(0,90,0.1),pars$muF2[500],sqrt(pars$sigmaF2[500]))
kd <- data.frame(Score=obdata$S[,3])
pd <- data.frame(score=seq(0,90,0.1),density=fT+fF1+fF2)

ggplot(kd,aes(x=Score))+geom_density(aes(color="Kernel density estimator"),size=1) +
  geom_line(data=pd,aes(x=score,y=density,color="Parametric density estimator"),size=1,linetype="dashed") +
  scale_colour_manual(NULL,values = c("black", "violetred")) +
  theme(legend.background = element_rect(size=2, linetype="solid", colour ="darkred"),
    legend.position = c(0.6,0.6), legend.text=element_text(size=30),
    axis.text.x = element_text(size=20), axis.text.y = element_text(size=20),
    axis.title.x = element_text(size=20), axis.title.y = element_text(,size=20))


##figure 8##############################################
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("OCplus", version = "3.8")
library(OCplus)

dat <- read.table("omija.txt",header=FALSE)
colnames(dat) <- rep(c("Chinese","Korean"),c(27,30))
tdat <- dat[21:(nrow(dat)-20),] #cutting top and bottom
ctdat <- log(1+tdat)
for(i in (1:57)){                           ## Normalizing
  ctdat[,i]=ctdat[,i]/sum(ctdat[,i])     ## Log version
}               


# 2d-fdr
fdr <- fdr2d(ctdat,colnames(ctdat),nperm = 500)
summary(fdr)

Tornadoplot(fdr,main="Tornado plot",label=T,cex = 2, cex.main = 2, col="gray50", lcol = "darkred")
Volcanoplot(fdr,55,main="Volcano plot",label=T,cex = 2, cex.main = 2, col="gray50", lcol = "darkred")
