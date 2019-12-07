################################################################
###################  FECUNDIDAD ################################
################################################################

rm(list = ls())
library(tidyverse)
library(readxl)
library(forecast)

fx5<-read.csv("fecu.csv", header = T)
fx5.proy<-read.csv("fx5.proy.csv", header = T)


TGF <- 5*colSums(fx5[,-1])
Vx <- fx5[,-1]
FxF <- Vx
for (i in 2:47) {
  Vx[,i-1] <- 5*cumsum(fx5[,i])
  FxF[,i-1] <- Vx[,i-1]/TGF[i-1]
}

x5 <- seq(17.5,47.5,5)

Yx.fec <- log(-log(FxF))
Yx.fec.lm <- list()
for (i in 1:46) {
  Yx.fec.lm[[i]] <- lm (Yx.fec[-7,i] ~ x5[-5])
}

a <- vector(length = 46)
b <- vector(length = 46)
for (i in 1:46) {
  a[i] <- Yx.fec.lm[[i]]$coefficients[1]
  b[i] <- Yx.fec.lm[[i]]$coefficients[2]
}

A <- exp(-exp(a))
B <- exp(b)
plot(A, type = "l")
plot(B, type = "l")
k1<-1
k2<-max(TGF)-k1+0.2
Yx<-log(k2/(TGF-k1)-1)
Yx.fit<-auto.arima(Yx,trace=T, d=1)
Yx.for<-forecast(Yx.fit, h=35, c(95))


x<-c(15:50)
FxF.Ajustada<-matrix(0,36,35)
TGF.Ajustada<-k1+k2/(1+exp(Yx.for$mean))

for (i in 1:35) {
  FxF.Ajustada[,i]<-TGF.Ajustada[i]*A[36]^(B[36]^(x))
}

fx<-matrix(0,36,35)
fx[1,]<-FxF.Ajustada[1,]
fx[2:36,]<-FxF.Ajustada[2:36,]-FxF.Ajustada[1:35,]
fx<-as.data.frame(fx)

matplot(fx, type = "l", main="Tasas especificas de fecundidad", xlab="Edad", ylab = "Fx")

write.csv(fx, file = "fx.csv")

#### Proyectando tasa global de fecundidad con series de tiempo
Yx.fit<-auto.arima(Yx,trace=TRUE, d=1)
TGF.fit<-k1+k2/(1+exp(Yx.fit$fitted))
Yx.for1<-forecast(Yx.fit, h=35, c(70))
TGF.for<-k1+k2/(1+exp(Yx.for1$mean))
TGF.for.up <- k1 + k2/(1+exp(Yx.for1$upper))
TGF.for.low <- k1 + k2/(1+exp(Yx.for1$lower))

TGF.dat <- data.frame(
  year=c(1970:2050),
  mean=c(TGF.fit, TGF.for),
  up=c(rep(NA,46), TGF.for.up),
  low=c(rep(NA,46), TGF.for.low)
)

library("ggplot2")

ggplot(TGF.dat, aes(x=year, y=mean))+
  geom_ribbon(aes(ymin=up, ymax=low), alpha=0.4, fill="blue")+
  geom_line(aes(y=mean), col="black")+
  ggtitle("Tasa Global de Fecundidad (TGF)")
  

###############################################################################
################################ MORTALIDAD ###################################
###############################################################################
library(tidyverse)

tabmort.hist <- readxl::read_xlsx("oaxaca.xlsx",skip = 3)
tabmort.hist <- tabmort.hist %>%
  select(c("Año y edad", "m(x)...3",  "m(x)...11"))

tabmort.hist <-
  tabmort.hist[is.na(tabmort.hist$"m(x)...3") == FALSE,]
names(tabmort.hist) <- c("edad", "mxH", "mxM")
tabmort.hist$yr <- rep(c(1970:2015), each = 110)

tabmort.df <- data.frame(matrix(0,220,46))

for(i in 1:46){
  tabmort.df[1:110,i] <-
    tabmort.hist[tabmort.hist$yr == (1969+i), "mxH"]
  tabmort.df[111:220,i] <-
    tabmort.hist[tabmort.hist$yr == (1969+i), "mxM"]
}

names(tabmort.df) <- c(1970:2015)
matplot(tabmort.df,main="Tasas específicas de mortalidad (Hombres-Mujeres)", xlab = "Edad", ylab = "Tasa", type="l", log="y")

lmx <- log(tabmort.df)
ax <- rowMeans(lmx)
plot(ax, lty=1, type="l")

lmx_ax <- lmx - ax
matplot(lmx_ax, type="l")

svdLC <- svd(lmx_ax)
options(scipen=999)
length(100*svdLC$d/sum(svdLC$d))

bx <- -svdLC$u[,1]
plot(bx, type="l")

D <- matrix(0,66,46)
diag(D) <- svdLC$d

kt <- -(D%*%t(svdLC$v))[1,]
plot(kt,type="l")

kt.fit <- forecast::auto.arima(kt, trace = T, d=1)
kt.for <- forecast::forecast(kt.fit, h = 35, c(95))

matplot(cbind(c(kt.fit$fitted, kt.for$mean),
              c(rep(NA,46), kt.for$upper),
              c(rep(NA,46), kt.for$lower)),type="l", lty = 1)

mx.fit <- exp(ax + bx%*%t(kt.fit$fitted))
mx.formean <- exp(ax + bx%*%t(kt.for$mean))
mx.forup <- exp(ax + bx%*%t(kt.for$upper))
mx.forlow <- exp(ax + bx%*%t(kt.for$lower))

matplot(cbind(mx.formean[,35],
              mx.forlow[,35],
              mx.forup[,35]),type="l",
        log="y")
tabmort <- function(m,edades,sex){
  
  mx <- m
  
  nax <- matrix(0.5,dim(mx)[1],dim(mx)[2])
  ## 1 MUJERES 2 HOMBRES
  if(sex==1){
    for(i in 1:dim(mx)[2]){
      if(mx[1,i]<0.01724){
        nax[1,i] <- 0.14903-2.05527*mx[1,i]
      }else if(mx[1,i]>=0.01724 & mx[1,i]<0.06891){
        nax[1,i] <- 0.04667+3.88089*mx[1,i]
      }else{nax[1,i] <- 0.31411}
    }
  }else{
    for(i in 1:dim(mx)[2]){
      if(mx[1,i]<0.023){
        nax[1,i] <- 0.14929-1.99545*mx[1,i]
      }else if(mx[1,i]>=0.023 & mx[1,i]<0.08307){
        nax[1,i] <- 0.02832+3.26021*mx[1,i]
      }else{nax[1,i] <- 0.29915}
    }
  }
  
  
  nax[edades,] <- 1/mx[edades,]
  
  qx<-matrix(1,dim(mx)[1],dim(mx)[2])
  
  for(i in 1:(dim(mx)[1])){
    qx[i,]<-mx[i,]/(1+(1-nax[i,])*mx[i,])
  }
  
  px <- 1-qx
  
  lx<-matrix(1,dim(mx)[1],dim(mx)[2])
  
  for(i in 2:dim(mx)[1]){
    lx[i,] <- lx[i-1,]*px[i-1,]
  }
  
  dx <- matrix(0,dim(mx)[1],dim(mx)[2])
  dx[dim(mx)[1],] <- lx[dim(mx)[1],]
  for(i in 1:(dim(mx)[1]-1)){
    dx[i,]<-lx[i,]-lx[i+1,]
  }
  
  
  Lx<-matrix(0,dim(mx)[1],dim(mx)[2])
  Lx[1,] <- dx[1,]/mx[1,]
  Lx[edades,] <- dx[edades,]/mx[edades,]
  for(i in 2:(edades-1)){
    Lx[i,]<-(lx[i,]+lx[i+1,])/2
  }
  
  Tx<-matrix(0,dim(mx)[1],dim(mx)[2])
  Tx[edades,]<-Lx[edades,]
  for(i in (edades-1):1){
    Tx[i,]<-Lx[i,]+Tx[i+1,]
  }
  
  ex <- Tx/lx
  
  Sx<-matrix(NA,(dim(mx)[1]+1),dim(mx)[2])
  Sx[1,]<-Lx[1,]/lx[1,]
  Sx[(edades+1),] <- Tx[edades,]/Tx[(edades-1),]
  for(i in 2:edades){
    Sx[i,]<-Lx[i,]/Lx[i-1,]
  }
  
  tabmort <- list(Edad=c(0:(edades-1)),mx=mx, nax=nax, qx=qx,
                  px=px, lx=lx, dx=dx, Lx=Lx, Tx=Tx, ex=ex, Sx=Sx)
}

tmH.fit <- tabmort(mx.fit[1:110,],
                   edades = 110, sex = 2)
tmM.fit <- tabmort(mx.fit[111:220,],
                   edades = 110, sex = 1)
tmH.mean <- tabmort(mx.formean[1:110,],
                    edades = 110, sex = 2)
tmH.up <- tabmort(mx.forup[1:110,],
                  edades = 110, sex = 2)
tmH.low <- tabmort(mx.forlow[1:110,],
                   edades = 110, sex = 2)
tmM.mean <- tabmort(mx.formean[111:220,],
                    edades = 110, sex = 1)
tmM.up <- tabmort(mx.forup[111:220,],
                  edades = 110, sex = 1)
tmM.low <- tabmort(mx.forlow[111:220,],
                   edades = 110, sex = 1)

tmH.mean$nax

e0.for <- data.frame(year = c(1970:2050),
                     e0Hmean = c(tmH.fit$ex[1,], tmH.mean$ex[1,]),
                     e0Mmean = c(tmM.fit$ex[1,], tmM.mean$ex[1,]),
                     e0Hup = c(rep(NA,46), tmH.up$ex[1,]),
                     e0Hlow = c(rep(NA,46), tmH.low$ex[1,]),
                     e0Mup = c(rep(NA,46), tmM.up$ex[1,]),
                     e0Mlow = c(rep(NA,46), tmM.low$ex[1,]))

Sx.for <- data.frame(year = c(1970:2050),
                     edad = c(1:110),
                     SxH = c(tmH.fit$Sx[1:110,], tmH.mean$Sx[1:110,]),
                     SxM = c(tmM.fit$Sx[1:110,], tmM.mean$Sx[1:110,]))

write.csv(Sx.for, file = "Sx.csv")

library(ggplot2)
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
e0H.plot <- ggplot(e0.for,aes(x=year,y=e0Hmean))+
  geom_ribbon(aes(ymin=e0Hlow,ymax=e0Hup),alpha=0.4,fill="green")+
  geom_line(aes(y=e0Hmean),col="black")+
  scale_x_continuous('Año')+
  scale_y_continuous('Esperanza de vida al nacimiento')+
  ggtitle("Hombres")+
  coord_cartesian(ylim=c(40, 90))+
  theme_bw()
e0M.plot <- ggplot(e0.for,aes(x=year,y=e0Mmean))+
  geom_ribbon(aes(ymin=e0Mlow,ymax=e0Mup),alpha=0.5,fill="pink")+
  geom_line(aes(y=e0Mmean),col="black")+
  scale_x_continuous('Año')+
  scale_y_continuous('Esperanza de vida al nacimiento')+
  ggtitle("Mujeres")+
  coord_cartesian(ylim=c(40, 90))+
  theme_bw()

multiplot(e0H.plot,e0M.plot,cols=2)

###########################################################################################
#################  MIGRACIÓN ##############################################################
##########################################################################################

Inmx<-read.csv("Inmg.csv", header = T)
Emigx<-read.csv("Emig.csv", header = T)
Nx<-read.csv("Nx.csv", header = T)


ixt.F<-as.matrix(Inmx[Inmx$Sexo=="Mujeres",-c(1:2)])/as.matrix(Nx[Nx$Sexo=="Mujeres", -c(1:2)])
ixt.M<-as.matrix(Inmx[Inmx$Sexo=="Hombres",-c(1:2)])/as.matrix(Nx[Nx$Sexo=="Hombres", -c(1:2)])
ext.F<-as.matrix(Emigx[Emigx$Sexo=="Mujeres",-c(1:2)])/as.matrix(Nx[Nx$Sexo=="Mujeres", -c(1:2)])
ext.M<-as.matrix(Emigx[Emigx$Sexo=="Hombres",-c(1:2)])/as.matrix(Nx[Nx$Sexo=="Hombres", -c(1:2)])
ixt.F<-ixt.F[1:(dim(ixt.F)[1]-14),]
ixt.M<-ixt.M[1:(dim(ixt.M)[1]-14),]
ext.F<-ext.F[1:(dim(ext.F)[1]-14),]
ext.M<-ext.M[1:(dim(ext.M)[1]-14),]


layout(matrix(c(1,2,3,4),2,2, byrow = T))
matplot(ixt.F, main="Tasas específicas de Inmigración de las Mujeres", xlab = "Edad", ylab = "Tasa", type="l")
matplot(ixt.M, main="Tasas específicas de Inmigración de los Hombres", xlab = "Edad", ylab = "Tasa", type="l")
matplot(ext.F, main="Tasas específicas de Emigración de las Mujeres", xlab = "Edad", ylab = "Tasa", type="l")
matplot(ext.M, main="Tasas específicas de Emigración de las Hombres", xlab = "Edad", ylab = "Tasa", type="l")



######## Usamos el modelo de Lee Carter 
Inmx<-read.csv("Inmg.csv", header = T)
Emigx<-read.csv("Emig.csv", header = T)
Nx<-read.csv("Nx.csv", header = T)

ixt.F<-as.matrix(Inmx[Inmx$Sexo=="Mujeres",-c(1:2)])/as.matrix(Nx[Nx$Sexo=="Mujeres", -c(1:2)])
ixt.M<-as.matrix(Inmx[Inmx$Sexo=="Hombres",-c(1:2)])/as.matrix(Nx[Nx$Sexo=="Hombres", -c(1:2)])
ext.F<-as.matrix(Emigx[Emigx$Sexo=="Mujeres",-c(1:2)])/as.matrix(Nx[Nx$Sexo=="Mujeres", -c(1:2)])
ext.M<-as.matrix(Emigx[Emigx$Sexo=="Hombres",-c(1:2)])/as.matrix(Nx[Nx$Sexo=="Hombres", -c(1:2)])
ixt.F<-ixt.F[1:(dim(ixt.F)[1]-19),]
ixt.M<-ixt.M[1:(dim(ixt.M)[1]-19),]
ext.F<-ext.F[1:(dim(ext.F)[1]-19),]
ext.M<-ext.M[1:(dim(ext.M)[1]-19),]

lnix<-log(ixt.F)
ax<-rowMeans(lnix)
matplot(lnix, type = "l")
lines(ax, lwd=2)
lines(lnix[,41], lwd=2, col="orange")

lnix_ax<-lnix-ax
matplot(lnix_ax, type="l")
100*cumsum(svd(lnix_ax)$d/sum(svd(lnix_ax)$d))
bx<-svd(lnix_ax)$u[,1]
D<-matrix(0,length(svd(lnix_ax)$d), length(svd(lnix_ax)$d))
diag(D)<-svd(lnix_ax)$d
it<-(D%*%t(svd(lnix_ax)$v))[1,]

matplot(bx, type = "l")
matplot(it, type = "l")

library(forecast)
it.fit<-auto.arima(it, trace=T, d=1)
it.for<-forecast(it.fit, h=35, c(95))

ixt.for<-exp(ax+bx%*%t(it.for$mean))
ixt.for.u<-exp(ax+bx%*%t(it.for$upper))
ixt.for.l<-exp(ax+bx%*%t(it.for$lower))

matplot(cbind(ixt.for.l[,1],
              ixt.for[,1],
              ixt.for.u[,1]),
        type = "l")

####Lee Carter

Inmx<-read.csv("Inmg.csv", header = T)
Emigx<-read.csv("Emig.csv", header = T)
Nx<-read.csv("Nx.csv", header = T)

ixt.F<-as.matrix(Inmx[Inmx$Sexo=="Mujeres",-c(1:2)])/as.matrix(Nx[Nx$Sexo=="Mujeres", -c(1:2)])
ixt.M<-as.matrix(Inmx[Inmx$Sexo=="Hombres",-c(1:2)])/as.matrix(Nx[Nx$Sexo=="Hombres", -c(1:2)])
ext.F<-as.matrix(Emigx[Emigx$Sexo=="Mujeres",-c(1:2)])/as.matrix(Nx[Nx$Sexo=="Mujeres", -c(1:2)])
ext.M<-as.matrix(Emigx[Emigx$Sexo=="Hombres",-c(1:2)])/as.matrix(Nx[Nx$Sexo=="Hombres", -c(1:2)])
ixt.F<-ixt.F[1:(dim(ixt.F)[1]-19),]
ixt.M<-ixt.M[1:(dim(ixt.M)[1]-19),]
ext.F<-ext.F[1:(dim(ext.F)[1]-19),]
ext.M<-ext.M[1:(dim(ext.M)[1]-19),]

lc.mig<-function(mig,h){
  require(forecast)
  lnix<-log(mig)
  ax<-rowMeans(lnix)
  lnix_ax<-lnix-ax
  bx<-svd(lnix_ax)$u[,1:5]
  D<-matrix(0,length(svd(lnix_ax)$d), length(svd(lnix_ax)$d))
  diag(D)<-svd(lnix_ax)$d
  it<-(D%*%t(svd(lnix_ax)$v))[1:5,]
  it.fit<-list()
  it.for<-list()
  for (i in 1:5) {
    it.fit[[i]]<-
      auto.arima(it[i,], trace=T)
    it.for[[i]]<-
      forecast(it.fit[[i]], h=h, c(95))
    
  }
  
  ixt.for<-exp(ax+bx[,1]%*%t(it.for[[1]]$mean)+
                 bx[,2]%*%t(it.for[[2]]$mean)+
                 bx[,3]%*%t(it.for[[3]]$mean)+
                 bx[,4]%*%t(it.for[[4]]$mean)+
                 bx[,5]%*%t(it.for[[5]]$mean))
  
  ixt.for.u<-exp(ax+bx[,1]%*%t(it.for[[1]]$upper)+
                   bx[,2]%*%t(it.for[[2]]$upper)+
                   bx[,3]%*%t(it.for[[3]]$upper)+
                   bx[,4]%*%t(it.for[[4]]$upper)+
                   bx[,5]%*%t(it.for[[5]]$upper))
  
  ixt.for.l<-exp(ax+bx[,1]%*%t(it.for[[1]]$lower)+
                   bx[,2]%*%t(it.for[[2]]$lower)+
                   bx[,3]%*%t(it.for[[3]]$lower)+
                   bx[,4]%*%t(it.for[[4]]$lower)+
                   bx[,5]%*%t(it.for[[5]]$lower))
  
  mig.list<-list(modelo=it.fit,
                 media=ixt.for,
                 inf=ixt.for.l,
                 sup=ixt.for.u)
  return(mig.list)
}

inm.M<-lc.mig(ixt.M, h=35)
inm.F<-lc.mig(ixt.F, h=35)
em.M<-lc.mig(ext.M, h=35)
em.F<-lc.mig(ext.F, h=35)

Mx.M<-exp(inm.M$media[,]-em.M$media[,])
Mx.F<-exp(inm.F$media[,]-em.F$media[,])

layout(matrix(c(1,2),2,2, byrow = T))
matplot(Mx.F, main="Saldo Neto Migratorio 2015-2050 (Mujeres)", xlab = "Edad", ylab = "Tasa", type="l")
matplot(Mx.M ,main="Saldo Neto Migratorio 2015-2050 (Hombres)", xlab = "Edad", ylab = "Tasa", type="l")





#################################################################
################# AGREGACION FINAL ##############################
#################################################################
#Fecundidad
SxF<-read.csv("SxF.csv")
SxM<-read.csv("SxH.csv")

#Migracion
Mx<-read.csv("MxT.csv")

#Población a inicio de year
Px<-read.csv("Px.csv", header = T)

#Matriz de leslie

A<- array(0, dim=c(110, 110, 35))
B<- array(0, dim=c(110, 110, 35))

#A<-matrix(0,110,110)
for (i in 1:35) {
  A[1,16:51,i]<-fx[,i]
  
  diag(A[-1,,i])<-SxF[SxF$year==2015 +i &
                        SxF$edad!=109,
                      "S.x.."]*mx[mx$sexo == "Mujeres" & mx$year ==2015+i & mx$edad!=61, "valor"]          
  
  A[110,110,i]<- SxF[SxF$year==2015+i &
                       SxF$edad==109,
                     "S.x.."] *mx[mx$sexo == "Mujeres" & mx$year ==2015+i & mx$edad==61, "valor"]
  
  diag(B[-1,,i])<-SxM[SxM$year==2015+i &
                        SxM$edad!=109,
                      "S.x.."]*mx[mx$sexo == "Hombres" & mx$year ==2015+i & mx$edad!=61, "valor"]          
  
  B[110,110,i]<- SxM[SxM$year==2015+i &
                       SxM$edad==109,
                     "S.x.."]*mx[mx$sexo == "Hombres" & mx$year ==2015+i & mx$edad==61, "valor"]
  
  
}

for (i in 1:35) {
  A[1,16:51,i]<-fx[,i]
  
  diag(A[-1,,i])<-SxF[SxF$year==2015 +i &
                        SxF$edad!=109,
                      "S.x.."]        
  
  A[110,110,i]<- SxF[SxF$year==2015+i &
                       SxF$edad==109,
                     "S.x.."]
  
  diag(B[-1,,i])<-SxM[SxM$year==2015+i &
                        SxM$edad!=109,
                      "S.x.."]        
  
  B[110,110,i]<- SxM[SxM$year==2015+i &
                       SxM$edad==109,
                     "S.x.."]
  
  
}


#Proyeccion de poblacionPx$X1970.00
n0F<-Px[Px$Sexo=="Mujeres" , 48]
n0M<-Px[Px$Sexo=="Hombres", 48]

nF.proy <- data.frame(matrix(0,110,36))
row.names(nF.proy) <- c(0:109)
names(nF.proy)<- c(2015:2050)
nF.proy[,1] <- n0F

nM.proy <- data.frame(matrix(0,110,36))
row.names(nM.proy) <- c(0:109)
names(nM.proy)<- c(2015:2050)
nM.proy[,1] <- n0M

for (i in 2:36) {
  nF.proy[,i] <- A[,,i-1]%*%nF.proy[,i-1]  
  nF.proy[1,i] <- A[1,,i-1]%*%nF.proy[,i-1]*(1/2.05)
  
  nM.proy[,i] <- B[,,i-1]%*%nM.proy[,i-1]  
  nM.proy[1,i] <- A[1,,i-1]%*%nF.proy[,i-1]*(1.05/2.05)
}

matplot(nF.proy, type = "l")
colSums(nF.proy)

matplot(nM.proy, type = "l")
colSums(nM.proy)

pob.tot.oaxaca<- colSums(nM.proy)+ colSums(nF.proy)

pob.conapo<-read.csv("ind_dem_proyecciones.csv", header = T)
Pob.Conapo<-pob.conapo[1686:1721,"POB_MIT_AÑO"]

layout(matrix(c(1,2),2,2, byrow = T))
matplot(pob.tot.oaxaca, main="Población para el estado de Oaxaca (2015-2015)", xlab = "year", ylab = "Población", type="l")
matplot(Pob.Conapo, main="Población para el estado de Oaxaca de Conapo (2015-2015)", xlab = "year", ylab = "Población", type="l")

###Proyecciones probabilisticas
NxM<-read.csv("NxM.csv", header = T)
NxF<-read.csv("NxF.csv", header = T)
PobTotM<-as.matrix(t(colSums(NxM)))
PobTotH<-as.matrix(t(colSums(NxF)))

PobTotM.fit<-auto.arima(t(PobTotM), d=1)
PobTotM.for<-forecast(PobTotM.fit, h=35, c(95))

PobTotH.fit<-auto.arima(t(PobTotH), d=1)
PobTotH.for<-forecast(PobTotH.fit, h=35, c(95))

layout(matrix(c(1,2),2,2, byrow = T))
plot(PobTotH.for, main = "Población Hombres-Oaxaca 1970 a 2050", xlab = "year", ylab = "Población")+abline(v=42, col="blue")
plot(PobTotM.for, main = "Población Mujeres-Oaxaca 1970 a 2050", xlab = "year", ylab = "Población")+abline(v=42, col="blue")

write.csv(nM.proy, file = "nM.proy.csv")
write.csv(nF.proy, file = "nF.proy.csv")
