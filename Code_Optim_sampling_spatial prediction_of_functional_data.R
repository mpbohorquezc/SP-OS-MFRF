########################################################################
#Optimal sampling for spatial prediction of functional data
#########################################################################

rm(list=ls())

library(fda)
library(fda.usc)
library(geoR)
library(gstat)
library(splines)
library(colorspace)
library(GenSA)
##Con Fourier

########################################################################
#Read discrete observations, build and plot functional data.  Fourier
#########################################################################

PM10=read.table("PartMat10.txt",head=T,dec=",")
MPM10=as.matrix(PM10,nrow=8761,ncol=10,dimnames=c(rownames(PM10),colnames=colnames(PM10)))
ubicaciones=c("0->Bosque","1->IDRD","2->CarvajalSony","3->Guaymaral","4->SubaCorpas","5->Fontibon","6->PteAranda","7->MAVDT","8->Kennedy","9->Tunal")
estaciones=c("Bosque","IDRD","Carvajal_Sony","Guaymaral","Suba_Corpas","Fontibon","PteAranda","MAVDT","Kennedy","Tunal")
#TMPM10=t(MPM10)
datosf=fdata(MPM10,argvals=1:nrow(MPM10))
nbasis <-91
hourange <- c(1,nrow(MPM10))
hourbasis <- create.fourier.basis(hourange,period=8761/365,nbasis)
lambda=0.00002
harmaccelLfd <- vec2Lfd(c(1,10), hourange)
PM10_fdPar_Fou <- fdPar(fdobj=hourbasis,Lfdobj=harmaccelLfd,lambda)
PM10_fd <- smooth.basis(MPM10,argvals=1:nrow(MPM10),PM10_fdPar_Fou)
PM10_fd_Fou= PM10_fd$fd
plot(PM10_fd_Fou)
lines(PM10_fd_Fou,col=1:10,lwd=3,lty=1)
legend(32,1.3,legend=c(expression(tilde(mu)(t)),expression(tilde(mu)[linf]*"(t),"*tilde(mu)[lsup](t)),expression(mu(t)==sin(2*pi*t))),lwd=3,col=1:10)
x11()
#plotfit.fd(MPM10, argvals=1:nrow(MPM10),PM10_fd_Fou,lwd=1,ylab=" ")

########################################################################
#FPCA. Select eigenfunctions and scores
#########################################################################

PM10_FPC=pca.fd(PM10_fd_Fou,centerfns=T)
plot.pca.fd(PM10_FPC)
puntaje=PM10_FPC$scores
colnames(puntaje)=c("f1","f2")
rownames(puntaje)=estaciones
puntajes=as.data.frame(puntaje)
datos=read.table("Coordenadas_planas_Estaciones_Dama_2014.txt",dec=",",sep="\t",header=T)

################################################################################################
#FPCA. Convert to spatial object and modeled cross variograms, predict and select variance
#################################################################################################
coordenadas=datos[,-1]
coordinates(puntajes)=coordenadas
f1.vgm = variogram(f1~1,puntajes)
f2.vgm = variogram(f2~1,puntajes)
#cross variogram
g = gstat(NULL,"f1",f1~1,puntajes)
g = gstat(g,"f2",f2~1,puntajes)
v = variogram(g)
g = gstat(g, model = vgm(1, "Exp", 300, 1), fill.all = TRUE)      ##por defecto exponencial
g.fit = fit.lmc(v, g)
g.fit
plot(v,g.fit)
##################################################################################################
###con geoR
##################################################################################################

puntaj=data.frame(datos[,-1],puntaje)
puntaj1_s=as.geodata(puntaj, coords.col = 1:2, data.col = 3)
puntaj2_s=as.geodata(puntaj, coords.col = 1:2, data.col = 4)
f1var=variog(puntaj1_s)
f2var=variog(puntaj2_s)
x11()
plot(f1var,main="variog f1")
plot(f2var,main="variog f2")
ini=eyefit(f1var)
plot(f1var,main="variog f1")
v1.fit=variofit(f1var,ini=c(2494708,15000), cov.model="gneiting.matern",fix.kappa = TRUE,kappa=0.95,weights="npairs")
##non_sampled locations
loci=datos[1,2:3]
kc <- krige.conv(puntaj1_s,loc=loci,krige=krige.control(,cov.model="gneiting.matern",cov.pars=v1.fit$cov.pars))
##             gneiting.matern
kc$krige.var

gneiting.maternModel=function(sigma2,h,phi,k1,k2)
{
sigma2-sigma2*(1+8*(10/47)*(h/(phi*k2))+25*(10/47)^2*(h/(phi*k2))^2+32*(10/47)^3*(h/(phi*k2))^3)*(1-(10/47)*h/(phi*k2))^8*matern(h,phi,k1)
}
################################################################################################
#1 week Fourier
################################################################################################

PM10_1semana=PM10[1:168,]
MPM10_1semana=as.matrix(PM10_1semana,nrow=8761,ncol=10,dimnames=c(rownames(PM10),colnames=colnames(PM10)))
datosf_1semana=fdata(MPM10_1semana,argvals=1:nrow(MPM10_1semana))
nbasis_1semana <-168
hourange_1semana <- c(1,nrow(MPM10_1semana))
hourbasis_1semana <- create.fourier.basis(hourange_1semana,period=8761/365,nbasis_1semana)
lambda=0.0000001
harmaccelLfd <- vec2Lfd(c(1,10), hourange_1semana)
PM10_fdPar_Fou_1semana <- fdPar(fdobj=hourbasis_1semana,Lfdobj=harmaccelLfd,lambda)
PM10_fd_1semana <- smooth.basis(MPM10_1semana,argvals=1:nrow(MPM10_1semana),PM10_fdPar_Fou_1semana)
PM10_Fou_1semana_datos= PM10_fd_1semana$fd
plot(PM10_Fou_1semana_datos)
lines(PM10_Fou_1semana_datos,col=rainbow(10),lwd=2,lty=1)
legend(32,1.3,legend=c(expression(tilde(mu)(t)),expression(tilde(mu)[linf]*"(t),"*tilde(mu)[lsup](t)),expression(mu(t)==sin(2*pi*t))),lwd=3,col=1:10)
x11()
#plotfit.fd(MPM10_1semana, argvals=1:nrow(MPM10_1semana),PM10_Fou_1semana_datos,lwd=1,ylab=" ")
x11()
PM10_FPC_Fou=pca.fd(PM10_Fou_1semana_datos,centerfns=T)
plot.pca.fd(PM10_FPC_Fou)

################################################################################################
#With Bsplines best fit
################################################################################################

PM10=read.table("PartMat10.txt",head=T,dec=",")
MPM10=as.matrix(PM10,nrow=8761,ncol=10,dimnames=c(rownames(PM10),colnames=colnames(PM10)))
#TMPM10=t(MPM10)
datosf=fdata(MPM10,argvals=1:nrow(MPM10))
nbasis <-191
hourange <- c(1,nrow(MPM10))
lambda=0.000001
harmaccelLfd <- vec2Lfd(c(1,10), hourange)
PM10_fdPar_Bspline<-fdPar(fdobj=hourbasis_Bsplines,Lfdobj=harmaccelLfd,lambda)
hourbasis_Bsplines <- create.bspline.basis(hourange,nbasis)
PM10_fd_Bspline <- smooth.basis(argvals=1:nrow(MPM10),MPM10,PM10_fdPar_Bspline)
PM10_fd_Bspl=PM10_fd_Bspline$fd
plot(PM10_fd_Bspl)
lines(PM10_fd_Bspl,col=rainbow(10),lwd=2,lty=1)
plotfit.fd(MPM10, argvals=1:nrow(MPM10),PM10_fd_Bspl,lwd=1,ylab=" ")
x11()
PM10_FPC_bspl=pca.fd(PM10_fd_Bspl,centerfns=T)
plot.pca.fd(PM10_FPC_bspl)
puntaje=PM10_FPC_bspl$scores
colnames(puntaje)=c("f1","f2")
rownames(puntaje)=estaciones
puntajes=as.data.frame(puntaje)
datos=read.table("Coordenadas_planas_Estaciones_Dama_2014.txt",dec=",",sep="\t",header=T)
coordenadas=datos[,-1]
coordinates(puntajes)=coordenadas
f1.vgm = variogram(f1~1,puntajes)
f2.vgm = variogram(f2~1,puntajes)
#cross variogram
g = gstat(NULL,"f1",f1~1,puntajes)
g = gstat(g,"f2",f2~1,puntajes)
v = variogram(g)
g = gstat(g, model = vgm(1, "Exp", 300, 1), fill.all = TRUE)
g.fit = fit.lmc(v, g)
g.fit
plot(v,g.fit)

#1semanaBspl
PM10_1semana=PM10[1:168,]
MPM10_1semana=as.matrix(PM10_1semana,nrow=8761,ncol=10,dimnames=c(rownames(PM10),colnames=colnames(PM10)))
datosf_1semana=fdata(MPM10_1semana,argvals=1:nrow(MPM10_1semana))
nbasis <-31
hourange_1semana <- c(1,nrow(MPM10_1semana))
lambda=0.000001
harmaccelLfd <- vec2Lfd(c(1,10), hourange_1semana)
hourbasis_Bsplines <- create.bspline.basis(hourange_1semana,nbasis)
PM10_fdPar_Bspline<-fdPar(fdobj=hourbasis_Bsplines,Lfdobj=harmaccelLfd,lambda)
PM10_fd_Bspline <- smooth.basis(argvals=1:nrow(MPM10_1semana),MPM10_1semana,PM10_fdPar_Bspline)
PM10_fd_Bspl=PM10_fd_Bspline$fd
plot(PM10_fd_Bspl,col="white",xlab="hour",ylab="PM10 (ppm)")
lines(PM10_fd_Bspline,lwd=2,lty=1:10,col="black")
legend(128,144.5,legend=ubicaciones,lwd=2,cex=0.6,lty=1:10)
x11()
plot(PM10_fd_Bspl,col="white",xlab="hour",ylab=expression("PM10"*"  ("*mu*"g"*"/m"^3*")"),xlim=c(15,160),cex.lab=0.8)
lines(PM10_fd_Bspl,col=rainbow(10),lwd=2,lty=1)
legend(126.5,128.5,legend=ubicaciones,lwd=2,cex=0.55,col=rainbow(10))
plotfit.fd(MPM10_1semana, argvals=1:nrow(MPM10_1semana),PM10_fd_Bspl,lwd=1,ylab=" ",col=rainbow(10))
plotfit.fd(MPM10_1semana,PM10_fd_Bspl,lwd=1,ylab=" ",col=rainbow(10))

#Bogotá Map
require(MASS)
require(akima)
require(gstat)
require(geoR)
require(lattice)
require(maptools)
require(rgdal)
require(ape)
require(vegan)
trellis.par.set(sp.theme())
poligonos=readShapePoly("Bogota.shp")
plot(poligonos,axes=T)
#muestra=spsample(poligonos,type="regular",n=100)
#muestra1=as.data.frame(muestra)
#text(muestra1,pch=3)
names(muestra1)=c("x","y")
gridded(muestra1)=c("x","y")
datos1=datos
datos$Z=rep(1,11)
rownames(datos)=datos[,1]
datos=datos[,-1]
xy=SpatialPoints(datos[c("X","Y")])
coordinates(datos)=c("X","Y")
li=list("sp.polygons",poligonos,lwd=2,lty=1,type="s")
pts=list("sp.points",xy,c("0","1","2","3","4","5","6","7","8","9","M"),col=rainbow(10),cex=2)
d1=list("SpatialPolygonsRescale",layout.north.arrow(),offset=c(90000,123500),scale=2000)
d2=list("SpatialPolygonsRescale",layout.scale.bar(),offset=c(85000,120000),scale=8000,fill=c("transparent","black"))
d3=list("sp.text",c(85000,122000),"0")
d4=list("sp.text",c(93500,122000),"5000m")
d5=list("sp.text",c(93000,105000),expression(bold(s)[0]^1))
d6=list("sp.text",c(100000,115000),expression(bold(s)[0]^2))
pru=idw(Z~1,datos,newdata=muestra1,idp=1.5)
spplot(pru,c("var1.pred"),as.table=F,main="",scales=list(draw=T),sp.layout=list(li,pts,d1,d2,d3,d4,d5,d6),contour=F,labels=T,pretty=F,col="red",
col.regions="transparent",colorkey = FALSE,xlim=c(83000,111000),ylim=c(85000,127000))

################################################################################################
###Optimal sampling
################################################################################################

################################################################################################
###If it is not possible to evaluate all options, we use simulated annealing in package GenSA
################################################################################################
#Optimizing the prediction with kriging, simulated annealing
#x <- available East coordinates
#y <- available North coordinates
#z=c(x,y)  As vec with x and y
##

rm(list=ls())
library(GenSA)
#x <- available East coordinates 
#y <- available North coordinates
#z=c(x,y)  As vec with x and y
## 
##coordenadas=data.frame("este"=z[1:n],"norte"=z[n1:N])

#Covariance functions
#cov.spatial1<-function(h,modelo,p){modelo(h,)}

#modelos
#Any model of cov.spatial or
exponencial=function(h,cs,a){cs*exp(-h/a)}
gaussiano=function(h,cs,a){cs*exp(-h^2/a^2)}

varSK<-function(z)
{
fn.call <<- fn.call + 1
N=length(z)
n=0.5*N
n1=n+1
coordenadas=data.frame("East"=z[1:n],"North"=z[n1:N])
   MatDistances=as.matrix(dist(coordenadas))
   covMatrix=exponencial(MatDistances,cs=1,a=2.5)
   coordenadas_and_s0=rbind(coordenadas,c(3.5,4.3))
   DistSample_s0=as.matrix(dist(coordenadas_and_s0))[nrow(coordenadas)+1,-1*(nrow(coordenadas)+1)]
   covSample_s0=exponencial(DistSample_s0,cs=1,a=2.5)
   vari=t(covSample_s0)%*%covMatrix%*%covSample_s0
   varOK=as.numeric(vari)
}
lower <- #array with minimal values for Este  Ex. rep(0,n)
upper <- #array with minimal values for Norte Ex. rep(1,n)
p0<- #Current position of network (initial values)
fn.call<-0
expected.val <- 0
absTol <- 1e-13
fn.call<-0
out.GenSA <- GenSA(par = NULL, lower = lower, upper = upper, fn = varSK)

################################################################################################
#If there are only a few candidates location, it is possible exhaustively
################################################################################################

##p: parameter vector
##cs0

varSK<-function(N,n,s0,Ds,modelo,p,kappa)
{
combinaciones=combn(N,n)
x=matrix(,nrow=n,ncol=2)
colnames(x)=c("este","norte")
tablaVar=data.frame()
for(i in 1:ncol(combinaciones))
      {
      x=Ds[combinaciones[,i],]
      MatDistances=as.matrix(dist(x))
        covMatrix=matrix(cov.spatial(MatDistances,cov.model=modelo,cov.pars=p,kappa =kappa),nrow=n,ncol=n)
         x_s0=rbind(x,s0)
         DistSample_s0=as.matrix(dist(x_s0))[nrow(x)+1,-1*(nrow(x)+1)]
         covSample_s0=matrix(cov.spatial(DistSample_s0,cov.model=modelo,cov.pars=p,kappa =kappa),nrow=1,ncol=n)
         vari=covSample_s0%*%covMatrix%*%t(covSample_s0)
         varOK=as.numeric(vari)
         tablaVar[i,]=c(i,varOK)
}
TablaVar=as.data.frame(tablaVar)
colnames(tablaVar)=c("combinación","varianzaKO")
names(tablaVar)=tablaVar[,1]
TablaVar=as.data.frame(tablaVar)
attach(TablaVar)
TablaVarOrder=TablaVar[order(varianzaKO),]
return(list(tablaVar,TablaVarOrder[1,]))
}

################################################################################################
################################################################################################
L <- function(z)
{
fn.call <<- fn.call + 1
N=length(z)
n=0.5*N
n1=n+1
coordenadas=data.frame(z[1:n],z[n1:N])
sum(solve(exponencial(as.matrix(dist(data.frame("East"=z[1:n],"North"=z[n1:N]))),1,2.3)))
}
L(z=z)
lower <- #array with minimal values for Este
upper <- #array with minimal values for Norte
p0<- #Current position of network (initial values)
fn.call<-0
expected.val <- 0
absTol <- 1e-13
fn.call<-0
out.GenSA <- GenSA(par = p0, lower = lower, upper = upper, fn = L)
                                                                                              
################################################################################################
##preliminar computations; auxiliar and alternative functions
################################################################################################
##file grilla: 2 columns with all possible locations to move stations
grilla=read.table("grillaBog.txt",head=T)
#plot(grilla$x,grilla$y)
#attach(grilla)

###Cross covariance matrix
matriz<-function(distancias,n2,n,cova1,cova12,cova21,cova2,phi1,phi2,phi12,phi21,sigma1,sigma2,sigma12,sigma21)
{
mc<-distancias
for(i in 1:n2)
  for(j in 1:n2)  
  if(i<=n && j<=n) mc[i,j]<-cova1(sigma1,phi1,distancias[i,j]) else 
  {if(i<=n && j>n)
    mc[i,j]<-cova12(sigma12,phi12,distancias[i,j]) else
    {
      if(i>n && j<=n) mc[i,j]<-cova21(sigma21,phi21,distancias[i,j])
     else mc[i,j]<-cova2(sigma2,phi2,distancias[i,j])}
  }
      return(mc)
}

################################################################################################
#If there are only a few candidates location, it is possible exhaustively
################################################################################################

##p: parameter vector
##cs0
muestras=combn(red,n)
todos=rbind(red,s0)
distancias=as.matrix(dist(todos))
cs0=distancias[-n,n] 

optlocation<-function(fcova,p,red,ubicaciones,s0)
{ 
  uno<-rep(1,nrow(red)+1)
  vars<-c()
    for(i in 1:nrow(ubicaciones))
    {
    red_i<-rbind(red,ubicaciones[i,]) 
    distancias_i<-as.matrix(dist(red_i))
    cova_i<-fcova(p,distancias_i)
    red_0i<-rbind(red_i,s0)
    red_0i<-red_0i[-nrow(red_0i),nrow(red_0i)]
    distancias_0i<-as.matrix(dist(red_0i))
    cova_0i<-fcova(p,distancias_0i)
    var_0_i<-(sum(cova_0i)*cova_i*sum(cova_0i))^(-1)
    vars<-c(vars,var_0_i)
     }
  resumen<-cbind(ubicaciones,vars)
  resumen_ord<-resumen[order(vars),]
  redNueva<-rbind(red, resumen_ord[1,1:2])
  best_location<-list(optloc=resumen_ord[1,],redNueva=redNueva,Ubicaciones_vars=resumen_ord)
  return(best_location)
    }
###variance   
#####pseudocode  it could be also in c, c++, fortran
###any covariance model, p: parameter vector
cov.spatial1<-function(h,modelo,p){modelo(h,p)}
variance=function(grilla,muestra,n,s0,cova){
sum(cova_0i)*cova(distancias)*sum(cova_0i)
muestra=combn(grilla,n)                                                         
nueva_grilla=grilla[muestra[,i],]
distancia_i=as.matrix(dist(nueva_grilla))
nueva_grilla_mas_s0=rbind(muestra[,i],s0)
distancia_0i=as.matrix(dist(nueva_grilla_mas_s0))
vector_cs0i=distancia_0i[-n+1,n+1]
cs0i=sum(cov.spatial1(vector_cs0,modelo,p))
omega_i=cov.spatial1(distancia_i,modelo,p)
variance_i=cs0i*omega_i*cs0i}
expected.val <- 0
absTol <- 1e-13
out.GenSA <- GenSA(par = NULL, lower = lower, upper = upper, fn = variance,
+ control = list(threshold.stop = expected.val + absTol))
out.GenSA[c("value", "par", "counts")]
##upper: sum of eigenvalues
###Cross validation process, plot fit and residual
for(i in 1:ncol(PM10))
PM10_i=PM10[,-1]
MPM10_i=as.matrix(PM10_i,nrow=8761,ncol=9,dimnames=c(rownames(PM10),colnames=colnames(PM10)))
PM10_fd_Bspline_i<- smooth.basis(argvals=1:nrow(MPM10_i),MPM10_i,PM10_fdPar_Bspline)
PM10_fd_Bspl_i=PM10_fd_Bspline_i$fd
PM10_FPC_bspl_i=pca.fd(PM10_fd_Bspl_i,centerfns=T)
#plot.pca.fd(PM10_FPC_bspl_1)
puntaje_1=PM10_FPC_bspl_i$scores
datos=read.table("Coordenadas_planas_Estaciones_Dama_2014.txt",dec=",",sep="\t",header=T)
coordenadas=datos[-i,-i]
puntaj_1=data.frame(coordenadas,puntaje_1)
puntaj1_s=as.geodata(puntaj_1, coords.col = 1:2, data.col = 3)
f1var=variog(puntaj1_s)
#ini=eyefit(f1var)
plot(f1var,main="variog f1")
v1.fit=variofit(f1var,ini=c(5500000,10000), cov.model="...",fix.kappa = TRUE,kappa=0.5,weights="npairs")
loci=datos[i,2:3]
kc <- krige.conv(puntaj1_s,loc=loci,krige=krige.control(,cov.model="...",cov.pars=v1.fit$cov.pars))
lines((kc$predict)*PM10_FPC_bspl_i$harmonics[1]+mean(PM10_fd_Bspl_i),cex=5,col="#FF0000FF",lwd=2)
plot(PM10_fd_Bspl[i]-((kc$predict)*PM10_FPC_bspl_i$harmonics[1]+mean(PM10_fd_Bspl_i)),ylab="Residuals",cex=5,col=4,lwd=2)
