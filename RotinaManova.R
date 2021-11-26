remove(list=ls())
library(MASS)
setwd("D:/Videos/_Analise multivariada/5.0 MANOVA")
D=read.table("Dados.txt",h=T)

#############################################################
MatrizIncidencia=function(Dados){
  Dados=as.numeric(as.factor(Dados))
  NumTrat=max(Dados);
  NumParc=length(Dados);
  Matrizout=matrix(0,ncol=NumTrat,nrow=NumParc)
  for (i in 1:NumParc){
    for (j in 1:NumTrat){
      if (Dados[i]==j){Matrizout[i,j]=1; }
    }}
  return(Matrizout)
}
##########################################################
Im=MatrizIncidencia(rep(1,nrow(D)))
Itrat=MatrizIncidencia(D$Trat)
Ibloc=MatrizIncidencia(D$Bloco)

Y=as.matrix(D[,3:5])

########## rm
X=Im
B=solve(t(X)%*%X)%*%t(X)%*%Y
Yp=X%*%B
rm=t(Yp)%*%Yp
rm

######## r modelo
X=cbind(Im,Itrat,Ibloc)
B=ginv(t(X)%*%X)%*%t(X)%*%Y
Yp=X%*%B
rmodelo=t(Yp)%*%Yp-rm
rmodelo

######## r total
X=diag(nrow(D))
B=ginv(t(X)%*%X)%*%t(X)%*%Y
Yp=X%*%B
rtotal=t(Yp)%*%Yp-rm
rtotal

############ SQ residuo
SQPresiduo=rtotal-rmodelo


########## r(m,bloco)
X=cbind(Im,Ibloc)
B=ginv(t(X)%*%X)%*%t(X)%*%Y
Yp=X%*%B
rm.bloco=t(Yp)%*%Yp-rm
rm.bloco


########## r(m,trat)
X=cbind(Im,Itrat)
B=ginv(t(X)%*%X)%*%t(X)%*%Y
Yp=X%*%B
rm.trat=t(Yp)%*%Yp-rm
rm.trat


rmodelo-rm.bloco
rmodelo-rm.trat



library(ExpDes.pt)
dbc(D$Trat,D$Bloco,D$PMV)
dbc(D$Trat,D$Bloco,D$PMS)
dbc(D$Trat,D$Bloco,D$PRC)


library(MultivariateAnalysis)
m=MANOVA(D)
Dist=Distancia(m$Med,Metodo = 7,m$CovarianciaResidual)
Dendrograma(Dist,Metodo = 1)
Dendrograma(Dist,Metodo = 2)
Dendrograma(Dist,Metodo = 3)

Tocher(Dist,Metodo = "original")
Tocher(Distancia(m$Med),Metodo = "original")


library(ExpDes.pt)
