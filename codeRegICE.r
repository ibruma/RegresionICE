## ----echo=FALSE---------------------------------------------------------------
options(width=80,digits=5)
rm(list=ls(all=TRUE))
library(xtable)
library(car)


## ----echo=TRUE----------------------------------------------------------------
library(readxl)
Purines <- read_xlsx("Purines_Metano2.xlsx")
Purines <- as.data.frame(Purines)


## ----echo=TRUE----------------------------------------------------------------
round(cor(Purines),2)


## ----echo=FALSE---------------------------------------------------------------
## Las 5 primeras observaciones del marco de datos Purines
head(Purines,n=5)
datos <- Purines


## ----echo=FALSE---------------------------------------------------------------
ModReg1<-lm(CH4~EE+FND+LAD,data=datos)
X<-model.matrix(ModReg1)
y<-datos$CH4
XX<-t(X)%*%X
Xy<-t(X)%*%y
iXX<-solve(XX)

b <- round(coef(ModReg1),2)
signo <- ifelse(sign(b)==1,"+","-")
if (signo[1]=="+") signo[1] <-""
names(b) <- c("",names(b)[-1])
formula <- paste(paste(signo,paste(abs(b),names(b),sep=""),sep=" "),collapse=" ")
#formula=paste(ModReg1$terms[[2]],formula,sep=" = ")


## ----echo=TRUE,eval=FALSE-----------------------------------------------------
## library(rgl)
## library(nlme)
## library(mgcv)
## scatter3d(CH4~EE+FND, data=Purines, surface=TRUE, residuals=TRUE, bg="white",
##   axis.scales=TRUE, grid=TRUE, ellipsoid=FALSE)


## ----echo=FALSE---------------------------------------------------------------
library(car)
modEE <- lm(CH4 ~ EE ,data=Purines)
modFND <- lm(CH4 ~ FND ,data=Purines)
modLAD <- lm(CH4 ~ LAD ,data=Purines)
sal <- compareCoefs(ModReg1,modEE,modFND,modLAD,print=FALSE,se=FALSE)
sal <- sal[c("EE","FND","LAD"),]
tabla <- t(sal)
row.names(tabla) <- c("RLM con EE, FND y LAD","RLS con EE","RLS con FND","RLS con LAD")


## ----echo=FALSE,results='asis'------------------------------------------------
xtable(tabla,caption="Estimas de los coeficientes de regresión de cada una de las variables explicativas en el MRLM con las tres variables explicativas y en los MRLS con cada una de las variables explicativas\\label{tab:tab3.1}")


## ----echo=FALSE---------------------------------------------------------------
tb<-summary(ModReg1)$coefficients


## ----echo=FALSE---------------------------------------------------------------
sm<-summary(ModReg1)


## ----echo=TRUE,eval=FALSE-----------------------------------------------------
## summary(powerTransform(ModReg1, family="bcPower"))


## ----echo=TRUE----------------------------------------------------------------
Purines$tCH4 <- sqrt(Purines$CH4)
ModReg2t <- lm(tCH4 ~ EE +FND+LAD,data=Purines)


## ----echo=TRUE,eval=FALSE-----------------------------------------------------
## crPlots(ModReg2t, smooth=list(span=0.5))


## ----echo=TRUE----------------------------------------------------------------
Purines$tEE <- log(Purines$EE)
ModReg3t <- lm(tCH4 ~ tEE +FND + LAD, data=Purines)


## ----echo=TRUE----------------------------------------------------------------
Purines$FND2 <- Purines$FND**2
ModReg4t <- lm(tCH4 ~ tEE +FND+ FND2 + LAD,data=Purines)


## ----echo=FALSE---------------------------------------------------------------
tabla <- signif(summary(ModReg4t)$coefficients,4)
R2adj <- summary(ModReg4t)$adj.r.squared


## ----echo=FALSE---------------------------------------------------------------
medFND <- signif(mean(Purines$FND),3)
medLAD <- signif(mean(Purines$LAD),3)
k0 <- tabla[1,1]+ tabla[3,1]*medFND +
tabla[4,1]* medFND**2 +tabla[5,1]* medLAD
k0 <- signif(k0,3)


## ----echo=FALSE---------------------------------------------------------------
medtEE <- signif(mean(Purines$tEE),3)
k0 <- tabla[1,1]+ tabla[2,1]*medtEE +tabla[5,1]* medLAD
k0 <- signif(k0,3)


## ----echo=FALSE---------------------------------------------------------------
k0 <- tabla[1,1]+ tabla[2,1]*medtEE +tabla[3,1]* medFND+tabla[4,1]*medFND**2
k0 <- signif(k0,3)


## ----echo=FALSE,results='asis'------------------------------------------------
x <- seq(2.9,25,by=5)
tx <- log(x)
ty <- k0+tabla["tEE",1]*tx
y <- ty**2
tabla2 <- signif(cbind(x,tx,ty,y),3)
colnames(tabla2) <- c("EE","tEE","tCH4","CH4")
print(xtable(tabla2,caption="Valores medios estimados de la variable \\texttt{tCH4} para distinto valores de la variable \\texttt{tEE} y sus valores en la escala original en el modelo de regresión estimado \\texttt{ModReg4t} cuando las variables \\texttt{FND} y \\texttt{LAD} están fijadas en su media\\label{tab:tab3.2}"),include.rownames = FALSE)


## ----echo=TRUE----------------------------------------------------------------
ModReg5t <- lm(sqrt(CH4)~log(EE)+poly(FND, degree = 2, raw = TRUE) + LAD,data=Purines)


## ----echo=TRUE,eval=FALSE-----------------------------------------------------
## plot(predictorEffects(ModReg5t))


## ----echo=TRUE,tidy.opts = list(blank = FALSE, width.cutoff = 50)-------------
Purines$tCH4 <- sqrt(Purines$CH4)
Purines$tEE <- log(Purines$EE)
ind <- sample(c(1,2,3), replace=TRUE,size=79)
 niter<-100
Solmod1<-matrix(0,ncol=3,nrow=niter)
Solmod2<-matrix(0,ncol=3,nrow=niter)
for(i in 1:niter){
  ind<-sample(c(1,2,3),replace=TRUE,size=79)
  for(j in 1:3){
    datf<-Purines[ind==j,]
    datp<-Purines[ind!=j,]
    fit1<-lm( tCH4 ~ LAD+poly(FND,2) + tEE, data = Purines)
    fit2<-lm( tCH4 ~ LAD+poly(FND,3) + tEE, data = Purines)
    pre1<-predict(fit1,datp)
    pre2<-predict(fit2,datp)
    Solmod1[i,j]<-sqrt(sum((datp$tCH4-pre1)^2)/length(pre1))
    Solmod2[i,j]<-sqrt(sum((datp$tCH4-pre2)^2)/length(pre2))
  }
}
modFND2<-round(mean(apply(Solmod1,1,mean)),3)
modFND3<-round(mean(apply(Solmod2,1,mean)),3)
modFND2
modFND3



## ----echo=TRUE----------------------------------------------------------------
#library(car)
library(MASS)
data(fgl)
fgl2<-fgl[,-10]
datos<-fgl2


## ----echo=FALSE---------------------------------------------------------------
EST<-function(X){
   s<-round(c(summary(X),SD=sd(X),n=length(X)),2)
   return(s)
}
apply(fgl2,2,EST)


## ----echo=TRUE----------------------------------------------------------------
scatterplotMatrix(fgl2,diagonal=list(method="boxplot"),smooth=FALSE)


## ----echo=TRUE,figure=TRUE,out.width="75%"------------------------------------
library(corrplot)
k<-cor(fgl2)
corrplot(k, type="upper", tl.col="black", tl.srt=45)


## ----echo=TRUE----------------------------------------------------------------
vary<-colnames(datos)[-1]
sol<-matrix(0,nrow=length(vary),ncol=4)
colnames(sol)<-c("lower","upper","p.beta","p.beta.adjust")
rownames(sol)<-vary
#Bucle que ajusta modelos variable por variable y recoge IC y p valor
for(i in 1:length(vary)){
  m<-lm(formula(paste("RI~",vary[i],sep="")),data=datos)
  sol[i,1:2]<-confint(m)[2,]
  sol[i,3]<-summary(m)$coefficients[2,4]
}
sol[,4]<-p.adjust(sol[,3])
round(sol,4)
#Ahora selecciono las variables que tienen un p valor corregido menor de 0.05
rownames(sol)[sol[,4]<0.05]
var.f<-rownames(sol)[sol[,4]<0.05]
mf<-lm(formula(paste("RI~",paste(var.f,collapse="+",sep=""))),data=datos)


## ----echo=FALSE,results='asis'------------------------------------------------
#print(xtable(mf,capion="Estimas del modelo ajustado"),digits=3)


## ----echo=TRUE----------------------------------------------------------------
set.seed(54321)
nvar<-100
nobs<-100
dat<-matrix(rnorm(101*100),ncol=nvar+1,nrow=nobs)
colnames(dat)<-c("y",paste("X",1:100,sep=""))
dat<-data.frame(dat)
#Usaremos el código anterior para lanzar 100 modelos
vary<-colnames(dat)[-1]
sol<-matrix(0,nrow=length(vary),ncol=3)
colnames(sol)<-c("lower","upper","p.beta")
rownames(sol)<-vary
#Bucle que ajusta modelos variable por variable y recoge IC y p valor
for(i in 1:length(vary)){
  m<-lm(formula(paste("y~",vary[i],sep="")),data=dat)
  sol[i,1:2]<-confint(m)[2,]
  sol[i,3]<-summary(m)$coefficients[2,4]
}
rownames(sol)[sol[,3]<0.05]


## ----echo=TRUE----------------------------------------------------------------

modT<-lm(RI~.,data=datos)
k<-step(modT,trace=0)# Provad sin trace=0
#El mejor modelo esta guardado en k$terms
modTf<-lm(formula(k$terms),data=datos)
summary(modTf)
ModReg1<-lm(RI~Na+Al+Si+K+Ca,data=datos)


## ----echo=TRUE----------------------------------------------------------------
ind <- sample(c(1,2,3), replace=TRUE,size=214)
 niter<-100
Solmod1<-matrix(0,ncol=3,nrow=niter)
Solmod2<-matrix(0,ncol=3,nrow=niter)
for(i in 1:niter){
  ind<-sample(c(1,2,3),replace=TRUE,size=79)
  for(j in 1:3){
    datf<-datos[ind==j,]
    datp<-datos[ind!=j,]
    fit1<-lm( RI~Na+Al+Si+K+Ca, data = datos)
    fit2<-lm( RI~Na+Mg+K+Ca+Ba, data = datos)
    pre1<-predict(fit1,datp)
    pre2<-predict(fit2,datp)
    Solmod1[i,j]<-sqrt(sum((datp$RI-pre1)^2)/length(pre1))
    Solmod2[i,j]<-sqrt(sum((datp$RI-pre2)^2)/length(pre2))
  }
}
mod.bruto<-round(mean(apply(Solmod1,1,mean)),3)
mod.AIC<-round(mean(apply(Solmod2,1,mean)),3)
mod.bruto
mod.AIC


## ----echo=FALSE---------------------------------------------------------------
par(mfrow=c(2,2))
plot(modTf)

############################
genera_datos <- function(dni=12355){
  
  set.seed(dni)   ## Fija la semilla
  n<-sample(150:250,1)
  n_covariables<-sample(3:8,1)
  n_factores<-sample(3:5,1)
  nt<-n_covariables+n_factores
  # Paso 1: Generación de las covariables (simuladas de una distribución normal)
  covariables <- matrix(rnorm(n * n_covariables), nrow = n, ncol = n_covariables)
  colnames(covariables) <- paste("X", 1:n_covariables, sep = "")
  covariables <- round(covariables, 2)
  # Paso 2: Generación de los factores (simulados como factores con 2 niveles)
  factores <- matrix(sample(0:1, n * n_factores, replace = TRUE), nrow = n, ncol = n_factores)
  colnames(factores) <- paste(LETTERS[1:n_factores], 1:n_factores, sep = "")
  
  # Paso 3: Crear una matriz de covariables y factores combinados
  X <- cbind(covariables, factores)
  
  # Paso 4: Definir los coeficientes para las covariables y los factores
  # Coeficientes de las covariables (se asume que afectan el logaritmo del tiempo de supervivencia)
  coef_covariables <- sample(seq(0,1.5*n_covariables),n_covariables,replace=TRUE)
  beta0<- n_covariables
  # Coeficientes de los factores
  coef_factores <- sample(seq(0,1.5*n_factores),n_factores,replace=TRUE)
  
  # Paso 5: Calcular y basado en las covariables y los factores
  # y = β0 + β * X + βf * Z donde Z_j son los factores
  y <- beta0+covariables%*%coef_covariables + factores %*% coef_factores
  y <- y+round(rnorm(n,mean=0,sd=(n_factores)),3)
  # Paso 6: Crear el conjunto de datos con tiempos y censura
  datos <- data.frame(y,X)
  
  # Paso 7: Agregar los factores como tal, al conjunto de datos
  for(i in 1:n_factores){
    datos[,(1+n_covariables+i)]<-factor(datos[,(1+n_covariables+i)])
  }
  
  sal<-list(datos=datos,coef.var=coef_covariables,coef.fact=coef_factores,
            error=(n_factores))
  
}
##################
##################
modelo<-function(df,y,x){
  # Ejemplo syntaxis
  # modelo("datosAustriacaza.xlsx","Hares",c("Sunflowers","Foxes"))
  library(readxl)
  library(corrplot)
  datos <- read_excel(df)
  datos<-data.frame(datos)
  k<-dim(datos)
  
  cat(paste(y,"~",paste(x,collapse="+")),sep="\n")
  cat("Variables explicativas",length(x),"\n")
  cat("observaciones:",k[[1]],"\n")
  cat("variables disponibles:",k[[2]],"\n")
  f<-formula(paste(y,"~",paste(x,collapse="+")))
  ft<-formula(paste(y,"~."))
  m<-lm(f,data=datos)
  mb<-step(m,trace=0)
  mb <-lm(formula(mb),data=datos)
  if(k[[1]]>k[[2]]) mt <-lm(ft,data=datos)
  sol<-list()
  sol[[1]]<-round(summary(m)$coefficients,3)
  sol[[2]]<-round(summary(m)$adj.r.squared,3)
  sol[[3]]<-round(summary(mb)$coefficients,3)
  sol[[4]]<-round(summary(mb)$adj.r.squared,3)
  C<-datos[,c(y,x)]
  M<-cor(C)
  corrplot(M, method = 'ellipse', order = 'AOE', type = 'upper')
  plot(mb)
  return(sol)
}
