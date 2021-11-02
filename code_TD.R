setwd('/Users/damianocerasuolo/Desktop/UBRC/2021_UBRC_Cours/M2-MODEL/TD')
source("functions_models.R")

# FUNCTION SIR  
# ELLE SERA DEFINI PAR LES PARAMETRES QUI SUIVENT DANS CE MODELE
SIR<-function(t,x,para){
        S<-x[1]
        In<-x[2]
        R<- x[3]
        N<-sum(x)
        
        b<-para["b"]
        g<-para["g"]
        
        dS<--b*S*In/N
        dIn<-b*S*In/N-g*In # g = ??? 
        dR<-g*In
        res <-c(dS,dIn,dR)
        list(res)
}

library('deSolve')

# TAILLE DE LA POPULATION
N<-1000000000
# TAILLE DU COMPARTIMENT DES INFECTIEUX
In0<-10
# BETA = TAUX DE CONTACTS EFFICACES DANS L'UNITE DE TEMPS 
# RAPPEL ; BETA = R0/ND
# OU D = DUREE DE LA CONTAGIOSITE (PAS DE L'EPIDEMIE !)
b<-2
# DUREE EXPRIMEE EN JOURS
Duree<-5

# VECTEUR DU 'POINT DE DEPART' DES DONNES
param<-c(b=b, g=1/Duree) 
# S = N - In0 : TAILLE DE LA POPULATION SUSCEPTIBLE AU TEMPS ZERO (IL Y A CLAIREMENT DEJA QUELQU'UN MALADE !)
# I = In0 : TAILLE DE LA POPULATION AU TEMPS ZERO
# R = R EFFECTIF AU TEMPS ZERO
x_init<-c(S=N-In0,I=In0,R=0)

# CREATION D'UN VECTEUR DE 0 A 40 AVEC UN PAS DE 0.1
times <- seq(0, 40, by=0.1) 
# REAFFFICHAGE DE LA FENETRE GRAPHIQUE
# INTEGRATION DU SYSTEME
# 'INIT' EST LA FONCTION INITIALE, 'PARAM' EST LES VECTEUR DES PARAMETRES DU MODELE
windows()
dynamique<-data.frame(lsoda(x_init, times, SIR, param)) 

valBeta<-round(b,1)
valGamma<-round(1/Duree,1)

plot(dynamique$time,dynamique$I,type="l",col=2,lwd=2,ylim=range(dynamique[,c("S","I","R")]),xlab="Temps",ylab="Individus",main=bquote(paste(beta==.(valBeta), "/jour, ",gamma==.(valGamma),"/jour")))
lines(dynamique$time,dynamique$S,col=4,lwd=2)
lines(dynamique$time,dynamique$R,col=1,lwd=2)
legend("right",c("S","I","R"),col=c(4,2,1),lwd=2)

###########################################################################################################################################

# INSTALLATION DU PACKAGE 'deSolve'
# (SI VOUS L'AVEZ DEJA FAIT >>> ETAPE SUIVANTE)
install.packages('deSolve')

# RAPPEL DU PACKAGE
library('deSolve')

sirmodel=function(t, y, parms){
     S=y[1]
     I=y[2]
     R=y[3]
     
     beta=parms['beta']
     mu=parms['mu']
     gamma=parms['gamma']
     N=parms['N']
     
     dS=mu * (N-S) - beta *S * I/N
     dI=beta * S * I/N - (mu+gamma) * I
     dR=gamma * I-mu * R
     res=c(dS,dI,dR)
     # GRADIENTS
     list(res)
}

# NOUS ALLONS UTILISER LA FONCTION 'ode'
# CETTE FONCTION PERMET DE RESOUDRE LES EQUATIONS DIFFERENTIELLES 
# POUR CETTE FONCTION, NOUS ALLONS DEFINIR
#   times = LE TEMPS, EASY
#   parms = LES PARAMETRES
#   start = LA SITATION QUI SE PRESENTE AU DEBUT DE L'EPIDEMIE

# NOUS ALLONS MODELISER 
times=seq(0, 25, by=1/10)
parms=c(mu=0, N=1, beta=2, gamma=1/2)
start=c(S=0.999, I=0.001, R=0)

# NOUS ALLONS TRANSFORMER LE RESULTAT DE NOTRE MODELISATION DANS UN DATA FRAME (PLUS LISIBLE)
out=ode(y=start, times=times, func=sirmodel, parms=parms)
out=as.data.frame(out)
# names(out)
# COMPARTIMENTS AUX DIFFERENTS TEMPS
# round = ARRONDI
# out = NOS OUTPUTS
# 3 = CHIFFRES DE L'ARRONDI
head(round(out, 3))

# VISULISATION GRAPHIQUE DU MODELE OBTENU
plot(x=out$time, y=out$S, ylab='fraction', xlab='temps', type="l")
# LE COMPARTIMENT DES SUSCEPTIBLES
lines(x=out$time, y=out$S, col='red')
# LE COMPARTIMENT DES RETIRED
lines(x=out$time, y=out$R, col='green')
# LE COMPARTIMENT DES INFECTIEUX
lines(x=out$time, y=out$I, col='blue')
legend("right", legend=c("S", "I", "R"), 
        lty=c(1,1,1),col=c("black", "red", "green"))

# COMMENT CALCULER LE R0 ?
# R0 EST DEFINI COMME LE NOMBRE DE INFECTIONS SECONDAIRES A PARTIR D'UN SEUL SUJET MALADE
# DANS UNE POPULATION CONSTITUE EXCLUSIVEMENT DE SUJETS SUSCEPTIBLES
# NOUS DEFINISSONS 
#   beta = TAUX DE CONTACTS EFFICACES DANS L'UNITE DE TEMPS 
#   gamma = SUJETS QUE DEVIENNENT 'RECOVERED'
#   *NEW ENTRY!* mu = DEATH RATE
R0 = parms['beta']/(parms['gamma'] + parms['mu'])

# AJOUT DU RE (LE R EFFECTIF) DE L'EPIDEMIE 
par(mar=c(5,5,2,5))
# VARIABLES D'ETAT
plot(x=out$time, y=out$S, ylab="fraction", xlab="temps", type="l")
lines(x=out$time, y=out$I, col="red")
lines(x=out$time, y=out$R, col="green")
xx=out$time[which.max(out$I)]
lines(c(xx,xx), c(1/R0,max(out$I)), lty=3)

# DEUXIEME PLOT (RE PLOT !)
par(new=TRUE)
plot(x=out$time, y=R0*out$S, type="l", lty=2, lwd=2,
     col="black", axes=FALSE, xlab=NA, ylab=NA, ylim=c(-.5,4.5))
lines(c(xx, 25, c(1,1), lty=3))
axis(side=4)
mtext(side=4, line=4, expression(R[E]))
legend("right", legend=c("S", "I", "R", expression(R[E])), lty=c(1,1,1,2), 
       col=c("black", "red", "green", "black"))

# CALCULER LA TAILLE DE L'EPIDEMIE 
install.packages('rootSolve')
library('rootSolve')

equil = runsteady(y=c(S=1-1E-5, I=1E-5, R=0), times=c(0,1E5), func=sirmodel, parms=parms)
round(equil$y,10)

# runsteady = RESOLUTION DE LA CONDITION DES EQUATIONS DIFFERENTIELLES ORDINAIRES ('ode') STEADY-STATE
# SOUS LA FORME : dy/dt = f(t,y)

# VOUS OBTENEZ LE RESULTAT SUIVANT
# COMMENT PEUT-ON L'INTERPRETER ? 
#S            I            R 
#0.0198298498 0.0000000241 0.9801701261 

# COMMENT PEUT-ON EXPLORER LA DEPENDANCE DE LA TAILLE DE L'EPIDEMIE DE SON R0 ?
# DANS LE MODELE ETUDIE R0 = BETA / GAMMA+MU
# MAIS MU = 0 CAR ON EST EN TRAIN D'ETUDIER UNE POPULATION CLOSE 
# ON PREND UNE PERIODE D'INFECTION DE 2 SEMAINES, GAMMA = 1/2
# LA VARIATION DE BETA EST ENTRE 0.1 ET 0.5

R0 = seq(0.1, 5, length=50)
betas = R0 * 1/2 # R0 = beta / gamma+mu = beta = R0 * gamma 
f = rep(NA, 50)
for(i in seq(from=1, to=50, by=1)){
     equil=runsteady(y=c(S=1-1E-5, I=1E-5, R=0), times=c(0,1E5), func=sirmodel, 
                     parms=c(mu=0, N=1, beta=betas[i], gamma=1/2))
     f[i]=equil$y["R"]
}
plot(R0, f, type='l', xlab=expression(R[0]))
curve(1-exp(-x), from=1, to=5, add=TRUE, col='red')

# LA TAILLE FINALE DE L'EPIDEMIE EST UNE FONCTION DE R0
# LA LIGNE NOIR EST UNE SOLUTION DE L'EQUATION DIFFERENTIELLE
# LA LIGNE ROUGE EST UNE APPROXIMATION
# ON EN DEDUIT QUE 
#   POUR R0 > 2.5 L'ESTIMATION DE R0 EST BONNE
#   POUR R0 < 2.5 ON SUR-ESTIME LA TAILLE DE L'EPIDEMIE

# POUR UNE EPIDEMIE CLOSE, ON PEUT CALCULER LA TAILLE EXACTE DE LA POPULATION QUI ECHAPPE 
# L'INFECTION (1-f)
# CETTE EQUATION EST DONNEE PAR exp(-R0(1-f))-f=0

fn = function(x, R0){
     exp(-(R0*(1-x)))-x
}

# POUR UN R0=2 LA TAILLE DE L'EPIDEMIE EST DONNEE
1-uniroot(fn, lower=0, upper=1-1E-9,
          tol=1E-9, R0=2)$root
# APPROX
exp(-2)-uniroot(fn, lower=0, upper=1-1E-9, tol=1e-9, R0=2)$root

###########################################################################################################################################

# ETUDE DU R0 - DONNEE ROUGEOLE
# L'ESTIMATION DU R0 SERA FAITE EN PRENANT EN COMPTE LA CROISSANCE EXPONENTIELLE ('r') ET
# L'INTERVALLE DE LA SERIE ('SERIAL INTERVAL', 'V')
install.packages('epimdr')
library('epimdr')

data(niamey)
head(niamey) # REMARQUEZ SURTOUT LA DERNIERE VARIABLE. ON VA L'UTILISER APRES
head(niamey[,1:5])

# VISUALISATION DE LA PERIODE DE CROISSANCE EXPONENTIELLE
par(mar=c(5, 5, 2, 5))
plot(niamey$absweek, niamey$tot_cases, type='b', xlab='semaine', ylab='incidence')
par(new=T)
plot(niamey$absweek, niamey$cum_cases, type='l', col='red', axes=FALSE, xlab=NA, ylab=NA, log='y')
axis(side = 4)
mtext(side = 4, line = 4, 'incidence cumulee')
legend('topleft', legend=c('cases', 'cumulee'), lty=c(1, 1), pch=c(1,NA), col=c('black', 'red'))

# L'INCIDENCE CUMULEE SEMBLE ETRE LOG-LINEAIRE SEULEMENT DANS LA PREMIERE PHASE DE L'EPIDEMIE 
# PUIS, LE CHAOS...
# 'SERIAL INTERVAL' = TEMPS MOYEN ENTRE INFECTION (1ER CAS) ET REINFECTION (2EME CAS)
# COMMENT ON CALCULE LE R0 ? 
# PAR LA REGRESSION (LINAIRE) DU LOG DE L'INCIDENCE CUMULEE (ELLE EST EXPONENTIELLE... ON NE POURRAIT
# PAS LA REGRESSER PAR UN MODELE LINEIARE SINON !)
# LE SERIAL INTERVAL DE LA ROUGEOLE EST DE 10-12 JOURS
# DONC LE SERIAL INTERVAL ('V') SERA DE 1.5 ET 1.8 SEMAINES

fit = lm(log(cum_cases) ~ absweek, subset=absweek<7, data = niamey)
r=fit$coef['absweek']
V=c(1.5, 1.8)

R0_rougeole = V * r + 1
# [1] 1.694233 1.833080
# VOICI LE RAPPORT DE REPRODUCTION DE CETTE EPIDEMIE (ON LE LIRA PLUTOT 1.5 ET 2 POUR LE GRAND PUBLIC ?)

# Lipsitch et al.
V = c(1.5, 1.8)
f = (5/7)/V
V * r + 1 + f * (1 - f) $ (V * r)^2

# DONNEES FLU
# LA DUREE DE LA MALADIE EST A CONSIDERER COMME COMPRISE ENTRE 5 ET 7 JOURS

data(flu)
head(flu)
plot(flu$day, flu$cases, type='b', xlab="jours", ylab='au lit')

# LES DONNES SONT LOG LINEAIRES ? SI OUI, SUR QUELLE PERIODE ? 
# VOUS POUVEZ UTILISER CETTE INFORMATION POUR VOTRE ESTIMATION DE R0
# ESTIMEZ RAPIDEMENT LE R0 DE LA MALADIE 

fit_flu = lm(log(cases) ~ day, subset=day<=5, data = flu)
r_flu = fit_flu$coef['day']
V = c(2,4)
R0_flu = V * r_flu + 1

# DONNES EBOLA
data(ebola)
head(ebola)
mean(ebola$cases)

par(mar=c(5, 5, 2, 5))
plot(ebola$day, ebola$cases, type='b', xlab='semaine', ylab='incidence')
par(new=T)
plot(ebola$day, ebola$cum_cases, type='l', col='red', axes=FALSE, xlab=NA, ylab=NA, log='y')
axis(side = 4)
mtext(side = 4, line = 4, 'incidence cumulee')
legend('topleft', legend=c('cases', 'cumulee'), lty=c(1, 1), pch=c(1,NA), col=c('black', 'red'))

V = 15
f = 0.5
fit=lm(log(cases)~day, subset = day<=14, data = ebola)
lambda=fit$coef['day']
V * lambda + 1 + f * (1-f) * (V*lambda)^2

cases=sapply(split(ebola$cases, floor((ebola$day-.1)/14)),sum)
sum(cases)
install.packages('bbmle')
library('bbmle')
fit=mle2(llik.cb, start=list(S0=20000, beta=2), 
         method='Nelder-Mead', data=list(I=cases))
summary(fit)
confint(fit, std.err=c(100,0.1))

###########################################################################################################################################

install.packages('splines')
install.packages('fields')

library('splines')
library('fields')

# FoI = FORCE DE L'INFECTION 
# ELLE CORRESPOND AU TAUX PER CAPITA AUQUEL LES SUSCEPTIBLES SONT EXPOSES A L'INFECTION
# DANS LE MODEL PLUS SIMPLE FoI EST INDEPENDENATE DE L'AGE E DU TEMPS 
# DANS CE CAS, LA PROBABILITE D'ETRE INFECTE A L'AGE a EST 
# 1-exp(-phi(a)
# SI ON CONNAIT LE NOMBRE DES SUJETS INFECTES EN FONCTION DE L'AGE NOUS POUVONS ESTIMER LA FORCE DE
# L'INFECTION PAR UN MODELE LINEAIRE GENERALISE (glm)

# LES glm ONT DEUX COMPOSTANTES 
# UNE DISTRIBUTION DE L'ERREUR (QUI PEUT ETRE 'BIONMIAL', 'POISSON', 'NEGATIVE BINOMIAL', 'NORMAL')
# ET UNE FONCTION QUI SPECIFIE COMMENT LES VALEURS ATTENDUES (PREDITES) SONT LIEES AU VALEURS PREDITES 
# LINEAIRES x = a + b1x1 + b2x2 + ... 

# glm(cbind(inf, noninf) ~ offset(log(a)), family=binomial(link="cloglog"))

# LA FONCTION POUR UNE POPULATION nA SUJETS, D'AGE a
# LA DISTRIBUTION (SIMPLE) QU'ON CHOISIT EST BIEN LA BINOMIALE 

# ON REPREND L'EXEMPLE DE LA ROUGEOLE MAIS SUR UNE AUTRE BASE DE DONNEES

data(black)
black
# f EST LA SEROPREVALENCE 
# ELLE EST TRES ELEVEE CHEZ LES PETITS - GRACE AUX ANTICORPS DE LA MERE
# ET REDEVIENT MAXIMALE AUTOURS DE 20 ANS

# CREATION DES GROUPES EN FONCITION DE L'AGE
b2=black[-c(1,8,9)]
fit = glm(cbind(pos, neg) ~ offset(log(mid)), family = binomial(link = 'cloglog'), data = b2)

phi = exp(coef(fit))
curve(1-exp(-phi*x), from=0, to=60, ylab='seroprev', xlab='age')
points(black$mid, black$f, pch='*', col='black')
points(b2$mid, b2$f, pch=8)
exp(fit$coef)

# FoI
# (Intercept) 
# 0.1116066 

###########################################################################################################################################

# LA RUGEOLE EST UNE MALADIE QUI PRESENTE UNE PROGRESSION PLUTOT BENIGNE
# SAUF PENDANT LA GROSSESSE - SYNDROME DE LA RUGEOLE CONGENITALE
# UN OBJECTIF IMPORTANT DE SANTE PUBLIQUE EST DONC DE MINIMISER LA FORCE OF INFECTION CHEZ LES FEMMES
# EN AGE DE PROCREER 

# NOUS CONSIDERONS DES DONNES PRESENTANT L'EVOLUTION D'UNE EPIDEMIE DE RUGEOLE 'AGE INTENSITY'
# AU PERU, ENTRE 1997 ET 2009, 24116 CAS DE ROUGEOLE ONT ETE DECLARES 
# LES DONNEES SONT CHARGEES DANS 'PERU'

install.packages('splines')
library('splines')

data(peru)
peru
head(peru)

# CALCUL DE L'INCIDENCE CUMULEE
peru$cumulative=cumsum(peru$incidence)
peru$cumulative
# DEFINITION DU DENOMINATEUR
peru$n=sum(peru$incidence)

par(mar = c(5, 5, 2, 5))
plot(peru$incidence ~ peru$age, type = 'b', ylab = 'incidence', xlab = 'age')
par(new = T)
plot(peru$cumulative ~ peru$age, type = 'l', col = 'red', axes = FALSE, xlab = 'NA', ylab = 'NA')
axis(side = 4)
mtext(side = 4, line = 4, 'cumulative')
legend('right', legend = c('incidence', 'cumumative'), lty = c(1,1), col = c('black', 'red'))

# CUT OFF SUPERIEURS
up=c(1:20, 30, 40, 50, 60, 70, 100)
para=rep(.1, length(up))
# LOG-VRAISAMBLANCE 
est2=optim(par=log(para),fn=llik.pc, age=peru$age, num=peru$cumulative, denom=peru$n, up=up, 
           method='Nelder-Mead', control = list(trace=2, maxit=2000))
# PLOT 
x=c(0,up)
y=exp(c(est2$par, est2$par[26]))
plot(x, y, ylab='FoI relative', xlab='age', type='l', ylim=c(0,0.25), xlim=c(0,80))

#####################################

# CONSEQUENCES DE LA VACCINATION DANS LA TRANCHE D'AGE 0-45 ANS 

data3 = peru[peru$age <45, ]
head(data3)
df = 5
para = rep(0.1, df +1)

spline_fn1 = function(x,dl){
     x=predict(dl, newdata=data.frame(age=x))
     exp(x)
}

dl=lm(cumulative ~ bs(age,df), data=data3)

spline_fn2 = function(par, data, df){
     dl=lm(cumulative ~ bs(age, df), data=data)
     dl$coefficients=par
     ll=0
     for(a in 1:length(data$age)){
          p=((1-exp(-integrate(spline_fn1, 0, data$age[a], dl=dl)$value)))
          ll=ll+dbinom(data$cumulative[a], data$n[a], p, log=T)
     }
     return(-ll)
}

spline_model = optim(par=log(para), fn=spline_fn2,
                     data=data3, df=df, method='Nelder-Mead',
                     control = list(trace=4, maxit=5000))

dl$coefficients=spline_model$par
plot(exp(predict(dl))~data3$age, xlab='age', ylab='FoI relative', type='l')

###########################################################################################################################################

# VARIATION DE R0
R0<-c(0.1,0.5,1,1.5,2,3)
LD<-5 # LD = D
N<-1000
I0<-1
par(mfrow=c(2,3))
for (i in 1:length(R0))
{
        Plot.Several.Gillespie(10,R0[i]/LD,1/LD,10,I0,N-I0,N)
}


# VARIATION DE I0
R0<-2
LD<-5
N<-1000
I0<-c(1,10,20,100)
par(mfrow=c(2,2))
for (i in 1:length(I0))
{
        Plot.Several.Gillespie(100,R0/LD,1/LD,10,I0[i],N-I0[i],N)
}

# VARIATION DE N
R0<-2
LD<-5
N<-c(10,100,1000,10000)
I0<-10
par(mfrow=c(2,2))
for (i in 1:length(N))
{
        Plot.Several.Gillespie(10,R0/LD,1/LD,10,I0,N[i]-I0,N[i])
}


##on compare au prediction du modele déterministe
#windows()
par(mfrow=c(1,2))
for (R0 in c(0.9,2.1))
{
        LD<-5
        N<-10000
        I0<-10
        nbjours<-30
        
        Plot.Several.Gillespie(30,R0/LD,1/LD,nbjours,I0,N-I0,N)
        ##on trace la prediction du modele déterministe
        
        In0<-I0
        b<-R0/LD
        Duree<-LD
        param<-c(b=b,g=1/Duree)
        x_init<-c(N-In0,In0,0)
        names(x_init)<-c("S","I","R")
        times<-seq(0,nbjours,by=0.1)
        ##Intégration du système
        dynamique<-data.frame(lsoda(x_init,times,SIR,param))
        lines(dynamique$time,dynamique$I,lwd=4,col=1)
}