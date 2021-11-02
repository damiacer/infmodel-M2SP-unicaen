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