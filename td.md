# modélisation maladies infectieuses

Le package *deSolve* [^1] est utilisé pour résoudre les équations *ODE* (*ordinary differential equations*).\
L'installation d'un package sur *R* se fait de la façon suivante :

[^1]: Document de référence : [deSolve - CRAN](https://cran.r-project.org/web/packages/deSolve/deSolve.pdf)

```
install.packages("deSolve")
library("deSolve")
```

## SIR MODEL

### Construction du modèle 

Bien que les modèles puissent être différents entre eux, la démarche à suivre est proche.

#### étape 1 : la fonction

Nous définissons l'équation (la fonction gradient) du système d'équations. Le package *deSolve* nécessite la présence des paramètres :

* le temps, *t* 
* *y*, un vecteur avec les valeurs des variables d'état (*S*, *I*, *R*)
* les valeurs des paramètres ($\beta$, $\mu$, $\gamma$ et $N$), *parms*

```
sirmod = function(t, y, parms){
        # variables d'état pour le vecteur y
        S = y[1]
        I = y[2]
        R = y[3]
       
        # paramètres
        beta = parms["beta"]
        mu = parms["mu"]
        gamma = parms["gamma"]
        N = parms["N"]
        
        # définition des équations
        dS = mu * (N-S) - beta * S * I/N
        dI = beta * S * I/N - (mu + gamma) * I
        dR = gamma * I - mu * R
        res = c(dS, dI, dR)
        
        # liste des gradients
        list(res)
}
```

Les fonctions *ODE* résolvent les équations numériquement.

#### étapes 2-4 : time points

Nous définissons à quel intervalle de temps les *ODE* enregistrent l'état du système. Dans l'exemple, nous utilisons le temps de 26 semaines, avec un pas de 10. La valeur des paramètres est définie par le vecteur *parms* et les conditions au début de l'épidémie.\
Pour cet exemple, nous considérons : 

* $N=1$, la fraction de sujets dans chaque classe
* une période infective de 2 semaines, $\gamma = 1/2$
* l'absence de décès ou naissances, $\mu = 0$
* un taux de transmission $\beta = 2$

````
times = seq(0, 26, by = 1/10)
parms = c(mu = 0, N = 1, beta = 2, gamma = 1/2)
start = c(S = 0.999, I = 0.001, R = 0)
````

#### étape 5 : input des valeurs

Cette étape permet de transformer le résultat dans un dataframe lisible, que nous pourrons visualiser par la fonction **head**. Cette fonction permet de montrer les 5 premières lignes du dataframe. **round (3)** nous permettra d'arrondir les valeur à 3 décimales.

````
out = ode(y = start, times = times, func = sirmod, parms = parms)
out = as.data.frame(out)
head(round(out, 3))
````

Nous pouvons donner aux résultats une représentation graphique. Il faut remarquer la présence des différents compartiments : *plot( ... y = out\$S)*, pour le compartiment des susceptibles, auquel nous ajoutons *lines( ... y = out\$I)*, pour le compartiment des contagieux, et *lines( ... y = out\$I)*, pour le compartiment des retirés.

````
plot(x = out$time, y = out$S, ylab = "fraction", xlab = "time", type = "l")
lines(x = out$time, y = out$I, col = "red")
lines(x = out$time, y = out$R, col = "green")
legend("right", legend=c("S", "I", "R"), 
        lty=c(1,1,1),col=c("black", "red", "green"))

````

Nous pouvons ajouter à cette image la droite correspondante à $R_E$, le taux de reproduction effectif. Cette quantité est donnée par : 

$R_E = R_0 s = 1$

où $s$ est la proportion des susceptibles.

Si la vaccination permet de tenir la proportion susceptible en dessous du seuil de $p_c = 1-1/R_0$, la diffusion du pathogène est dissipée et il ne sera pas capable de diffuser dans la population. Il est possible d'ajouter cette information à notre graphique.

On calcule $R_0$ : 
````
R0 = parms["beta"] / (parms["gamma"] + parms["mu"])
````
Le graphique : 

````
# définition des marges
par(mar = c(5,5,2,5))
# plot
plot(x = out$time, y = out$S, ylab = "fraction", xlab = "time", type = "l")
lines(x = out$time, y = out$I, col = "red")
lines(x = out$time, y = out$R, col = "green")
# ligne verticale au turnover
xx = out$time[which.max(out$I)]
lines(c(xx, xx), c(1/R0, max(out$I)), lty = 3)
````
Auquel nous pouvons ajouter le $R_E$ [^2] :

[^2]: R_E et **transmission** : la fréquence d'une maladie infectieuse est souvent définie par la *prévalence* (la proportion de sujets atteints par la maladie au sein d'une population, à tout temps *t*) et l'*incidence* (risque ou taux de nouveau cas per unité de temps per 1000 individus à risque ou susceptibles). Ces informations, bien que indispensables, ne nous permettent pas de définir à quel point une infection est transmissible.\ 
Pour définir la transmissibilité d'une infection, nous utilisons le *secondary attack rate* et le *reproduction number*, $R_0$, le nombre (moyen) de transmission efficaces par personne. La valeur de $R_0$ est élevée quand une personne contagieuse est introduite dans une population naïve. Néanmoins, la présence de l'infection au sein de notre population déterminera une diminution de la proportion des susceptibles (immunisation). Le nouveau *reproduction number* sera représenté par le $R$ effectif ou *net* : $R_E = R_0 \times s$.

````
par(new = TRUE) # superposition du plot 
plot(x = out$time, y = R0*out$S, type = "l", lty = 2, lwd = 2,
  col = "black", axes = FALSE, xlab = NA, ylab = NA,
  ylim = c(-.5, 4.5))
lines(c(xx, 26), c(1,1) tly = 3)

# légende
legend("right", legend = c("S", "I", "R",
  expression(R[E])), lty = c(1,1,1,2),
  col = c("black", "red", "green", "black"))
````

### La taille finale de l'épidémie 

Le modèle SIR fermé a deux situations d'équilibre : 

* ${S=1, I=0, R=0}$, instable quand $R_0 >1$, et 
* ${S*, I*, R*}$, qui reflet l'équilibre de l'épidémie à sa fin[^3]

[^3]: Où : I*=0 si l'épidémie est en train de s'éteindre en absence d'un recrutement efficace de sujets susceptibles ; $S*$ est la taille de la population susceptible qui échappe à l'infection et $R*$ est la taille finale de l'épidémie.

Le package **rootSolve** [^4] cherche de trouver l'équilibre des systèmes de équations différentielles par l'intégration. Installation :  

[^4]: Document de référence : [rootSolve - CRAN](https://cran.r-project.org/web/packages/rootSolve/rootSolve.pdf)

````
install.packages("rootSolve")
library("rootSolve")
````

La fonction **runsteady** est une fonction *"wrapper"* autour des ODE, qui intègre jusqu'au moment où elle trouve un *steady-state* (si cela existe). 

````
equil = runsteady(y = c(S=1-(1e-5), I=(1e-5), R=0),
  times = c(0, 1e5), func = sirmod, parms = parms)
round(equil$y, 3)
````

Le résultat de runsteady (qui utilise les paramètres, *parms*, et la fonction, *sirmod*, précédemment écrite) nous permet de connaître le nombre de sujets qui échapperont à l'infection ou qui seront infectés (et *retired* du système).\


Dans notre épidémie, $\mu = 0$ (il n'y a pas de décès). Dans l'exemple, la période de contagiosité était égale à $2$ ($\gamma = 1/2$). Cela nous permet de faire varier $\beta$ afin que $R_0$ soit compris entre 0.1 et 5. Pour des $R_0$ de grandeur modérée, $f = 1-exp(-R_0)$ [^5].

[^5]: *f = final epidemic size*. Anderson and May, 1982

````
# valeurs de R0 et des betas
R0 = seq(0.1, 5, length = 50)
betas = R0 * 1/2

# vecteur de données manquantes à remplacer par des valeurs 
f = rep(NA, 50)

for(i in seq(from = 1, to = 50, by = 1)){
    equil = runsteady(y = c(S=1-(1E-5), I=(1E-5), R=0),
    times = c(0, 1E5), func = sirmod, 
    parms = c(mu = 0, N = 1, beta = betas[i], gamma = 1/2))
    f[i] = equil$y["R"]
}
````

Le graphique: 

````
plot(R0, f, type = "l", xlab = expression(R[0]))
curve(1-exp(-x), from = 1, to = 5, add = TRUE, col = "red")

````

Interprétation de l'image : l'approximation est valable pour $R_0>2.5$ mais surestime la taille finale de l'épidémie pour des valeurs de $R_0<2.5$ (en particulier pour $R0<1$).

***

## R0

Pour les maladies transmissibles directement, $R_0$ est défini comme le nombre de cas secondaires qui dérivent d'un cas *index* dans une population susceptible.
$R_0$ joue un rôle fondamental dans différents aspects de la dynamique de la maladie. Pour des pathogènes responsables d'une immunisation parfaite : 

* Le seuil nécessaire afin que le pathogène puisse déclencher la maladie. Quand $R_0>1$, un pathogène peut infecter. Quand $R_0<1$, la chaine de transmission ralenti et interrompue [^6]. Ce seuil peut être associé au *critical host density* ;
* Le seuil nécessaire afin que le vaccin puisse induire la *herd immunity*. Ce seuil est calculé comme $p_c = 1-1/R_0$. Dans le cas de la *rougeole*, $R_0=20$ : il est indispensable vacciner 95% de la population pour éliminer la maladie. 

[^6]: Lloyd-Smith, 2009

### Estimation de R0 pour une épidémie simple

Une large variété de méthodes a été proposée pour estimes $R_0$ (ou le *effective reproductive ratio* $R_E$[^7][^8]). Certaines méthodes sont fondées exclusivement sur des modèles.

L'hypothèse la plus simple est d'estimer que, pendant la première phase de diffusion d'un pathogène, la réduction de la fraction susceptible est négligeable et que, conséquemment, la croissance de l'épidémie puisse être considérée **exponentielle**. Le taux de croissance exponentielle peut être traduit par $r=log(R_0)/G$, où $G$ représente le temps-génération. Cela nous permet de calculer $R_0=exp(rG)$.
La croissance exponentielle d'une population peut être calculée par $N(t) = N(0)exp(rt)$. Le temps nécessaire pour doubler sa taille est $log(2)/r$.

Nous pouvons appliquer ces formules à la phase précoce d'une épidémie.
$N$ représente la *prévalence* et $G$ l'intérvalle sérial ($V$), le temps moyen entre infection et reinfection. Le taux de croissance de l'*incidence* (plus souvent represéntée dans les données ou les résultat) est toujour exponentiel. Pour calculer $R_0$, nous pouvons régresser *log(incidence cumulée) sur le temps* (calcul du taux de croissance exponentielle ($r$)) et, ensuite, $R_0 = Vr +1$ [^9].  

[^7]: $R_E$ le nombre effectif de cas dans une population partiellement immunisé :  $R_E = s \times R_0$.
[^8]: $R_E = R_N$ 
[^9]: Anderson et May, 1991

Pour illustrer ces concepts, nous allons utiliser le package *epimdr*[^10] et les *données Niamey* [^11].

[^10]: Document de référence : [epimdr - CRAN](https://cran.r-project.org/web/packages/epimdr/epimdr.pdf)
[^11]: *Weekly measles incidence from 2003-04 in Niamey, Niger*, page 24 de [epimdr](https://cran.r-project.org/web/packages/epimdr/epimdr.pdf)

````
install.packages("epimdr")
library("epimdr")
data(niamey)
````

Pour visualiser les données, la fonction *head* :

````
head(niamey[, 1:5])
names(niamey)
````

Nous étudions graphiquement l'*incidence cumulée* :

```
par(mar = c(5,5,2,5))
plot(niamey$absweek, niamey$tot_cases, type = "b", xlab = "week", ylab = "incidence")
# où absweek = week since beginning of outbreak
# et tot_cases = weekly incidence for the whole city

par(new = T)
plot(niamey$absweek, niamey$cum_cases, type = "l", col = "red", axes = "FALSE", 
  xlab = NA, ylab =NA, log = "y")
# où cum_cases = weekly cumulative incidence for the whole city
axis(side = 4)
mtext(side = 4, line = 4, "Cumulative incidence")
legend("topleft", legend = c("Cases", "Cumulative"), 
  lty = c(1,1), pch = c(1,NA), col = c("black", "red"))
````

L'incidence cumulée semble être log-linéaire pendant la première phase de l'épidémie.
L'intervalle serial pour la rougeole est de 10-12 jours et, étant les données renseignées chaque semaine, $V=1.5-1.8$ semaines. Nous pouvons calculer le $R_0$[^12] en utilisant l'intervalle de 1.5-1.8 : 

[^12]: Dans ce calcul, nous allons appliquer la transformation logarithmique, car l'incidence cumulée est estimée être exponentielle - et, donc, le modèle de régression logistique ne peut pas être appliqué. 

````
fit = lm(log(cum_cases) ~ absweek, subset = absweek < 7, data = niamey)
r = fit$coef["absweek"]
V = c(1.5, 1.8)
R0.rugeole = V * r + 1
```` 
En tenant en compte le fait que le pays (Niger) a mis en place plusieurs campagnes de vaccination, le $R_0$ retrouvé corréspond à une estimation du taux de reproduction *effectif*, $R_E$. 

Lipsitch[^13], en utilisant les données SARS, a montré qu'il faut distinger les périodes de contagiosité et latente : $R = Vr +1 + f(1-f)(Vr)^2$.

[^13]: Lipsitch et al., 2003

````
V = c(1.5, 1.8)
f = (5/7)/V # rapport entre période de contagiosité et serial interval
R0.rugeole2 = (V * r) + 1 + (f * (1 - f)) * ((V * r)^2)
````
### R0 de l'épidémie A/H1N1, 1977

Nous utiliserons les données *flu*[^14] du package *epimdr* : 

[^14]: *flu, Boarding school influenza data*, page 13 de [epimdr](https://cran.r-project.org/web/packages/epimdr/epimdr.pdf)

`````
data(flu)
names(flu)
head(flu)
`````

La durée typique de la maladie est de 5-7 jours.
Nous pouvons tracer le plot de ces données. Il faut remarque que les cas sont définis comme sujet en isolement par jour : il s'agit donc d'une mesure de *prévalence*, et pas de l'incidence. 

````
plot(flu$day, flu$case, type = "b", xlab = "jour", ylab = "en isolement", log = "y")
tail(flu)
````
Notre estimation "rapide" du $R_0$ est :

```
fit = lm(log(cases) ~ day, subset = day<=5, data = flu)
lambda = fit$coef["day"]
V = c(2,3)
R0.FLU = V * lambda + 1
````
### Ebola en Sierra Leone

Nous utiliserons les données *ebola*[^15] du package *epimdr* : 

[^15]: *ebola, Sierra-Leone Ebola 2015 data.*, page 52 de [epimdr](https://cran.r-project.org/web/packages/epimdr/epimdr.pdf)


````
data(ebola)
head(ebola)
names(ebola)
tail(ebola)
````
L'intervalle serial pour Ebola est estimé à 15 jours avec une période d'incubation d'environ 11 jours. Le délai moyen d'hospitalisation (en 5 jours) est de 5 jours et le délai moyen de décès est de 11 jours.

````
par(mar=c(5, 5, 2, 5))
plot(ebola$day, ebola$cases, type = "b", xlab = "semaine", ylab = "incidence")
par(new=T)
plot(ebola$day, ebola$cum_cases, type = "l", col = "red", axes = FALSE, xlab = NA, ylab = NA, log = "y")
axis(side = 4)
mtext(side = 4, line = 4, "incidence cumulee")
legend("topleft", legend = c("cases", "cumulee"), lty=c(1, 1), pch=c(1,NA), col=c("black", "red"))
````

Nous pouvons calculer le $R_0$ en utilisant la correction proposée par Lipsitch (note 13), incluant la période de contagiosité et la période latente : 


````
V = 15
f = 0.5
fit = lm(log(cases) ~ day, subset = day<=14, data = ebola)
lambda = fit$coef["day"]
R0.ebola = (V * lambda) + 1 + (f * (1-f)) * ((V*lambda)^2)
````

***

## FoI et incidence dépendante de l'âge

### *Burden of disease*

Les termes *infection* et *maladie* sont souvent utilisés comme synonymes. Néanmoins, en modélisation des maladies infectieuses, il est nécessaire de bien garder en tête le fait que la maladie correspond exclusivement à la symptomatologie déterminée par la colonisation du pathogène. Il est également nécessaire distinguer entre : 

* la période de latence : la période entre colonisation et moment où le sujet infecté devient contagieux (et peut transmettre l'infection) ;
* l'incubation : la période entre la colonisation et l'apparition des symptômes.

Nous connaissons la différence entre *VIH* et *SIDA*. Pour la "grippe", les périodes peuvent ne pas être simplement identifiables, à cause de la permanence des symtômes. 

La séverité des maladies dépend souvent de l'âge de la personne atteinte. 
L'effet de l'âge peut avoir un impact majeur sur le *burden* de la maladie et sur la dynamique de l'infection.

### **FoI**, force de l'infection
La force de l'infection est le taux *per capita* auquel les sujets susceptibles sont exposés à l'infection.
Dans le modèle à compartiment S(E)IR, la force de 'l'infection est $\phi = \beta I/N$. 

$I/N$ représente la fraction des contacts que chaque sujet susceptible peut avoir avec les autres membres de la même population. $\beta$ est le taux de contacts efficaces (transmission de l'infection). 

*FoI* est en taux, donc dans une population sans différences d'âge (ou dans une population où l'âge "ne compte pas") et dans le cas de contacts aléatoires (*randomly mixing population*), le temps d'infection est $1/\phi$.
Pour une infection qui produit une immunisation permanente, dans une population constante, $R_0 \approx 1+L/\bar{a}$. 

$L$ est l’espérance de vie du sujet atteint. La moyenne d'âge de l'infection est $\bar{a} \approx L/(R_0-1)$. [^16]

[^16]: La population n'est pas fermée : $\bar{a} \approx L/(\mu(R_0-1))$, où $\mu$ est le taux de naissance. Dietz et Schenzle, 1985.

Le taux auquel chaque sujet d'une population peut être infecté est $\phi(a,t)$ et dépend donc de l'âge, $a$, et du temps, $t$. La probabilité d'être infecté à un âge donné est : $p(a) = 1 - e^{-\int_{0}^{a} \phi(a)da}$. Cette formulation de la propbabilité d'infection corréspond au modèle catalytique[^17], qui, dans sa formulation plus simple (*FoI* indépendenate de l'âge et du temps) peut être écrite comme $1-exp(-\phi a)$.

Nous pouvons estimer le *FoI* par un modèle linéaire généralisé (*glm*) si nous disposons du nombre de sujets infectés par classe d'âge. Ce modèle est codé : 

````
glm(cbind(inf, noninf) ~ offset(log(a)), family=binomial(link="cloglog"))
````

[^17]: Muench 1959 et Hens 2010.

Pour expliquer le modèle catalytique, nous allons utiliser les données issues de la rougeole Il s'agit des données *black*[^18] du package *epimdr*.

[^18]: *black, Black’s measles seroprevalence data.*, page 4 de [epimdr](https://cran.r-project.org/web/packages/epimdr/epimdr.pdf)

````
data(black)
black
````

$f$ dans la base de données représente la séroprévalence de l'infection. Elle est très élevée chez les plus petits (d'âge inférieur à 1 an) en raison des anticorps maternels. Cette séroprévalence élevée est suivie par une diminution rapide et une complète protection autour de l'âge de 20 ans. Il peut avoir une perte d'immunité chez les plus âgés. \
Nous allons utiliser le modèle de régression binomiale, après avoir réalisé le *subset* des données en fonction de l'âge :

````
b2 = black[-c(1,8,9)]
fit = glm(cbind(pos, neg) ~ offset(log(mid)), family = binomial(link = "cloglog"), data = b2)
````
Réalisation du graphique : 

````
phi = exp(coef(fit)) # extraction des coefficients
curve(1-exp(-phi*x), from = 0, to = 60, ylab= "seroprev", xlab = "age")
points(black$mid, black$f, pch="*", col = "black")
points(b2$mid, b2$f, pch=8)
exp(fit$coef)
````
On retrouve une *FoI$ de 0.16/year, avec une moyenne d'âge prédite de 6 ans. 

***

## Rubéole

La rubéole est une maladie préentant une progression globalemtn bénigne, à l'éxception de la femme enceinte, produisant un syndrome congénita (*congenital rubella syndrome, CRS*). Le premier objectif de santé publique serait donc de réduire la *FoI* chez les femmes en âge de procréer. 
L'importance de cet objectif a été soulignée par une épidémie de CRS en Grèce, pendant les années 1990. 
Nous allons utiliser les données *peru*[^19] (oui, nous ne sommes plus en Grèce !), issue du package *epimdr*. Nous aurons besoin de quelques packages de plus : 

[^19]: *peru, Rubella in Peru data.*, page 28 de [epimdr](https://cran.r-project.org/web/packages/epimdr/epimdr.pdf)

````
install.packages('splines')
library('splines')

data(peru)
peru
head(peru)
# names(peru)
````
Nous allons calculer l'incidence cumulée : 
````
peru$cumulative = cumsum(peru$incidence)
peru$cumulative
peru$n = sum(peru$incidence)
````
Représentation graphique 
````
par(mar = c(5, 5, 2, 5)) # création des axes et du graphique
plot(peru$incidence ~ peru$age, type = "b", ylab = "incidence", xlab = "age")
par(new = T)
plot(peru$cumulative ~ peru$age, type = "l", col = "red", axes = FALSE, xlab = "NA", ylab = "NA")
axis(side = 4)
mtext(side = 4, line = 4, "cumulative")
legend("right", legend = c("incidence", "cumumative"), lty = c(1,1), col = c("black", "red"))
````
Nous allons appliquer le modèle en faisant l'hypothèse que les différentes tranches d'âge présentent des *FoI* différents. La représentation graphique est la suivante : 
````
up = c(1:20, 30, 40, 50, 60, 70, 100) # cut-off par âge
para = rep(.1, length(up))

# éstimation de la log-vraisamblance 
est2 = optim(par = log(para),fn = llik.pc, age = peru$age, num = peru$cumulative, 
        denom = peru$n, up = up, method = "Nelder-Mead", control = list(trace = 2, maxit = 2000))

x = c(0, up)
y = exp(c(est2$par, est2$par[26]))
plot(x, y, ylab = "FoI relative", xlab = "age", type = "l", ylim = c(0,0.25), xlim = c(0,80))
````
Nous pouvons (ou nous devrons) remarquer une *FoI* élevé chez les sujets d'âge compris entre 8 et 10 ans. Ce schéma est cohérent avec la biologie de la rubéole. Le Peru a une espérance de vie de 75 ans et le $R_0$ de la rubéole est estimé dans l'interval 4-10. 
La moyenne d'âge de l'infection peut être calculée par $\bar{a} \approx L/(R_0-1)$. Elle devrait être autour de 10 ans. 

##### Conséquences de la vaccination dans la tranche d'âge 0-45 ans
Nous allons réaliser un *subset* de nos données et les stocker dans le dataframe *data3*. Ensuite nous allons créer une fonction. 
````
data3 = peru[peru$age <45, ]
head(data3)
df = 5
para = rep(0.1, df +1)

tmpfn = function(x,dl){
  x = predict(dl, newdata=data.frame(age = x))
  exp(x)
}
````
Nous avons utilisé une log-transformation afin de pouvoir constraidre la *FoI* afin qu'il soit positif. \
Pour cela nous utilisons un objet *dummy*.

````
# dummy lm
dl=lm(cumulative ~ bs(age,df), data = data3)

# log-likelihood function
tmpfn2 = function(par, data, df){
  dl = lm(cumulative ~ bs(age, df), data = data)
  dl$coefficients=par
  ll = 0
  for(a in 1:length(data$age)){
    p = ((1-exp(-integrate(tmpfn, 0, data$age[a], dl = dl)$value)))
    ll = ll+dbinom(data$cumulative[a], data$n[a], p, log = T)
  }
  return(-ll)
}
````

Nous devrions obtenir une courbe avec deux peak. Le peak dominant a lieu autour des 10 ans, alors qu'un deuxième peak a lieu autour de 35 ans. L'hypothèse que nous pouvons faire est que la plupart des sujets est infecté pendant l'école, mais que, ceux qui échappent à cette modalité dominante d'infection, sont ensuite infectés lors qu'ils sont des enfants.

````
spline_model = optim(par = log(para), fn = tmpfn2,
                     data = data3, df = df, method = "Nelder-Mead",
                     control = list(trace = 4, maxit = 5000))

dl$coefficients = spline_model$par
plot(exp(predict(dl)) ~ data3$age, xlab = "age", ylab = "FoI relative", type = "l")
````

La proportion de cases qui sont censés avoir lieu dans notre intervalle d'âge est représenté par la probabilité jointe. \
Cette valeur est donnée par le modèle spline ; 

````
(exp(-integrate(tmpfn, 0, 15, dl = dl)$value))*(1-
  exp(-integrate(tmpfn, 15, 40, dl = dl)$value))
`````

Nous pouvons donc conclure que 9% des cas ont lieu dans le groupe à risque. 

***

## Modèlisation stochatsique 

````
rm(list=ls())

Gillespie<-function(B,G,Tmax,I0,S0,N)
{
	Time=NULL
	S=NULL
	I=NULL
	R=NULL

	# Initial conditions

	S[1]=S0
	I[1]=I0
	R[1]=N-S0-I0
	Time[1]=0
	t=0 
	idx=0    # indice, ça permet d'inférer mise-à-jour du modèle
	NbInf=I[idx+1]

while(t<Tmax & NbInf>0)	 
	{
		idx = idx+1
	
		# drawing when next event will happen

		p1 = B*S[idx]*I[idx]/N 
		
		#idx indique la ligne dont on se sert pour calculer les intensité de transition, donc la dernière ligne que nous avons rempli
		
		p2 = G*I[idx]

		waiting.time=rexp(1, rate = p1+p2 ) 
		#fonction qui permet de tirer au sort dans une loi exp. "1" une valeur, "rate" somme de ces transitions
		
		dt = waiting.time

		# drawing which event will happen next
		
		# Quelle transition va avoir lieu? -> Probabilité d'une infection et probabilité d'une guérison
		# Nous allons estimer "r" afin de pouvoir inférer qu'il s'agit bien d'une guérison ou d'une infection
		# Nous allons décrire les probabilités sur une droite 0-1 où P-infection (prob infect) est représentée
		# Cette probabilité : 
		# 0 ------ p infection ------ 1 
		# si r est entre 0 et p-infection ==> infection
		# si r est entre p-infection et 1 ==> guérison
		# donc nous sommes en face d'une grandeur qui suit une loi normale - centrée reduite 

		p.Infection = p1/(p1+p2)
		p.Recovery = p2/(p1+p2)

		r = runif(1,min=0,max=1)

		if (r<p.Infection) 
		{
			S[idx+1]=S[idx]-1; I[idx+1]=I[idx]+1 ; R[idx+1]=R[idx] # infection
		}else{
			S[idx+1]=S[idx]; I[idx+1]=I[idx]-1 ; R[idx+1]=R[idx]+1 # recovery
		}

		t=t+dt

		Time[idx+1]=t
		NbInf=I[idx+1]
	}

	tab=as.data.frame(cbind(Time,S,I,R)) #la sortie est représentée par un tableau avec quatre colonnes
	return(tab)
}
````

Définition des paramètres :

* *while* permet de répeter plusieurs fois l'instruction, sans savoir à l'avance quand on peut arrêter la répetion.
* *Nbinf* à un temps $t$ est $I_t$ 
* $I$ est le vecteur dans lequel nous stockons les valeurs succesives retrouvées 
* En effet, à chaque tour de boucle, on remplit les vecteurs Time, $S$, $I$, $R$ (mise à jour de *idx*)

\

La fonction *simulSIR* peut être utilisée pour simuler une fois la fonction épidémique et la représenter

````
SimulSIR.Gillespie<-function(B,G,Tmax,I0,S0,N,plot.bool)
{
	tab=Gillespie(B,G,Tmax,I0,S0,N)
	Time=tab$Time
	
	if(plot.bool==TRUE)
	{
		x11()
		par(las=1,mar=c(5,10,5,10),cex.main=1.5,cex.axis=1.5)
		plot(Time,tab$S[1:length(Time)],type="l",col="green",main=paste("Simulation of a stochastic SIR Model with Gillespie algorithm\n",expression(beta),"=",B,", ",expression(gamma),"=",G,sep=""),xlab="Time",ylab="",ylim=c(0,N),bty="n")
		lines(Time,tab$I[1:length(Time)],type="l",col="red")
		lines(Time,tab$R[1:length(Time)],type="l",col="black")
		legend("topright",c("S","I","R"),col=c("green","red","black"),lty=1,bty="n")
	}
	
	return(tab)

}
````

Focntion Gillespie \
\
À noter :

* Si $n = 10$, la courbe décrite est répetée 10 fois
* n = nombre de courbes épidemiques

````
Plot.Several.Gillespie<-function(n,B,G,Tmax,I0,S0,N)
{
  par(las=1,mar=c(5,10,5,1),cex.main=1.5,cex.axis=1.5)
  Time<-vector("list",length=n)
  In<-vector("list",length=n)
	for (i in 1:n)
	{
		tab=Gillespie(B,G,Tmax,I0,S0,N)
    Time[[i]]<-tab$Time
    In[[i]]<-tab$I
  }
  
  #for c'est pour une boucle
  #dans l'exemple, seulement les cas contagieux sont représentés
  
 	for (i in 1:n)
	{  
      R0=B/G
      if (i==1){plot(Time[[i]],In[[i]],type="s",col=1+i,main=paste(N="N", R0="R0", I0="I0", sep=""), 
        xlab = "Time", ylab = "I", xlim = c(0,Tmax), ylim = c(0,range(In)[2]+1), bty="n")}
		  if(i>1){lines(Time[[i]],In[[i]],col=1+i, type = "s")}
	}

}
````

***

## Remerciements et références

* Ottar N. Bjørnstad, "Epidemics"
* Emilia Vynnycky and Richard G. White, "An Introduction to Infectious Disease Modelling"
* Dr Angelo D'Ambrosio, European Centre for Disease Prevention and Control - Stockholm
* Dr Judith Legrand, AgroParisTech - Paris
