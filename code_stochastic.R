rm(list=ls())

Gillespie<-function(B,G,Tmax,I0,S0,N)
{
  Time=NULL
  S=NULL
  I=NULL
  R=NULL
  
# SITUATION AU DEBUT DE L'EPIDEMIE
  S[1]=S0
  I[1]=I0
  R[1]=N-S0-I0
  Time[1]=0
# TEMPS COURANT - IL SERA MIS A JOUR DANS LE PROCESSUS
  t=0 
# INDIATION DE LA MISE A JOUR DU MODELE (A QUEL POINT IL FAUT L'UPDATER ?)
  idx=0 
  NbInf=I[idx+1]
  
# COMMANDE 'while'
# CETTE COMMANDE PERMET DE REPETER PLUSIEURS FOIS L'INSTRUCTION, SANS QUE LE MOMENT D'ARRET DE LA 
# REPETITION SOIT CONNU
# AFIN D'ARRETER LA REPETITION
# CETTE INSTRUCTION EST INDISPENSABLE POUR ARRETER LE SYSTEME

# 'NbInf' AU TEMPS t EST = I
# 'I' EST UN VECTEUR DE STOCKAGE DES VALEURES RETROUVEES A CHAQUE BOUCLE 
# LE PROCESSUS CONSISTE A REPLIR LE VECTEURS 'Time', 'S', 'I' ET 'R'
# REMIS A JOUR A CHAQUE BOUCLE 'idx'

  while(t<Tmax & NbInf>0)	 
  {
    idx=idx+1
  
# ESTIMER QUAND LE PROCHAIN EVENEMENT AURA LIEU
# 'idx' INDIQUE LA LIGNE DONT ON SE SER POUR LE CALCUL DE L'INTENSITE DE TRANSITION
# LA DERNIERE LIGNE QU'ON A REMPLIT
    p1=B*S[idx]*I[idx]/N  
    p2=G*I[idx]
# LA FONCTION QUI SUIT PERMET DE TIRER AU SORT DANS UNE LOI EXP
# '1' = VALEUR
# 'RATE' = SOMME DES TRANSITIONS
    waiting.time=rexp(1, rate = p1+p2 )
    dt=waiting.time
    
# PREVOIR QUEL EVENEMENT AURA LIEU
# QUELLE TRANSITION VA AVOIR LIEU ? QUELLE EST LA PROBABILITE D'UNE INFECTION ET CELLE D'UNE GUERISON ?
# ET COMMENT SAVOIR DE QUELLE TRANSITION IL S'AGIT ? 
    #quelle transition va avoir lieu? probabilité d'une infection et probabilité d'une guérison
    #comment j'estime qu'il s'agit d'une infection ou d'un guérison ? 
    
# J'ESTIME 'r'
# LES PROBABILITES D'INFECTION SONT REPRESENTEES PAR UNE DROITE ENTRE 0 ET 1
# 0 ------- P INFECTION ------- 1
# SI r EST COMPRIS ENTRE 0 ET P INFECTION ==> INFECTION
# SI r EST COMPRIS ENTRE P INFECTION ET 1 ==> GUERISON
# LOI NORMALE N(0,1)
    
    p.Infection=p1/(p1+p2)
    p.Recovery=p2/(p1+p2)
    
    r=runif(1,min=0,max=1)
    
    if (r<p.Infection) 
    {
      S[idx+1]=S[idx]-1; I[idx+1]=I[idx]+1 ; R[idx+1]=R[idx] # INFECTION
    }else{
      S[idx+1]=S[idx]; I[idx+1]=I[idx]-1 ; R[idx+1]=R[idx]+1 # GUERISON
    }
    
    t=t+dt
    
    Time[idx+1]=t
    NbInf=I[idx+1]
  }
  
  tab=as.data.frame(cbind(Time,S,I,R)) 
  return(tab)
}

#########################################################################

# simulSIR PERMET DE SIMULER ET DE REPRESENTER LA FONCTION EPIDEMIQUE

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

# LA FONCTION GILLESPIE EST FAITE TOURNER PLUSIEURS FOIS
# SI n=10, ON REPETE 10 FOIS LA COURBE DECRITE
# 'n' EST DONC LE NOMBRE DE COURBES EPIDEMIQUES
# LES AUTRES PARAMETRES : B=beta, G=gamma, I0=INFECTION, S0=SUSCEPTIBLES, N=NOMBRE DE PERSONNES DANS LA POP

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

# FOR ==> CREATION DE LA BOUBLE 'i in 1:n'
# REPRESENTATION DES CAS CONTAGIEUX SEULEMENT 
  
  for (i in 1:n)
  {  
    R0=B/G
    if (i==1){plot(Time[[i]],In[[i]],type="s",col=1+i,main=paste("N=",N,", R0=",R0,", I0=",I0,sep=""),xlab="Time",ylab="I",xlim=c(0,Tmax),ylim=c(0,range(In)[2]+1),bty="n")}
    if(i>1){lines(Time[[i]],In[[i]],col=1+i,type="s")}
  }
  
  
}