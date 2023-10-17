
############################
#  Figura para densidades Beta
###########################

par(mfrow=c(1,3))

t<-seq(from=0,to=1,length=1000)

densidade<-NULL

alfa<-c(0.2,0.7,3,10)
for(a in alfa){densidade<-cbind(densidade,dbeta(t,a,3*a))}

matplot(t,densidade,type='l',ylim=c(0,6),xlab='',lwd=2)
mtext(expression(theta),side=1,at=0.5,line=2)
title(main='(a)')
legend(locator(1),legend=c('Beta(0.2,0.6)','Beta(0.7,2.1)','Beta(3,9)','Beta(10,30)'),
       lty=1:4,col=1:length(alfa),bty='n',cex=1,lwd=2)

densidade<-NULL

alfa<-c(0.2,1,3,10)
for(a in alfa){densidade<-cbind(densidade,dbeta(t,a,a))}

matplot(t,densidade,type='l',ylim=c(0,6),xlab='',lwd=2)
mtext(expression(theta),side=1,at=0.5,line=2)
title(main='(b)')
legend(locator(1),legend=c('Beta(0.2,0.2)','Beta(1,1.1) (Unif)',
                           'Beta(3,3)','Beta(10,10)'),
       lty=1:4,col=1:length(alfa),bty='n',cex=1,lwd=2)

densidade<-NULL

beta<-c(0.2,0.7,3,10)
for(b in beta){densidade<-cbind(densidade,dbeta(t,3*b,b))}

matplot(t,densidade,type='l',ylim=c(0,6),xlab='',lwd=2)
mtext(expression(theta),side=1,at=0.5,line=2)
title(main='(c)')
legend(locator(1),legend=c('Beta(0.6,0.2)','Beta(2.1,0.7)','Beta(9,3)','Beta(30,10)'),
       lty=1:4,col=1:length(alfa),bty='n',cex=1,lwd=2)

# Continuacao

p<-0.90
t.I<-0.40
t.S<-0.60
alfa<-30

pbeta(t.S,alfa,alfa)-pbeta(t.I,alfa,alfa)


prob<-function(a){pbeta(t.S,a,a)-pbeta(t.I,a,a)-p}
alfa<-uniroot(f=prob,lower=10,upper=100)$root
pbeta(t.S,alfa,alfa)-pbeta(t.I,alfa,alfa)
alfa

alfa<-33.39
beta<-33.39
s<-93
n<-200

s/n    # estimador de MV de theta
(alfa+s)/(alfa+beta+n)    # média a posteriori
pbeta(0.6,alfa+s,beta+n-s)-pbeta(0.4,alfa+s,beta+n-s)  # nova probabilidade

qbeta(0.05,alfa+s,beta+n-s)  # limites de intervalo crível 90%
qbeta(0.95,alfa+s,beta+n-s)

t<-seq(from=0.3,to=.7,length=1000)
priori<-dbeta(t,alfa,beta)
vero<-dbeta(t,s+1,n-s+1)
post<-dbeta(t,alfa+s,beta+n-s)

par(mfrow=c(1,1))

matplot(t,cbind(priori,vero,post),type='l',ylab='densidade',
        xlab=expression(theta),lwd=3)
legend(locator(1),legend=c('priori','Verossimilhanca',
                           'posteriori'),
       lty=1:3,col=1:3,bty='n',cex=1,lwd=3)

# dist. preditiva (beta-binomial)

a<-1
b<-1
n<-30
m<-15
s.1n<-6

pbb.log<-rep(NA,m+1)
j<-0:15
pbb.log<-lgamma(m+1)-lgamma(j+1)-lgamma(m-j+1)+
  lbeta(a+s.1n+j,b+n+m-s.1n-j)-lbeta(a+s.1n,b+n-s.1n)

pbb<-exp(pbb.log)
sum(pbb[9:16])

barplot(height=pbb,xlab='',ylab='',)
mtext('j',side=1,line=2,cex=2)
mtext('prob.',side=2,line=2,cex=2)

#######################################
exemplo Normal-Normal
#######################################

n<-10
x.barra<-3.53
sigma<-0.2

y.barra<-3.25
s.y<-0.5

y.barra-1.96*s.y/sqrt(50);y.barra+1.96*s.y/sqrt(50)

mu.0<-y.barra
tau<-s.y

(mu.1<-(mu.0/tau^2+x.barra*n/sigma^2)/(1/tau^2+n/sigma^2))

(tau.1<-1/sqrt(1/tau^2+n/sigma^2))

mu.1-qnorm(0.975)*tau.1;mu.1+qnorm(0.975)*tau.1

mu<-seq(from=2.5,to=4.5,length=10000)
dens.priori<-dnorm(mu,mu.0,tau)
dens.post<-dnorm(mu,mu.1,tau.1)
vero<-dnorm(mu,x.barra,sigma)

matplot(mu,cbind(dens.priori,vero,dens.post),type='l',
        xlab=expression(mu),ylab='densidade',lwd=3,cex=2)
legend(locator(1),legend=c('priori','Verossimilhanca',
                           'posteriori'),
       lty=1:3,col=1:3,bty='n',cex=1,lwd=3)

######################

sigma<-10

t<-seq(from=0.01,to=0.99,length=1000)
eta<-log(t/(1-t))
priori<-dnorm(eta,0,sigma)*exp(eta)/((1+exp(eta))^2)
post<-priori*t^9*(1-t)^3
post<-post/sum(post)

posts<-cbind(1000*post,dbeta(t,10,4))
matplot(t,posts,type='l')
sum(post)
post

#############################

# Regiões críveis

nivel<-0.90
alfa<-13.1
beta<-10.1

qgamma(0.01,alfa,beta);qgamma(0.91,alfa,beta)
qgamma(0.09,alfa,beta);qgamma(0.99,alfa,beta)
qgamma(0.45,alfa,beta);qgamma(0.55,alfa,beta)

qgamma(0.01,alfa,beta)-qgamma(0.91,alfa,beta)
qgamma(0.09,alfa,beta)-qgamma(0.99,alfa,beta)

# intervalo HPD

# Primeira forma (minimizando o comprimento do intervalo)

comp<-function(a.1){
  a.2<-1-nivel-a.1
  comprimento<-qgamma(1-a.2,alfa,beta)-qgamma(a.1,alfa,beta)
  return(comprimento)
}

(temp<-optimize(f=comp,lower=0,upper=1-nivel))

# intervalo
(inferior<-qgamma(temp$minimum,alfa,beta))
(superior<-qgamma(nivel+temp$minimum,alfa,beta))

# Checando:
pgamma(superior,alfa,beta)-pgamma(inferior,alfa,beta)
dgamma(inferior,alfa,beta)
dgamma(superior,alfa,beta)

#  Outra forma (fazendo a densidade igual nos extremos)

dif<-function(a.1){
  a.2<-1-nivel-a.1
  diferenca<-dgamma(qgamma(1-a.2,alfa,beta),alfa,beta)-
    dgamma(qgamma(a.1,alfa,beta),alfa,beta)
  return(diferenca)
}

(temp<-uniroot(f=dif,lower=0,upper=1-nivel))

# checando
(inferior<-qgamma(temp$root,alfa,beta))
(superior<-qgamma(nivel+temp$root,alfa,beta))

# Checando:
pgamma(superior,alfa,beta)-pgamma(inferior,alfa,beta)
dgamma(inferior,alfa,beta)
dgamma(superior,alfa,beta)

# gráfico para RC:

alfa<-13.1
beta<-10.1
theta<-seq(from=0,to=3,length=1000)
post.dens<-dgamma(theta,alfa,beta)
plot(theta,post.dens,type='l',lwd=2,xlab=expression(theta))

nivel<-0.50
(temp<-optimize(f=comp,lower=0,upper=1-nivel))
(inferior<-qgamma(temp$minimum,alfa,beta))
(superior<-qgamma(nivel+temp$minimum,alfa,beta))

limites<-c(inferior,superior)
extremos<-c(0.4,2.5)

lines(extremos,dgamma(limites,alfa,beta),lwd=0.75,lty=6,col='purple')
lines(limites,dgamma(limites,alfa,beta),lwd=1.5,lty=6,col='purple')
lines(limites,dgamma(limites,alfa,beta),lwd=1.5,lty=6,col='purple')


lines(rep(inferior,2),c(0,dgamma(inferior,alfa,beta)),lwd=1.5,lty=6,col='purple')
lines(rep(superior,2),c(0,dgamma(inferior,alfa,beta)),lwd=1.5,lty=6,col='purple')

nivel<-0.80
(temp<-optimize(f=comp,lower=0,upper=1-nivel))
(inferior<-qgamma(temp$minimum,alfa,beta))
(superior<-qgamma(nivel+temp$minimum,alfa,beta))

limites<-c(inferior,superior)
extremos<-c(0.4,2.5)

lines(extremos,dgamma(limites,alfa,beta),lwd=0.75,lty=5,col='orange')
lines(limites,dgamma(limites,alfa,beta),lwd=1.5,lty=5,col='orange')
lines(rep(inferior,2),c(0,dgamma(inferior,alfa,beta)),lwd=1.5,lty=5,col='orange')
lines(rep(superior,2),c(0,dgamma(inferior,alfa,beta)),lwd=1.5,lty=5,col='orange')

nivel<-0.90
(temp<-optimize(f=comp,lower=0,upper=1-nivel))
(inferior<-qgamma(temp$minimum,alfa,beta))
(superior<-qgamma(nivel+temp$minimum,alfa,beta))

limites<-c(inferior,superior)
extremos<-c(0.4,2.5)

lines(extremos,dgamma(limites,alfa,beta),lwd=0.75,lty=4,col='green')
lines(limites,dgamma(limites,alfa,beta),lwd=1.5,lty=4,col='green')
lines(rep(inferior,2),c(0,dgamma(inferior,alfa,beta)),lwd=1.5,lty=4,col='green')
lines(rep(superior,2),c(0,dgamma(inferior,alfa,beta)),lwd=1.5,lty=4,col='green')

nivel<-0.95
(temp<-optimize(f=comp,lower=0,upper=1-nivel))
(inferior<-qgamma(temp$minimum,alfa,beta))
(superior<-qgamma(nivel+temp$minimum,alfa,beta))

limites<-c(inferior,superior)
extremos<-c(0.4,2.5)

lines(extremos,dgamma(limites,alfa,beta),lwd=0.75,lty=3,col='blue')
lines(limites,dgamma(limites,alfa,beta),lwd=1.5,lty=3,col='blue')
lines(rep(inferior,2),c(0,dgamma(inferior,alfa,beta)),lwd=1.5,lty=3,col='blue')
lines(rep(superior,2),c(0,dgamma(inferior,alfa,beta)),lwd=1.5,lty=3,col='blue')



nivel<-0.99
(temp<-optimize(f=comp,lower=0,upper=1-nivel))
(inferior<-qgamma(temp$minimum,alfa,beta))
(superior<-qgamma(nivel+temp$minimum,alfa,beta))

limites<-c(inferior,superior)
extremos<-c(0.4,2.5)

lines(extremos,dgamma(limites,alfa,beta),lwd=0.75,lty=2,col='red')
lines(limites,dgamma(limites,alfa,beta),lwd=1.5,lty=2,col='red')
lines(rep(inferior,2),c(0,dgamma(inferior,alfa,beta)),lwd=1.5,lty=2,col='red')
lines(rep(superior,2),c(0,dgamma(inferior,alfa,beta)),lwd=1.5,lty=2,col='red')

legend(locator(1),legend=c('50%','80%','90%','95%','99%'),lty=6:2,
       col=c('purple','orange','green','blue','red'),bty='n')

############################

pgamma(1.798,13.1,10.1)-pgamma(0.611,13.1,10.1)
pgamma(2.273,13.1,10.1)-pgamma(0.848,13.1,10.1)
pgamma(1.220,13.1,10.1)+1-pgamma(1.309,13.1,10.1)


###########################
par(mfrow=c(1,1))

m<-1000
t<-seq(from=0,to=1,length=m)
vero<-dcauchy(-2,t,1)*dcauchy(-1,t,1)*dcauchy(0,t,1)*
  dcauchy(1.5,t,1)*dcauchy(2.5,t,1)
post<-dunif(t,0,1)*vero
post<-post/sum(post)

(teta.esp<-sum(t*post))
(tetaquad.esp<-sum(post*t^2))

tetaquad.esp-teta.esp^2

plot(t,post,type='l')

qbeta(0.5,10,4)
10/14
9/12
10*11/(14*15)


M<-50000
alfa<-10
beta<-4

t<-rbeta(M,alfa,beta)
quantile(t^2,probs=0.5)

quantile(rbeta(50000,10,4)^2,probs=0.5)
