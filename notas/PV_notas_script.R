# Exemplo 1: Dados:  S S S F S F S S S S S F

x<-9
n<-12
k<-n-x
theta.0<-0.5

# Inferencia Classica, delineamento (parada) binomial

dbinom(0:n,n,theta.0) # probabilidades sob H_0
sum(dbinom(x:n,n,theta.0)) # p-valor
1-pbinom(x-1,n,theta.0)    # mesmo p-valor

# Infer?ncia Classica, delineamento (parada) binomial negativo

dnbinom(0:(2*x),k,1-theta.0) # probabilidades sob H_0
1-sum(dnbinom(0:(x-1),k,1-theta.0)) # p-valor
1-pnbinom(x-1,k,1-theta.0)   # mesmo p-valor

#############################
#     Fig. 1 (duas verossimilhanças para exemplo do PV)
#############################

par(mfrow=c(1,2))
theta<-seq(from=0.01,to=0.99,length=1000)

plot(theta,220*theta^9*(1-theta)^3,type='l',col='red',
     xlab=expression(theta),ylab=expression(P[theta]),lty=1,lwd=2)
lines(theta,55*theta^9*(1-theta)^3,col='blue',lty=2,lwd=2)
title(main='(a)')

theta<-seq(from=0.1,to=0.9,length=1000)
plot(theta,log(220)+9*log(theta)+3*log(1-theta),type='l',col='red',
     xlab=expression(theta),ylab=expression(paste('log',' ',P[theta])),lty=1,lwd=2)
lines(theta,log(55)+9*log(theta)+3*log(1-theta),col='blue',
      ylim=c(log(55)+9*log(0.1)+3*log(0.9),log(220)+9*log(0.9)+3*log(0.1)),lty=2,lwd=2)
title(main='(b)')

te<-c(0.2,0.4,0.6,0.8)
for(t in te){
  y.l<-log(55)+9*log(t)+3*log(1-t)
  lines(c(t,t),c(y.l,log(220)+y.l-log(55)))
}

#   Exemplo 2 (PC, misturas)

sigma.a<-0.1
sigma.b<-0.9
nivel<-0.95
z.nivel<-qnorm((1+nivel)/2)

x<-3.7

x-z.nivel*sigma.a;x+z.nivel*sigma.a   # IC, Resposta 2.1

# Calculo dos percentis da mistura

F.mistura<-function(z){0.5*pnorm(z,0,sigma.a)+0.5*pnorm(z,0,sigma.b)-(1+nivel)/2}
mist.upper<-uniroot(f=F.mistura,interval=c(0,3))$root

F.mistura<-function(z){0.5*pnorm(z,0,sigma.a)+0.5*pnorm(z,0,sigma.b)-(1-nivel)/2}
mist.lower<-uniroot(f=F.mistura,interval=c(-3,0))$root

mist.lower;mist.upper

#############################
#     Fig. 2 (mistura de va's para o PC)
#############################

par(mfrow=c(1,2))

mu<-3
z<-seq(from=-1.8,to=1.8,length=10000)
z<-mu+z
f.a<-dnorm(z,mu,sigma.a)
f.b<-dnorm(z,mu,sigma.b)
f.mist<-0.5*f.a+0.5*f.b

matplot(z,cbind(f.a,f.b,f.mist),type='l',
        ylab='f',xlab='z',xaxt='n',lwd=3,col=c('black','red','blue'))
mtext(expression(mu),side=1,line=1,at=mu)
title(main='(a)')

z<-mu-z
plot(z,f.mist,type='l',
     ylab='f',xlab=expression(paste('z','-',mu)),xaxt='n',lwd=2)
mtext('0',side=1,at=0,line=1)
title(main='(b)')

mtext(as.character(round(mist.lower,3)),side=1,line=1,at=mist.lower)
mtext(as.character(round(mist.upper,3)),side=1,line=1,at=mist.upper)

z.dots<-seq(from=mist.lower,to=mist.upper,length=100)
for(z in z.dots){
  lines(c(z,z),c(0,0.5*dnorm(z,0,0.1)+0.5*dnorm(z,0,1)),lty=1,lwd=0.3)
}

# Exemplo 1: inferencia sob probabilidade inversa

alfa<-x+1
beta<-n-x+1

1/beta(alfa,beta)  # calculo do reciproco da integral do denominador

1-pbeta(theta.0,alfa,beta)   #  probabilidade da hipótese alternativa
pbeta(theta.0,alfa,beta)     #  probabilidade da hipótese nula

(1-pbeta(theta.0,alfa,beta))/pbeta(theta.0,alfa,beta) # chances relativas

qbeta(0.025,alfa,beta); qbeta(0.975,alfa,beta)  # intervalo de credibilidade para theta, 95%

alfa/(alfa + beta)  # média da distribuição
qbeta(0.5,alfa,beta)   # mediana

#############################
#     Fig. 3 (dist. Beta(10,4)
#############################

par(mfrow=c(1,1))
theta<-seq(from=0,to=1,length=1000)

plot(theta,dbeta(theta,10,4),type='l',xlim=c(0,1.01),lwd=2,
     xlab=expression(theta),ylab='densidade')
lines(c(theta.0,theta.0),c(0,dbeta(theta.0,alfa,beta)),lwd=1,col='red')

for(t in seq(theta.0,1,0.01)) {
  lines(c(t,t),c(0,dbeta(t,alfa,beta)),lty=2,col='red',lwd=0.5)
  mtext(as.character(theta.0),at=theta.0,side=1,line=1,col='red')
}