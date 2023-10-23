#### 1
n2 = 27
t_barra = 3
s1 = t_barra*n2


beta=1/10^2
alpha=beta

beta1 = beta
alfa1 = alpha

# 1.1
(beta1+s1)/(alfa1+n2-1)

# 1.2
var1 <- ((beta1+s1)^2)/(((alfa1+n2-1)^2)*(alfa1+n2-2))
sqrt(var1)

# 1.3
# paciente sobreviver pelo menos y anos
y <- 5
if (!require('Pareto')) install.packages('Pareto'); library('Pareto')

alfa2 <- (alfa1+n2) # shape
k <- (beta1+s1) # location

1-pPareto(k+y, k, alfa2)
##### 2
n1 = 10
sigma2 = 2^3
mu0 = 7
tau2 = 1
media = 10

valor1 = 7.99
valor2 = 9.93

# 2.1 valor esperado de mu a posteriori
(mu<-round((1/tau2*mu0+n1/sigma2*media)/(1/tau2+n1/sigma2),2))

# 2.2
var<-(1)/(1/tau2+n1/sigma2) 
(dp<-round(sqrt(var),2))

# 2.3
(distribuicao<-c(paste0("Normal~(",mu,",",dp,")")))
(round(pnorm(valor2,mu,dp)-pnorm(valor1,mu,dp),2))

##### 3 

# 2.1
alfa <- 4
s <- 10
betaa <- 4
n <- 10

(esperado <- (alfa+s)/(betaa+n))

# 2.2
(desvio <- sqrt((alfa+s)/((betaa+n)^2)))

# 2.3 probabilidade de theta ser maior que x
x <- 1
1-pgamma(x, s+alfa, n+alfa)