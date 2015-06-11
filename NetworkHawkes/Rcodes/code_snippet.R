M <- 51
T <- 23214
dt <- 92856 / M
l <- 0.15
cov  <- matrix(0, M, M)
for(i in 1:M){
	for(j in 1:M){
		cov[i, j] <- exp(-2 * (sin((i-j)*dt/T*pi)^2) / l^2)
	}
}
solve(cov)

##
## test LGCP mean and variance
##
m0 = 0
m1 = 26
m = 13
d = max(m1 - m, m - m0)  

A <- exp(1) - 1
B <- exp(2) - 1
delta <- max(0, (A+B)*d^2 - 4*A*B*m^2)
s = (2 * A * m - sqrt( delta))/(2 * A + 2)
t = m - s

A*s^2+B*t^2 - d^2/4

mu0 = log(s) - .5
mu1 = log(max(t,1e-4)) - 1
mu0
mu1

sim = 10000
y = exp(rnorm(sim, mean = mu0, sd = 1)) + exp(rnorm(sim, mean = mu1, sd = 1)) * exp(rnorm(sim)) 
summary(y)
length(which(y < m0)) /sim
length(which(y < m1)) /sim
mean(y) - 2 * sd(y)
mean(y) + 2 * sd(y)

mean(y)
s + t

var(y)
A * s^2 + B * t^2

##
## test LGCP mean and variance
##
a <- 3.46
b <- 1.22
z = exp(rnorm(sim, mean = a, sd = 0.1)) + exp(rnorm(sim, mean = 0, sd = 1)) * exp(rnorm(sim, mean = b, sd = 0.1))
mean(z)
summary(z)

exp(b+1) + exp(a + .5)

A <- exp(1) - 1
B <- exp(2) - 1
var(z)
A * exp(2*a + 1) + B * exp(2*b+2) 


