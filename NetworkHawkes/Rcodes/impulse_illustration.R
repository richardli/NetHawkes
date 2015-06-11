tmax <- 10
toplot <- seq(0.01, 9.99, len = 10000)
g <- function(t, tau, mu){(tmax / (t * (tmax - t))) * sqrt(tau / 2 / pi) * exp(-tau/ 2 *(log(t / (tmax - t)) - mu)^2)}

pdf("../figures/Impulse.pdf", width = 7, height = 7)
tau <- 10
mu <- -1
plot(toplot, g(toplot, tau, mu), type = "l", ylim = c(0, 1.3), lwd = 2, 
	ylab = "impulse function", xlab = "time", 
	xaxt = 'n', col = 8)
axis(1, at = c(0, 5, 10), labels = c("0", 
	expression(T[max]/2), expression(T[max])))	
tau <- 5
mu <- -1
lines(toplot, g(toplot, tau, mu), col = 2, lwd = 2)	
tau <- 1
mu <- -1
lines(toplot, g(toplot, tau, mu), col = 3, lwd = 2)	
tau <- 1
mu <- -0.5
lines(toplot, g(toplot, tau, mu), col = 4, lwd = 2)	
tau <- 1
mu <- -2
lines(toplot, g(toplot, tau, mu), col = 5, lwd = 2)	
tau <- 5
mu <- -2
lines(toplot, g(toplot, tau, mu), col = 6, lwd = 2)	
tau <- 5
mu <- -5
lines(toplot, g(toplot, tau, mu), col = 1, lwd = 2)
legend("topright", col = c(8,2,3,4,5,6,1), lty = rep(1, 7), lwd = rep(2, 7), 
	c("tau = 10, mu = -1",
	  "tau = 5,   mu = -1",
	  "tau = 1,   mu = -1",
	  "tau = 1,   mu = -1/2",
	  "tau = 1,   mu = -2",
	  "tau = 5,   mu = -2",
	  "tau = 5,   mu = -5"))
dev.off()