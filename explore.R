# autocorrelated t
T = 53
sig = .1
mu = log(3)
x = rep(0, T)
x[1] = mu
rho = .9
for(t in 2:T) {
    x[t] = rt(1, df = 5)*sig + mu + rho*(x[t-1]-mu)
}

plot(1:T, exp(x), col = 'red', ylim = c(0,8))
lines(1:T, exp(x))
    
