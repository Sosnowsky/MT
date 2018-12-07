#E = sum J_1S_i S_j + sum D_2 S_i times S_j

library(latex2exp)

J1 <- 1
D2 <- 0.01

theta <- seq(from = -pi, by = 0.01, to = pi)
E <- -J1*cos(theta) + D2*sin(2*theta)

solp <- asin((J1/D2 + sqrt((J1/D2)^2+32))/8)
Ep <- -J1*cos(solp) + D2*sin(2*solp)
soln <- asin((J1/D2 - sqrt((J1/D2)^2+32))/8)
En <- -J1*cos(soln) + D2*sin(2*soln)

plot(x = theta, y = E, type = 'l')
lines(x = solp, y = Ep, lwd = 5, type = 'p', col = 'green')
lines(x = soln, y = En, lwd = 5, type = 'p', col = 'green')
