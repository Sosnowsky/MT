#units: hbar=1 (Energy*Time)  t1=1(Energy),
#electric field is measured in units of t1/ea

library(latex2exp)
t1 <- 1 #NN hopping integral
t2 <- 0.1 #NNN hoping integral
delta <- 0.5 #NNN intrinsic SOI coupling

U <- 10 #Interaction energy
w <- 6 #Laser frecuency
bessel_trunc_order <- 10
lattice_cte <- 1
field_energies <- seq(from=0, by=0.01, to=20)

#J1 = 2*t1^2/U*J1_factor
J1_factor <- vector(mode="numeric", length=length(field_energies))

for (n in -bessel_trunc_order:bessel_trunc_order) {
  J1_factor <- J1_factor + besselJ(field_energies/sqrt(2), nu = n)^2 / (1 + n * w / U)
}

#D2 = 4*t2*delta/U*D2_factor
D2_factor <- vector(mode="numeric", length=length(field_energies))

for (n in -bessel_trunc_order:bessel_trunc_order) {
  D2_factor <- D2_factor + besselJ(sqrt(3)*field_energies/sqrt(2), nu = n)^2 / (1 + n * w / U)
}

ratioJD0 <- t1^2 / ( 2 * t2 * delta)
ratioJD <-  ratioJD0*J1_factor/D2_factor

jpeg('rplot.jpg')
plot(x,y)


plot(x=field_energies, y=J1_factor, type='l', col = "green", xlab = TeX('$E$'), ylab = '')
lines(x=field_energies, y = D2_factor, type = 'l', col = "red")
legend("topright", legend = c(TeX('$J_1/J_1^0$'), TeX('$D_2/D_2^0$')), fill = c("green", "red"))

dev.off()
jpeg('rplot.jpg')
plot(x=field_energies, y = ratioJD, ylim = c(-2*ratioJD0,2*ratioJD0), type = 'l', col = "blue", xlab = TeX('$E$'), ylab = TeX('$J/D$'))
dev.off()
