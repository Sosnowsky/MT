#units: hbar=1 (Energy*Time)  t1=1(Energy),
#electric field is measured in units of t1/ea

setwd("~/Desktop/MT/Code")
library(latex2exp)
library(tikzDevice)
options(tikzMetricPackages = c("\\usepackage[utf8]{inputenc}",
                               "\\usepackage[T1]{fontenc}", "\\usetikzlibrary{calc}",
                               "\\usepackage{amssymb}"))

plot_disconnect <- function(x, y, ylim, ...) {
  Y_DISC <- 10
  plot(x = range(x), y = ylim, type = 'n', ...)
  discontinuities <- which(abs(diff(y)) > Y_DISC)
  last_point <- 1
  for (current_point in discontinuities) {
    if (current_point > last_point + 1) {
      cat('plotting from ', x[last_point], ' to ', x[current_point], ' with avg: ', range(y[last_point:current_point]) ,'\n')
      lines(x = x[last_point:current_point], y = y[last_point:current_point], type = 'l', ...)
    }
    last_point <- current_point + 1
  }
  if (last_point < length(x)) {
    cat('plotting from ', x[last_point], ' to ', x[length(x)],'\n')
    lines(x = x[last_point:length(x)], y = y[last_point:length(x)], type = 'l', ...)
  }
}

t1 <- 1 #NN hopping integral
t2 <- 0.1 #NNN hoping integral
delta <- 0.5 #NNN intrinsic SOI coupling

U <- 10 #Interaction energy
w <- 6 #Laser frecuency
bessel_trunc_order <- 10
lattice_cte <- 1
field_energies <- seq(from=0, by=0.001, to=5)

#J1 = 2*t1^2/U*J1_factor
J1_factor <- vector(mode="numeric", length=length(field_energies))

for (n in -bessel_trunc_order:bessel_trunc_order) {
  J1_factor <- J1_factor + besselJ(field_energies/sqrt(2), nu = n)^2 / (1 + n * w / U)
}

#D2 = 4*t2*delta/U*D2_factor
D2_factor <- vector(mode="numeric", length=length(field_energies))

for (n in -bessel_trunc_order:bessel_trunc_order) {
  D2_factor <- D2_factor + besselJ(sqrt(2)*field_energies, nu = n)^2 / (1 + n * w / U)
}

ratioJD0 <- t1^2 / ( 2 * t2 * delta)
ratioJD <-  ratioJD0*J1_factor/D2_factor

current_file <- "NNvsNNN"
tikz(paste0(current_file, '.tex'), width = 4, height = 4, standAlone = TRUE,
     packages = c("\\usepackage{tikz}",
                  "\\usepackage[active,tightpage,psfixbb]{preview}",
                  "\\PreviewEnvironment{pgfpicture}",
                  "\\setlength\\PreviewBorder{0pt}",
                  "\\usepackage{amssymb}"))

plot(x=field_energies, y=J1_factor, type='l', col = "green", xlab = "$\\mathcal{E}$", ylab = '', lwd = 2)
lines(x=field_energies, y = D2_factor, type = 'l', col = "red", lwd = 2)
legend("topright", legend = c(TeX('$J_1/J_1^0$'), TeX('$D_2/D_2^0$')), fill = c("green", "red"))

dev.off()
tools::texi2pdf(paste0(current_file, '.tex'))
system(paste0("mv ", current_file, ".pdf ../Chapters/", current_file, ".pdf"))
#system(paste(getOption("pdfviewer"), "formula.pdf"))

current_file <- "ratio"
tikz(paste0(current_file, '.tex'), width = 4, height = 4, standAlone = TRUE,
     packages = c("\\usepackage{tikz}",
                  "\\usepackage[active,tightpage,psfixbb]{preview}",
                  "\\PreviewEnvironment{pgfpicture}",
                  "\\setlength\\PreviewBorder{0pt}",
                  "\\usepackage{amssymb}"))

plot_disconnect(x=field_energies, y = ratioJD, ylim = c(-3*ratioJD0,2*ratioJD0), 
                col = "blue", lwd = 2, xlab = "$\\mathcal{E}$", ylab = '$\\frac{J_1}{D_2}$')

dev.off()
tools::texi2pdf(paste0(current_file, '.tex'))
system(paste0("mv ", current_file, ".pdf ../Chapters/", current_file, ".pdf"))

EMax <- which(D2_factor <= 0)[1]-1
res_ratio <- ratioJD[1:EMax]
res_field_energies <- field_energies[1:EMax]
theta <- asin((res_ratio-sqrt(res_ratio^2+32))/8)

current_file <- "theta"
tikz(paste0(current_file, '.tex'), width = 4, height = 4, standAlone = TRUE,
     packages = c("\\usepackage{tikz}",
                  "\\usepackage[active,tightpage,psfixbb]{preview}",
                  "\\PreviewEnvironment{pgfpicture}",
                  "\\setlength\\PreviewBorder{0pt}",
                  "\\usepackage{amssymb}"))

plot(x = res_field_energies, y = theta, type = 'l', lwd =2,
     col = "purple", xlab = "$\\mathcal{E}$", ylab="$\\theta$")

dev.off()
tools::texi2pdf(paste0(current_file, '.tex'))
system(paste0("mv ", current_file, ".pdf ../Chapters/", current_file, ".pdf"))

