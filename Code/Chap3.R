#units: hbar=1 (Energy*Time)  t1=1(Energy),
#electric field is measured in units of t1/ea

setwd("~/Desktop/MT/Code")
library(latex2exp)
library(tikzDevice)
options(tikzMetricPackages = c("\\usepackage[utf8]{inputenc}",
                               "\\usepackage[T1]{fontenc}", "\\usetikzlibrary{calc}",
                               "\\usepackage{amssymb}"))

plot_disconnect <- function(x, y, ...) {
  Y_DISC <- 10
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
frequencies <- list(6,16) #Laser frecuency
bessel_trunc_order <- 10
lattice_cte <- 1
field_energies <- seq(from=0, by=0.001, to=5)

#J1 = 2*t1^2/U*J1_factor
J1_factor <- lapply(frequencies, function(x) vector(mode="numeric", length=length(field_energies)))

#D2 = 4*t2*delta/U*D2_factor
D2_factor <- lapply(frequencies, function(x) vector(mode="numeric", length=length(field_energies)))

ratioJD0 <- t1^2 / ( 2 * t2 * delta)
ratioJD <-  lapply(frequencies, function(x) vector(mode="numeric", length=length(field_energies)))

for (i in 1:length(frequencies)) {
  for (n in -bessel_trunc_order:bessel_trunc_order) {
    J1_factor[[i]] <- J1_factor[[i]] + besselJ(field_energies/sqrt(2), nu = n)^2 / (1 + n * frequencies[[i]] / U)
    D2_factor[[i]] <- D2_factor[[i]] + besselJ(sqrt(2)*field_energies, nu = n)^2 / (1 + n * frequencies[[i]] / U)
  }
  ratioJD[[i]] <- ratioJD0*J1_factor[[i]]/D2_factor[[i]]
}

current_file <- "NNvsNNN"
tikz(paste0(current_file, '.tex'), width = 4, height = 4, standAlone = TRUE,
     packages = c("\\usepackage{tikz}",
                  "\\usepackage[active,tightpage,psfixbb]{preview}",
                  "\\PreviewEnvironment{pgfpicture}",
                  "\\setlength\\PreviewBorder{0pt}",
                  "\\usepackage{amssymb}"))

plot(x=range(field_energies), y = range(J1_factor), type='n', xlab = "$\\mathcal{E}$", ylab = '')
for (i in 1:length(frequencies)) {
  lines(x = field_energies, y = J1_factor[[i]], type = 'l', col = "green", lwd = 2, lty = i)
  lines(x = field_energies, y = D2_factor[[i]], type = 'l', col = "red", lwd = 2, lty = i)
}
legend("topright", legend = c("$\\frac{J_1}{J_1^0}$", "$\\frac{D_2}{D_2^0}$"), fill = c("green", "red"))

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

ylim <- c(-3*ratioJD0,2*ratioJD0)
plot(x = range(field_energies), y = ylim, type = 'n', 
     xlab = "$\\mathcal{E}$", ylab = '$\\frac{J_1}{D_2}$')

for (i in 1:length(frequencies)) {
  plot_disconnect(x=field_energies, y = ratioJD[[i]], col = "blue", lwd = 2, lty = i, ylim = ylim)
}

dev.off()
tools::texi2pdf(paste0(current_file, '.tex'))
system(paste0("mv ", current_file, ".pdf ../Chapters/", current_file, ".pdf"))

EMax <- lapply(1:length(frequencies), function(i) which(D2_factor[[i]] <= 0)[1]-1 )
res_ratio <- lapply(1:length(frequencies), function(i)  ratioJD[[i]][1:EMax[[i]]] )
res_field_energies <- lapply(1:length(frequencies), function(i)  field_energies[1:EMax[[i]]])
theta <- lapply(1:length(frequencies), function(i) asin((res_ratio[[i]]-sqrt(res_ratio[[i]]^2+32))/8) )  

current_file <- "theta"
tikz(paste0(current_file, '.tex'), width = 4, height = 4, standAlone = TRUE,
     packages = c("\\usepackage{tikz}",
                  "\\usepackage[active,tightpage,psfixbb]{preview}",
                  "\\PreviewEnvironment{pgfpicture}",
                  "\\setlength\\PreviewBorder{0pt}",
                  "\\usepackage{amssymb}"))

plot(x = range(sapply(res_field_energies, range)), y = range(sapply(theta, range)), type = 'n',
     xlab = "$\\mathcal{E}$", ylab="$\\theta$" )
for (i in 1:length(frequencies)) {
  lines(x = res_field_energies[[i]], y = theta[[i]], type = 'l', lwd =2, lty = i,
         col = "purple")
}

dev.off()
tools::texi2pdf(paste0(current_file, '.tex'))
system(paste0("mv ", current_file, ".pdf ../Chapters/", current_file, ".pdf"))

