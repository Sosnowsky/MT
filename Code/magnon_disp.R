library(plotly)
setwd("~/Desktop/MT/Code")
library(latex2exp)
library(tikzDevice)
options(tikzMetricPackages = c("\\usepackage[utf8]{inputenc}",
                               "\\usepackage[T1]{fontenc}", "\\usetikzlibrary{calc}",
                               "\\usepackage{amssymb}"))
gamma <- function(k) {
  delta1=0.5*c(1,sqrt(3))
  delta2=0.5*c(1,-sqrt(3))
  delta3=c(-1,0)
  result <- exp(complex(real = 0, imaginary = 1)*(k%*%delta1))+
    exp(complex(real = 0, imaginary = 1)*(k%*%delta2))+
    exp(complex(real = 0, imaginary = 1)*(k%*%delta3))
  return(result)
}

delta <- function(k) {
  delta1=0.5*c(3,sqrt(3))
  delta2=c(0,-sqrt(3))
  delta3=0.5*c(-3,sqrt(3))
  result <- -2*D/J*
    (sin(k%*%delta1)+sin(k%*%delta2)+sin(k%*%delta3))
  return(result)
}

energy_p <- function(k) {
  return(2*(delta(k)+sqrt(9-abs(gamma(k)^2))))
}

energy_m <- function(k) {
  return(2*(-delta(k)+sqrt(9-abs(gamma(k)^2))))
}

kmax = 3
res = 0.01
J=1
D=0.1
S=2.5
kx_vec<-seq(from=-kmax, to=kmax, by=res)
samples=length(kx_vec)
energyp <- matrix(0, length(kx_vec), length(kx_vec))
energym <- matrix(0, length(kx_vec), length(kx_vec))

for(kx in 1:samples){
  for(ky in 1:samples){
    energyp[kx,ky]<-energy_p(c(kx_vec[kx],kx_vec[ky]))[1][1]
    energym[kx,ky]<-energy_m(c(kx_vec[kx],kx_vec[ky]))[1][1]
  }
}
ky=0

plot_for_ky <- function(ky) {
  energyp<-sapply(kx_vec, function(x) energy_p(c(x,ky))[1][1])
  energym<-sapply(kx_vec, function(x) energy_m(c(x,ky))[1][1])
  p<-plot_ly(x= ~kx_vec, y=~energyp,type='scatter',mode='lines',name=paste0("e+, ky=",ky)) %>%
    add_trace(x=~kx_vec, y=~energym, name=paste0("e-, ky=",ky),type='scatter',mode='lines')
  return(p)
}
p<-plot_for_ky(0)
plot_for_ky(-1)
subplot(plot_for_ky(-1),plot_for_ky(0),plot_for_ky(1), nrows = 3)
?subplot

p <- plot_ly(z = ~energyp, x=kx_vec,y=kx_vec) %>% 
  layout(
    title = "",
    scene = list(
      xaxis = list(title = ""),
      yaxis = list(title = ""),
      zaxis = list(title = "")
    )) %>% add_surface() %>% add_surface(z = ~energym)
p
p <- plot_ly(z = energyp, x=kx_vec,y=kx_vec, type = "heatmap") %>% layout(scene=list(      xaxis = list(title = "x"),
                                                                                           yaxis = list(title = "y"),
                                                                                           zaxis = list(title = "")))
m <- plot_ly(z = energym, x=kx_vec,y=kx_vec, type = "heatmap")
p <- subplot(p, m)
p

chart_link = api_create(p, filename="surface-1")
chart_link
