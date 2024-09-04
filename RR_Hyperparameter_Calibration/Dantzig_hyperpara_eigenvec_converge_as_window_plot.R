rm(list = ls())

setwd("~/Documents/GitHub/Network-Portfolio/RR_Hyperparameter_Calibration")

window_size = 125

window_size = 250

window_size = 375

window_size = 500

load(paste0("Dantzig_lambda_rolling_window_windowsize",window_size,".RData"))

a = unlist(lmd.EC.Dantzig.list)

b = cumsum(a) / (1:length(a))

c = sapply(1:length(a), function(i) median(a[1:i]))

x = 1:length(a)

Colors_Ret<-rainbow(12)

pngname<-paste0(getwd(),"/Dantzig_hyper_eigenvector_converge_rollingwindow_250_375_500.png")
png(file = pngname, width=500, height=400, bg = "transparent")
plot(x, a, col=Colors_Ret[c(1)], ylab="", xlab="")
# plot(x, b, type = 'l', col=Colors_Ret[c(1)], ylab="", xlab="")
lines(x, c, col = Colors_Ret[c(2)])
legend("topright",legend=c("rho^E (WS250)", "median(rho^E) (WS250)"), col=Colors_Ret[c(1,2)], cex=0.8)
dev.off()


#####  plot 4 window sizes in the same figure ######

window_size = 125
load(paste0("Dantzig_lambda_rolling_window_windowsize",window_size,".RData"))
a = unlist(lmd.EC.Dantzig.list)
x_125 = 1:length(a)
median_rho_125 = sapply(1:length(a), function(i) median(a[1:i]))

window_size = 250
load(paste0("Dantzig_lambda_rolling_window_windowsize",window_size,".RData"))
a = unlist(lmd.EC.Dantzig.list)
x_250 = 1:length(a)
median_rho_250 = sapply(1:length(a), function(i) median(a[1:i],n))

window_size = 375
load(paste0("Dantzig_lambda_rolling_window_windowsize",window_size,".RData"))
a = unlist(lmd.EC.Dantzig.list)
x_375 = 1:length(a)
median_rho_375 = sapply(1:length(a), function(i) median(a[1:i]))

window_size = 500
load(paste0("Dantzig_lambda_rolling_window_windowsize",window_size,".RData"))
a = unlist(lmd.EC.Dantzig.list)
x_500 = 1:length(a)
median_rho_500 = sapply(1:length(a), function(i) median(a[1:i]))

Colors_Ret<-rainbow(12)

pngname<-paste0(getwd(),"/Dantzig_hyper_eigencentr_median_rollingwindow_125_250_375_500.png")
png(file = pngname, width=500, height=400, bg = "transparent")
plot(x_500, median_rho_500, type='l', col=Colors_Ret[c(1)], 
     ylab = "Median values", xlab="")
# plot(x_125, median_rho_125, type='l', col=Colors_Ret[c(1)], 
#      ylab = "Median values", xlab="")
lines(x_250, median_rho_250, col = Colors_Ret[c(2)])
lines(x_375, median_rho_375, col = Colors_Ret[c(3)])
lines(x_500, median_rho_500, col = Colors_Ret[c(4)])
legend("topright",legend=c("WS125", "WS250", "WS375", "WS500"), col=Colors_Ret[c(1,2,3,4)], cex=0.8, lty=1)
dev.off()

#####  plot 2 window sizes in the same figure ######

window_size = 375
load(paste0("Dantzig_lambda_rolling_window_windowsize",window_size,".RData"))
a = unlist(lmd.EC.Dantzig.list)
x_375 = 1:length(a)
median_rho_375 = sapply(1:length(a), function(i) median(a[1:i]))
median_rho_375

window_size = 500
load(paste0("Dantzig_lambda_rolling_window_windowsize",window_size,".RData"))
a = unlist(lmd.EC.Dantzig.list)
x_500 = 1:length(a)
median_rho_500 = sapply(1:length(a), function(i) median(a[1:i]))

Colors_Ret<-rainbow(12)

pngname<-paste0(getwd(),"/Dantzig_hyper_eigencentr_median_rollingwindow_375_500.png")
png(file = pngname, width=500, height=400, bg = "transparent")
plot(x_375, median_rho_375, type='l', col=Colors_Ret[c(1)], 
     ylab = "Median values of rho^E", xlab="")
lines(x_500, median_rho_500, col = Colors_Ret[c(8)])
legend("topright",legend=c("WS375", "WS500"), col=Colors_Ret[c(1,8)], cex=0.8, lty=1)
dev.off()

