rm(list = ls())
library(latex2exp)
library(painter)
library(splines)

path_prefix = "../../result/ConjugatePairs/"
material = c("ngk", "zhou", "ogden")
postfix = ""
mtp = 3
point_sz = 0.5; line_sz = 4.
lineclr = adjustcolor('#000000', alpha.f = 0.9)
with_lines = TRUE

# Function to plot color bar
color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
  scale = (length(lut)-1)/(max-min)
  
  #dev.new(width=1.75, height=5)
  par(fig=c(0.8,0.95,0.05,0.7))
  plot(c(0,5), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
  axis(2, ticks, las=1)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=lut[i], border=NA)
  }
  par(fig=c(0,1,0,1))
}

for (mtp in 1:3){
dat = read.csv(paste(path_prefix, "distribution_", material[[mtp]], postfix, ".csv", sep=""), header = TRUE)
take = function() { dat[["ps"]]>=0.0 }
z <- dat[["ps"]][take()] #dat[["R"]][take()]**2
colors = adjustcolor(ColorBy(z,Palette("Blue", "Red", 100)), alpha.f = 0.5) 

plot.new()
x <- dat[["delta"]][take()]; y <- dat[["pi"]][take()];
plot(x, y, type="p",
     ylab= TeX("$\\pi$"), 
     xlab = TeX("$\\delta$"), 
     main =  TeX(paste(material[[mtp]], ": $\\delta(\\pi)$")),
     pch=16, cex = point_sz,
     #xlim = c(1, 100000),
     #ylim = c(1e-3, 5),
     col=colors, 
     panel.first=grid())
if (with_lines) {
  par(new=TRUE)
  fit <- lm( y~poly(x,3) )
  d = max(x) - min(x)
  xx <- seq(min(x)-d*0.05,max(x)+d*0.05, length.out=250)
  lines(xx, predict(fit, data.frame(x=xx)), col=lineclr, lwd=line_sz)
}
par(new=TRUE)
color.bar(Palette("Blue", "Red", 100), 0, max(z))#, title=TeX("$p^*$"))
dev.print(pdf, paste(path_prefix, material[[mtp]], "_dp.pdf"))
dev.print(svg, paste(path_prefix, material[[mtp]], "_dp.svg"))

plot.new()
x <- dat[["epsilon"]][take()]; y <- dat[["sigma"]][take()];
plot(x, y, type="p",
     ylab= TeX("$\\epsilon$"),
     xlab = TeX("$\\sigma$"), 
     main =  TeX(paste(material[[mtp]], ": $\\epsilon(\\sigma)$")),
     pch=16, cex = point_sz,
     #xlim = c(1, 100000),
     #ylim = c(1e-3, 5),
     col=colors, 
     panel.first=grid())
if (with_lines) {
  par(new=TRUE)
  fit <- lm( y~poly(x,3) )
  d = max(x) - min(x)
  xx <- seq(min(x)-d*0.05,max(x)+d*0.05, length.out=250)
  lines(xx, predict(fit, data.frame(x=xx)), col=lineclr, lwd=line_sz)
}
par(new=TRUE)
color.bar(Palette("Blue", "Red", 100), 0, max(z))
dev.print(pdf, paste(path_prefix, material[[mtp]], "_es.pdf"))
dev.print(svg, paste(path_prefix, material[[mtp]], "_es.svg"))

plot.new()
x <- dat[["gamma"]][take()]; y <- dat[["tau_s"]][take()];
plot(x, y, type="p",
     ylab= TeX("$\\gamma$"), 
     xlab = TeX("$\\tau_*$"), 
     main =  TeX(paste(material[[mtp]], ": $\\gamma(\\tau_*)$")),
     pch=16, cex = point_sz,
     #xlim = c(1, 100000),
     #ylim = c(1e-3, 5),
     col=colors, 
     panel.first=grid())
if (with_lines) {
  par(new=TRUE)
  fit <- lm( y~poly(x,3) )
  d = max(x) - min(x)
  xx <- seq(min(x)-d*0.05,max(x)+d*0.05, length.out=250)
  lines(xx, predict(fit, data.frame(x=xx)), col=lineclr, lwd=line_sz)
}
par(new=TRUE)
color.bar(Palette("Blue", "Red", 100), 0, max(z))
dev.print(pdf, paste(path_prefix, material[[mtp]], "_gt.pdf"))
dev.print(svg, paste(path_prefix, material[[mtp]], "_gt.svg"))
}