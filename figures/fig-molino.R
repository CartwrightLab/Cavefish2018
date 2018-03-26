library(RColorBrewer)
library(stringr)
library(data.table)

scriptname = commandArgs() %>% str_subset("^--file=") %>% str_replace("^--file=","")
pdfname = scriptname %>% str_replace("\\.[^.]*$", ".pdf")

bcz = brewer.pal(10,"RdYlBu")
keylevels = seq(0,1,length.out=length(bcz)+1)
keytext = keylevels

a = fread("zstdcat ../molino/results.csv.zst")
b = a[, .(qbar = mean(q)), by=list(m,s,g)]

xaxs = c(0, 2e-6, 4e-6, 6e-6, 8e-6, 10e-6)
yaxs = c(1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1)

pdf(pdfname, width=5.2, height=5.2)
par(mai=c(48,48,12,0)/72,mgp=c(2,0.6,0),las=2)
bb = b[g == 91000+71000, ]
z = dcast(bb, m~s, value.var="qbar")
z = as.matrix(z[-1,])
filled.contour(z,col=bcz,levels=keylevels,
    plot.title = {
        title(xlab=expression(m*" ("%*%10^-6*")"),
            ylab="s",cex.lab=1.2^1)
    },
    key.title = title(xlab=expression(bar(q)),line=1,cex.lab=1.2^2),
    plot.axes = {
        box()
        axis(1, at=seq(0,1,length.out=6), seq(0,10,2))
        axis(2, at=seq(0,1,length.out=6), expression(10^-5, 10^-4, 10^-3, 0.01, 0.1, 1))
        text(0,0.95, "Blindness",col="white",pos=4,font=2)
        text(1,0.05, "Sightedness",col="white",pos=2,font=2)
    }
)
invisible(dev.off())
embedFonts(pdfname, options="-DPDFSETTINGS=/prepress")

gg = sort(unique(b$g))
unlink("molino*.png")
png("molino%03d.png", height=600,width=600)
par(mai=c(48,48,48,48)/72,mgp=c(2,0.6,0),las=2)
for(gen in gg) {
    bb = b[g == gen, ]
    z = dcast(bb, m~s, value.var="qbar")
    z = as.matrix(z[-1,])
    filled.contour(z,col=bcz,levels=keylevels,
        plot.title = {
            title(line=1,main="Molino Cave",cex.main=1.2^3)
            title(xlab="m", ylab="s",cex.lab=1.2^2)
        },
        key.title = title(xlab=expression(bar(q)),line=1,cex.lab=1.2^3),
        plot.axes = {
            box()
            axis(1, at=seq(0,1,length.out=6), xaxs)
            axis(2, at=seq(0,1,length.out=6), yaxs)
        }
    )
    txt = sprintf("Gen %d", gen)
    text(0.80,0.95, txt,col="white",pos=2,cex=1.2^2,font=2)
    
}
invisible(dev.off())

system("convert -delay 10 molino*.png molino.mp4")
unlink("molino*.png")


