## script for plotting inserts length
source('common.R')

dat = readSingleColumnFile('insert.wig', 'a-zA-Z')
dat = dat[dat != 0]
m = mean(dat)
s = sd(dat)
dat = abs(dat - m)

plot(dat, cex=0.3, col="blue", pch=19, title="",
    xlab="insert number", ylab="insert length")
abline(h = s * 2, col = "orange", cex=1.5)
abline(h = s * 3, col = "red", cex=1.5)
dev.copy(png,'inserts.png')
dev.off()
