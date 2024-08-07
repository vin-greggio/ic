library(wavethresh)

myspec <- cns(1024)

myspec <- putD(myspec, level=2, sin(seq(from=0, to=4*pi, length=1024))^2)
myspec <- putD(myspec, level=4, cos(seq(from=0, to=4*pi, length=1024))^2)

burst <- c(rep(0,200), rep(1,100), rep(0,724))

myspec <- putD(myspec, level=9, v=burst)
myLSWproc <- LSWsim(myspec)

ts.plot(myLSWproc, main="My Spectrum")