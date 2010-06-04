\name{faz}
\alias{faz}
\docType{data}
\title{Frankfurter frequencies}
\description{
Frequencies of references to previous years in issues of the Frankfurter 
Allgemeine Zeiting published in 1994.
}
\usage{data(faz)}
\format{
  A data frame with 800 observations on the following 2 variables.
  \describe{
    \item{\code{Year}}{a numeric vector coding years referenced in articles
      published in 1994.}
    \item{\code{Frequency}}{a numeric vector for the frequencies with which
      years are referenced.}
  }
}
\references{
Pollman, T. and Baayen, R. H. (2001)
Computing historical consciousness.  A quantitative inquiry
into the presence of the past in newspaper texts, 
\emph{Computers and the Humanities}, 35, 237-253.
}
\examples{
\dontrun{
data(faz)
faz$Distance = 1:nrow(faz)

# ---- visualization

plot(log(faz$Distance), log(faz$Frequency + 1), 
xlab = "log Distance", ylab = "log Frequency")
abline(v = log(49), lty=1, col="red")   # 1945
abline(v = log(54), lty=1, col="red")   # 1940
abline(v = log(76), lty=2, col="blue")  # 1918
abline(v = log(80), lty=2, col="blue")  # 1914

# ---- breakpoint analysis

deviances = rep(0, nrow(faz)-1)
faz$LogFrequency = log(faz$Frequency + 1)
faz$LogDistance = log(faz$Distance)
for (pos in 1 : (nrow(faz)-1)) {                             # be patient
  breakpoint = log(pos)
  faz$ShiftedLogDistance = faz$LogDistance - breakpoint
  faz$PastBreakPoint = as.factor(faz$ShiftedLogDistance > 0)
  faz.both = lm(LogFrequency~ShiftedLogDistance:PastBreakPoint, data = faz)
  deviances[pos] = deviance(faz.both)
}
breakpoint = log(which(deviances == min(deviances)))

# ---- refit and plot

faz$ShiftedLogDistance = faz$LogDistance - breakpoint
faz$PastBreakPoint = as.factor(faz$ShiftedLogDistance > 0)
faz.both = lm(LogFrequency ~ ShiftedLogDistance:PastBreakPoint, data = faz)

plot(faz$LogDistance, faz$LogFrequency, 
xlab = "log Distance", ylab = "log Frequency", col = "darkgrey")
lines(faz$LogDistance, fitted(faz.both))
}

}
\keyword{datasets}
