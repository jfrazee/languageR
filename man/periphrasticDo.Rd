\name{periphrasticDo}
\alias{periphrasticDo}
\docType{data}
\title{The development of periphrastic do in English}
\description{
  The development of periphrastic \emph{do} in English: Ellegard's counts
  for the use of \emph{do} across four sentence types in 11 consecutive
  time periods between 1390 and 1710.
}
\usage{data(periphrasticDo)}
\format{
  A data frame with 44 observations on the following 5 variables.
  \describe{
    \item{\code{begin}}{a numeric vector with beginnings of the time periods
    used by Ellegard.}
    \item{\code{end}}{a numeric vector with ends of these time periods.}
    \item{\code{type}}{a factor for sentence type, with levels 
       \code{affdecl} (affirmative declarative),
       \code{affquest} (affirmative question), 
       \code{negdecl} (negative declarative) and 
       \code{negquest} (negative question).}
    \item{\code{do}}{a numeric vector with the count of sentences with
       \emph{do}.}
    \item{\code{other}}{a numeric vector with the count of sentences without
       \emph{do}.}
  }
}
\source{
Ellegard, A. (1953) \emph{ The auxiliary do: The establishment and regulation
of its use in English}, Stockholm: Almquist & Wiksell.
}
\references{
Vulanovic, R. and Baayen, R. H. (2006) Fitting the development of periphrastic
do in all sentence types, in Grzybek, P. and Koehler, R. (eds.), Festschrift
fuer Gabriel Altmann, Berlin: Walter de Gruyter, p. 679-688.
}
\examples{
data(periphrasticDo)

# add midpoints of time periods

periphrasticDo$year = periphrasticDo$begin +
  (periphrasticDo$end-periphrasticDo$begin)/2

# and ad an indicator variable distinguishing the first three time periods
# from the others

periphrasticDo$Indicator = rep(c(rep(0, 3), rep(1, 8)), 4)

# fit a logistic regression model

periphrasticDo.glm = glm(cbind(do, other) ~
(year + I(year^2) + I(year^3)) * type + Indicator * type + 
Indicator * year, data = periphrasticDo, family = "binomial")

anova(periphrasticDo.glm, test = "F")

# visualization of data and model predictions

periphrasticDo$predict = predict(periphrasticDo.glm, type = "response")
par(mfrow=c(2, 2))
for (i in 1:nlevels(periphrasticDo$type)) {
  subset = periphrasticDo[periphrasticDo$type == 
    levels(periphrasticDo$type)[i], ]
  plot(subset$year,
    subset$do/(subset$do + subset$other), 
    type = "p", ylab = "proportion", xlab = "year", 
    ylim = c(0, 1), xlim = c(1400, 1700))
  mtext(levels(periphrasticDo$type)[i], line = 2)
  lines(subset$year, subset$predict, lty = 1)
}
par(mfrow=c(1, 1))

}
\keyword{datasets}
