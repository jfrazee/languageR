\name{compare.richness.fnc}
\alias{compare.richness.fnc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compare Lexical Richness of Two Texts}
\description{
Comparisons of lexical richness between two texts are carried out on the basis
of the vocabulary size (number of types) and on the basis of the vocabulary
growth rate.  Variances of the number of types and of the number of hapax
legomena required for the tests are estimated with the help of LNRE models.
}
\usage{
compare.richness.fnc(text1, text2, digits = 5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{text1}{First text in the comparison.}
  \item{text2}{Second text in the comparison.}
  \item{digits}{Number of decimal digits required for the growth rate.}
}
\details{
  The comparison for the vocabulary size is carried out with the test statistic

  \deqn{Z = \frac{E[V_1] - E[V_2]}{\sqrt{\sigma(V_1)^2 + \sigma(V_2)^2}}}{Z = (E[V_1] - E[V_2])/sqrt(VAR[V_1] + VAR[V_2])}

  and the comparison of the growth rates with the test statistic

  \deqn{Z = \frac{\frac{1}{N_1}E[V_1(1)] - \frac{1}{N_2}E[V_2]}{\sqrt{\frac{1}{N_1^2}\sigma(V_1(1))^2 + \frac{1}{N_2^2}\sigma(V_2(1))^2}}}{Z = (E[V_1(1)]/N_1 - E[V_2(1)]/N_2)/sqrt(VAR[V_1(1)]/N_1^2 + VAR[V_2(1)]/N_2^2)}

  where \eqn{N} denotes the sample size in tokens, \eqn{V} the vocabulary size,
  and \eqn{V(1)} the number of hapax legomena.
}
\value{
  A summary listing the Chi-Squared measure of goodness of fit for the
  LNRE models (available in the zipfR package) used to estimate variances, 
  a table listing tokens, types, hapax legomena and the vocabulary growth rate, 
  and two-tailed tests for differences in the vocabulary sizes and growth rates with
  Z-score and p-value.
}
\references{
Baayen, R. H. (2001) \emph{Word Frequency Distributions},
Kluwer Academic Publishers, Dordrecht.
}
\author{
R. Harald Baayen 
  Radboud University Nijmegen and 
  Max Planck Institute for Psycholinguistics, Nijmegen, The Netherlands. 
  baayen@mpi.nl
}
\note{
  It is probably unwise to attempt to apply this function to texts comprising
  more than 500,000 words.
}
%\seealso{See Also \code{\link{growth.fnc}}, \code{\link{growth2vgc.fnc}}, 
%\code{\link{text2spc.fnc}} and the \code{\link{zipfR}} package.
%}
\examples{
data(alice, through, oz)
compare.richness.fnc(tolower(alice), tolower(through[1:length(alice)]))
compare.richness.fnc(tolower(alice), tolower(oz[1:25942]))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
