\name{affixProductivity}
\alias{affixProductivity}
\docType{data}
\title{Affix productivity}
\description{
Affix productivity, gauged by the P* productivity measure, 
for 27 English affixes in 44 texts.
}
\usage{data(affixProductivity)}
\format{
  A data frame with 44 observations on the following 30 variables.
  \describe{
    \item{\code{semi}}{a numeric vector of P*-values}
    \item{\code{anti}}{a numeric vector of P*-values}
    \item{\code{ee}}{a numeric vector of P*-values}
    \item{\code{ism}}{a numeric vector of P*-values}
    \item{\code{ian}}{a numeric vector of P*-values}
    \item{\code{ful}}{a numeric vector of P*-values}
    \item{\code{y}}{a numeric vector of P*-values}
    \item{\code{ness}}{a numeric vector of P*-values}
    \item{\code{able}}{a numeric vector of P*-values}
    \item{\code{ly}}{a numeric vector of P*-values}
    \item{\code{unV}}{a numeric vector of P*-values}
    \item{\code{unA}}{a numeric vector of P*-values}
    \item{\code{ize}}{a numeric vector of P*-values}
    \item{\code{less}}{a numeric vector of P*-values}
    \item{\code{erA}}{a numeric vector of P*-values}
    \item{\code{erC}}{a numeric vector of P*-values}
    \item{\code{ity}}{a numeric vector of P*-values}
    \item{\code{super}}{a numeric vector of P*-values}
    \item{\code{est}}{a numeric vector of P*-values}
    \item{\code{ment}}{a numeric vector of P*-values}
    \item{\code{ify}}{a numeric vector of P*-values}
    \item{\code{re}}{a numeric vector of P*-values}
    \item{\code{ation}}{a numeric vector of P*-values}
    \item{\code{in.}}{a numeric vector of P*-values}
    \item{\code{ex}}{a numeric vector of P*-values}
    \item{\code{en}}{a numeric vector of P*-values}
    \item{\code{be}}{a numeric vector of P*-values}
    \item{\code{AuthorCodes}}{a factor with levels 
      \describe{
        \item{\code{BLu}}{(King James Version: Luke-Acts)}
        \item{\code{BMo}}{(Book of Mormon) }
        \item{\code{CAs}}{(Aesop's fables, translation by Townsend)}
        \item{\code{CBo}}{(Baum, The Marvelous Land of Oz)}
        \item{\code{CBp}}{(Barrie, Peter Pan and Wendy)}
        \item{\code{CBw}}{(Baum, The Wonderful Wizard of Oz)}
        \item{\code{CCa}}{(Carroll, Alice's Adventures in Wonderland)}
        \item{\code{CCt}}{(Carroll, Through the Looking Glass and what Alice Found There)}
        \item{\code{CGr}}{(Grimm Fairy Tales, translations)}
        \item{\code{CKj}}{(Kipling, The Jungle Book)}
        \item{\code{LAp}}{(Austen, Pride and Prejudice)}
        \item{\code{LBp}}{(Burroughs, A Princess of Mars)  }
        \item{\code{LBw}}{(Bronte, Wuthering Heights) }
        \item{\code{LCl}}{(Conrad, Lord Jim)}
        \item{\code{LCn}}{(Conrad, Nigger of the Narcissus) }
        \item{\code{LDb}}{(Doyle, The Casebook of Sherlock Holmes)}
        \item{\code{LDc}}{(Dickens, The Chimes: a Goblin Story) }
        \item{\code{LDC}}{(Dickens, A Christmas Carol)     }
        \item{\code{LDh}}{(Doyle, The Hound of the Baskervilles)}
        \item{\code{LDv}}{(Doyle, The Valley of Fear)  }
        \item{\code{LJc}}{(James, Confidence)}
        \item{\code{LJe}}{(James, The Europeans)}
        \item{\code{LLc}}{(London, The Call of the Wild)}
        \item{\code{LLs}}{(London, The Sea Wolf)}
        \item{\code{LMa}}{(Montgomery, Anne of Avonlea)}
        \item{\code{LMm}}{(Melville, Moby Dick)}
        \item{\code{LMn}}{(Morris, News from Nowhere)}
        \item{\code{LMp}}{(Milton, Paradise Lost)}
        \item{\code{LOs}}{(Orczy, The Scarlet Pimpernel)}
        \item{\code{LSd}}{(Stoker, Dracula)}
        \item{\code{LSs}}{(Chu, More than a Chance Meeting (Startrek))}
        \item{\code{LTa}}{(Trollope, Ayala's Angel)}
        \item{\code{LTe}}{(Trollope, The Eustace Diamonds)}
        \item{\code{LTf}}{(Trollope, Can you Forgive her?)}
        \item{\code{LTy}}{(Twain, A Connecticut Yankee in King Arthur's Court)}
        \item{\code{LWi}}{(Wells, The Invisible Man)}
        \item{\code{LWt}}{(Wells, The Time Machine)}
        \item{\code{LWw}}{(Wells, The War of the Worlds)}
        \item{\code{OAf}}{(The Federalist Papers)}
        \item{\code{OCh}}{(Texts sampled from Congress Hearings) }
        \item{\code{OCl}}{(Texts sampled from Clinton's Election Speeches)}
        \item{\code{ODo}}{(Darwin, On the Origin of the Species)}
        \item{\code{OGa}}{(Selected Texts from the Government Accounting Office)}
        \item{\code{OJe}}{(James, Essays in Radical Empiricism)}
      }
    }
    \item{\code{Registers}}{a factor with levels \code{B} (Biblical texts) 
    \code{C} (Children's books) \code{L} (Literary texts) \code{O} (other)}
    \item{\code{Birth}}{a numeric vector for the author's year of birth 
         (where available)}
  }
}
\source{
 Most texts were obtained from the Gutenberg Project 
 (\url{http://www.gutenberg.org/wiki/Main_Page}) and the Oxford Text
 Archive (http://ota.ahds.ac.uk/).
}
\references{
  Baayen, R. H. (1994) Derivational Productivity and Text Typology,
  \emph{Journal of Quantitative Linguistics}, 1, 16-34.
}
\examples{
\dontrun{
data(affixProductivity)
affixes.pr = prcomp(affixProductivity[,1:(ncol(affixProductivity)-3)], 
center = TRUE, scale. = TRUE)
library(lattice)
trellis.device()
super.sym = trellis.par.get("superpose.symbol")
splom(data.frame(affixes.pr$x[,1:3]), 
groups = affixProductivity$Registers, 
panel = panel.superpose,
key = list(title  = "texts in productivity space",
text   = list(c("Religious", "Children", "Literary", "Other")),
points = list(pch = super.sym$pch[1:4], col = super.sym$col[1:4])))
}
}
\keyword{datasets}
