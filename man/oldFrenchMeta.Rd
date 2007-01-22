\name{oldFrenchMeta}
\alias{oldFrenchMeta}
\docType{data}
\title{Meta data for the oldFrench data}
\description{
Meta data for the oldFrench data, a matrix of frequencies for texts (rows) by tag trigrams (columns).
The meta data provide information on the texts, manuscript variants, their authors, 
their region and approximate date of origin, their general topic, and their genre.
}
\usage{data(oldFrenchMeta)}
\format{
  A data frame with 342 observations on the following 7 variables.
  \describe{
    \item{\code{Textlabels}}{a factor with texts coded as follows:
      \describe{
         \item{\code{Abe}}{ J. de Meun, Traduction de la premiere epitre de P. Abelard, 1--821} 
         \item{\code{Hyl1}}{ Anon, La vie de saint Hylaire} 
         \item{\code{Art}}{ J. de Meun, L'art de chevalerie} 
         \item{\code{Bar}}{ Anon, L'histoire de Barlaam et Josaphat} 
         \item{\code{Cathy}}{ Anon, La vie de sainte Catherine d'Alexandrie} 
         \item{\code{Hyl2}}{ Anon, La vie de saint Hylaire} 
         \item{\code{Ch1}}{ Chretien de Troyes, Le Chevalier au lion} 
         \item{\code{Ch2}}{ Chretien de Troyes, Le chevalier au lion} 
         \item{\code{Clari}}{ Robert de Clari, La conquete de Constantinople} 
         \item{\code{Marie}}{ Rutebeuf, Sainte Marie l'Egyptienne} 
         \item{\code{Fab4c}}{ Anon, Fabliau nr 4 ms C} 
         \item{\code{Fab4e}}{ Anon, Fabliau nr 4 ms E} 
         \item{\code{Fab4f}}{ Anon, Fabliau 4f} 
         \item{\code{Faba}}{ Anon, Fabliaux nrs 1,2,4,23 et 29 du ms A} 
         \item{\code{Fabb}}{ Anon, Fabliaux nrs 2 et 4 du ms B} 
         \item{\code{Fabd.COD}}{ Anon, Fabliaux nrs 2 et 4 du ms D} 
         \item{\code{Hyl3}}{ Anon,  La vie de saint Hylaire} 
         \item{\code{Jacobi}}{ Pierre de Beauvais, The Liber Sancti Jacobi} 
         \item{\code{Louis}}{ J. de Joinville, La vie de saint Louis} 
         \item{\code{Cathy1}}{ Anon, La passion saynte Katherine} 
         \item{\code{Lancelot}}{ Anon, Lancelot do Lac, p. 1.1--20.13} 
         \item{\code{Merlin1}}{ Merlin, Robert  de Boron} 
         \item{\code{Marga}}{ Anon, La vie de Sainte Marguerite de Wace} 
         \item{\code{Martin}}{ Anon, Leben und Wunderthaten des heiligen Martin} 
         \item{\code{Merlin2}}{ Anon, Merlin, p.1--29 (ms.  Huth)} 
         \item{\code{RoseA}}{ J. de Meun, Le Roman de la Rose} 
         \item{\code{Arthur}}{ Anon, La mort le roi Artu, par.1--35}
         \item{\code{NimAf}}{ Anon, charroi de Nimes, ms. A, fragment} 
         \item{\code{NimB1}}{ Anon, Le charroi de Nimes, ms B1} 
         \item{\code{NimB2}}{ Anon, Le charroi de Nimes, ms B2} 
         \item{\code{Nouvel}}{ Jacquemart Gielee, Renart le Nouvel} 
         \item{\code{Jehan}}{ Anon, La vie de saint Jehan Bouche d'Or}
         \item{\code{Per0}}{ Chretien de Troyes, Perceval} 
         \item{\code{perL}}{ Chretien de Troyes, Perceval}
         \item{\code{PerQ}}{ Chretien de Troyes, Perceval} 
         \item{\code{PerS}}{ Chretien de Troyes, Perceval}
         \item{\code{PerU}}{ Chretien de Troyes, Perceval} 
         \item{\code{Queste}}{ Anon, La queste del saint Graal, p.1.1--41.17} 
         \item{\code{Rob}}{ Anon, Robert le Diable, v.1--808} 
         \item{\code{RomB}}{ Anon, Le Roman de Renart, br.VI, ms B} 
         \item{\code{RomD}}{ Anon, Roman de Renard, br.VI, ms D}
         \item{\code{RomL}}{ Anon, Roman de Renard, br. VI, ms L} 
         \item{\code{RomO}}{ Anon, Le Roman de Renart, br. VI, ms 0} 
         \item{\code{RoseB}}{ Guillaume de Lorris, Le roman de la rose} 
         \item{\code{Sapi}}{ Anon, Sermo de sapientia, dans: Li dialoge Gregoire lo pape} 
         \item{\code{Troi}}{ Anon, Le roman de Troie en prose, par.1--19} 
         \item{\code{Conqueste}}{ Josfroi de Vileharduyn, La conqueste de Costentinoble} 
         \item{\code{YvA}}{ Chretien de Troyes, Le chevalier au lion, v.1--1000} 
         \item{\code{YvP}}{ Chretien de Troyes, Le chevalier au lion, v.1--1000} 
         \item{\code{YvS}}{ Chretien de Troyes, Le chevalier au lion, v.1--1000} 
         \item{\code{YvV}}{ Chretien de Troyes, Le chevalier au lion, v 1--1000}
      }
    }
    \item{\code{Codes}}{a factor with manuscript variants, indicated by extensions to the text codes.}
    \item{\code{Author}}{a factor with levels \code{Anon}, \code{ChretienDeTroyes}, \code{GuillaumeDeLorris}, 
         \code{Joinville}, \code{Meun}, \code{NouvelRenart}, \code{PierreDeBeauvais}, \code{RobertDeBoron},
         \code{RobertDeClari}, \code{RobertLeDiable}, \code{Rutebeuf}, and \code{Villeharduyn}.}
    \item{\code{Topic}}{a factor with levels \code{Knight}, \code{Other}, and \code{Saint}.}
    \item{\code{Genre}}{a factor with levels \code{poetry} and \code{prose}.}
    \item{\code{Region}}{a factor with levels \code{R1} (Picardie), \code{R2} (Champenois), and
          \code{R3} (Nievre-Allier).}
    \item{\code{Year}}{a numeric vector indicating approximate year of origin.}
  }
}
\source{
Data from Nouveau Corpus d'Amsterdam, \url{http://www.uni-stuttgart.de/lingrom/stein/corpus/}.
}
\references{
Ernestus, M., van Mulken, M. and Baayen, R. H. (2007)
De syntax van Oud-Franse ridders en heiligen in ruimte en tijd
To appear in \emph{Onze Taal}.
}
\examples{
data(oldFrench)
data(oldFrenchMeta)

oldFrench.ca = corres.fnc(oldFrench)

plot(oldFrench.ca, rlabels = oldFrenchMeta$Genre, 
rcol = as.numeric(oldFrenchMeta$Genre), rcex = 0.5, 
extreme = 0.1, ccol = "blue")
}
\keyword{datasets}
