.First.lib = function(lib, pkg) {

packageStartupMessage("attaching languageR\n")

#setMethod("show", "corres", function(object) print.corres(object))
#setMethod("show", "growth", function(object) show.growth(object))

options(show.signif.stars=FALSE)
}
