.First.lib = function(lib, pkg) {

cat("attaching languageR\n")

require(methods, quietly = TRUE)

#setMethod("show", "corres", function(object) print.corres(object))
#setMethod("show", "growth", function(object) show.growth(object))

options(show.signif.stars=FALSE)
}
