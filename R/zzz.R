.First.lib = function(lib, pkg) {

cat("attaching languageR\n")

require(methods, quietly = TRUE)

setClass("corres", representation(data="list"))
setMethod("show", "corres", function(object) print.corres(object))

setClass("growth", representation(data="list"))
setMethod("show", "growth", function(object) show.growth(object))

options(show.signif.stars=FALSE)
}
