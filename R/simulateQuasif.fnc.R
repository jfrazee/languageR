simulateQuasif.fnc =
  function(dat,  with = TRUE, nruns = 100, nsub = NA, nitem = NA, ...) {

  require("MASS", quietly = TRUE, character = TRUE)

  data(simulateQuasiffnc)
  return(simulateQuasif.fnc(dat=dat, with=with, nruns=nruns, nsub=nsub, nitem=nitem,...))
}
