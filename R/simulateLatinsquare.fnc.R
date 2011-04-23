simulateLatinsquare.fnc = 
  function(dat, with = TRUE, trial = 0, nruns = 100, nsub = NA, nitem = NA, ...) {

  # to avoid namespace conflicts with lme4, we redefine the function
  data(simulateLatinsquarefnc)
  return(simulateLatinsquare.fnc(dat=dat, with=with, trial=trial, nruns=nruns, nsub=nsub, nitem=nitem, ...))
  

}

