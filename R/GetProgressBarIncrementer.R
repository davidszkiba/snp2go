## This function is used to show the progress 

GetProgressBarIncrementer <- function(verbose, environment, min, max, style) {
  if(verbose && max > min ){
    progEnv <- environment
    
    if(!exists("pb", envir=progEnv)) {
      assign("pb", txtProgressBar(min=min, max=max, style=style), envir=progEnv)
    } 
    else { 
      print("pb exists!") 
    }
    if(!exists("pbCounter", envir=progEnv)) {
      assign("pbCounter", min, envir=progEnv)
    }
    else { 
      print("pbCounter exists!") 
    }
    if(!exists("pbMax", envir=progEnv)) {
      assign("pbMax", max, envir=progEnv)
    }
    else { 
      print("pbMax exists!") 
    }	

    ## Return a function
    function(progEnv){
      with(progEnv, pbCounter <- pbCounter + 1)
      with(progEnv, setTxtProgressBar(pb, pbCounter))
      if(with(progEnv, pbCounter >= pbMax)) {
        with(progEnv, close(pb))
        rm("pb", envir=progEnv)
        rm("pbCounter", envir=progEnv)
        rm("pbMax", envir=progEnv)
      }
    }
  }
  else {
    function(foo){} # return empty function
  }
}
