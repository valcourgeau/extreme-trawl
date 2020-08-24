ExponentialTrawl <- new.env()

ExponentialTrawl$TrawlB1 <- function(param, h){
  stopifnot(length(param) == 1)
  return((1.0- exp(-param*h))/param)
}

ExponentialTrawl$TrawlB2 <- function(param, h){
  stopifnot(length(param) == 1)
  return(exp(-param*h)/param)
}

ExponentialTrawl$TrawlB3 <- function(param, h){
  return(ExponentialTrawl$TrawlB1(param, h))
}

ExponentialTrawl$Config <- function(){
  return(list(n_params=1, lower=0.05, upper=0.99))
}

