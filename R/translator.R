ParametrisationTranslator <- function(params, parametrisation, target='noven', target_alpha=1.0){
  testit::assert(length(params) == 3)
  testit::assert(parametrisation %in% c('standard', 'noven'))

  # from parametrisation to noven
  params_target <- params

  if(parametrisation == target){
    return(params)
  }

  if(parametrisation == 'standard' & target == 'noven'){
    params_target[1] <- 1/params[1]
    params_target[2] <- params[2]/abs(params[1]) - params[3]
  }else{
    if(parametrisation == 'noven' & target == 'standard'){
      print('o')
      params_target[1] <- 1/params[1]
      params_target[2] <- (params[2] + params[3])/abs(params[1])
      print(params)
      print(params_target)
      print('o')
    }else{
      if(parametrisation == 'standard' & target == 'transform'){
        params_target[1] <- 1.0/target_alpha
        params_target[2] <- (1.0+params[3])*abs(target_alpha)
      }else{
        if(parametrisation == 'noven' & target == 'transform'){
          params_target[1] <- 1.0/target_alpha
          params_target[2] <- (1.0+params[3])/abs(target_alpha)
        }
      }
    }
  }

  return(params_target)
}
