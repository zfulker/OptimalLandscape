##########################################
######## Ininital Sampling Functs ########
##########################################

randomRows <- function(n_init, sobol, sample_datagrid){
  #Use sobol to speed up parameter space search, but random to model humans
  if(sobol=="Sobol"){
    return(getSobol(sample_datagrid, n_init))
  }else{
    samp <- sample_datagrid[sample(nrow(sample_datagrid), n_init),]
    while(any(duplicated(samp)==T)){
      samp <- sample_datagrid[sample(nrow(sample_datagrid), n_init),]    
    }
    return(samp)
  }
}

getSobol <- function(sample_datagrid, n_init){
  #NOT CURRENTLY SETUP!!!!!!!!!!!
  # Generates a random Sobol Sequence... this is a quasi-random sequence, with a random seed
  df <- rbind(c(0,0), sobol(n_init, 2, scrambling=3, seed=sample(c(1:1000), 1)))
  colnames(df) <-c('X', 'Y')
  df <- data.frame(df)
  sampy <- min(sample_datagrid$Y) + df$Y * (max(sample_datagrid$Y) - min(sample_datagrid$Y))
  sampx <- min(sample_datagrid$X) + df$X * (max(sample_datagrid$X) - min(sample_datagrid$X))
  out <- data.frame("X"=as.integer(sampx),"Y"=as.integer(sampy))
  out <- out[-1,]
  rownames(out) <- 1:nrow(out)
  return(out)    
}

getInfoNumber <- function(x, y, map, datagrid_dict) { 
  x <- datagrid_dict[[as.character(round(x,3))]]
  y <- datagrid_dict[[as.character(round(y,3))]]
  return(map[x,y])
}

initialize_process <- function(datagrid, datagrid_dict, SearchType, n_init, map, sobol){
  first_digs <- randomRows(n_init, sobol, datagrid)
  history <- data.frame(X=first_digs$X, Y=first_digs$Y, Round=1:n_init, Oil=0)
  for(i in history$Round){
    history[i,]$Oil <- getInfoNumber(history[i,]$X, history[i,]$Y, map, datagrid_dict)
  }
  history$SearchType <- SearchType        
  history$Stop <- 0
  history <- calc_regret(history, map)
  GPredict_curr <- run_gp(history, datagrid)
  return(list(history, GPredict_curr))
}

##############
### Regret ###
##############

calc_regret <- function(history, map){
  global_max <- max(map)
  for(i in history$Round){
    history_new <- history[1:i,]
    current_max <- max(history_new$Oil)
    current_regret <- global_max-current_max
    history$Regret[i] <- current_regret
    history$CumulRegret[i] <- sum(history$Regret)/i
  }
  history$CumulRegret <- history$Regret/history$Round
  return(history)
}

get_avg_results <- function(history, col_name){
  for(i in sort(unique(history$Iteration))){
    iteration_history <- subset(history, Iteration==i)
    search_history <- iteration_history[,col_name, with=FALSE]
    #if(length(search_history) < 15){ # add obs for early termination THIS IS NOT WORKING BC LENGTH = 1!!!!!!!!!!!
    #last_val = search_history[length(search_history)]
    #for(j in (length(search_history)+1):15){search_history[j]<-last_val}
    #}
    if(i == 1){
      avg_cumreg <- search_history
    } else {
      avg_cumreg <- cbind(avg_cumreg, search_history)
    }
  }
  return(rowMeans(avg_cumreg))
}

#########################
### Algorithm Process ###
#########################

normalize_col <- function(x){ 
  a <- min(x) 
  b <- max(x) 
  (x - a)/(b - a) 
} 

run_gp <- function(history, datagrid){
  evidence <- unique(history[,c("X", "Y", "Oil")])
  gp_model <- GP_fit(evidence[,c("X", "Y")], evidence$Oil) #evidence[,c("X", "Y")]
  GPredict <- predict(gp_model, datagrid)
  return(GPredict)
}

sample_next <- function(history, datagrid, datagrid_dict, current_round, map, GP_Prev){
  method <- history$SearchType[1]
  if(method=="GPUCB-Seeking"){
    return_list <- add_gpucb_seeking_point(history, datagrid, datagrid_dict, current_round, map, GP_Prev)
  }
  if(method=="GPUCB-Averse"){
    return_list <- add_gpucb_averse_point(history, datagrid, datagrid_dict, current_round, map, GP_Prev)
  }
  return(return_list)
}

run_algorithm_once <- function(datagrid, map, datagrid_dict, searchType, num_digs, n_init=4, rho0=0.001){
  current_round <- n_init
  return_list <- initialize_process(datagrid, datagrid_dict, searchType, n_init, map, sobol = 'Random')
  history <- data.frame(return_list[1])
  GP_Prev <- return_list[2]
  tick <- 1
  while(current_round < num_digs){
    return_list <- sample_next(history, datagrid, datagrid_dict, current_round, map, GP_Prev)
    history <- data.frame(return_list[1])
    GP_Prev <- return_list[2]
    current_round <- current_round + tick
    #Check stoping criteria
    #if(sum(1-tail(history,L)$Stop)<4*rho0){
    #return(history)
    #}		
  }
  return(history)
}

run_many_times <- function(datagrid, map, datagrid_dict, num_samp_per_model, num_digs, searchType){
  out <- foreach(i = 1:num_samp_per_model, .export=c('run_algorithm_once','initialize_process','randomRows','getInfoNumber','calc_regret','sample_next','run_gp','GP_fit','stopping_criteria','discounted_rank_dissimilarity','get_rank_distances','get_max_distance','add_gpucb_averse_point','add_gpucb_seeking_point')) %dopar% {
    history <- run_algorithm_once(datagrid, map, datagrid_dict, searchType, num_digs)
    history$Iteration <- i
    history
  }
  out <- rbindlist(out)
  return(out)
}

#########################
### Models of digging ###
#########################

# Risk seeking (97.5% of pop) GPUCB Digging: choose coordinate with maximum (Mean)^1.55 + Upper MSE
add_gpucb_seeking_point <- function(history, datagrid, datagrid_dict, current_round, map, GP_Prev){ 
  GPredict_prev <- GP_Prev[[1]]
  
  # find upper confidence bound
  y_plus_sigma  <- (GPredict_prev$Y_hat)**1.55 + 2*(GPredict_prev$MSE)**.5
  UCB_point <- GPredict_prev$complete_data[which.max(y_plus_sigma),]	
  coords <- as.numeric(UCB_point[1:2])
  new_point <- data.frame(X=coords[1], Y=coords[2], Round= current_round +1)	
  new_point$Oil <- getInfoNumber(coords[1], coords[2], map, datagrid_dict)
  new_point$SearchType <- history$SearchType[1]
  new_point$Stop <- 2
  new_point$Regret <- 1
  new_point$CumulRegret <- 1
  history_curr <- rbind(history, new_point)
  history_curr <- calc_regret(history_curr, map)	
  GPredict_curr <- run_gp(history_curr, datagrid)	
  stop <- stopping_criteria(GPredict_prev, GPredict_curr)
  history <- calc_regret(history_curr, map)
  history$Stop[max(history$Round)] <- stop
  return(list(history, GPredict_curr))
}

# Risk averse (2.5% of pop) GPUCB Digging: choose coordinate with maximum (Mean)^0.75 + Upper MSE
add_gpucb_averse_point <- function(history, datagrid, datagrid_dict, current_round, map, GP_Prev){ 
  GPredict_prev <- GP_Prev[[1]]
  
  # find upper confidence bound
  y_plus_sigma  <- (GPredict_prev$Y_hat)**0.75 + 2*(GPredict_prev$MSE)**.5
  UCB_point <- GPredict_prev$complete_data[which.max(y_plus_sigma),]	
  coords <- as.numeric(UCB_point[1:2])
  new_point <- data.frame(X=coords[1], Y=coords[2], Round= current_round +1)	
  new_point$Oil <- getInfoNumber(coords[1], coords[2], map, datagrid_dict)
  new_point$SearchType <- history$SearchType[1]
  new_point$Stop <- 2
  new_point$Regret <- 1
  new_point$CumulRegret <- 1
  history_curr <- rbind(history, new_point)
  history_curr <- calc_regret(history_curr, map)	
  GPredict_curr <- run_gp(history_curr, datagrid)	
  stop <- stopping_criteria(GPredict_prev, GPredict_curr)
  history <- calc_regret(history_curr, map)
  history$Stop[max(history$Round)] <- stop
  return(list(history, GPredict_curr))
}


######################
### Stopping Rules ###
######################
stopping_criteria <- function(GPredict_prev, GPredict_curr){
  # stop when the procedure ceases to learn about the landscape (page 10, Tsunami)
  # this means, measure the global changes in mu between two successive iterations
  # and focus more on the highest values; measure the correlation between mu and mu+1		
  # first, n_v is the size of the validation dataset, X_V
  # and G_nv is the set of all permutations of size n_v
  # pi_t is a ranking function with respect to pi_t_plus1
  prev_df <- data.frame(GPredict_prev$complete_data)
  curr_df <- data.frame(GPredict_curr$complete_data)
  prev_df$num <- c(1:nrow(prev_df)) 
  curr_df$num <- c(1:nrow(prev_df))
  pi_prev_sorted <- prev_df[order(prev_df$Y_hat),]
  pi_curr_sorted <- curr_df[order(curr_df$Y_hat),]
  pi_prev_rank <- pi_prev_sorted$num
  pi_curr_rank <- pi_curr_sorted$num
  numerator <- discounted_rank_dissimilarity(pi_curr_rank, pi_prev_rank)
  denominator <- get_max_distance(length(pi_curr_rank))
  rhoXv <- 1 - numerator/denominator
  return(rhoXv)
}	

get_rank_distances <- function(pi_t1, pi_t0){
  return(c(pi_t1 - pi_t0)**2)
}

discounted_rank_dissimilarity <- function(pi_t1, pi_t0){
  numerator <- get_rank_distances(pi_t1, pi_t0)
  denominator <- pi_t1**2
  dists <- numerator/denominator
  d <- sum(dists)
  return(d)
}

get_max_distance <- function(nv){
  curr_max <- nv
  for(i in 1:nv){
    normal <- c(1:i)
    revers <- c(i:1)
    d <- discounted_rank_dissimilarity(normal,revers)
    if(d>curr_max){
      curr_max <- d
    }
  }
  return(curr_max)
}

########################
### Optimal Design   ###
########################

get_divergance <-  function(GPUCB_Averse_oil, GPUCB_Seeking_oil){
  # averse as baseline 
  averse_div <- KL.plugin(GPUCB_Averse_oil, GPUCB_Seeking_oil)
  
  # seeking as baseline 
  seeking_div <- KL.plugin(GPUCB_Seeking_oil, GPUCB_Averse_oil)
  
  overall_div <- (averse_div +  seeking_div)/2
  return(overall_div)
}

get_landscape_search_results <- function(x, y, z, datagrid, datagrid_dict, optimal_datagrid_dict, num_samp_per_model, num_digs){
  x <- optimal_datagrid_dict[[as.character(round(x,3))]]
  y <- optimal_datagrid_dict[[as.character(round(y,3))]]
  z <- optimal_datagrid_dict[[as.character(round(z,3))]]
  
  map <- getMap(n=x, m=y,exponent=z, dimx=length(datagrid_dict), dimy=length(datagrid_dict))
  
  GPUCB_Averse_history <- run_many_times(datagrid, map, datagrid_dict, num_samp_per_model, num_digs, searchType='GPUCB-Averse')
  GPUCB_Averse_oil <- get_avg_results(GPUCB_Averse_history, 'Oil')
  
  GPUCB_Seeking_history <- run_many_times(datagrid, map, datagrid_dict, num_samp_per_model, num_digs, searchType='GPUCB-Seeking')
  GPUCB_Seeking_oil <- get_avg_results(GPUCB_Seeking_history, 'Oil')
  
  overall_div <- get_divergance(GPUCB_Averse_oil, GPUCB_Seeking_oil)
  return(overall_div)
}

run_gp_optimal <- function(history, datagrid){
  evidence <- history[,c("X", "Y", "Z", "Div")]
  evidence <- evidence[!duplicated(evidence[1:3]),] #place holder until figure out how to do stochastic functions
  gp_model <- GP_fit(evidence[,c("X", "Y", "Z")], evidence$Div) 
  GPredict <- predict(gp_model, datagrid)
  return(GPredict)
}

add_gpucb_pe_point_optimal <- function(history, datagrid, datagrid_dict, optimal_datagrid, current_round, num_samp_per_model, num_digs, k=1){
  ############
  ### UCB step
  ############
  GPredict_prev <- run_gp_optimal(history, optimal_datagrid)
  
  # find upper confidence bound
  y_plus_sigma  <- GPredict_prev$Y_hat + 2*(GPredict_prev$MSE)**.5
  UCB_point <- GPredict_prev$complete_data[which.max(y_plus_sigma),]
  coords <- as.numeric(UCB_point[1:3])
  new_point <- data.frame(X=coords[1], Y=coords[2], Z=coords[3], Round= current_round +1)	
  new_point$Div <- get_landscape_search_results(coords[1], coords[2], coords[3], datagrid, datagrid_dict, optimal_datagrid_dict, num_samp_per_model, num_digs)
  new_point$SearchType <- history$SearchType[1]
  new_point$Stop <- 2
  history <- rbind(history, new_point)
  GPredict_curr <- run_gp_optimal(history, optimal_datagrid)	
  stop <- stopping_criteria(GPredict_prev, GPredict_curr)
  history$Stop[max(history$Round)] <- stop
  
  current_round <- current_round + 1
  
  ###########
  ### PE step
  ###########
  y_minus_sigma <- GPredict_prev$Y_hat - 2*(GPredict_prev$MSE)**.5
  region <- which(y_plus_sigma > max(y_minus_sigma))
  if(length(region) > 1){
    for(i in 1:k){
      space_to_sample <- GPredict_prev$complete_data[region,]	
      PE_point <- space_to_sample[which.max(space_to_sample[,4]),]
      coords <- as.numeric(PE_point[1:3])
      new_point <- data.frame(X=coords[1], Y=coords[2], Z=coords[3], Round= current_round +1)
      new_point$Div <- get_landscape_search_results(coords[1], coords[2], coords[3], datagrid, datagrid_dict, optimal_datagrid_dict, num_samp_per_model, num_digs)
      new_point$SearchType <- history$SearchType[1]
      new_point$Stop <- 2
      history <- rbind(history, new_point)
      GPredict_curr <- run_gp_optimal(history, optimal_datagrid)	
      stop <- stopping_criteria(GPredict_prev, GPredict_curr)
      history$Stop[max(history$Round)] <- stop		
      current_round <- current_round + 1		
    }
    return(history)
  }
  else {
    PE_point <- history[which.max(history$Div),]
    coords <- as.numeric(PE_point[1:3])
    new_point <- data.frame(X=coords[1], Y=coords[2], Z=coords[3], Round= current_round +1)
    new_point$Div <- get_landscape_search_results(coords[1], coords[2], coords[3], datagrid, datagrid_dict, optimal_datagrid_dict, num_samp_per_model, num_digs)
    new_point$SearchType <- history$SearchType[1]
    new_point$Stop <- 2
    history <- rbind(history, new_point)
    GPredict_curr <- run_gp_optimal(history, optimal_datagrid)	
    stop <- stopping_criteria(GPredict_prev, GPredict_curr)
    history$Stop[max(history$Round)] <- stop					
    return(history)		
  }
}

initialize_process_optimal <- function(datagrid, datagrid_dict, optimal_datagrid, optimal_datagrid_dict, SearchType, n_init, num_samp_per_model, num_digs, sobol){
  first_digs <- randomRows(n_init, sobol, optimal_datagrid)
  history <- data.frame(X=first_digs$X, Y=first_digs$Y, Z=first_digs$Z, Round=1:n_init, Div=0)
  for(i in history$Round){
    history[i,]$Div <- get_landscape_search_results(history[i,]$X, history[i,]$Y, history[i,]$Z, datagrid, datagrid_dict, optimal_datagrid_dict, num_samp_per_model, num_digs)
  }
  history$SearchType <- SearchType        
  history$Stop <- 0
  return(history)
}

run_optimal_search <- function(datagrid, datagrid_dict, optimal_datagrid, optimal_datagrid_dict, num_samp_per_model, num_digs, failsafe = 30, rho0=0.001){
  current_round <- 4
  history <- initialize_process_optimal(datagrid, datagrid_dict, optimal_datagrid, optimal_datagrid_dict, 'GPUCB-PE', 4, num_samp_per_model, num_digs, sobol = 'Random')
  tick <- 2
  L <- 4
  for(i in 1:L){
    history <- add_gpucb_pe_point_optimal(history, datagrid, datagrid_dict, optimal_datagrid, current_round, num_samp_per_model, num_digs)
    current_round <- current_round + tick
  }
  while(current_round < failsafe){
    history <- add_gpucb_pe_point_optimal(history, datagrid, datagrid_dict, optimal_datagrid, current_round, num_samp_per_model, num_digs)
    current_round <- current_round + tick
    if(sum(1-tail(history,L)$Stop)<4*rho0){
      return(history)
    }		
  }
  return(history)
}