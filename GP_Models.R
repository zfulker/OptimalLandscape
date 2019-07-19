##########################################
######## Ininital Sampling Functs ########
##########################################

randomRows <- function(n_init, sobol, sample_datagrid){
  #range=c(1,length(datagrid_dict))
  #avoid searching the corner, reflects human behavior (removed for now)
  #sample_datagrid <- expand.grid(X=seq(range[1], range[2], 1), Y=seq(range[1], range[2], 1)) #+20, -20
  #sample_datagrid <- as.data.frame(apply(sample_datagrid,2,normalize_col))
  if(sobol=="Sobol"){
    return(getSobol(sample_datagrid, n_init))
  } else {
    samp <- sample_datagrid[sample(nrow(sample_datagrid), n_init),]
    while(any(duplicated(samp)==T)){
      samp <- sample_datagrid[sample(nrow(sample_datagrid), n_init),]    
    }
    return(samp)
  }
}

getSobol <- function(sample_datagrid, n_init){
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
  history <- data.frame(X=first_digs$X, Y=first_digs$Y,
                        Round=1:n_init, Oil=0)
  for(i in history$Round){
    history[i,]$Oil <- getInfoNumber(history[i,]$X, history[i,]$Y, map, datagrid_dict)
  }
  history$SearchType <- SearchType        
  history$Stop <- 0
  history <- calc_regret(history, map)
  return(history)
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
    search_history <- iteration_history[,col_name]
    if(length(search_history) < 15){ # add obs for early termination
      last_val = search_history[length(search_history)]
      for(j in (length(search_history)+1):15){search_history[j]<-last_val}
    }
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

sample_next <- function(history, datagrid, datagrid_dict, current_round, map){
  method <- history$SearchType[1]
  if(method=="Random"){
    history <- add_random_sampled_point(history, datagrid, datagrid_dict, current_round, map)
  }
  if(method=="GPUCB"){
    history <- add_gpucb_point(history, datagrid, datagrid_dict, current_round, map)
  }
  if(method=="GPUCB_PE"){
    history <- add_gpucb_pe_point(history, datagrid, datagrid_dict, current_round, map)
  }
  if(method=="Exploration"){
    history <- add_exploration_point(history, datagrid, datagrid_dict, current_round, map)
  }
  if(method=="Exploitation"){
    history <- add_exploitation_point(history, datagrid, datagrid_dict, current_round, map)
  }
  return(history)
}

run_algorithm_once <- function(datagrid, map, datagrid_dict, searchType, n_init=4, rho0=0.001, failsafe = 15, L = 4){
  current_round <- n_init
  history <- initialize_process(datagrid, datagrid_dict, searchType, n_init, map, sobol = 'Random')
  tick <- 1
  if(searchType == "GPUCB_PE"){
    tick <- 2
    L <- round(L/2)
  }
  for(i in 1:L){
    history <- sample_next(history, datagrid, datagrid_dict, current_round, map)
    current_round <- current_round + tick
  }
  while(current_round < failsafe){
    history <- sample_next(history, datagrid, datagrid_dict, current_round, map)
    current_round <- current_round + tick
    #print(paste0("Current Round is: ", current_round," and stop == ",sum(1-tail(history,L)$Stop)))
    if(sum(1-tail(history,L)$Stop)<4*rho0){
      return(history)
    }		
  }
  return(history)
}

run_many_times <- function(datagrid, map, datagrid_dict, searchType, n_init=4, rho0=0.001, failsafe = 15, L = 4, n_times=25){
  out <- data.frame()
  i <- 1
  while(i <= n_times){
    history <- run_algorithm_once(datagrid, map, datagrid_dict, searchType)
    #		history <- try(run_algorithm_once(total_grid,search,n_init,rho0,failsafe,L))
    #		if(class(history)!="try-error") {
    history$Iteration <- i
    out <- rbind(out, history)
    #print(paste("Finished iteration: ", i, sep=""))
    i <- i + 1
    #		} 
  }
  return(out)
}

#########################
### Models of digging ###
#########################

# Random Digging
add_random_sampled_point <- function(history, datagrid, datagrid_dict, current_round, map){
  GPredict_prev <- run_gp(history, datagrid)
  coords <- as.numeric(datagrid[sample(nrow(datagrid), 1),])
  while(all(coords %in% history[,c("X", "Y")])){ 
    coords <- as.numeric(datagrid[sample(nrow(datagrid), 1),])
  }
  new_point <- data.frame(X=coords[1], Y=coords[2], Round=current_round+1)
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
  return(history)
}

# Exploit Digging: choose coordinate with maximum mean
add_exploitation_point <- function(history, datagrid, datagrid_dict, current_round, map){ 
  GPredict_prev <- run_gp(history, datagrid)
  
  exploit_point <- GPredict_prev$complete_data[which.max(GPredict_prev$Y_hat),]	
  coords <- as.numeric(exploit_point[1:2])
  new_point <- data.frame(X=coords[1], Y=coords[2], Round= current_round+1)
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
  return(history)
}


# Explore Digging: choose coordinate with maximum Mean Squared Error
add_exploration_point <- function(history, datagrid, datagrid_dict, current_round, map){
  GPredict_prev <- run_gp(history, datagrid)
  
  explore_point <- GPredict_prev$complete_data[which.max(GPredict_prev$MSE),]	
  coords <- as.numeric(explore_point[1:2])
  new_point <- data.frame(X=coords[1], Y=coords[2], Round= current_round+1)
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
  return(history)
}

# GPUCB Digging: choose coordinate with maximum Mean + Upper MSE
add_gpucb_point <- function(history, datagrid, datagrid_dict, current_round, map){ 
  GPredict_prev <- run_gp(history, datagrid)
  
  # find upper confidence bound
  y_plus_sigma  <- GPredict_prev$Y_hat + 2*(GPredict_prev$MSE)**.5
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
  return(history)
}

# GPUCB Digging: choose coordinate with maximum Mean + Upper MSE, AND a point
# with the highest *lower bound* within a certain region
add_gpucb_pe_point <- function(history, datagrid, datagrid_dict, current_round, map, k=1){
  ############
  ### UCB step
  ############
  GPredict_prev <- run_gp(history, datagrid)
  
  # find upper confidence bound
  y_plus_sigma  <- GPredict_prev$Y_hat + 2*(GPredict_prev$MSE)**.5
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
      coords <- as.numeric(PE_point[1:2])
      new_point <- data.frame(X=coords[1], Y=coords[2], Round= current_round+1)
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
      current_round <- current_round + 1		
    }
    return(history)
  }
  else {
    PE_point <- history[which.max(history$Oil),]
    coords <- as.numeric(PE_point[1:2])
    new_point <- data.frame(X=coords[1], Y=coords[2], Round= current_round+1)
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
    return(history)		
  }
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

get_divergance <-  function(rand_oil, exploit_oil, explore_oil, GPUCB_oil, GPUCB_PE_oil){
  # rand as baseline 
  averaged_against_rand <- rowMeans(cbind(exploit_oil, explore_oil, GPUCB_oil, GPUCB_PE_oil))
  rand_div <- KL.plugin(rand_oil, averaged_against_rand)
  
  # exploit as baseline
  averaged_against_exploit <- rowMeans(cbind(rand_oil, explore_oil, GPUCB_oil, GPUCB_PE_oil))
  exploit_div <- KL.plugin(exploit_oil, averaged_against_exploit)
  
  # explore as baseline
  averaged_against_explore <- rowMeans(cbind(rand_oil, exploit_oil, GPUCB_oil, GPUCB_PE_oil))
  explore_div <- KL.plugin(explore_oil, averaged_against_explore)
  
  # GPUCB as baseline
  averaged_against_GPUCB <- rowMeans(cbind(rand_oil, exploit_oil, explore_oil, GPUCB_PE_oil))
  GPUCB_div <- KL.plugin(GPUCB_oil, averaged_against_GPUCB)
  
  # GPUCB as baseline
  averaged_against_GPUCB_PE <- rowMeans(cbind(rand_oil, exploit_oil, explore_oil, GPUCB_oil))
  GPUCB_PE_div <- KL.plugin(GPUCB_PE_oil, averaged_against_GPUCB_PE)
  
  overall_div <- (rand_div + exploit_div + explore_div + GPUCB_div + GPUCB_PE_div)/5
  return(overall_div)
}

get_landscape_search_results <- function(x, y, z, datagrid, datagrid_dict, optimal_datagrid_dict){
  x <- optimal_datagrid_dict[[as.character(round(x,3))]]
  y <- optimal_datagrid_dict[[as.character(round(y,3))]]
  z <- optimal_datagrid_dict[[as.character(round(z,3))]]
  
  map <- getMap(n=x, m=y,exponent=z, dimx=length(datagrid_dict), dimy=length(datagrid_dict))
  
  rand_history <- run_many_times(datagrid, map, datagrid_dict, searchType='Random')
  rand_oil <- get_avg_results(rand_history, 'Oil')
  
  explore_history <- run_many_times(datagrid, map, datagrid_dict, searchType='Exploration')
  explore_oil <- get_avg_results(explore_history, 'Oil')
  
  exploit_history <- run_many_times(datagrid, map, datagrid_dict, searchType='Exploitation')
  exploit_oil <- get_avg_results(exploit_history, 'Oil')
  
  GPUCB_history <- run_many_times(datagrid, map, datagrid_dict, searchType='GPUCB')
  GPUCB_oil <- get_avg_results(GPUCB_history, 'Oil')
  
  GPUCB_PE_history <- run_many_times(datagrid, map, datagrid_dict, searchType='GPUCB_PE')
  GPUCB_PE_oil <- get_avg_results(GPUCB_PE_history, 'Oil')
  
  overall_div <- get_divergance(rand_oil, exploit_oil, explore_oil, GPUCB_oil, GPUCB_PE_oil)
  return(overall_div)
}

run_gp_optimal <- function(history, datagrid){
  evidence <- unique(history[,c("X", "Y", "Z", "Div")])
  gp_model <- GP_fit(evidence[,c("X", "Y", "Z")], evidence$Div) #evidence[,c("X", "Y")]
  GPredict <- predict(gp_model, datagrid)
  return(GPredict)
}

add_gpucb_pe_point_optimal <- function(history, datagrid, datagrid_dict, optimal_datagrid, current_round, k=1){
  ############
  ### UCB step
  ############
  GPredict_prev <- run_gp_optimal(history, optimal_datagrid)
  
  # find upper confidence bound
  y_plus_sigma  <- GPredict_prev$Y_hat + 2*(GPredict_prev$MSE)**.5
  UCB_point <- GPredict_prev$complete_data[which.max(y_plus_sigma),]
  coords <- as.numeric(UCB_point[1:3])
  new_point <- data.frame(X=coords[1], Y=coords[2], Z=coords[3], Round= current_round +1)	
  new_point$Div <- get_landscape_search_results(coords[1], coords[2], coords[3], datagrid, datagrid_dict, optimal_datagrid_dict)
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
      new_point$Div <- get_landscape_search_results(coords[1], coords[2], coords[3], datagrid, datagrid_dict, optimal_datagrid_dict)
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
    new_point$Div <- get_landscape_search_results(coords[1], coords[2], coords[3], datagrid, datagrid_dict, optimal_datagrid_dict)
    new_point$SearchType <- history$SearchType[1]
    new_point$Stop <- 2
    history <- rbind(history, new_point)
    GPredict_curr <- run_gp_optimal(history, optimal_datagrid)	
    stop <- stopping_criteria(GPredict_prev, GPredict_curr)
    history$Stop[max(history$Round)] <- stop					
    return(history)		
  }
}

initialize_process_optimal <- function(datagrid, datagrid_dict, optimal_datagrid, optimal_datagrid_dict, SearchType, n_init, sobol){
  first_digs <- randomRows(n_init, sobol, optimal_datagrid)
  history <- data.frame(X=first_digs$X, Y=first_digs$Y, Z=first_digs$Z, Round=1:n_init, Div=0)
  for(i in history$Round){
    history[i,]$Div <- get_landscape_search_results(history[i,]$X, history[i,]$Y, history[i,]$Z, datagrid, datagrid_dict, optimal_datagrid_dict)
  }
  history$SearchType <- SearchType        
  history$Stop <- 0
  return(history)
}

run_optimal_search <- function(datagrid, datagrid_dict, optimal_datagrid, optimal_datagrid_dict, failsafe = 30, rho0=0.001){
  current_round <- 4
  history <- initialize_process_optimal(datagrid, datagrid_dict, optimal_datagrid, optimal_datagrid_dict, 'GPUCB-PE', 4 ,sobol = 'Random')
  tick <- 2
  L <- 4
  for(i in 1:L){
    history <- add_gpucb_pe_point_optimal(history, datagrid, datagrid_dict, optimal_datagrid, current_round, k=1)
    current_round <- current_round + tick
  }
  while(current_round < failsafe){
    history <- add_gpucb_pe_point_optimal(history, datagrid, datagrid_dict, optimal_datagrid, current_round, k=1)
    current_round <- current_round + tick
    #print(paste0("Current Round is: ", current_round," and stop == ",sum(1-tail(history,L)$Stop)))
    if(sum(1-tail(history,L)$Stop)<4*rho0){
      return(history)
    }		
  }
  return(history)
}

########################
### Plotting process ###
########################
plot_full_map <- function(datagrid, datagrid_dict, map){
  #fullmap <- total_grid[,c("x1", "x2", "Round")]
  #fullmap$y <- apply(fullmap, 1, getInfoNumber)
  #fullmap$y[fullmap$x1==-4 & fullmap$x2==-4] <- 1
  for(i in 1:length(datagrid_dict)){
    datagrid$Z <- getInfoNumber(datagrid$X[i], datagrid$Y[i], map, datagrid_dict)
  }
  full<-ggplot(datagrid, aes(x=x1,y=x2,color=y,guides=F))+
    geom_point(alpha=.7, size = 2)+xlim(c(-4,4))+ylim(c(-4,4))+
    scale_color_gradientn(colors=jet.colors)+
    ggtitle("True Function")+theme_minimal()+xlab("x1")+ylab("x2")+guides(size = F, alpha=F) +
    theme(plot.title = element_text(hjust = 0.5, size = 11))
  return(list(full))
}  

plot_method_results <- function(rand_data, exploit_data, explore_data, GPUCB_data, GPUCBPE_data, y_lab){
  #x <- c(1:25)
  data <- rbind(rand_data, exploit_data, explore_data, GPUCB_data, GPUCBPE_data)
  #plot(x,rand_data,type="l",col="red", xlab="Round", ylab=y_lab, lwd=3,ylim=c(0, 1))
  #lines(x,exploit_cumregret,col="green4", lwd=3)
  #lines(x,explore_data,col="orange", lwd=3)
  #lines(x,GPUCB_data,col="yellow", lwd=3)
  #lines(x,GPUCBPE_data,col="blue", lwd=3)
  #legend(1, 1, legend=c("Random", "Exploit", "Explore", "GPUCB", "GPUCB-PE"), col=c("red", "green4", "orange", "yellow", "blue"), lty=1, cex=0.8, lwd=3)
  write.csv(data, file = paste(y_lab, "newbest.csv", sep=""))
}