library(GPfit)
library(lhs)
library(entropy)
library(foreach)
library(doParallel)
library(data.table)#/scratch/fulker.z/MLProject/
source("/scratch/fulker.z/MLProject/Risk_GP_Models.R") #C:/Users/Owner/Dropbox/Northeastern/Lab/Search/Code/OptimalLandscape/GP_Models.R 
source("/scratch/fulker.z/MLProject/Perlin_Landscape_Generator.R") #C:/Users/Owner/Dropbox/Northeastern/Lab/Search/Code/OptimalLandscape/Perlin_Landscape_Generator.R
cores <- detectCores()
registerDoParallel(cores=cores)

########################################
############### Inputs ################# 
########################################
landscape_side_length <- 50
num_samp_per_model <- 64
num_digs <- 30
x <- 15
y <- 12
z <- 2
########################################

# datagrid for search landscape
range<-c(1,landscape_side_length)
datagrid <- expand.grid(X=seq(range[1], range[2], 1), Y=seq(range[1], range[2], 1))
datagrid <- as.data.frame(apply(datagrid,2,normalize_col))
datagrid_dict <- c(1:range[2])
names(datagrid_dict) <- round(datagrid$X[1:range[2]], 3)

euc_dist <- function(x1, x2){sqrt(sum((x1 - x2) ^ 2))}

history_to_distances <- function(history, name){
  search_dist_list <- c()
  round_list <- c()
  for(i in sort(unique(history$Iteration))){
    iteration_history <- subset(history, Iteration==i)
    for(j in iteration_history$Round){
      if(j != max(iteration_history$Round)){
        search_dist <- euc_dist(c(iteration_history[j]$X, iteration_history[j]$Y), c(iteration_history[j+1]$X, iteration_history[j+1]$Y))
        search_dist_list <- c(search_dist_list, search_dist)
        round_list <- c(round_list, j)
      }
    }
  }
  search_dist_results <- data.frame(Distances=search_dist_list, Rounds=round_list)
  write.csv(search_dist_results, file = paste('dist_results_', name,'_', landscape_side_length**2, '_', x, '_', y, '_', z,'_RA'))
}

plot_landscape_search_results <- function(x, y, z, datagrid, datagrid_dict, optimal_datagrid_dict, num_samp_per_model, num_digs){
  map <- getMap(n=x, m=y,exponent=z, dimx=landscape_side_length, dimy=landscape_side_length)
  
  GPUCB_Averse_history <- run_many_times(datagrid, map, datagrid_dict, num_samp_per_model, num_digs, searchType='GPUCB-Averse')
  GPUCB_Averse_oil <- get_avg_results(GPUCB_Averse_history, 'Oil')
  write.csv(GPUCB_Averse_oil, file = paste('GPUCB_Averse_oil_', landscape_side_length**2, '_', x, '_', y, '_', z, '_RA'))
  history_to_distances(GPUCB_Averse_history, 'Averse')
  
  GPUCB_Seeking_history <- run_many_times(datagrid, map, datagrid_dict, num_samp_per_model, num_digs, searchType='GPUCB-Seeking')
  GPUCB_Seeking_oil <- get_avg_results(GPUCB_Seeking_history, 'Oil')
  write.csv(GPUCB_Seeking_oil, file = paste('GPUCB_Seeking_oil_', landscape_side_length**2, '_', x, '_', y, '_', z, '_RA'))
  history_to_distances(GPUCB_Seeking_history, 'Seeking')
  
  overall_div <- get_divergance(GPUCB_Averse_oil, GPUCB_Seeking_oil)
  print(overall_div)
}

# Driver
plot_landscape_search_results(x, y, z, datagrid, datagrid_dict, optimal_datagrid_dict, num_samp_per_model, num_digs)
stopImplicitCluster()