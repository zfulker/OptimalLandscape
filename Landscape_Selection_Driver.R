library(GPfit)
library(lhs)
library(entropy)
#library(ggplot2)
#library(scales)
#library(RColorBrewer)
#library(gridExtra)
#library(sensitivity)

source("C:/Users/Owner/Dropbox/GaussianProcessPerlinLandscapes/Tsunami/MLProject/Perlin_landscape_generator.R") #C:/Users/Owner/Dropbox/GaussianProcessPerlinLandscapes/Tsunami/MLProject/
source("C:/Users/Owner/Dropbox/GaussianProcessPerlinLandscapes/Tsunami/MLProject/Optimal_Perlin_helper.R") #C:/Users/Owner/Dropbox/GaussianProcessPerlinLandscapes/Tsunami/MLProject/

# datagrid for search landscape
range=c(1,20)
datagrid <- expand.grid(X=seq(range[1], range[2], 1), Y=seq(range[1], range[2], 1))
datagrid <- as.data.frame(apply(datagrid,2,normalize_col))
datagrid_dict <- c(1:range[2])
names(datagrid_dict) <- round(datagrid$X[1:range[2]], 3)

# datagrid for optimal search
optimal_range=c(2,7)
optimal_datagrid <- expand.grid(X=seq(optimal_range[1], optimal_range[2], 1), Y=seq(optimal_range[1], optimal_range[2], 1), Z=seq(optimal_range[1], optimal_range[2], 1))
optimal_datagrid <- as.data.frame(apply(optimal_datagrid,2,normalize_col))
optimal_datagrid_dict <- c(2:optimal_range[2])
names(optimal_datagrid_dict) <- round(optimal_datagrid$X[1:(optimal_range[2]-1)], 3)

# history <- run_optimal_search(datagrid, datagrid_dict, optimal_datagrid, optimal_datagrid_dict)
# 
history <- read.csv('C:/Users/Owner/Dropbox/GaussianProcessPerlinLandscapes/Tsunami/MLProject/results_2ndtry.csv')
GPredict_final <- run_gp_optimal(history, optimal_datagrid)
experiment_point <- as.data.frame(GPredict_final$complete_data[order(-GPredict_final$Y_hat),])$xnew.3[] #which.max(GPredict_final$Y_hat),
worst_point <- GPredict_final$complete_data[which.min(GPredict_final$Y_hat),]
print(experiment_point)
print(worst_point)
# 
# write.csv(history, file = 'results_longrun.csv')

# plot rand land 
#my_palette <- colorRampPalette(c("aquamarine", "yellow", 'red'))(n = 299)
# map1 <- getMap(n=3, m=2,exponent=6, dimx=range[2], dimy=range[2])
#map1 <- getMap(n=7, m=2,exponent=4, dimx=range[2], dimy=range[2])
#heatmap(map1, Rowv=NA, Colv = NA, labRow = NA, labCol = NA, col = my_palette)
# 
# # currently assuming single given landscape
# map1 <- getMap(n=3, m=2,exponent=6, dimx=range[2], dimy=range[2])
# 
# rand_history <- run_many_times(datagrid, map1, datagrid_dict, searchType='Random')
# rand_oil <- get_avg_results(rand_history, 'Oil')
# rand_regret <- get_avg_results(rand_history, 'Regret')
# 
# explore_history <- run_many_times(datagrid, map1, datagrid_dict, searchType='Exploration')
# explore_oil <- get_avg_results(explore_history, 'Oil')
# explore_regret <- get_avg_results(explore_history, 'Regret')
# 
# exploit_history <- run_many_times(datagrid, map1, datagrid_dict, searchType='Exploitation')
# exploit_oil <- get_avg_results(exploit_history, 'Oil')
# exploit_regret <- get_avg_results(exploit_history, 'Regret')
# 
# GPUCB_history <- run_many_times(datagrid, map1, datagrid_dict, searchType='GPUCB')
# GPUCB_oil <- get_avg_results(GPUCB_history, 'Oil')
# GPUCB_regret <- get_avg_results(GPUCB_history, 'Regret')
# 
# GPUCB_PE_history <- run_many_times(datagrid, map1, datagrid_dict, searchType='GPUCB_PE')
# GPUCB_PE_oil <- get_avg_results(GPUCB_PE_history, 'Oil')
# GPUCB_PE_regret <- get_avg_results(GPUCB_PE_history, 'Regret')
# 
# plot_method_results(rand_oil, exploit_oil, explore_oil, GPUCB_oil, GPUCB_PE_oil, 'Oil')
# 
# plot_method_results(rand_regret, exploit_regret, explore_regret, GPUCB_regret, GPUCB_PE_regret, 'Regret') 

#oil_results <- read.csv('C:/Users/Owner/Dropbox/GaussianProcessPerlinLandscapes/Tsunami/MLProject/Oilnewbest.csv')
#get_divergance(oil_results[1,][2:16], oil_results[2,][2:16], oil_results[3,][2:16], oil_results[4,][2:16], oil_results[5,][2:16])