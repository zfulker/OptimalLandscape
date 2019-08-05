library(GPfit)
library(lhs)
library(entropy)
library(foreach)
library(doParallel)
#library(data.table)
source("C:/Users/Owner/Dropbox/Northeastern/Lab/Search/Code/OptimalLandscape/GP_Models.R")
source("C:/Users/Owner/Dropbox/Northeastern/Lab/Search/Code/OptimalLandscape/Perlin_Landscape_Generator.R")
cores <- detectCores()
registerDoParallel(cores=cores)

########################################
############### Inputs ################# 
########################################
landscape_side_length <- 20
lower_search_para <- 2
upper_search_para <- 7
########################################

# datagrid for search landscape
range<-c(1,landscape_side_length)
datagrid <- expand.grid(X=seq(range[1], range[2], 1), Y=seq(range[1], range[2], 1))
datagrid <- as.data.frame(apply(datagrid,2,normalize_col))
datagrid_dict <- c(1:range[2])
names(datagrid_dict) <- round(datagrid$X[1:range[2]], 3)

# datagrid for optimal landscape parameters
optimal_range=c(lower_search_para,upper_search_para)
optimal_datagrid <- expand.grid(X=seq(optimal_range[1], optimal_range[2], 1), Y=seq(optimal_range[1], optimal_range[2], 1), Z=seq(optimal_range[1], optimal_range[2], 1))
optimal_datagrid <- as.data.frame(apply(optimal_datagrid,2,normalize_col))
optimal_datagrid_dict <- c(optimal_range[1]:optimal_range[2])
names(optimal_datagrid_dict) <- round(optimal_datagrid$X[1:(optimal_range[2]-(optimal_range[1]-1))], 3)

# run model
ptm <- proc.time()
history <- run_optimal_search(datagrid, datagrid_dict, optimal_datagrid, optimal_datagrid_dict)
stopImplicitCluster()
print(proc.time() - ptm)


# results
GPredict_final <- run_gp_optimal(history, optimal_datagrid)
optimal_point <- GPredict_final$complete_data[which.max(GPredict_final$Y_hat),]
worst_point <- GPredict_final$complete_data[which.min(GPredict_final$Y_hat),]
print(paste('Optimal Experiment Point: ',optimal_point[1],optimal_point[2],optimal_point[3]))
print(paste('Least Optimal Experiment Point: ',worst_point[1],worst_point[2],worst_point[3]))

# visualize results
my_palette <- colorRampPalette(c("aquamarine", "yellow", 'red'))(n = 299)
optimal_n = optimal_datagrid_dict[[as.character(round(optimal_point[1],3))]]
optimal_m = optimal_datagrid_dict[[as.character(round(optimal_point[2],3))]]
optimal_e = optimal_datagrid_dict[[as.character(round(optimal_point[3],3))]]
optimal_map <- getMap(n=optimal_n, m=optimal_m,exponent=optimal_e, dimx=range[2], dimy=range[2])
worst_n = optimal_datagrid_dict[[as.character(round(worst_point[1],3))]]
worst_m = optimal_datagrid_dict[[as.character(round(worst_point[2],3))]]
worst_e = optimal_datagrid_dict[[as.character(round(worst_point[3],3))]]
worst_map <- getMap(n=worst_n, m=worst_m,exponent=worst_e, dimx=range[2], dimy=range[2])
heatmap(optimal_map, Rowv=NA, Colv = NA, labRow = NA, labCol = NA, col = my_palette, main = 'Optimal Landscape')
heatmap(worst_map, Rowv=NA, Colv = NA, labRow = NA, labCol = NA, col = my_palette, main = 'Worst Landscape')

# save results 
write.csv(history, file = paste('results_', landscape_side_length**2, '_', lower_search_para, '_', upper_search_para, '_RA'))
