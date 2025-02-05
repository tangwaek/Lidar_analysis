#=========================LIDAR PROCESSING =====================================
# R version 4.3.1
# Authors:  
# Date Jan 2025
#===============================================================================

# Load required libraries
library(here)
library(sf)
library(terra)       # for raster processing
library(parallel)    # check number of core 
library(future)      # parallel processing 
library(lidR)        # processing point clouds
library(ggplot2)     # visualization 
library(dplyr)       # data frame processing

detectCores()
core_num <- detectCores()-8 # change according to file size
plan(multisession, workers = core_num) 

# Step store the lasfile in your R project folder as native R (rds) files

# lasfiles <- here("../../") # go 2 steps back in the file directory
# ctg<- readLAScatalog(paste0(lasfiles,"/las_files"))
# saveRDS(ctg, paste0(lasfiles, "/GIT_Work/Lidar_analysis/Data/lascatalog.rds"))

lasfiles <-here()

field_data <- read.csv(paste0(lasfiles,"/Data/Vol_cal_trees_2005_2015B_ver2.csv")) %>% 
  dplyr::select(c(2,6,8,19,20,52:54,33,35,36))%>% 
  st_as_sf(coords = c("X_tree", "Y_tree"), crs = 2180)

abnormal = function(las){
  if (is(las, "LAS")){
    las = filter_poi(las, Z>=0, Z<=45)
    return(las)
  }
  if (is(las, "LAScatalog")) {
    options = list(
      need_output_file = TRUE,    # Throw an error if no output template is provided
      need_buffer = FALSE)        
    res <- catalog_map(las, abnormal, .options = options)
    return(res)
  }
}

#=========== create a CHM (2D) from 3D point clouds ============================
# x is supposed to be a folder with las files 

ctg_to_chm <- function (ctg){
  opt_select(ctg) <- "xyzir"
  opt_chunk_size(ctg) <- 250
  opt_chunk_buffer(ctg) <- 20
  
  #  Create .lax files for indexing and speedy processing
  lidR:::catalog_laxindex(ctg) 
  
 # save outputs 2 dirctories back beginning with the home directory 
  ras_dir <-  paste0(lasfiles,"../../output/ras_dir") # go 2 step up the root directory and change 
  if(dir.exists(ras_dir)){
    unlink(ras_dir,recursive = TRUE)
  } 
  dir.create(ras_dir, recursive = TRUE) 
  
  opt_output_files(ctg) <- paste0(lasfiles,"../../output/ras_dir/dtm_{XLEFT}_{YBOTTOM}")
  dtm <- rasterize_terrain(ctg, algorithm = tin(), res = 1, overwrite = TRUE)
  
  # normalize point clouds (without DEM) and write out put to folder external drive  
  opt_output_files(ctg) <- paste0(lasfiles,"../../output/ctg_dir/nc_{XLEFT}_{YBOTTOM}")
  
  nc_ctg <-normalize_height(ctg, algorith =tin(), 
                            na.rm =TRUE, res =1)  
  
  # filter and remove abnormal height values 
  opt_output_files(nc_ctg) = paste0(lasfiles, "../../output/filtered/f_{XLEFT}_{YBOTTOM}") 
  filtered_ctg = abnormal(nc_ctg) 
  
  # Create CHM  from normalized point clouds  
  opt_output_files(filtered_ctg) <- paste0(lasfiles,"../../output/ras_dir/chm_{XLEFT}_{YBOTTOM}")
  chm <- rasterize_canopy(filtered_ctg, res = 1,
                          algorithm = p2r(subcircle = 0.2), pkg = "terra")
  
  # writeRaster(chm, paste0(lasfiles,"../output/chm_raw.tif"), overwrite = TRUE)
    
  return(list(filtered_ctg, chm, dtm))
  #  gc() # gabbage collection
}

out_chm <- ctg_to_chm(readRDS(paste0(lasfiles,"/Data/lascatalog.rds")))

col <-height.colors(25)
plot(out_chm[[2]], col = col)
plot(out_chm[[3]])

#st_crs(out_chm[[1]]) <- st_crs(2180) # just in case point cloud is not projected

saveRDS(out_chm, paste0(lasfiles, "/Output/pdt_ctgprocessing.rds")) # save as a light weight R object

#===============================================================================
#           Fill CHM NAs / pits (POST PROCESSING)  
#===============================================================================
# provide a matrix 3x3 matrix that checks if the value in the 5th position is NA
# it then fills NA with the mean value of neighboring cells in the 3x3 matrix. 

fill.na <- function(x, i=5) { 
  if (is.na(x)[i]) {
    pix_val<- mean(x, na.rm = TRUE)
    return(pix_val)
  } 
  else {
    return(x[i])
  }
}

w <- matrix(1, 3, 3) # consider 3 x3  window

filled_chm <- terra::focal(out_chm[[2]], w, fun = fill.na)

writeRaster(filled_chm, file.path(lasfiles,
                                  "/Output/filled_chm.tif"),overwrite = TRUE)

# Smooth CHM using mean filter (3x3)
# to reduce noise due to small variation and enhance visibility of large trees
smoothed <- terra::focal(filled_chm, w, fun = mean, na.rm = TRUE)

# Define a Gaussian filter kernel
gauss_filter <- matrix(c(1, 2, 1, 2, 4, 2, 1, 2, 1), nrow = 3, ncol = 3) / 16

# Apply the Gaussian filter using focal
chm_gauss <- terra::focal(filled_chm, w = gauss_filter, fun = mean)

writeRaster(chm_gauss, paste0(lasfiles,
                              "/Output/chm_gauss.tif"),overwrite = TRUE)

#===============================================================================
#           TREE DETECTION and CROWN SEGMENTATION
#===============================================================================
# determined appropriate window size which defines the search radius for different trees heights. 
# the custom fxn for window mmay also be deter if infor about crown width is known. 

app_ws <- function(h){
  calc_ws <- 2.1 * (-(exp(-0.06*(h-2)) - 1)) + 2.2 
  calc_ws[h <= 7] <- 2.5
  calc_ws[h >= 32] <- 5
  return(calc_ws)
}
#  visualize window size variability with height 
h_range <- seq(0, 43, 5) 
ws <- app_ws(h_range)
plot(h_range, ws, type = "l", ylim = c(0,7))

# ============= Identify individual tree tops [locate_tree()] ======================

# uniqueness () merge similar tree id due to parallel processing, removing duplicates 
# t_tops <- locate_trees(out_chm[[1]], lmf(app_ws), uniqueness = "bitmerge")

t_tops_chm <- locate_trees(chm_gauss, lmf(app_ws, hmin = 7.0)) 

chm_Ras <- as(chm_gauss, "Raster")

# identical(st_crs(t_tops_chm), st_crs(chm_gauss))

plot(chm_gauss, col =col)
plot(t_tops_chm, col = "black", add = TRUE, cex = 0.5)

t_tops_prj <- st_transform(t_tops_chm, st_crs(chm_Ras))

opt_output_files(out_chm[[1]]) <- paste0(lasfiles, 
                                         "../../output/seg_dir/tseg_{XLEFT}_{YBOTTOM}")
trees_watSeg <- segment_trees(out_chm[[1]],
                              algorithm = watershed(chm = out_chm[[2]]))

# check segmentation for one of the catalogiue files  
segmented_las <- readLAS(paste0(lasfiles,"../../output/seg_dir/tseg_378500_401750.las"))
names(segmented_las)
#length(unique(segmented_las$treeID) |> na.omit())
plot(segmented_las, color = "treeID")


metrics_sf  <- list()
#  match field_data to  sf object with matching CRS
field_data <- st_as_sf(field_data, coords = c("X", "Y"), crs = st_crs(trees_watSeg))

#st_write(field_data, paste0(lasfiles,"/Data/trees.shp"))
# Define the path to your LAS files

las_path <- "C:/Users/ElvisTangwa/Documents/GIT_Work/output/seg_dir"
output_path <- paste0(lasfiles, "/Data/matched_trees.shp")

# List all .las files in the directory
las_files <- list.files(path = las_path, pattern = "*.las$", full.names = TRUE)

  


# Save as a shapefile for visualization
st_write(matched_trees, dsn =output_path,delete_layer = TRUE,
         layer_options = "SHPT=POINTZ")

#===============================================================================


# Save as a shapefile for visualization
output_path2b <- paste0(lasfiles, "../../Lidar_analysis/Data/height_matchB.shp")
st_write(matched_trees2b, dsn =output_path2b,delete_layer = TRUE,
         layer_options = "SHPT=POINTZ")

# APPROACH 1====================================================
  library(sf)
  library(dplyr)

match_als_to_reference <- function(las_files, field_data, search_radius = 3, height_threshold = 0.1) {
  
  matched_trees <- list()
  
  for (f in seq_along(las_files)) {
    # Read LAS file safely
    las <- readLAS(las_files[f])
    if (is.null(las)) next  
    
    # Compute tree metrics & convert to sf object
    als_trees <- tree_metrics(las, .stdtreemetrics) %>%
      st_as_sf(coords = c("X", "Y"), crs = st_crs(las))
    
    # Find candidate trees within the specified search radius
    candidate_matches <- st_is_within_distance(als_trees, field_data, dist = search_radius)
    
    matched_list <- list()
    
    for (i in seq_along(candidate_matches)) {
      if (length(candidate_matches[[i]]) > 0) {
        als_tree <- als_trees[i, ]
        ref_candidates <- field_data[candidate_matches[[i]], ]
        
        # Compute height differences
        ref_candidates$height_diff <- abs(als_tree$Z - ref_candidates$height_2007)
        
        # Filter based on height threshold
        ref_candidates <- ref_candidates %>%
          filter((height_diff / height_2007) <= height_threshold)
        
        if (nrow(ref_candidates) > 0) {
          # Select the closest reference tree
          distances <- st_distance(als_tree, ref_candidates)
          closest_ref <- ref_candidates[which.min(distances), ]
          
          matched_list[[i]] <- data.frame(
            matched_ref_id = closest_ref$plot_nr_2015,
            height_ref = closest_ref$height_2007,
            height_diff = closest_ref$height_diff
          )
        }
      }
    }
    
    matched_trees[[f]] <- do.call(rbind, matched_list)
  }
  
  # Combine all matched trees into a single data frame
  results <- do.call(rbind, matched_trees)
  
  return(list(rlts = results, matched_lst = matched_list, dst= distances,
              cls_ref =closest_ref, ref_candi = ref_candidates, candi_mat = candidate_matches))
}

canditate_tree <- match_als_to_reference(las_files, field_data, search_radius = 3, height_threshold =0.1)


# APPROACH 2  (recommended) =========================================
# Load libraries
library(sf)      # For handling spatial data
library(FNN)     # Fast Nearest Neighbor search
library(dplyr)   # Data manipulation

find_matching_trees <- function(las_files, field_data) {
  # Initialize an empty data frame to store results
  results <- data.frame()
  
  # Tree Position Extraction
  for (f in seq_along(las_files)) {
    # Read LAS file safely, skipping invalid files
    las <- readLAS(las_files[f])
    if (is.null(las)) next
    
    # Compute tree metrics & convert to sf object
    als_trees <- tree_metrics(las, .stdtreemetrics) %>%
      st_as_sf(coords = c("X", "Y"), crs = st_crs(las))
    
    # Convert to a matrix for nearest neighbor search
    als_sf <- als_trees %>% st_coordinates()
    reference_trees <- field_data
    
    # Convert to a matrix for nearest neighbor search
    ref_sf <- reference_trees %>% st_coordinates() %>% cbind(field_data$height_2007)
    colnames(ref_sf) <- c("X", "Y", "Z")
    
    # Perform Nearest Neighbor Search within 3 meters
    # "kd_tree" algorithm: Uses a k-d tree for efficient nearest neighbor search in low-dimensional spaces.
    nn_result <- get.knnx(data = ref_sf, query = als_sf, k = 1, algorithm = "kd_tree")
    
    # Filter matches within 3 meters
    valid_matches <- nn_result$nn.dist <= 3
    
    # Extract matched reference tree IDs
    als_trees$matched_ref_id <- ifelse(valid_matches, 
                                       field_data$plot_nr_2015[nn_result$nn.index], NA)
    
    # Compute height differences
    als_trees$height_ref <- ifelse(valid_matches,
                                   reference_trees$height_2007[nn_result$nn.index], NA)
    als_trees$height_diff <- ifelse(valid_matches,
                                    abs(als_trees$Z - als_trees$height_ref), NA)
    
    # Create classification column
    als_trees$classification <- ifelse(valid_matches & (als_trees$height_diff / als_trees$height_ref <= 0.1), 1, 0)
    
    # Select relevant columns for the output
    result <- als_trees %>%
      select(matched_ref_id, height_ref,
             tree_id = treeID, height_diff, classification)
    
    # Append to results
    results <- rbind(results, result)
  }
  
  return(results)
}
candidates <- find_matching_trees(las_files, field_data)

#%>% group_by(tree_id) %>% summarise() %>% ungroup()

# Accuracy Assessment: Errors
# by comparing the number of matched trees against the total number of reference trees. 
# Count matched trees
num_matched <- nrow(ALS_trees_matched)

# Total reference trees
num_reference <- nrow(Reference_trees_sf)

# Detection Rate
detection_rate <- num_matched / num_reference

# Print results
print(paste("Matched Trees:", num_matched))
print(paste("Total Reference Trees:", num_reference))
print(paste("Detection Rate:", round(detection_rate * 100, 2), "%"))

# Omissions (reference trees with no matching ALS tree)

# Commissions (ALS trees without a corresponding reference tree)








# Calculate RMSE and Bias v
rmse_value <- rmse(field_heights, als_heights)
bias_value <- mean(als_heights - field_heights)

# Print accuracy metrics
print(paste("RMSE:", rmse_value))
print(paste("Bias:", bias_value))




#---------METRICS-----------------------------------------
# height of each tree
# crown area of each tree
# competition
# gap size

h_stats <- crown_metrics(las = trees_watSeg, 
                         func = ~list(n = length(Z),
                                      z_max = max(Z),
                                      z_mean = mean(Z)))
head(h_stats)
plot(h_stats["n"], cex=0.5)

std_met <- crown_metrics(las=trees_watSeg  , func = .stdtreemetrics)
std_met    
plot(x = std_met["Z"], cex = 0.5)

crn_area <- function(x, y) {
  #Get x-y coords for each tree and computes subset of points within the convex hull
  ch <- chull(coords)
  # Close coordinates
  ch <- c(ch, ch[1])
  ch_coords <- coords[ch, ]
  # Generate polygon and calculate area
  p <- sf::st_polygon(list(ch_coords))
  area <- sf::st_area(p)
  return(list(A = area))
}
metrics <- crown_metrics(las = trees_watSeg, func = ~crn_area(X, Y))

#====================MATCHING TREES=============================
#   suggestion workflow from copilot 

# Extract tree metrics such as height and crown area
metrics1 <- tree_metrics(las1, .stdtreemetrics)
metrics2 <- tree_metrics(las2, .stdtreemetrics)

# Convert tree metrics to spatial points
coords1 <- SpatialPointsDataFrame(coords = metrics1[, c("X", "Y")], data = metrics1)
coords2 <- SpatialPointsDataFrame(coords = metrics2[, c("X", "Y")], data = metrics2)

# Perform spatial join to match trees based on proximity
matched_trees <- sp::over(coords1, coords2, returnList = TRUE)

# Match trees based on spatial proximity (within a certain radius)
matched_trees <- sp::over(als_data, field_data, returnList = TRUE)

# Extract matched tree attributes
als_heights <- als_data$Height[!is.na(matched_trees)]
field_heights <- field_data$Height[unlist(matched_trees)]

# Calculate RMSE and Bias
rmse_value <- rmse(field_heights, als_heights)
bias_value <- mean(als_heights - field_heights)

# Print accuracy metrics
print(paste("RMSE:", rmse_value))
print(paste("Bias:", bias_value))



# Define a height similarity threshold (e.g., 10% of tree height)
height_threshold <- 0.1

# Function to match trees based on height similarity

match_trees <- function(ref_data, als_data, height_threshold) {
  matches <- list()
  for (i in 1:nrow(ref_data)) {
    ref_tree <- ref_data[i, ]
    potential_matches <- als_data[abs(als_data$Height - ref_tree$Height) / ref_tree$Height <= height_threshold, ]
    if (nrow(potential_matches) > 0) {
      closest_match <- potential_matches[which.min(spDistsN1(as.matrix(potential_matches[, c("X", "Y")]), as.matrix(ref_tree[, c("X", "Y")])), longlat = FALSE), ]
      matches[[i]] <- closest_match
    }
  }
  return(matches)
}

# Apply the matching function
matched_trees <- match_trees(ref_data, als_metrics, height_threshold)







#===========================================
# try using the catalog_apply()
funName <- function(){
  if (is(las,"LASCatalog")) {
    list =(option)
    #identified_trees = 
  }
  return (identified_trees)
}


# 
# # crown_map <- as(crown_map, "SpatRaster")
# 
# identical(st_crs(t_tops_chm), st_crs(out_chm[[1]]))
# 
# out_chm[[1]] <- st_transform(out_chm[[1]], st_crs(chm_gauss))

#--------------------- Alternative Watershed Algorithm  ---------------------------------------

# the function mcws() implements the watershed algorithm

# crowns_ras <- mcws(treetops = t_tops_chm, 
#                    CHM = chm_gauss, minHeight = 5) # previously i used 4
# 
# writeRaster(crowns_ras,"Output/ras_dir/watershed_crowns.tif", overwrite = TRUE)
# # open in QGIS "Raster >> conversion >> polygonize(ras to vec)
# 
# lidrWS_crown <- watershed(chm_gauss, th_tree = 2, tol = 1, ext = 1)()




#  -------------------TREE LEVEL METRICS  -----------------------------------

rfile <- readRDS("C:/Users/ElvisTangwa/Documents/las_files/test_ctg_chm.RData")
ttops <- locate_trees(rfile[[1]], lmf(app_ws, hmin = 5.0), uniqueness = "bitmerge")

opt_output_files(out_chm[[1]]) <- "Output/watershed_seg_algo/seg_chm{*}"
#seg_algo <- segment_trees(out_chm[[1]], watershed(chm_gauss))
t_tops <- locate_trees(out_chm[[1]], lmf(app_ws), uniqueness = "bitmerge")


# segment_trees(las, algorithm, attribute = "treeID", uniqueness = "incremental")


#==============================================================
# interactive viewing and selection of lidar tile from a lascatalog
library(mapview)
library(mapedit)

new_ctg <- catalog_select(lasArea1)

a <- "ras"
if(dir.exists(a)){
  unlink(a,recursive = TRUE)
} 


