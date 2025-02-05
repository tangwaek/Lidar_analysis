#=========================LIDAR PROCESSING =====================================
# R version 4.3.1
# Authors:  
# Date Dec 2024 
#===============================================================================

# Load required libraries
library(here)
library(sf)
library(terra)       # for raster processing
library(parallel)    # check number of core 
library(future)      # parallel processing 
library(lidR)        # processing point clouds
library(ForestTools) # watershed segmentation
library(ggplot2)
library(dplyr)

detectCores()
core_num <- detectCores()-8 # change according to file file size
plan(multisession, workers = core_num) # process 2 chuncks at a  time. #8L

lasfiles <- here("../../") # this goes 2 steps back in the file directory
#------------------------------------------------------------------------------
#++++++ test with a single las file  --------------------------

las <- readLAS(paste0(lasfiles,"/las_files/M-33-23-D-a-1-3-1-1.las"))

dtm <- rasterize_terrain(las, res = 1, algorithm = tin())

las<- filter_poi(las, Z > 45)

nlas <- las - dtm

# fill NAs  by interpolation 
chm <- rasterize_canopy(nlas, res = 1, algorithm = p2r(0.2, na.fill = tin()))
col <- height.colors(25)
plot(chm, col=col)

# function to fill NAS  based on moving window.
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
schm <- terra::focal(chm, w, fun = fill.na)

ttops <- locate_trees(las = schm, algorithm = lmf(ws = 2.5, hmin = 7.0))

# t_tops_chm <- locate_trees(chm_gauss, lmf(app_ws, hmin = 5.0)) 
ttops

plot(schm, col = col)
plot(ttops, col = "black", add = TRUE, cex = 0.5)

#++++++++++++++++ POINT CLOUD SEGMENTATION++++++++++++++++++++++++++++++++++++
# dalponte segmentation 
las_seg <- segment_trees(las = las, algorithm = dalponte2016(chm = schm, treetops = ttops))

# Count number of trees detected and segmented
length(unique(las$treeID) |> na.omit())

# Visualize all trees
plot(las_seg, color = "treeID")

# Select trees by ID
tree525 <- filter_poi(las = las_seg, treeID == 525)
tree2000 <- filter_poi(las = las_seg, treeID == 2000)

plot(tree2000, size = 4)

# watershed segmentation --- BETTER than Dalponte2016
# noticed that unlike watershed does not need treetop info as input  
trees_watSeg <- segment_trees(las, algorithm = watershed(chm = schm))
trees_watSeg
plot(trees_watSeg, color = "treeID")

tree2000b <- filter_poi(las = trees_watSeg, treeID == 2000)
plot(tree2000b, size = 4)

# compare segemented algorithms ....watershed is better than Dalponte2016
crowns_dalponte <- crown_metrics(las_seg, func = NULL, attribute = "treeID", geom = "concave")
crowns_watSeg <- crown_metrics(trees_watSeg, func = NULL, attribute = "treeID", geom = "concave")

par(mfrow=c(1,2),mar=rep(0,4))
#extract shape of crown allowing the possibility add my layers to the plot 
plot(sf::st_geometry(crowns_dalponte),main="Dalponte2016 segementation", reset = FALSE)
plot(sf::st_geometry(crowns_watSeg),main = "Watershed segementation",  reset = FALSE)    

# --------Raster based segementation ---------------------------------------
# output are raster  layers which not so help for calculating tree metrics. 

chm2 <- as(chm, "Raster") # Dalponte won't work with SpatRaster 

ttops <- st_transform(ttops, st_crs(chm))  
identical(st_crs(chm), st_crs(ttops)) 
trees <- dalponte2016(chm = chm, treetops = ttops)() # chm must be of data type == "RasterLayer"
trees
col1 <- pastel.colors(900)
plot(trees, col=col1)
trees2 <- as(trees, "SpatRaster") # better viz
plot(trees2, col=col1)

t1
#------------Metrics----------------------------------------------------
#  crown analytics
# notice that a function is need because metrics are actually being calculated.

metrics <- crown_metrics(las = trees_watSeg, func = ~list(n = length(Z)))
head(metrics)
plot(metrics["n"], cex=0.5)

# User defined function for area calculation
f <- function(x, y) {
  # Get xy for tree
  coords <- cbind(x, y)
  
  # Convex hull
  ch <- chull(coords)
  
  # Close coordinates
  ch <- c(ch, ch[1])
  ch_coords <- coords[ch, ]
  
  # Generate polygon
  p <- sf::st_polygon(list(ch_coords))
  
  #calculate area
  area <- sf::st_area(p)
  
  return(list(A = area))
}

# Apply user-defined function for convex hull area calculation
metrics <- crown_metrics(las = trees_watSeg, func = ~f(X, Y))
metrics
plot(metrics["A"], cex = 0.5)

metrics <- crown_metrics(las=trees_watSeg  , func = .stdtreemetrics)
metrics  #  height (Z) values are not correct.  
plot(x = metrics["Z"], cex = 0.5)

# Delineating crowns
cvx_hulls <- crown_metrics(las = trees_watSeg, func = .stdtreemetrics, geom = 'convex')
cvx_hulls
plot(cvx_hulls)

z.max_mean <- crown_metrics(trees_watSeg, ~list(z_max = max(Z), z_mean = mean(Z)))
head(z.max_mean) # height needs to be normalized to be correct


int_met <- crown_metrics(trees_watSeg, ~list(imean = mean(Intensity))) # calculate tree intensity metrics
int_met <- int_met[int_met$imean > 80,] # filter intensity

# check for memebership ie treeIDs in trees_watSeg present int_met$treeID
subset <- filter_poi(trees_watSeg, treeID %in% int_met$treeID)
x <- plot(trees_watSeg, bg = "white", size = 4)
plot(subset, add = x + c(-100, 0), size = 5) # some plotting
#=======================================================================

library(ForestGapR)
library(raster)

# ----- Gaps (size and density) influence light availability ------- 

gap_raster <- chm2
gap_raster[gap_raster > 5] <- NA
plot(gap_raster, main = "Canopy Gaps")


threshold <- 10  # Set height threshold for gap detection
size <- c(1, 500)  # Set size range for gaps in square meters
gaps <- getForestGaps(chm = chm2, threshold = threshold, size = size)
plot(gaps, col = "red", add = TRUE)

# set the height thresholds
nthresholds <- c(10, 15, 20, 25)
size <- c(1, 10^4) # m2

# creating an empy raster stack to store multiple gaps as RasterLayers
gaps_stack <- stack()
# Gap detection
for (i in nthresholds) {
  gaps_i <- getForestGaps(chm_layer = chm2, threshold = i, size = size)
  # iteration over all 4 thresholds however, only the last is displayed
  # which justifies stacking 
  names(gaps_i) <- paste0("gaps_", i, "m")
  gaps_stack <- stack(gaps_stack, gaps_i)
}

gap_stats <- GapStats(gap_layer = gaps, chm_layer = chm2)
print(gap_stats)

#  this can be applicable at singe tree level by 
#    gap --> gap_stats,  gap size, distance to nearest gap. 
#    try linking treeID with gapID if this not possible consider 
#   1.   buffer analysis and spatial join (st_intersect, st_within) 
#   2.   nearest neighbor join == library(FNN)

# ++++++++++++++++++++ COMPETITION METRICES +++++++++++++++++++++++++++++++
# crown overlap
# height-based competition
# neighbourhood  density 
# Extract tree positions and heights

# convert 3d data to 2d
tree_data <- ttops
#tree_data <- st_zm(ttops, drop = TRUE, what = "ZM")
coords <- st_coordinates(tree_data)
tree_data$x <- coords[, "X"]; tree_data$y <- coords[, "Y"]
head(tree_data)
# drop spatial component of data 
tree_data <-sf::st_drop_geometry(tree_data[,-1])
# reorder  and then rename columns
tree_data <- tree_data[, c(2,3,1)]
colnames(tree_data)[c(1,2,3)] <-  c("X","Y","Height")

# Neighborhood density (number of neighboring trees within a radius of 5 m distance)
# - 1 means exclude the main tree. 
tree_data$Neighbors <- sapply(1:nrow(tree_data), function(i) {
  sum(sqrt((tree_data$X - tree_data$X[i])^2 + (tree_data$Y - tree_data$Y[i])^2) < 5) - 1
})

# Basal area of neighbors (proxy for competition)
# finds trees within 5 m and sums the height  of trees
tree_data$BAL <- sapply(1:nrow(tree_data), function(i) {
  sum(tree_data$Height[sqrt((tree_data$X - tree_data$X[i])^2 + (tree_data$Y - tree_data$Y[i])^2) < 5])
})

# Plot competition metrics
ggplot(tree_data, aes(X, Y, color = Neighbors)) +
  geom_point(size = 2) +
  scale_color_viridis_c() +
  labs(title = "Tree-Level Competition - Neighborhood density", color = "Neighbors")

# Plot competition metrics
ggplot(tree_data, aes(X, Y, color = BAL)) +
  geom_point(size = 2) +
  scale_color_viridis_c() +
  labs(title = "Tree-Level Competition - Basal Area of Neighbors", color = "BAL")

#------------------------------------------------------
# Computing  crown overlap index (COI)
# crowns_watSeg <- crown_metrics(trees_watSeg, func = NULL, attribute = "treeID", geom = "concave")
# compute pairwise crown overlap

crown_sf <- st_as_sf(crowns_watSeg)
overlap_matrix <- st_intersects(crown_sf, crown_sf)

overlap_areas <- sapply(1:length(overlap_matrix), function(i) {
  sapply(overlap_matrix[[i]], function(j) {
    intersection <- st_intersection(crown_sf[i, ], crown_sf[j, ])
    if (!is.null(intersection)) st_area(intersection) else 0
  })
})

# Example: Converting a list of numeric vectors to a matrix
overlap_matrix <- do.call(rbind, overlap_areas)

dim(overlap_matrix)

crown_areas <- st_area(crown_sf)
crown_overlap_index <- rowSums(overlap_areas) / crown_areas

crown_overlap_index <- rowSums(overlap_matrix) / crown_areas

crown_sf$COI <- crown_overlap_index

# Convert crown areas and overlap values to numeric if needed
crown_sf$COI <- as.numeric(crown_sf$COI)

# Alternatively, if the COI was computed with units:
crown_areas <- as.numeric(crown_areas)  # Convert crown areas to numeric
crown_sf$COI <- rowSums(overlap_matrix) / crown_areas  # Recompute COI as numeric

ggplot() +
  geom_sf(data = crown_sf, aes(fill = COI)) +
  scale_fill_viridis_c() +
  labs(title = "Crown Overlap Index (COI)", fill = "COI")

#----------------------------------------------------------------
# additional metrics 
library(lidRmetrics)

# LASfile <- system.file("extdata", "Megaplot.laz", package="lidR")
# las <- readLAS(LASfile, select = "*", filter = "-keep_random_fraction 0.5")

# you can run any metrics_* function with cloud_metrics()
m1 <- cloud_metrics(trees_watSeg, ~metrics_basic(Z))

# or you can run one of the metric sets in pixel_metrics()

# m2 <- pixel_metrics(trees_watSeg, ~metrics_set2(Z, ReturnNumber, NumberOfReturns), res = 20)

# each metrics_* function has a convenient shortcut to run it with default parameters: 
m3 <- pixel_metrics(las, .metrics_set3, res = 20)



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
#          SCALE UP WORKFLOW
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

abnormal = function(las)
{
  if (is(las, "LAS"))
  {
    las = filter_poi(las, Z>=0, Z<=45)
    return(las)
  }
  
  if (is(las, "LAScatalog"))
  {
    options = list(
      need_output_file = TRUE,    # Throw an error if no output template is provided
      need_buffer = FALSE)         # Throw an error if buffer is 0
    res <- catalog_map(las, abnormal, .options = options)
    return(res)
  }
}

#=========== create a CHM (2D) from 3D point clouds ============================
# x is supposed to be a folder with las files 

ctg_to_chm <- function (x){
  ctg <- readLAScatalog(x, progress = TRUE,
                        select = "xyzir",
                        # filter = "Z>=0, Z<=45",
                        chunk_size = 200,     # 0
                        chunk_buffer = 20)    # 10
  
  #  Create .lax files for indexing and speedy processing
  lidR:::catalog_laxindex(ctg) 
  # st_crs(ctg) <- 2180
  
  # normalize point clouds and write out put to folder external drive  
  opt_output_files(ctg) <- "Output/ctg_dir/nc_{XLEFT}_{YBOTTOM}"
  nc_ctg <-normalize_height(ctg, algorith =tin(), 
                            na.rm =TRUE, res =1)  
  
  # filter and remove abnormal height values 
  opt_output_files(nc_ctg) = "Output/filtered/f_{XLEFT}_{YBOTTOM}" 
  filtered_ctg = abnormal(nc_ctg) 
  
  # Create CHM  from normalized point clouds  
  opt_output_files(filtered_ctg) <- "Output/ras_dir/chm_{XLEFT}_{YBOTTOM}"
  chm <- rasterize_canopy(filtered_ctg, res = 1, 
                          algorithm = p2r(subcircle = 0.2), pkg = "terra", overwrite=TRUE)   
  
  writeRaster(chm, "Output/ras_dir/chm_raw.tif", overwrite = TRUE)
  
  return(list(filtered_ctg, chm))
}

out_chm <- ctg_to_chm("C:/Users/ElvisTangwa/Documents/las_files")

#st_crs(out_chm[[1]]) <- st_crs(2180) # just in case point cloud is not projected

saveRDS(out_chm, file = "test_ctg_chm.RData") # save as a light weight R object

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

writeRaster(filled_chm, "Output/ras_dir/filled_chm.tif", overwrite = TRUE)

# Smooth CHM using mean filter (3x3)
# to reduce noise due to small variation and enhance visibility of large trees
smoothed <- terra::focal(filled_chm, w, fun = mean, na.rm = TRUE)

# Define a Gaussian filter kernel
gauss_filter <- matrix(c(1, 2, 1, 2, 4, 2, 1, 2, 1), nrow = 3, ncol = 3) / 16

# Apply the Gaussian filter using focal
chm_gauss <- terra::focal(filled_chm, w = gauss_filter, fun = mean)

writeRaster(chm_gauss,"Output/ras_dir/chm_gauss.tif", overwrite = TRUE)

#===============================================================================
#           TREE DETECTION and CROWN SEGMENTATION
#===============================================================================
# check out the algorithm of Sterenczak et al;2020
# determined appropriate window size which defines the search radius for different trees heights. 

# Scenario 1---------

app_ws <- function(h){
  calc_ws <- y <- 2.1 * (-(exp(-0.06*(h-2)) - 1)) + 2.2 
  calc_ws[h <= 7] <- 2.5
  calc_ws[h >= 32] <- 5
  
  return(calc_ws)
}
#  visualize window size variability with height 
h_range <- seq(0, 43, 5) 
ws <- app_ws(h_range)
plot(h_range, ws, type = "l", ylim = c(0,7))

# Scenario 2   ----------

f <- function(x) {x * 0.1 + 1.5}
heights <- seq(0,42,5)
ws <- f(heights)
plot(ws)

# ============= Identify individual trees [locate_tree()] ======================

# uniqueness () merge similar tree id due to parallel processing, removing duplicates 
# t_tops <- locate_trees(out_chm[[1]], lmf(app_ws), uniqueness = "bitmerge")

t_tops_chm <- locate_trees(chm_gauss, lmf(app_ws, hmin = 5.0)) 

chm_gauss <- as(chm_gauss, "Raster")

# identical(st_crs(t_tops_chm), st_crs(chm_gauss))

plot(out_chm[[2]])
plot(t_tops_chm, col = "red", add = TRUE, cex = 0.5)

t_tops_chm <- st_transform(t_tops_chm, st_crs(chm_gauss))

crown_map <- dalponte2016(chm_gauss, t_tops_chm, 
                          th_tree = 5, th_seed = 0.40, th_cr = 0.55,
                          max_cr = 15, ID = "treeID")()

writeRaster(crown_map, "crown_dalponte_s14.tif", overwrite=TRUE)

crown_segments <- as(crown_map, "SpatRaster")
#crown_segments <- as.polygons(crown_map)
crown_segments <- st_as_sf(crown_segments)

st_write(crown_segments, "crown_dalponte_s14.shp", append=FALSE)

# apply the segment_tree () CHM based segemntation 
# Segment trees using dalponte algorithm 


#opt_output_files(out_chm[[1]]) = "Output/seg/las_{XLEFT}_{YBOTTOM}" 
las <- segment_trees(out_chm[[2]], 
                     algorithm = dalponte2016(chm = chm_gauss, treetops = t_tops_chm))

#opt_output_files(out_chm[[1]]) = "Output/seg/las_{XLEFT}_{YBOTTOM}" 
las <- segment_trees(las = out_chm[[1]], algorithm = li2012())()

length(unique(las$treeID) |> na.omit())

# crown_map <- as(crown_map, "SpatRaster")

identical(st_crs(t_tops_chm), st_crs(out_chm[[1]]))

out_chm[[1]] <- st_transform(out_chm[[1]], st_crs(chm_gauss))

#--------------------- Alternative Watershed Algorithm  ---------------------------------------

# the function mcws() implements the watershed algorithm

crowns_ras <- mcws(treetops = t_tops_chm, 
                   CHM = chm_gauss, minHeight = 5) # previously i used 4

writeRaster(crowns_ras,"Output/ras_dir/watershed_crowns.tif", overwrite = TRUE)
# open in QGIS "Raster >> conversion >> polygonize(ras to vec)

lidrWS_crown <- watershed(chm_gauss, th_tree = 2, tol = 1, ext = 1)()

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

# Create .lax files for indexing
lidR:::catalog_laxindex(lasArea1)

# Example: Compute a canopy height model
chm <- rasterize_canopy(new_ctg, res = 1, p2r())
plot(chm)

dtm <- rasterize_terrain(new_ctg, algorithm = tin()) # Generate DTM

nc_las <- normalize_height(new_ctg, dtm)  


