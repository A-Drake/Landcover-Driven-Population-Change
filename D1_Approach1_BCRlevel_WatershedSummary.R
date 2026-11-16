--------------------------------------------------------------
# Summarizing landbird data - landcover change only
# BCR, sub-basin, FMU and tenure
# Anna Drake
# Sept 4, 2026
# Script to tally Landcover change only prediction rasters
#-------------------------------------------------------------

library(terra)
library(dplyr)
library(tidyr)
library(purrr)
library(sf)
library(data.table)

# Call Values ------
years <- seq(1985, 2020, by = 5)
nboot <- 32

#--------------------------------------------------------------------------------
# Function to summarize rasters across different land-divisions/conditions
#--------------------------------------------------------------------------------

process_species <- function(spp.i, bcr.i,bcrn.i,root2, wshed.i, tenure.i, watermask.i, FMU.i, distevent.i) {

    
    print(paste0("Running ", spp.i, " ", bcrn.i))

    results <- list()
    
    for (yr in years) {
      
      print(paste0("Year ", yr))
      f <- file.path(root2,paste0(bcrn.i, "_LandChangeOnly_", spp.i, "_", yr, ".tif"))
      
# Check if file exists, else skip (redundant with the year check) ---------------------
   if (!file.exists(f)) { stop(paste("Missing:", f))
   }
      
     # Import if does --------------- 
      r <- rast(f) # prediction outputs
    
# Check that there are 32 bootstrap estimates, else skip -------------
   if (nlyr(r) != nboot) {
   stop(paste0(spp.i, " ", yr, " has ", nlyr(r), " layers rather than ", nboot))
   }
      
      # Crop, mask, and remove water ----------------
      r <- r %>%crop(bcr.i)%>%mask(bcr.i) %>% mask(vect(watermask.i), inverse = TRUE)
      
      # Binary raster for sample size ----------------
      r_valid <- ifel(is.na(r[[1]]), 0, 1) %>% crop(bcr.i)%>%mask(bcr.i)
      
      # ---------------------------------------------------------
      # 1. WHOLE BCR
      # ---------------------------------------------------------
      
      # Sum abundance across all pixels for each bootstrap ---------
      bcr_sum <- global(r, fun = "sum", na.rm = TRUE)
      bcr_n <- global(r_valid, fun = "sum", na.rm = TRUE)
      
      bcr_df <- data.frame(
        spp = spp.i,
        YEAR = yr,
        boot = seq_len(nrow(bcr_sum)),
        spatial_level = "BCR",
        spatial_id = as.character(bcrn.i),
        total_pixels = bcr_n[, 1],
        abundance = bcr_sum[, 1]
      )
      
rm(bcr_sum,bcr_n)
gc()
      print("done BCR stats")
      
      # ---------------------------------------------------------
      # 2. SUBBASIN
      # ---------------------------------------------------------
      if (nrow(wshed.i) > 0) { # conditional on data existing -----
      w_vals <- terra::extract(r, wshed.i, fun = sum, na.rm = TRUE)    
      w_npix <- terra::extract(r_valid, wshed.i, fun = sum, na.rm = TRUE) # sum pixels with valid ----
      sub_map <- data.frame(ID = 1:nrow(wshed.i), spatial_id = as.character(wshed.i$HYBAS_ID))
       
     sub_df <- w_vals %>%
        left_join(sub_map, by = "ID") %>%
        select(-ID) %>%
        pivot_longer(cols = -spatial_id, names_to = "boot_name", values_to = "abundance") %>%
        left_join(
          data.frame(ID = 1:nrow(wshed.i), total_pixels = w_npix[, 2]) %>% 
          left_join(sub_map, by = "ID") %>% select(-ID),
          by = "spatial_id"
        ) %>%
        group_by(spatial_id) %>%
        mutate(spp = spp.i, YEAR = yr, spatial_level = "subbasin", boot = row_number()) %>%
        ungroup() %>%
        select(spp, YEAR, boot, spatial_level, spatial_id, abundance, total_pixels)
        
rm(w_vals,w_npix)
gc()
        
      } else {
  print("No subbasins found in this region. Skipping subbasin stats.")
  sub_df <- data.frame() # Empty placeholder
}
print("done Subbasin stats")

# ---------------------------------------------------------
# 2. Forest Management Unit - split by BCR so in patches
# ---------------------------------------------------------
      if (nrow(FMU.i) > 0) {
      e_vals <- terra::extract(r, FMU.i, fun = sum, na.rm = TRUE)
      e_npix <- terra::extract(r_valid, FMU.i, fun = sum, na.rm = TRUE)# sum pixels with valid ---- 
      fmu_map <- data.frame(ID = 1:nrow(FMU.i), spatial_id = as.character(FMU.i$FMU_Nam))
      
# Aggregate e_vals by spatial_id while keeping all bootstrap columns
val_summary <- e_vals %>%
    left_join(fmu_map, by = "ID") %>%
    select(-ID) %>%
    group_by(spatial_id) %>%
    summarise(across(everything(), \(x) sum(x, na.rm = TRUE)), .groups = "drop")
  
# Aggregate total pixels uniquely by spatial_id
  pixel_summary <- data.frame(ID = 1:nrow(FMU.i), total_pixels = e_npix[, 2]) %>%
    left_join(fmu_map, by = "ID") %>%
    group_by(spatial_id) %>%
    summarise(total_pixels = sum(total_pixels, na.rm = TRUE), .groups = "drop")
    
  FMU_df <- val_summary %>%
    pivot_longer(cols = -spatial_id, names_to = "boot_name", values_to = "abundance") %>%
    left_join(pixel_summary, by = "spatial_id") %>%
    group_by(spatial_id) %>%
    mutate(spp = spp.i, YEAR = yr, spatial_level = "FMU", boot = row_number()) %>%
    ungroup() %>%
    select(spp, YEAR, boot, spatial_level, spatial_id, abundance, total_pixels)

  rm(e_vals, e_npix, val_summary, pixel_summary)
  gc()

 } else {
  print("No FMUs found in this region. Skipping FMU stats.")
  FMU_df <- data.frame() # Empty placeholder
}


print("done Forest Management Unit stats")
      
      
      # ---------------------------------------------------------
      # 3. TENURE
      # ---------------------------------------------------------
      
      tenure_sum <- terra::zonal(r, tenure.i,fun = "sum",na.rm = TRUE)
      
      tenure_npix <- terra::zonal(r_valid, tenure.i,fun = "sum",na.rm = TRUE)%>%
        rename(spatial_id = 1,total_pixels = 2)%>% mutate(spatial_id = as.character(spatial_id))
      
      tenure_df <- tenure_sum %>%
        rename(spatial_id = 1) %>%
        pivot_longer(
          cols = -spatial_id,
          names_to = "boot",
          values_to = "abundance"
        ) %>%
        mutate(
          spp = spp.i,
          YEAR = yr,
          spatial_level = "Tenure",
          spatial_id = as.character(spatial_id),
          boot = rep(seq_len(nlyr(r)), times = nrow(tenure_sum))
        ) %>%
        left_join(
          tenure_npix,
          by = "spatial_id"
        ) %>%
        select(
          spp, 
          YEAR, 
          boot,
          spatial_level, 
          spatial_id,
          abundance, 
          total_pixels
        )
     
     rm(tenure_sum,tenure_npix)
     gc()
      
      print("done Tenure stats")
      
  
# --------------------------------------------------------- 
# Combine the four spatial summaries 
# ---------------------------------------------------------
     
        results[[as.character(yr)]] <- bind_rows(
          bcr_df,
          sub_df,
          FMU_df,
          tenure_df
        )
        
        print(paste0("done ",yr))
    }
    
    final_results <-bind_rows(results) 
    rm(bcr_df, sub_df,FMU_df,tenure_df)
    gc()
    
    # Explicitly return the data frame
    return(final_results)
}

########## end of function #####################3

#---------------------------------------------------
# Get summary data 
#---------------------------------------------------

##---------------------------------------------------------------------------------
#root <- "G:/Shared drives/BAM_NationalModels5"
#root2 <-"C:/Users/andrake/OneDrive - NRCan RNCan/Desktop/Partitioning Change2"
root <-"./data"
root2 <-"./LandChangeOutputs"
root3 <-"./LandChangeSummary"
##---------------------------------------------------------------------------------

# HPC Setup & Args --------------
args <- commandArgs(trailingOnly = TRUE)
task_id <- as.numeric(args[1]) 
#--------------------------------------------
# Extract job list
#-------------------------------------------

#files to process:

jobs<-read.csv(file.path("./jobtracking", "AllJobsLandcover.csv"))
jobs <- jobs %>%
  select(BCR, species) %>%
  distinct()
      
    
#files already processed:
filelist<-list.files(root3) %>% .[. %like% "SummaryStats_"]

tasksdone <-   data.frame(filename =filelist, stringsAsFactors = FALSE) %>%
  tidyr::extract(
    col = filename, 
    into = c("type","BCR","spp"), 
    regex = "^([^_]+)_([^_]+)_([^_]+)\\..*$", 
    remove = TRUE
  ) %>%
  select(BCR, spp) %>%
  distinct()

#files still to do:
remainingtasks <- anti_join(jobs, tasksdone, by = c("BCR","species"="spp"))

#remainingtasks <- read.csv(file.path("./ChangeFigures","MissingInd.csv"))
# Get task --------------
current_task <- remainingtasks[task_id, ]

# Get current species and BCR --------------
#bcrn.i<-current_task$BCR
#spp.i<-current_task$spp
bcrn.i<-current_task$bcr
spp.i<-current_task$species

#----------------------------------------------------------------------
# Check that all prediction years are there? If not, stop the process
#----------------------------------------------------------------------
 
      s <- list.files(root2)%>% .[. %like% "_LandChangeOnly_"] %>%
      .[. %like% spp.i]%>% 
      .[. %like% bcrn.i]
      
      # Check if file exists, else skip ---------------------
      if (length(s)<8) { stop(sprintf("Too few files (%d/8) found for spp: %s, bcr: %s. Aborting job.", 
      length(s), spp.i, bcrn.i))
      }
rm(s)
gc()
     
#--------------------------------------
# Bring in Spatial components
#--------------------------------------

# Get BCR extent --------------------
can <- read_sf(file.path(root,"CAN_adm0.shp")) |> #"Regions", "CAN_adm",
  st_transform(crs=5072) # Canadian boundary

bcr<- read_sf(file.path(root,"BAM_BCR_NationalModel.shp")) |> #, "Regions"
  st_transform(crs=5072)|> st_intersection(can)

bcr$subUnit<-paste0("can",bcr$subUnit)
bcr.i<-bcr%>%filter(subUnit==bcrn.i)

rm(can)
gc()

#-----------------------------------------
# Prepare all layers to this BCR extent 
#-----------------------------------------
# Watermask -----------
watermask.i<-read_sf(file.path(root,"water5072.shp"))%>% #2, "gis"
  st_intersection(bcr.i)%>% st_make_valid() %>% filter(!st_is_empty(.))

# Watersheds --------
wshed.i<-read_sf(file.path(root,"Level8AllCanadaWatersheds.shp"))%>% #2, "gis"
  st_intersection(bcr.i)%>% st_make_valid() %>% filter(!st_is_empty(.))

# FMU --------
FMU.i<-read_sf(file.path(root,"FMU_5072.shp"))%>% 
  st_crop(bcr.i)%>%st_intersection(bcr.i)%>% st_make_valid() %>% filter(!st_is_empty(.))

# Tenure layer ----------
tenure.i<-rast(file.path(root,"LandTenure1k.tif"))%>%crop(bcr.i)%>%mask(bcr.i) %>% mask(vect(watermask.i), inverse = TRUE) #2

# Events ----
distevent.i<-rast(file.path(root,"DisturbanceEvents.tif"))%>%crop(bcr.i)%>%mask(bcr.i) %>% mask(vect(watermask.i), inverse =TRUE) #2

#------------------------------------------
# Now run summary for spp.i and bcr.i:
#------------------------------------------

  out <- process_species(
    spp.i = spp.i,
    bcrn.i=bcrn.i,
    bcr.i = bcr.i,
    root2 = root2,
    wshed.i = wshed.i,
    tenure.i = tenure.i,
    watermask.i = watermask.i,
    FMU.i=FMU.i,
    distevent.i = distevent.i
  )%>% as.data.frame()
 
 # Add columns --------------
  
 out$bcr <- bcrn.i
 out$dens <-out$abundance/out$total_pixels
 
print("done summary...saving....")

#---------------------------------------
# Write out the outputs 
#---------------------------------------  
write.csv(out,file.path(root3,paste0("SummaryStats_", bcrn.i, "_", spp.i, ".csv")),row.names = FALSE)

########## End of Code ###############
