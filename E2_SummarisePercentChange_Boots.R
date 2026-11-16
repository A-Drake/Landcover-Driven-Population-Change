#-----------------------------
# Organize outputs of Landcover only change analysis
# assessing abundance change
# Anna Drake
# Sept 10,2026 - pulling from mosaics - will have to cut to BCR extents afterwards
#-----------------------------

library(sf)
library(terra)
require(data.table)
library(tidyverse) # includes dplyr, ggplot2, tidyr, and purrr (for map_df)

`%notin%`<-negate(`%in%`)

# --------------------------------------
# Directories
# ---------------------------------------

root <- "./data" 
root2 <- "./LandChangeSummary"
root3 <- "./LandChangeOutputs/TruncatedMosaics"
out<- "./LandChangeFigures"

#------------------------------------
# HPC Setup & Args 
# -----------------------------------

args <- commandArgs(trailingOnly = TRUE)
task_id <- as.numeric(args[1]) 

# Get BCR jobs --------------------
bcrlist<-c("can5","can9","can10","can12","can13","can14",
           "can40","can41","can42","can60","can61","can70",
           "can71","can72","can80","can81","can82","ALL")
           
# Validate index range ------------------
if (is.na(task_id) || task_id < 1 || task_id > length(bcrlist)) {
  stop(sprintf("Invalid task_id: %s. Must be between 1 and %d.", args[1], length(bcrlist)))
}



bcrn.i <- bcrlist[task_id]

#bcrn.i <- "can9"

# Get V5 models deemed to be too poor to use based on AUC -----------------
removeSp<-read.csv(file.path(root,"SpeciesRemoved.csv"))
removeSp<-removeSp%>%subset(.,.$region==bcrn.i)%>%pull(spp)

# Species where the indicator model was poor ---------

PosthocID<-data.frame(spp=c("BOWA","BTBW","MAWA","WBNU","WWCR","YBCU"), BCR=c("can5","can80","can40","can81","can5","can13"))
PosthocID<-PosthocID%>%subset(.,.$BCR==bcrn.i)%>%pull(spp)

# Species with extreme predictions within the 99% ("raster stacks summary) ------------
Extreme<-data.frame(spp=c("CLSW","LISP","GWWA"), BCR=c("can14","can14","can13"))
Extreme<-Extreme%>%subset(.,.$BCR==bcrn.i)%>%pull(spp)

# Filter out nomadic breeding species - driving odd results - 
irruptive<-c("WWCR","RECR","PISI","EVGR")

# Final removal list ---------
removeSp<-c(removeSp,irruptive,PosthocID,Extreme)%>%unique()

#----------------------------------
# Indicators for region
#----------------------------------
Indicators<-read.csv(file.path(root,"IndicatorSpecies_Drake_etal_CJFR.csv")) %>%
subset(.,.$Class!=".") %>%
subset(.,.$scale=="200m" & .$BCR==bcrn.i & .$S>0.4055)%>% 
subset(.,.$species %notin% removeSp)

if(bcrn.i!="ALL"){
# Age -----------------------
AgeI<-Indicators%>%subset(.,.$Group=="ForestAge")#.$S>0.32)lowering threshold from double to 1.4 to see

V_Old<-subset(AgeI,AgeI$Class %in% c("very old")) %>% pull(species) %>% unique()
Old<-subset(AgeI,AgeI$Class %in% c("very old","old","v-old/old")) %>% pull(species) %>% unique()
Mature<-subset(AgeI,AgeI$Class %in% c("mature")) %>% pull(species) %>% unique()
Young<-subset(AgeI,AgeI$Class %in% c("young")) %>% pull(species) %>% unique()

# Composition --------------
ConI<-Indicators%>%subset(.,.$Group=="PropConifer")#.$S>0.32)lowering threshold from double to 1.4 to see

Conifer<-subset(ConI,ConI$Class %in% c("conifer")) %>% pull(species) %>% unique()
Mixedwood<-subset(ConI,ConI$Class %in% c("mixedwood")) %>% pull(species) %>% unique()
Deciduous<-subset(ConI,ConI$Class %in% c("deciduous")) %>% pull(species) %>% unique()

# Cover -----------------------------
CovI<-Indicators%>%subset(.,.$Group=="CanopyCover")

Open<-subset(CovI,CovI$Class %in% c("nonforest")) %>% pull(species) %>% unique()
Forest<-subset(CovI,CovI$Class %in% c("forest")) %>% pull(species) %>% unique()

# Stack them all -------------
Stack<-list(V_Old,Old,Mature,Conifer,Mixedwood,Deciduous,Open,Forest)
names<-c("Very old forest","Old forest","Mature forest","Conifer forest","Mixedwood forest","Deciduous forest","Low Cover","All Forest")  #removed "Young forest", hardly any species in this group

print("Finished Stacks")


SaveOut <- list(V_Old = V_Old,Old = Old,Mature = Mature,Conifer = Conifer, Mixedwood = Mixedwood, Deciduous = Deciduous,Open = Open,Forest = Forest)

SaveOutdf<- as.data.frame(
  lapply(SaveOut, `length<-`, max(lengths(SaveOut)))
)
print("save indicators")
write.csv(SaveOutdf, file.path(out,paste0(bcrn.i,"_Indicator_Lists.csv")), row.names = FALSE, na = "")
}


if(bcrn.i=="ALL"){
# All species in outputs ----------
All<-list.files(root3)%>% substr(.,1,4) %>%unique() %>% .[. %notin% removeSp]
Stack<-list(All)
names<-c("AllGroupedSpp")
}


###################################################################################
for (j in c(1:length(Stack))){ 
    name1<-names[j]
    
    if (is.null(Stack[[j]]) || length(Stack[[j]]) < 3) {
      next(print(paste0("Skipping:", name1, "is empty or has <3 species")))
    }
    
    codes_pattern <- paste(Stack[[j]], collapse = "|")
    
    list1<-list.files(file.path(root2))%>%.[. %like% "SummaryStatsMosaic_"]%>%
      .[. %like% codes_pattern] 
    
if (length(Stack[[j]]) > length(list1)) {
      print(paste0("Skipping: ", name1, 
                   " | Expected species: ", length(Stack[[j]]), 
                   " | Found summary files: ", length(list1)))
      next
    }
    
    n<-length(list1) # number imported
    
    df<-map_df(file.path(root2, list1), ~ read_csv(.x, show_col_types = FALSE))
    print(paste0("imported ", name1))
    
    diff_df <- df %>%
      pivot_wider(
        names_from = YEAR, 
        values_from = dens, 
        id_cols = c(spp, boot, spatial_level, spatial_id,total_pixels),
        names_prefix = "yr_"
      ) %>% mutate(per_diff = yr_2020-yr_1990)
    
cIdfb <- diff_df %>%
  group_by(boot, spatial_level, spatial_id) %>%
  summarize(
    mds = sum(per_diff, na.rm = TRUE), # total bird change for species group
    y90 = sum(yr_1990, na.rm = TRUE),
    total_pixels = first(total_pixels), # 
    .groups = "drop"
  ) %>%
  group_by(spatial_level, spatial_id) %>%
  summarize(
    total_pixels = mean(total_pixels, na.rm = TRUE),
    m_sds = mean(mds, na.rm = TRUE),
    dsLCI = quantile(mds, 0.025, na.rm = TRUE),
    dsUCI = quantile(mds, 0.975, na.rm = TRUE),
    m_y90 = mean(y90, na.rm = TRUE),
    y90LCI = quantile(y90, 0.025, na.rm = TRUE),
    y90UCI = quantile(y90, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

    cIdfb <- cIdfb %>%
    mutate(
  #  mpCH = (exp(m_mds) - 1) * 100,
  #  mpLCI = (exp(dsLCI) - 1) * 100,
  #  mpUCI = (exp(dsUCI) - 1) * 100,
    sig  = if_else(dsLCI < 0 & dsUCI > 0, 1, 0),
    nspp = n
      )

    missed<-length(Stack[[j]])-length(list1)
    print(paste0("Missing", missed," species in ",name1))
  
    write.csv(cIdfb, file.path(out,paste0("IndicatorMeanDensityChange_",bcrn.i,"_",name1,".csv")))
    
    print(paste0("done", name1))
  } # end of groups
  
  # End of Job -----------------------------------
