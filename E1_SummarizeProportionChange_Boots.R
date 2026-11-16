#-----------------------------
# Calculate significance of summary stats
# Landcover only change analysis
# Anna Drake
# Sept 10,2026 - modified to work with mosaics
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
root3 <- "./LandChangeOutputs/TruncatedMosaic/"
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

#fix part of map cut off           
#bcrlist<-c("can9")

# Validate index range
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

write.csv(SaveOutdf, file.path(out,paste0(bcrn.i,"_Indicator_Lists.csv")), row.names = FALSE, na = "")
}


if(bcrn.i=="ALL"){
# All species in outputs ----------
All<-list.files(root3)%>% substr(.,1,4) %>%unique() %>% .[. %notin% removeSp]
Stack<-list(All)
names<-c("AllGroupedSpp")
}
print("save indicators")

#########################################################
for (j in c(1:length(Stack))){

name1<-names[j]
print(name1)

if (is.null(Stack[[j]]) || length(Stack[[j]]) < 3) {
next(print(paste0("Skipping:", name1, "is empty or has <3 species")))
 }

codes_pattern <- paste(Stack[[j]], collapse = "|")

list1<-list.files(file.path(root2))%>%.[. %like% "SummaryStatsMosaic_"]%>%
  .[. %like% codes_pattern]

print("imported sumlist")

if (length(Stack[[j]]) > length(list1)) {
      print(paste0("Skipping: ", name1, 
                   " | Expected species: ", length(Stack[[j]]), 
                   " | Found summary files: ", length(list1)))
      next
    }
     
df<-map_df(file.path(root2, list1), ~ read_csv(.x, show_col_types = FALSE))
print(paste0("imported ", name1))

# Get comparison years -----------------------

diff_df <- df %>%
  pivot_wider(
    names_from = YEAR, 
    values_from = dens, 
    id_cols = c(spp, boot, spatial_level, spatial_id,total_pixels),
    names_prefix = "yr_"
  ) %>% mutate(dens_diff = yr_2020 - yr_1990)

cIdfb<-diff_df%>%group_by(spp,spatial_level,spatial_id,total_pixels)%>%
  summarize(mds=mean(dens_diff),
            dsUCI=quantile(dens_diff, probs=0.975, na.rm=T),
            dsLCI=quantile(dens_diff, probs=0.025, na.rm=T)
  )%>%ungroup()
  
cIdfb <- cIdfb %>%
  mutate(
    inc = if_else(dsLCI > 0, 1, 0),
    dec = if_else(dsUCI < 0, 1, 0),
    ns  = if_else(dsLCI <= 0 & dsUCI >= 0, 1, 0)
  )
  
missed<-length(Stack[[j]])-length(list1)
print(paste0("Missing", missed," species in ",name1))
n<-length(unique(cIdfb$spp))

tally <- cIdfb %>% group_by(spatial_level, spatial_id, total_pixels)%>%
  summarize(inc=sum(inc, na.rm=T),dec=sum(dec,na.rm=T),ns=sum(ns, na.rm=T), CIr=mean(dsUCI-dsLCI, na.rm=T),dCI=mean(dsUCI-mds, na.rm=T))%>%
  mutate(
    net_score = inc - dec,
    diffInd = (inc + (0.5*ns))/(inc + ns + dec),
    nspp =n,
    trend_strength = case_when(
    
      # Neutral states ------------
      inc == 0 & dec == 0 & ns> 0  ~ "Stable",
      inc == 0 & dec == 0 & ns==0 ~ "NA",
      inc == dec ~ "Equal loss/gain",
      
      # Gain tiers (positive net)
      net_score >= n*0.5 ~ "Gain dominated (Strong)",
      net_score >= n*0.25 ~ "Gain dominated (Moderate)",
      net_score > 0  ~ "Gain dominated (Weak)",
      
      # Loss tiers (negative net)
      net_score <= n*-0.5 ~ "Loss dominated (Strong)",
      net_score <= n*-0.25 ~ "Loss dominated (Moderate)",
      net_score < 0  ~ "Loss dominated (Weak)"
    ),
    
    # Set a logical factor level order for mapping/legends
    trend_strength = factor(
      trend_strength, 
      levels = c( "Loss dominated (Strong)",  "Loss dominated (Moderate)", "Loss dominated (Weak)", 
                  "Stable", "Equal loss/gain", 
                  "Gain dominated (Weak)",  "Gain dominated (Moderate)",  "Gain dominated (Strong)")
    )
  )
gc()


write.csv(tally, file.path(out,paste0("IndicatorChange_",bcrn.i,"_",name1,".csv")))

print(paste0("done", name1))
} # end of groups

# End of Job -----------------------------------
